import numpy as np
from collections import namedtuple
from csbdeep.utils import _raise, load_json, save_json, move_image_axes, normalize
from pathlib import Path

from tifffile import imread, TiffFile
from csbdeep.io import save_tiff_imagej_compatible
import imreg_dft as ird

import keras.backend as K
from stardist import export_imagej_rois
from stardist.models import StarDist2D

try:
    from tqdm.notebook import tqdm as tqdm_notebook
except ModuleNotFoundError:
    from tqdm import tqdm_notebook
# from tqdm import tqdm
tqdm = tqdm_notebook



_config = namedtuple('Config',(
    'base_dir',
    'data_dir',
    'registered_dir',
    'results_dir',
    'model_dir',

    'channel_order',

    'membrane_model',
    'membrane_prob_thresh',
    'membrane_nms_thresh',
    'dna_model',
    'dna_prob_thresh',
    'dna_nms_thresh',

    'channels_segment',
    'channel_track',

    'channel_drift_correction',
    'do_curation',
    'export_xlsx_file',
    'frame_interval_seconds',
    'pixel_size_um',
));



class Config(_config):

    def save(self, path):
        save_json(self._asdict(), path)

    @staticmethod
    def load(path):
        return Config(**load_json(path))



class Starchaea:

    def __init__(self, config):
        if isinstance(config,str):
            self.config = Config.load(config)
        else:
            # isinstance(config, Config) or _raise(ValueError('not a config'))
            self.config = config
        self._check_config()


    def _check_config(self):
        # be as through as you want, or simply hope for the best
        c = self.config
        channels_allowed = set(('membrane','dna'))
        assert Path(c.base_dir).exists()
        assert 1 <= len(c.channel_order) <= 2 and channels_allowed.union(set(c.channel_order)) == channels_allowed
        assert c.channel_drift_correction is None or c.channel_drift_correction in channels_allowed
        assert 1 <= len(c.channels_segment) <= 2 and channels_allowed.union(set(c.channels_segment)) == channels_allowed
        assert c.channel_track in channels_allowed


    def init(self):
        c = self.config
        self.base_dir = Path(c.base_dir)
        self.raw_file_path = self.base_dir  / c.data_dir
        self.raw_files = sorted(self.raw_file_path.rglob('*.tif*'))
        print(f'Base folder is {self.base_dir}')

        n_raw_files = len(self.raw_files)
        print(f'There are {n_raw_files} datasets to analyse :-)')

        if c.channel_drift_correction is not None:
            self.registered_dir = self.base_dir / c.registered_dir
            self.registered_dir.mkdir(exist_ok=True, parents=True)

        self.stardist_dir = self.base_dir / c.results_dir
        self.stardist_dir.mkdir(exist_ok=True, parents=True)
        for name in c.channels_segment:
            (self.stardist_dir / name).mkdir(exist_ok=True, parents=True)

        print('Successfully created save directories, yay!')


    def load_models(self):
        c = self.config
        d = c._asdict()
        self.models = {}
        K.clear_session()
        for name in c.channels_segment:
            self.models[name] = dict (
                model       = StarDist2D(None, name=d[name+'_model'], basedir=c.model_dir),
                prob_thresh = d[name+'_prob_thresh'],
                nms_thresh  = d[name+'_nms_thresh'],
            )


    def load_timelapse(self, file):
        c = self.config
        print(f'Loading image from {file}')
        T = imread(str(file))

        if T.ndim==3:
            axes = 'TYX'
        elif T.ndim==4:
            axes = 'TCYX'
        else:
            raise ValueError('Image shape has incorrect number of dimensions')

        # normalize image to always have channel axis
        T = move_image_axes(T, axes, 'TCYX', adjust_singletons=True)

        print(f'Data has axes {axes} with shape {T.shape}')

        if c.channel_drift_correction is not None:
            return self.drift_correction(T, file)
        else:
            return T#, axes


    def register(self, T, reg_ch):

        if T.ndim==3:
            reg_ch = None

        def _reg(x):
            return x if reg_ch is None else x[reg_ch]

        R = [T[0]]

        print('Running drift correction...')

        for frame in tqdm(T[1:]):
            result = ird.translation(_reg(R[-1]), _reg(frame))
            if reg_ch is None:
                freg = ird.transform_img(frame, tvec=result["tvec"])
            else:
                freg = np.stack([ird.transform_img(c, tvec=result["tvec"]) for c in frame])
            R.append(freg)

        reg = np.stack(R)

        return reg


    def drift_correction(self, T, file):
        c = self.config
        assert c.channel_drift_correction is not None

        # check if file already exists?
        reg_file = self.registered_dir / ('DRIFTCORRECTED_' + file.name)

        reg_ch = c.channel_order.index(c.channel_drift_correction)
        T_reg = self.register(T, reg_ch)
        T_reg = T_reg.astype(T.dtype)

        # with TiffFile(str(file)) as _file:
        #     imagej_metadata = _file.imagej_metadata
        #     ome_metadata = _file.ome_metadata
        # save_tiff_imagej_compatible(str(reg_file), T_reg, axes=axes, metadata=imagej_metadata)

        save_tiff_imagej_compatible(str(reg_file), T_reg, axes='TCYX')

        return T_reg


    def _predict_stardist(self, model, file, T, channel, prob_thresh, nms_thresh, out_dir):

        axes = 'TCYX'
        # if T.ndim==3:
        #     timelapse = T
        if T.ndim==4:
            timelapse = T[:,channel]
        else:
            raise ValueError('Data has unexpected number of dimensions. Weird.')

        # normalise
        print(f'Normalizing each frame to run Stardist', flush=True)
        timelapse = np.stack([normalize(frame, 1,99.8) for frame in timelapse])
        print(f"Timelapse has axes {axes.replace('C','')} with shape {timelapse.shape}")

        polygons = [model.predict_instances(frame, nms_thresh=nms_thresh, prob_thresh=prob_thresh)[1] for frame in tqdm(timelapse)]

        if prob_thresh is None:
            prob_string = 'default'
        else:
            prob_string = f'{prob_thresh:.2f}'

        if nms_thresh is None:
            nms_string  = 'default'
        else:
            nms_string = f'{nms_thresh:.2f}'

        roi_path = out_dir / f"{file.stem}_prob={prob_string}_nms={nms_string}"
        roi_path.parent.mkdir(parents=True, exist_ok=True)
        rois_python = Path(str(roi_path)+'.npz')
        rois_imagej = Path(str(roi_path)+'.zip')

        print(f'Saving ImageJ ROIs to {rois_imagej}')
        export_imagej_rois(str(rois_imagej), [poly['coord'] for poly in polygons])

        print(f'Saving Python rois to {rois_python}')
        np.savez(str(rois_python),
            coord  = [p['coord']  for p in polygons],
            points = [p['points'] for p in polygons],
            prob   = [p['prob']   for p in polygons],
        )


    def predict_stardist(self, file, T):
        c = self.config
        for channel, model in self.models.items():
            stardist_model = model['model']
            channel_ind = c.channel_order.index(channel)
            prob_thresh = model['prob_thresh']
            nms_thresh = model['nms_thresh']
            out_dir = self.stardist_dir / channel
            print(f'\n~~ Running predictions on channel {channel} ~~')
            self._predict_stardist(stardist_model, file, T, channel_ind, prob_thresh, nms_thresh, out_dir)

