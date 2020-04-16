# adapted from https://forum.image.sc/t/trackmate-with-manual-detection-and-automatic-tracking-in-jython-script/29164/7

#@ String (visibility=MESSAGE, value="<html><b>Assumptions:</b><br/><ul><li>Channel order if multi-channel: Membrane, DNA.</ul><br></html>", required="false") doc

#@ CommandService command
#@ UIService ui
#@ RoiManager rm

#@ String (label="Track Channel", choices={"Membrane", "DNA"}, style="radioButtonHorizontal") tracking_channel
#@ File (label="Output Directory", style="directory") save_dir

##@ String (visibility=MESSAGE, value="<html><br></html>", required="false") spacer1

#@ ImagePlus imp

#@ String (visibility=MESSAGE, label="<html><br/><b>Membrane StarDist Model</b></html>", value="<html><br/><hr width='100'></html>", required="false") membrane_msg
#@ Boolean (label="Enable", value="true") stardist_membrane_enabled
#@ File (label="Model File", required="false") stardist_membrane
#@ Float (label="Probability Threshold", stepSize="0.05", min="0", max="1", style="slider", value="0.844602", persist="false") prob_thresh_membrane
#@ Float (label="Overlap Threshold", stepSize="0.05", min="0", max="1", style="slider", value="0.7", persist="false") nms_thresh_membrane

#@ String (visibility=MESSAGE, label="<html><br/><b>DNA StarDist Model</b></html>", value="<html><br/><hr width='100'></html>", required="false") dna_msg
#@ Boolean (label="Enable", value="true") stardist_dna_enabled
#@ File (label="Model File", required="false") stardist_dna
#@ Float (label="Probability Threshold", stepSize="0.05", min="0", max="1", style="slider", value="0.696284", persist="false") prob_thresh_dna
#@ Float (label="Overlap Threshold", stepSize="0.05", min="0", max="1", style="slider", value="0.7", persist="false") nms_thresh_dna

#@ String (visibility=MESSAGE, label="<html><br/><b>Tracking parameters</b></html>", value="<html><br/><hr width='100'></html>", required="false") tracking_msg
#@ Float (label="Frame to frame linking max. distance", stepSize="0.5", min="0", max="20", style="slider", value="10.0", persist="false") frame_link_dist
#@ Float (label="Gap closing max. distance", stepSize="0.5", min="0", max="20", style="slider", value="15.0", persist="false") gap_close_dist
#@ Float (label="Segment splitting max. distance", stepSize="0.5", min="0", max="20", style="slider", value="7.0", persist="false") seg_split_dist


import sys
from math import pi
from math import sqrt
from random import shuffle
import os

from java.awt import Color
from java.io import FileWriter
from com.google.gson import Gson

# from ij import WindowManager
from ij.measure import ResultsTable, Measurements
from ij.plugin import ChannelSplitter
# from ij.plugin.frame import RoiManager
from ij.plugin.filter import Analyzer

from fiji.plugin.trackmate import Logger
from fiji.plugin.trackmate import Model
from fiji.plugin.trackmate import SelectionModel
from fiji.plugin.trackmate import Settings
from fiji.plugin.trackmate import Spot
from fiji.plugin.trackmate import SpotCollection
from fiji.plugin.trackmate import TrackMate
from fiji.plugin.trackmate.detection import ManualDetectorFactory
from fiji.plugin.trackmate.tracking import LAPUtils
from fiji.plugin.trackmate.providers import SpotAnalyzerProvider
from fiji.plugin.trackmate.providers import EdgeAnalyzerProvider
from fiji.plugin.trackmate.providers import TrackAnalyzerProvider
from fiji.plugin.trackmate.tracking.sparselap import SparseLAPTrackerFactory
from fiji.plugin.trackmate.visualization.hyperstack import HyperStackDisplayer
from fiji.plugin.trackmate.gui import TrackMateGUIController
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter
from org.jfree.chart.renderer.InterpolatePaintScale import Jet

from org.scijava.ui import DialogPrompt
from de.csbdresden.stardist import StarDist2D



def spots_from_results_table( results_table, frame_interval ):
	"""
	Creates a spot collection from a results table in ImageJ.
	Requires the current results table, in which the results from
	particle analysis should be. We need at least the center
	of mass, the area and the slice to be specified there.
	We also query the frame interval to properly generate the
	POSITION_T spot feature.
	"""

	# hyperstack
	frames = results_table.getColumnAsDoubles( results_table.getColumnIndex( 'Frame' ) )
	# stack
	if frames is None:
		frames = results_table.getColumnAsDoubles( results_table.getColumnIndex( 'Slice' ) )

	xs = results_table.getColumnAsDoubles( results_table.getColumnIndex( 'XM' ) )
	ys = results_table.getColumnAsDoubles( results_table.getColumnIndex( 'YM' ) )
	z = 0.
	# Get radiuses from area.
	areas = results_table.getColumnAsDoubles( results_table.getColumnIndex( 'Area' ) )
	spots = SpotCollection()

	for i in range( len( xs ) ):
		x = xs[ i ]
		y = ys[ i ]
		frame = frames[ i ]
		area = areas[ i ]
		t = ( frame - 1 ) * frame_interval
		radius = sqrt( area / pi )
		quality = i # Store the line index, to later retrieve the ROI.
		spot = Spot( x, y, z, radius, quality )
		spot.putFeature( 'POSITION_T', t )
		spots.add( spot, int( frame - 1 ) )

	return spots


def create_trackmate( imp, results_table, frame_link_dist, gap_close_dist, seg_split_dist ):
	"""
	Creates a TrackMate instance configured to operated on the specified
	ImagePlus imp with cell analysis stored in the specified ResultsTable
	results_table.
	"""

	cal = imp.getCalibration()

	# TrackMate.

	# Model.
	model = Model()
	model.setLogger( Logger.IJ_LOGGER )
	model.setPhysicalUnits( cal.getUnit(), cal.getTimeUnit() )

	# Settings.
	settings = Settings()
	settings.setFrom( imp )

	# Create the TrackMate instance.
	trackmate = TrackMate( model, settings )

	# Add ALL the feature analyzers known to TrackMate, via
	# providers.
	# They offer automatic analyzer detection, so all the
	# available feature analyzers will be added.
	# Some won't make sense on the binary image (e.g. contrast)
	# but nevermind.

	spotAnalyzerProvider = SpotAnalyzerProvider()
	for key in spotAnalyzerProvider.getKeys():
		print( key )
		settings.addSpotAnalyzerFactory( spotAnalyzerProvider.getFactory( key ) )

	edgeAnalyzerProvider = EdgeAnalyzerProvider()
	for  key in edgeAnalyzerProvider.getKeys():
		print( key )
		settings.addEdgeAnalyzer( edgeAnalyzerProvider.getFactory( key ) )

	trackAnalyzerProvider = TrackAnalyzerProvider()
	for key in trackAnalyzerProvider.getKeys():
		print( key )
		settings.addTrackAnalyzer( trackAnalyzerProvider.getFactory( key ) )

	trackmate.getModel().getLogger().log( settings.toStringFeatureAnalyzersInfo() )
	trackmate.computeSpotFeatures( True )
	trackmate.computeEdgeFeatures( True )
	trackmate.computeTrackFeatures( True )

	# Skip detection and get spots from results table.
	spots = spots_from_results_table( results_table, cal.frameInterval )
	model.setSpots( spots, False )

	# Configure detector. We put nothing here, since we already have the spots
	# from previous step.
	settings.detectorFactory = ManualDetectorFactory()
	settings.detectorSettings = {}
	settings.detectorSettings[ 'RADIUS' ] = 1.

	# Configure tracker
	settings.trackerFactory = SparseLAPTrackerFactory()
	settings.trackerSettings = LAPUtils.getDefaultLAPSettingsMap()
	settings.trackerSettings[ 'LINKING_MAX_DISTANCE' ] 		= frame_link_dist
	settings.trackerSettings[ 'GAP_CLOSING_MAX_DISTANCE' ]	= gap_close_dist
	settings.trackerSettings[ 'MAX_FRAME_GAP' ]				= 3
	settings.trackerSettings[ 'ALLOW_TRACK_SPLITTING' ]		= True
	settings.trackerSettings[ 'SPLITTING_MAX_DISTANCE' ]	= seg_split_dist

	settings.trackerSettings

	settings.initialSpotFilterValue = -1.

	### print(model.getFeatureModel().getTrackFeatureNames())
	# TRACK_START: Track start,
	# TRACK_INDEX: Track index,
	# NUMBER_MERGES: Number of merge events,
	# TRACK_STD_SPEED: Velocity standard deviation,
	# TRACK_ID: Track ID,
	# TRACK_MEDIAN_QUALITY: Median quality,
	# TRACK_STD_QUALITY: Quality standard deviation,
	# TRACK_X_LOCATION: X Location (mean),
	# TRACK_MEDIAN_SPEED: Median velocity,
	# NUMBER_SPOTS: Number of spots in track,
	# TRACK_MIN_SPEED: Minimal velocity,
	# NUMBER_GAPS: Number of gaps,
	# TRACK_Z_LOCATION: Z Location (mean),
	# TRACK_STOP: Track stop,
	# TRACK_MEAN_SPEED: Mean velocity,
	# NUMBER_SPLITS: Number of split events,
	# TRACK_MAX_SPEED: Maximal velocity,
	# TRACK_Y_LOCATION: Y Location (mean),
	# TRACK_DISPLACEMENT: Track displacement,
	# NUMBER_COMPLEX: Complex points,
	# TRACK_MEAN_QUALITY: Mean quality,
	# TRACK_DURATION: Duration of track,
	# TRACK_MAX_QUALITY: Maximal quality,
	# LONGEST_GAP: Longest gap,
	# TRACK_MIN_QUALITY: Minimal quality

	settings.addTrackFilter(FeatureFilter('NUMBER_SPLITS', 0.9, True))

	return trackmate


def process( trackmate ):
	"""
	Execute the full process BUT for the detection step.
	"""
	# Check settings.
	ok = trackmate.checkInput()
	# Initial filtering
	print( 'Spot initial filtering' )
	ok = ok and trackmate.execInitialSpotFiltering()
	# Compute spot features.
	print( 'Computing spot features' )
	ok = ok and trackmate.computeSpotFeatures( True )
	# Filter spots.
	print( 'Filtering spots' )
	ok = ok and trackmate.execSpotFiltering( True )
	# Track spots.
	print( 'Tracking' )
	ok = ok and trackmate.execTracking()
	# Compute track features.
	print( 'Computing track features' )
	ok = ok and trackmate.computeTrackFeatures( True )
	# Filter tracks.
	print( 'Filtering tracks' )
	ok = ok and trackmate.execTrackFiltering( True )
	# Compute edge features.
	print( 'Computing link features' )
	ok = ok and trackmate.computeEdgeFeatures( True )

	return ok


def display_results_in_GUI( trackmate, imp ):
	"""
	Creates and show a TrackMate GUI to configure the display
	of the results.

	This might not always be desriable in e.g. batch mode, but
	this allows to save the data, export statistics in IJ tables then
	save them to CSV, export results to AVI etc...
	"""

	gui = TrackMateGUIController( trackmate )

	# Link displayer and GUI.

	model = trackmate.getModel()
	selectionModel = SelectionModel( model)
	displayer = HyperStackDisplayer( model, selectionModel, imp )
	gui.getGuimodel().addView( displayer )
	displaySettings = gui.getGuimodel().getDisplaySettings()

	for key in displaySettings.keySet():
		displayer.setDisplaySettings( key, displaySettings.get( key ) )
	displayer.render()
	displayer.refresh()

	gui.setGUIStateString( 'ConfigureViews' )


def color_and_export_rois_by_track( trackmate, rm, path ):
	"""
	Colors the ROIs stored in the specified ROIManager rm using a color
	determined by the track ID they have.

	We retrieve the IJ ROI that matches the TrackMate Spot because in the
	latter we stored the index of the spot in the quality feature. This
	is a hack of course. On top of that, it supposes that the index of the
	ROI in the ROIManager corresponds to the line in the ResultsTable it
	generated. So any changes to the ROIManager or the ResultsTable is
	likely to break things.
	"""
	model = trackmate.getModel()
	track_colors = {}
	track_indices = []
	track_rois = {}
	for i in model.getTrackModel().trackIDs( True ):
		track_indices.append( i )
	shuffle( track_indices )

	for i, track_id in enumerate(track_indices):
		color = Jet.getPaint( float(i) / max(1, len(track_indices)-1) )
		track_colors[ track_id ] = color

	spots = model.getSpots()
	for spot in spots.iterable( True ):
		q = spot.getFeature( 'QUALITY' ) # Stored the ROI id.
		roi_id = int( q )
		roi = rm.getRoi( roi_id )

		# Get track id.
		track_id = model.getTrackModel().trackIDOf( spot )
		if track_id is None or track_id not in track_indices:
			color = Color.GRAY
		else:
			color = track_colors[ track_id ]
			roi_name = rm.getName( roi_id )
			track_rois.setdefault(track_id,[]).append(roi_name)

		roi.setFillColor( color )

	with open(path, 'w') as f:
		f.writelines([', '.join(track_rois[track_id])+'\n' for track_id in track_rois.keys()])


def error(msg):
	ui.showDialog(msg, DialogPrompt.MessageType.ERROR_MESSAGE);


def rename_rois(rm, is_hyperstack):
	frame = (lambda roi: roi.getTPosition()) if is_hyperstack else (lambda roi: roi.getPosition())
	rois = rm.getRoisAsArray()
	k = 1
	last_frame = -1
	for i,roi in enumerate(rois):
		t = frame(roi)
		if last_frame != t:
			last_frame = t
			k = 1
		rm.rename( i, 't%03d-%05d' % (t,k) )
		k += 1


def export_rois(rm, is_hyperstack, path):
	frame = (lambda roi: roi.getTPosition()) if is_hyperstack else (lambda roi: roi.getPosition())
	rois = rm.getRoisAsArray()
	polys = {}
	for roi in rois:
		fp = roi.getFloatPolygon()
		polys[roi.getName()] = {'t': frame(roi), 'x': fp.xpoints, 'y': fp.ypoints}

	writer = FileWriter(path)
	Gson().toJson(polys, writer)
	writer.close() # important


def export_calibration(imp, path):
	cal = imp.getCalibration()
	meta = dict(w=cal.pixelWidth,    w_unit=cal.getXUnit(),
				h=cal.pixelHeight,   h_unit=cal.getYUnit(),
				t=cal.frameInterval, t_unit=cal.getTimeUnit())
	writer = FileWriter(path)
	Gson().toJson(meta, writer)
	writer.close() # important


def save_path(save_dir, imp_name, name):
	out_dir = os.path.join(save_dir.getAbsolutePath(), imp_name)
	if not os.path.exists(out_dir):
		os.makedirs(out_dir)
	return os.path.join(out_dir, name)


def main():

	#------------------------------
	# 			MAIN
	#------------------------------

	imp_name = imp.getTitle()
	imp_name, ext = os.path.splitext(imp_name)

	models, channel_names, prob_threshs, nms_threshs = [], [], [], []
	if stardist_membrane_enabled:
		models.append(stardist_membrane)
		channel_names.append('Membrane')
		prob_threshs.append(prob_thresh_membrane)
		nms_threshs.append(nms_thresh_membrane)
	if stardist_dna_enabled:
		models.append(stardist_dna)
		channel_names.append('DNA')
		prob_threshs.append(prob_thresh_dna)
		nms_threshs.append(nms_thresh_dna)

	if len(models) == 0:
		return error("no stardist model enabled")

	if tracking_channel not in channel_names:
		return error("channel %s cannot be tracked, must be one of %s" % (tracking_channel, channel_names))

	n_channels = imp.getNChannels()
	n_frames = imp.getNFrames()
	is_hyperstack = n_channels > 1
	if n_frames < 2:
		return error("input must be a timelapse")
	if n_channels != len(models):
		return error("input image has %d channels, but %d stardist model(s) enabled" % (n_channels, len(models)))

	export_calibration( imp, save_path(save_dir, imp_name, 'calibration.json') )

	channel_imps = ChannelSplitter.split(imp)

	args = zip(channel_names, channel_imps, models, prob_threshs, nms_threshs)

	if tracking_channel == 'Membrane':
		args = reversed(args) # tracking_channel must come last

	params = {}
	params['modelChoice'] = "Model (.zip) from File"
	params['outputType'] = "ROI Manager"
	# params['roiPosition'] = "Automatic" # doesn't work because single channels are fed to StarDist, but result may be displayed on hyperstack
	params['roiPosition'] = "Hyperstack" if n_channels > 1 else "Stack"


	print "\n===============================\n"
	for channel_name, channel, model, prob_thresh, nms_thresh in args:
		params['input'] = channel
		params['modelFile'] = model.getAbsolutePath()
		params['probThresh'] = prob_thresh
		params['nmsThresh'] = nms_thresh

		# print 'StarDist', channel_name, ':', params, '\n'
		command.run(StarDist2D, False, params).get()
		rename_rois( rm, is_hyperstack )
		rm.runCommand( "Save",          save_path(save_dir, imp_name, 'rois_%s.zip'  % channel_name.lower()) )
		export_rois( rm, is_hyperstack, save_path(save_dir, imp_name, 'rois_%s.json' % channel_name.lower()) )

	assert channel_name == tracking_channel

	# backup global user-chosen measurements
	measurements = Analyzer.getMeasurements()
	# set needed measurements
	Analyzer.setMeasurements(Measurements.AREA + Measurements.CENTER_OF_MASS + Measurements.STACK_POSITION)
	# create measurements table
	rm.runCommand(imp, "Measure");
	# restore global user-chosen measurements
	Analyzer.setMeasurements(measurements)

	# close/hide measurements table
	results_window = ResultsTable.getResultsWindow()
	results_window.close(False)
	# results_window.setVisible(False)

	# Remove overlay if any.
	imp.setOverlay( None )

	# Get results table.
	results_table = ResultsTable.getResultsTable()
	# print results_table

	# Create TrackMate instance.
	trackmate = create_trackmate( imp, results_table, frame_link_dist, gap_close_dist, seg_split_dist )

	#-----------------------
	# Process.
	#-----------------------

	ok = process( trackmate )
	if not ok:
		sys.exit(str(trackmate.getErrorMessage()))

	#-----------------------
	# Display results.
	#-----------------------

	# TODO: close trackmate gui?

	# Create the GUI and let it control display of results.
	display_results_in_GUI( trackmate, imp )

	color_and_export_rois_by_track( trackmate, rm, save_path(save_dir, imp_name, 'tracks_%s.csv' % tracking_channel.lower()) )


main()