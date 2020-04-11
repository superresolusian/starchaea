[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/superresolusian/starchaea/master)

# starchaea
notebooks for stardist x trackmate x archaea project

## Data and models
Find 'em here: https://www.dropbox.com/sh/jbczc4devj1cjby/AAD-LELLm6_v74o9B7CLXJUUa?dl=0

## Order to run notebooks/code things in. General notes and musings.
0. Model training - haven't included this. Don't know if it's worth making a separate notebook in this repo or just point people to Stardist training example? Will gather together some examples of annotated data for both channels either way.

1. `Collated_process_up_to_trackmate.ipynb` - I've tested the fuck out of this for all the example data, should be pretty stable.

2. `Tracking_helper.ijm` in Fiji (needs to have `my_tracking.py` in Fiji plugins folder). Probably not worth trying to call this from a notebook is it? I got a bit over excited when I realised that you can open Fiji from a jupyter notebook (`Probably_a_bad_idea.ipynb`). <font color=red> Maybe should have GUI options for settings inside my_tracking? E.g. gap lengths etc </font>

3. `Process_trackmate.ipynb` - I've tested this on 2 colour data but not single colour data.

4. `Curation_helper.ijm` in Fiji. This is optional, instructions for use are provided at end of `Process_trackmate.ipynb`.

5. `Measure_polygons.ipynb` - again, I've tested this on 2 colour data but not single colour.

## Soundtrack
https://www.youtube.com/watch?v=jyO-MyJ4R1g - I LOVE this song and also it's by Starcadian which is basically starchaea :)
