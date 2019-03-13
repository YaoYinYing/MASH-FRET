---
layout: default
title: Workflow
parent: /video-processing
nav_order: 2
---

# Video processing workflow
{: .no_toc }

In this section you will learn how to process single molecule videos to obtain single molecule coordinates and intensity-time traces. Exported data, in particular to the 
[mash project](../../output-files/mash-mash-project) file, can further be used for 
[data analysis](../tutorials/analyze-data).

The procedure includes four steps:

1. TOC
{:toc}

## Register experiment settings

In theory, the software is compatible with:
* an unlimited number of video channels, 
* an unlimited number of alternated lasers,
* FRET calculations for an unlimited number of donor-acceptor network.

To adapt the software functionalities to your own experiment setup, MASH must be informed about the particular settings.

To register your experiment settings:

{: .bg-grey-lt-000 .pt-3 .pb-2 .pl-7 .pr-4}
1. Set parameters:  
     
   [Number of alternated lasers](panels/panel-experiment-settings#number-of-alternated-lasers)  
   [Laser wavelengths](panels/panel-experiment-settings#laser-wavelengths)  
   [Number of video channels](panels/panel-experiment-settings#number-of-video-channels)  
   [Project options](panels/panel-experiment-settings#project-options)  
   [Exposure time](panels/panel-experiment-settings#exposure-time) (optional); see 
   [Remarks](#remarks) for more details.  
     
1. Modify experiment settings whenever the laser or channel order change, or when experimental conditions varies.


## Localize bright spots

Bright spots coordinates are different from single molecule coordinates to the extent that molecule coordinates are paired across all video channels and spots coordinates concern individual channels.
Spot detection in single molecule videos (SMV) is often made more complicated by the presence of competitive background noise and by the variation of single molecule brightness in time.

To increase the contrast between bright spots and background, and thus ease spot detection, SMVs can be treated with image filters.
MASH offers a set of basic or smFRET-specific image filters that can be used for such purpose.
Some image filters can also be used as background correction prior creating intensity-time traces; see 
[Remarks](#remarks) for more details.

*Noise SMV to filter SMV and .sira file*  
[Figure]

Image filter or not, it is recommended to export every video analyzed to a 
[.sira file](../output-files/sira-mash-video). 
This will eventually reduce the processing time.

To get rid of brightness variations in time and smooth the background noise, video frames are averaged into one average image.
Finally, spot detection is performed on the average image with the Spotfinder tool.

*SMV averaging and spot detection*  
[Figure]

To localize bright spots in the SMV:

{: .bg-grey-lt-000 .pt-3 .pb-2 .pl-7 .pr-4}
1. [Load](panels/area-visualization#load-videoimage-file) the SMV file by pressing 
   ![Load...](../assets/images/gui/VP-but-load.png "Load...")  
     
1. Set parameters in 
   [Filter settings](panels/panel-edit-video#filter-settings) for each channel and apply the filter by pressing 
   ![Add](../assets/images/gui/VP-but-add.png "Add"). 
   Several filters can be cumulated with the top-filter in the 
   [Filter list](panels/panel-edit-video#filter-list) being applied first.  
     
1. Set parameters in 
   [Frame range](panels/panel-edit-video#frame-range) and export the modified (or not) video to a 
   [.sira file](../output-files/sira-mash-video) by pressing 
   ![Export...](../assets/images/gui/VP-but-export.png "Export...")   
     
1. [Load](panels/area-visualization#load-videoimage-file) the newly exported SMV file by pressing 
   ![Load...](../assets/images/gui/VP-but-load.png "Load...")  
     
1. Set parameters in [Average image](panels/panel-molecule-coordinates#average-image) 
     
1. Calculate and export the average image to a 
   [_ave.* file](../output-files/ave-average-image) by pressing  
   ![Go](../assets/images/gui/VP-but-go.png "Go")  
     
1. Load the average image by pressing 
   ![...](../assets/images/gui/VP-but-3p.png "...") in
   [Average image](panels/panel-molecule-coordinates#average-image)  
     
1. Set parameters in 
   [Spotfinder settings](panels/panel-molecule-coordinates#spotfinder-settings) for each channel  
     
1. Start the detection procedure by pressing 
   ![Find](../assets/images/gui/VP-but-find.png "Find"). 
   If necessary, refine the spots set using the 
   [Exclusion rules](panels/panel-molecule-coordinates#exclusion-rules)  
     
1. Export spots coordinates to a 
   [.spots file](../output-files/spots-spots-coordinates) by pressing 
   ![Save](../assets/images/gui/VP-but-save.png "Save")


## Transform spots coordinates

To obtain single molecules coordinates, spots coordinates must be transformed into other channels.
This implies that the spatial transformation specific to your setup must be calculated.

The spatial transformation is calculated from an already co-localized set of coordinates called the reference coordinates.
They are mapped manually using a reference image where reference samples emits in all channels.

*Map reference to transformation*  
[Figure]

The transformation from one channel to another uses a combination of symmetry operations which is specific to the recording setup. 
A list of transformation types is available for such purpose, going from the most simple to the most complex combination.

At first, the proper transformation type is usually unknown.
A good practice is to use the most simple one and check the quality of the calculated transformation.
The quality of the transformation is judged by the user's eye from the superposition of the original (in red) and transformed (in green) reference images.
If the transformation is correct, reference emitters will appear as yellow dots on a dark background.
If the quality is poor, red and green dots will be visible and the operation must be repeated with increasing complexity of the combination.

*Iteration calculate-check*  
[Figure]

To calculate and export the spatial transformation:

{: .bg-grey-lt-000 .pt-3 .pb-2 .pl-7 .pr-4}
1. Load the reference image by pressing 
   ![Map](../assets/images/gui/VP-but-map.png "Map")  
     
1. [Use the mapping tool](functionalities/use-mapping-tool) to map single emitters positions in every channels   
     
1. Export reference coordinates to a 
   [.map file](../output-files/map-mapped-coordinates) by closing the mapping tool.  
     
1. Select the 
   [Transformation type](panels/panel-molecule-coordinates#transformation-type)
     
1. Calculate and export the transformation to a 
   [.mat file](../output-files/mat-transformation) by pressing 
   ![Calculate](../assets/images/gui/VP-but-calculate.png "Calculate")
     
1. Load the reference image by pressing 
   ![Check quality](../assets/images/gui/VP-but-check-quality.png "Check quality") and judge the transformation quality. 
   If the quality is not satisfying, return to step 4.  
     
1. Recalculate the transformation whenever the setup is realigned. 
   Otherwise, the same transformation [.mat file](../output-files/mat-transformation) can be re-used; see 
   [Coordinates transformation](panels/panel-molecule-coordinates#coordinates-transformation) for more information.
     
1. Transform the spots coordinates and export the single molecule coordinates to a 
   [.coord file](../output-files/coord-transformed-coordinates) by pressing 
   ![Transform](../assets/images/gui/VP-but-transform.png "Transform")


## Create and export intensity-time traces

Intensities are calculated by summing up the brightest pixels laying in a square area drawn around the molecule coordinates.
The operation is repeated on every video frame to build the single molecule intensity-time traces.

Intensity-time traces and project parameters are automatically saved to a 
[.mash project](../output-files/mash-mash-project) and statistics on intensities are automatically exported to a 
[.tbl file](../output-files/tbl-intensity-statistics), but intensity data can be also exported to other optional files.

To calculate and export the intensity-time traces:

{: .bg-grey-lt-000 .pt-3 .pb-2 .pl-7 .pr-4}
1. Load the unmodified SMV file by pressing 
   ![...](../assets/images/gui/VP-but-3p.png "...") in 
   [Input video](panels/panel-intensity-integration#input-video).
   Background-corrected video files can also be used here; see 
   [Remarks](#remarks) for more details.   
     
1. Load the single molecule coordinates 
   [.coord file](../output-files/coord-transformed-coordinates) by pressing 
   ![...](../assets/images/gui/VP-but-3p.png "...") in 
   [Input coordinates](panels/panel-intensity-integration#input-coordinates)  
     
1. Set parameters:  
     
   [Integration parameters](panels/panel-intensity-integration#integration-parameters)  
   [Export options](panels/panel-intensity-integration#export-options)  
     
1. Calculate and export intensity-time traces to a 
   [.mash project](../output-files/mash-mash-project) and to selected ASCII files by pressing 
   ![Create & Export...](../assets/images/gui/VP-but-export.png "Create & Export..."). 
   The process is rather slow; see 
   [Remarks](#remarks) for more details.
   
Et voilà!


## Remarks
{: .no_toc}

The exposure time is automatically updated from the video file. 
If the file does not contain such information, the exposure time must be manually set in 
[Exposure time](panels/panel-experiment-settings#exposure-time).

Even though background correction is more accurate when performed in module Trace processing, some image filters can be used as background correction; see 
[Filter settings](panels/panel-edit-video#filter-settings) for more information.
In that case, the background-corrected video file can be used to create intensity-time traces, but background correction must be deactivated in module Trace processing; see
[Background correction](../trace-processing/panels/panel-subimage-background-correction#background-correction) for more information.

Creating intensity-time traces is a rather slow process because the script was written for low-RAM computers.
For more information about how the script works, please refer to the respective function in the source code:

```
MASH-FRET/source/traces/creation/getIntTrace.m
```