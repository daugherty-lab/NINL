# Microscopy scripts (from Reck-Peterson Lab)
for Stevens and Beierschmitt et al., <b>Antiviral function and viral antagonism of the rapidly evolving dynein activating adapter NINL</b>. in review at <i>eLife</i>.

# Basic Description
## Folders
- "ProjectsFiles" folder contain ImageJ macro and Cellprofiler project
- Result folder contain processed results
    - Measurements folder contain csv files
    - Images folder contain threshold image of DAPI and STAT1 channel

## BatchChannelSplitter.ijm
This ImageJ macro will read original image and identify DAPI and STAT1 channel from metadata.
Save each image as tif into user specified destination folder
It also perform background subtraction

## DAPI_IN_OUT_RATIO_RobustBackground.cpproj
This Cellprofiler project will identify DAPI region and STAT1 region in 3D. Then measure total and mean intensity of STAT1 channel inside and outside of DAPI region.
It output image intensity measurements into csv file.
- STAT1_Images.csv - all STAT1 intensity measurement has been background subtracted
    - Image Number : Image number assigned by Cellprofiler
    - Math_STAT1_Mean : Mean intensity of entire STAT1 channel inside of threshold region
    - Math_STAT1_MeanDAPI : Mean intensity of STAT1 channel inside of DAPI region - masked with STAT1 and DAPI threshold image
    - Math_STAT1_MeanNonDAPI : Mean intensity of STAT1 channel outside of DAPI region - maskek with STAT1 and inverse of DAPI threshold image
    - Math_STAT1_MeanRatio : Ratio of Math_STAT1_MeanDAPI over Math_STAT1_MeanNonDAPI
    - Math_STAT1_Total : Totalintensity of entire STAT1 channel inside of threshold region
    - Math_STAT1_TotalDAPI : Totalintensity of STAT1 channel inside of DAPI region - masked with STAT1 and DAPI threshold image
    - Math_STAT1_TotalNonDAPI : Totalintensity of STAT1 channel outside of DAPI region - maskek with STAT1 and inverse of DAPI threshold image
    - Math_STAT1_TotalRatio : Ratio of Math_STAT1_TotalDAPI over Math_STAT1_TotalNonDAPI

# How to Use it

## Run BatchChannelSplitter.ijm ImageJ Macro
1. Open macro in ImageJ - use Fiji
2. Run Macro
3. Select Source Directory where image files located
    - Image file can be in subfolders
4. Select Destination Directory where splitted image will be saved
5. Dialogbox will ask following options
    - Option to process subfolders under Source folder to find images
        - If you have subfolders, check this option
    - Option to create same subfolder structure in destination folder
        - Strongly recommned not to use this option
6. Push OK and wait finishing processing
7. Destination folder should have Filename_DAPI.tif and Filename_STAT1.tif

## Run DAPI_IN_OUT_RATIO_RobustBackground.cpproj Cellprofiler Project
1. Open project in Cellprofiler - require Cellprofiler version 3.0
2. Set Default output folder
    - Push "View output settings" button at the bottom left area
    - Specify "Default Output Folder" to the folder you want to save result.
3. Select "Images" module - top of the module list and add image files
    - Drag and Drop image files and confirm image list show up
    - Drop both _DAPI and _STAT1 images
    - You can use top folder of image
4. Select "Metadata" module and update metadata
    - Push "Update" button and confirm metadata are displayed
5. Select "NamesAndTypes" module and update settings
    - Push "Update" button and confirm processing images are displayed
    - DAPI column has DAPI file name and STAT1 column has STAT1 file name
6. Run project by pressing "Analyze Images" at the left bottom area
7. Wait for finishing processing
    - Images folder will be created and contain threshold of DAPI and STAT1 images
    - Measurements folder will be created and contain STAT1_Images.csv

# Technical Description

## BatchChannelSplitter.ijm ImageJ Macro
Since original image channel assignment order is random, it require to read image meta data to figure out which channel contain DAPI and STAT1 channel.
Luckily Name #2, Name #3, Name #4 contain channel information "405", "561", "640" to identify excitation of the channel.
This macro read these metadata and identify DAPI and STAT1 channel. Then save these two images to specified destination folder.

## IntensityMeasurements.cpproj Cellprofiler Project
- All processing is in 3D
- Nucleus region will be identified by using threshold DAPI channel
    - Use Minimum cross entropy method to automatically determine threshold
        - Threshold smoothing scale : 4
    - Save threshold image into "images" folder under default output folder
- STAT1 region will be identified by using threshold STAT1 channel
    - Use Robust Background Method to determine threshold
        - Threshold smoothing scale : 4
    - Save threshold image into "images" folder under default output folder
- Mask Images for image intensity measurements
    - Entire STAT1 Region : Mask STAT1 image with threshold STAT1 image for STAT1 region measurements
    - STAT1 inside DAPI : Mask STAT1 region image created above with DAPI threshold image
    - STAT1 outsside DAPI : Mask STAT1 region image created above with inverse of DAPI threshold image
- Measure intensity of masked images
    - Use Total and Mean intensity
    - STAT1, STAT1 in DAPI, STAT1 ouside DAPI  
- Calibrate measured intensity into original image intensity
    - Multiply 65535
    - STAT1 total and Mean intensity
    - STAT1 in DAPI total and Mean intensity
    - STAT1 outside DAPI total and mean intensity
- Calculate ratio values
    - Ratio of Total STAT1 in DAPI over Total STAT1 ouside DAPI
    - Ratio of Mean STAT1 in DAPI over Mean STAT1 ouside DAPI
- Save calibrated values and calculated ratio into csv file.











