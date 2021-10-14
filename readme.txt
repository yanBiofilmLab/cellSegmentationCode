Description
This program segments 3d confocal images of dense biofilms into its single cell constituents. The algorithm takes in a deconvolved image (with background removed) and outputs a list of 3d coordinates corresponding to each bacterium. This program has been tested on MATLAB R2018a. Contact Jing Yan, jing.yan@yale.edu, MCDB, Yale University for further information.
Installation
This program is meant to be used as a stand-alone algorithm in MATLAB. To use, add the underlying functions to the path, set the location of the image of interest in runsegmentation.m and execute.
Example
An example biofilm image, corresponding to the segmented biofilm in Supplementary Figure 2, is included. Runtimes depend on the size of the image; segmenting the example image takes approximately 2-3 hours on a laptop computer. 

