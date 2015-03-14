_A library for sorting, processing, and analysis of massive point clouds generated with terrestrial laser scanning_

## Features ##
  * Exporting from multi-GB ASCII datasets (Cyclone PTS format).
  * CPU-parallel I/O and processing on memory-mapped files
  * Sorting into an xy grid
  * Centroid, zmin, zmax -thinning;
  * Per-cell statistics
  * Per-cell detrending against a ground level defined on centroid triangulation.
  * Statistical detrended-point-based measures of surface roughness.

  * C++ dynamic library compilation
  * Python script interface

## System requirements ##
  * 64-bit OS and 64-bit Python
  * 64 version of wx if you want to use GUI for file selection: [Python Extension Packages for Windows](http://www.lfd.uci.edu/~gohlke/pythonlibs/)

We will be providing Win64 builds; Linux 64 builds by request.

Building instructions will be provided later.

### Getting started ###
On a Win64 install Python 2.6 64-bit into "C:\Python26(64)". Unzip the pctools\_win64...zip into a test folder. Run the pctools.bat. Refer to the "pctools.py" for the work-flow and parameters description.


## Used in ##
_"Terrestrial Laser Scanning in Fluvial Geomorphology: Retrieving Morphological and Sedimentological Models of Gravel Bed Rivers"_, **EGU2009**, http://www.reesscan.com/Publications/EGU_2009a.pdf?attredirects=0&d=1

## Links ##
[CloudCompare - Open Source project](http://www.danielgm.net/cc/)

to surface reconstruction projects.
to the Lidar and other remote sensing projects.
to laser scanning groups in geophysics.
to point-based visualization in computer graphics.



## Changes ##
0.1.1
exporting decimated point cloud for populated above threshold

improving the pctools.py: now it checks for the size of the proposed grid to be below 4GB (you can change this limit)