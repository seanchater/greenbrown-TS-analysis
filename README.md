# greenbrown-TS-analysis
Workflow for identifying breakpoints in gridded raster data using the R package greenbrown.

------
## Files:
- collect_greenbrown.js => a GEE script for collecting a Landsat timeseries over multiple years with harminization between different sensors and clour masking.
- greenbrown_script.R => an R script for pre-processing (LOESS gapfilling, Savitzky-Golay filtering, and control differencing) as well as running differnt mehtods of greenbrown algorithm
------
## Other files:
- are testing files and contain duplicate code
