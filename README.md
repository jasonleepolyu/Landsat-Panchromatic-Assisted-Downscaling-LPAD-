# Landsat-Panchromatic-Assisted-Downscaling (LPAD)
Pansharpening of Landsat 30 m data to 15 m resolution and then reprojection of the Landsat 15 m data into Sentinel-2 20 m resolution

REQUIREMENT:
Matlab R1014a+ (8.3.0.532) on Windows or Unix (lower version Matlab could be run but not tested yet)
Matlab Mapping Toolbox (this is only for geotiff reading/writing)

Usage:
1 Use 'ProPADS30to15' to finish the step one of PADS, i.e., downscaling 30 m Landsat 8 OLI bands into 15 m using 15 m panchromatic band; 
Parameters need changes: 
	%% (1) 'dataDir': the data directory to store 15 m panchromatic and 30 m OLI files
    %% (2) 'PanFile': the pan file name
    %% (3) 'MSFile': the low resolution MS file name, keep it the original size (do not resample)
    
    %% the following parameters also need change (optional)
    %% (4) 'rgb_index':  which three bands are blue, green and reb bands, setting them correctly is extremely important 
    %% (5) 'sharped_index':  which bands are to be downscaled/sharpened 
    %% (6) 'scaleF':  the re-scale factor of the reflectance 
	
	Note RGB bands must be included and the rgb_index must be correctly specified.
	
2 Note that geotiff formats with correct UTM coordinates information is needed for 30 m and 15 m bands;

3 Both clipped (without resampling) L1T and original size L1T data can be used as input;

4 When clipped the data, make sure the geographic extent of 30 m images are larger than that of the 15 m panchromatic image, the software will automatically clip the MS image according to the panchromatic extent; 

5 If the size of the panchromatic has odd number, the last row or column will not be processed. This applied to the original L1T images. Since original L1T images have m*n 30 m pixels and (2m-1)*(2n-1) 15 m panchromatic pixels, the output will have (2m-2)*(2n-2) 15 m downscaled pixels.
Reference:

Aiazzi B.; Baronti S.; Selva M.; Alparone L., Bi-cubic interpolation for shift-free pan-sharpening. ISPRS Journal of Photogrammetry and Remote Sensing. 2013, 86, 65-76.

Zhang H.K., Roy D.P., Computationally inexpensive Landsat 8 operational land imager (OLI) pansharpening. Remote Sensing. 2016, 8(3), 180.

Li Z., Zhang H.K., Roy D.P., Yan L., Panchromatic Assisted Downscaling of Landsat-8 30 m Reflective Band Data to Sentinel-2 20 m Resolution. Remote Sensing. In review.

