# Landsat-Panchromatic-Assisted-Downscaling (LPAD)
Pansharpening of Landsat 30 m data to 15 m resolution and then reprojection of the Landsat 15 m data into Sentinel-2 20 m resolution

REQUIREMENT:

	Matlab R1014a+ (8.3.0.532) on Windows or Unix (lower version Matlab could be run but not tested yet)
	Matlab Mapping Toolbox (this is only for geotiff reading/writing)

Contact information:

	Dr. David P. Roy (david.roy@sdstate.edu), GSCE, SDSU, US
	Dr. Hankui K. Zhang (hankui.zhang@sdstate.edu), GSCE, SDSU, US
	Dr. Zhongbin Li (zhongbin.li@sdstate.edu), GSCE, SDSU, US

Usage:

First, 'ProLPAD30to15.m' is used to finish the first step of LPAD, i.e., downscaling 30 m Landsat 8 OLI bands into 15 m using 15 m panchromatic band. Parameters need changes: 

    	%% (1) 'dataDir': the data directory to store 15 m panchromatic and 30 m OLI files
    	%% (2) 'PanFile': the pan file name
    	%% (3) 'MSFile': the low resolution MS file name, keep it the original size (do not resample)
    
    	%% the following parameters also need change (optional)
   	%% (4) 'rgb_index':  which three bands are blue, green and reb bands, setting them correctly is extremely important 
    	%% (5) 'sharped_index':  which bands are to be downscaled/sharpened 
    	%% (6) 'scaleF':  the re-scale factor of the reflectance 
	
    	1) Note RGB bands must be included and the rgb_index must be correctly specified.
	
	2) Note that geotiff formats with correct UTM coordinates information are needed for 30 m and 15 m bands;

	3) Both cropped (without resampling) L1T and original L1T data can be used as input;

	4) When cropping the data, make sure the geographic extent of the 30 m images are larger than that of the 15 m panchromatic image, the software will automatically crop the MS image according to the panchromatic extent; 

	5) If the size of the panchromatic image has odd number, the last row or column will not be processed. Since original L1T images have m*n 30 m pixels and (2m-1)*(2n-1) 15 m panchromatic pixels, the output will have (2m-2)*(2n-2) 15 m downscaled pixels.


Then, 'ProLPAD15to20.m' is used to finish the second step of LPAD, i.e., reprojecting the 15 m (or 30 m) Landsat data into the Sentinel-2 20 m resolution.

Notes:

	1) the input data of 'ProLPAD15to20.m' is the output of 'ProLPAD30to15.m';
	
	2) both bilinear and cubic convolution resampling approaches are provided; set the variable 'resampler='bl' or 'cc'' to choose;
	
	3) the 'bl_cc_resampling_at.m' will be invoked, in which the affine transformation coefficients between Landsat 20 m and Sentinel-2 20 m data are provided.  


References:

Aiazzi B.; Baronti S.; Selva M.; Alparone L., Bi-cubic interpolation for shift-free pan-sharpening. ISPRS Journal of Photogrammetry and Remote Sensing. 2013, 86, 65-76.

Zhang H.K. and Roy D.P., Computationally inexpensive Landsat 8 operational land imager (OLI) pansharpening. Remote Sensing. 2016, 8(3), 180.

Yan, L., Roy, D. P., Zhang, H., Li, J., & Huang, H. (2016). An automated approach for sub-pixel registration of Landsat-8 Operational Land Imager (OLI) and Sentinel-2 Multi Spectral Instrument (MSI) imagery. Remote Sensing, 8(6), 520.

Li Z., Zhang H.K., Roy D.P., Yan L., and Li J. Landsat 15-m Panchromatic Assisted Downscaling (LPAD) of the 30 m Reflective Wavelength Bands to Sentinel-2 20-m Resolution (2017). Remote Sensing.
