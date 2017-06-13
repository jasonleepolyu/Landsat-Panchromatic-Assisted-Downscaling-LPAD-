% function ProPADS30to15

    %% Input data parameters need change:
	%% (1) 'dataDir': the data directory to store 15 m panchromatic and 30 m OLI files
    %% (2) 'PanFile': the pan file name
    %% (3) 'MSFile': the low resolution MS file name, keep it the original size (do not resample)
    
    %% the following parameters also need change (optional)
    %% (4) 'rgb_index':  which three bands are blue, green and reb bands, setting them correctly is extremely important 
    %% (5) 'sharped_index':  which bands are to be downscaled/sharpened 
    %% (6) 'scaleF':  the re-scale factor of the reflectance 

	%% START *********************************************************************************************************************************
	%% important input data parameters need change for your own data
        
    dataDir = 'L:\\david_roy\\Hankui\\dull\\PADS_OLI_delivered\\';
    PanFile = '15m_B8.tif';
    MSFile = '30m_B1_B7.tif';
    outprefix='15m_B1_B7';

    sharped_index = 1:7;
    rgb_index = 2:4;
    scaleF = 10000; 
	%% END   *********************************************************************************************************************************

    %% different spectral weights (set to combine_case = 'RGB' in PADS)
    combine_case = 'Regress';
    beta = 'Regress'; %% Regression based spectral weights determination, see Zhang and Roy 2016

    combine_case = 'RGB.E';
	beta = [1.0/3	1.0/3	1.0/3];%% Equal spectral weights, see Zhang and Roy 2016
    
    combine_case = 'RGB';
    beta = [0.08017	0.5177	0.4030];%% OLI spectral response function spectral weights, see Zhang and Roy 2016 and Li et al. (in review)

   
    %% different fusion methods (set to method = 'Brovey' in PADS)
    method = 'NONE';%EXP
    method = 'GS'; % GS method, see Zhang and Roy 2016 
    method = 'GSLocal'; % local adaptive GS method, see Zhang and Roy 2016 
    method = 'Brovey'; % Btovey method, see Zhang and Roy 2016 and Li et al. (under review)

    windowsize = 6; % only used by method = 'GSLocal' and not used by PADS (method = 'Brovey') 

    %% input parameters
    ResolutionRatio = 2; 
    fprintf('Input parameters: Method = %s\t RGB.combination = %3s\t datadir = %s\n',method,combine_case,dataDir);

    %% read the input files
    Pinfo = geotiffinfo(sprintf('%s%s', dataDir,PanFile));
    [P,R] = geotiffread(sprintf('%s%s', dataDir,PanFile));
    MSinfo = geotiffinfo(sprintf('%s%s', dataDir,MSFile));
    MLow = geotiffread(sprintf('%s%s', dataDir,MSFile));
    [mp,np] = size(P);
    if (mod(mp,ResolutionRatio)~=0 || mod(np,ResolutionRatio)~=0 )
        fprintf('\tNote: The size of 15 m panchromatic image is odd, the last column or row is not processed \n');
        mp = int16(mp/ResolutionRatio-0.1)*ResolutionRatio;
        np = int16(np/ResolutionRatio-0.1)*ResolutionRatio;
        P = P(1:mp, 1:np);
%         xi = 50; yi = 100;
%         P = P(1:yi, 1:xi);
        [xlimits1, ylimits1] = intrinsicToWorld(R, 0.5, 0.5);
        [xlimits2, ylimits2] = intrinsicToWorld(R, double(np)+0.5, double(mp)+0.5);
        subR = R;
        subR.RasterSize = size(P);
        subR.XLimWorld = sort([xlimits1,xlimits2]);
        subR.YLimWorld = sort([ylimits1,ylimits2]);
        R = subR;

        [mp,np] = size(P);
    end  
    [xpUL,ypUL] = pix2map(Pinfo.RefMatrix, 1, 1);
    [xmUL,ymUL] = pix2map(MSinfo.RefMatrix, 1, 1);
    [m,n,nb] = size(MLow);
    [xpLR,ypLR] = pix2map(Pinfo.RefMatrix, mp, np);
    [xmLR,ymLR] = pix2map(MSinfo.RefMatrix, m, n);
    
    if (xpUL<xmUL || ypUL>ymUL || xpLR>(xmLR+15) || ypLR<(ymLR-15))
        fprintf('\tWarning: The 15 m panchromatic image has samller extent than the 30 m OLI data, downscale cannot be proceeded\n');
        return;
    end

    xclip = 1:m;
    yclip = 1:n;
    if (xpUL~=xmUL || ypUL~=ymUL || xpLR~=(xmLR+15) || ypLR~=(ymLR-15))
        fprintf('\tThe 15 m panchromatic and 30 m OLI data have different extents.... so a clipping will be applied to the panchromatic extent\n');
        
        xclip = (1+int16((ymUL-ypUL)/30)):(m-int16((ypLR-ymLR+15)/30));
        yclip = (1+int16((xpUL-xmUL)/30)):(n-int16((xmLR+15-xpLR)/30));
    end
    
    %% covert the image into single/float
    MLow = single(MLow(xclip,yclip,:))*scaleF;
    P = single(P)*scaleF;
    [m,n,nb] = size(MLow);
    [mp,np] = size(P);
    if (mp~=ResolutionRatio*m || np~=ResolutionRatio*n)
        fprintf('\tWarning: the dimensions of 15 m panchromatic and 30 m OLI images do not match, please check the data and clip\n');
        return;
    end

    %% resampling the original low-resolution image
    resample='odd'; %% the 'odd' resampling should be used for PADS see Table 2 in Aiazzi et al. 2013
    if strcmp(resample,'even') 
        MLow_resampled = getHigh2(MLow,ResolutionRatio,6); % even without shift
    else
        MLow_resampled = getHigh2(MLow,ResolutionRatio,7); % odd with half-pixel shift towards lower right
    end
	
    
	x_MLow = MLow(:,:,sharped_index);
	x_MLow_resampled = MLow_resampled(:,:,sharped_index);
	RGB_MLow = MLow(:,:,rgb_index);
	RGB_MLow_resampled = MLow_resampled(:,:,rgb_index);
	
    PLow = getLow2(P,ResolutionRatio,6);
    % MLow = MLow_resampled;

%% starting to pan-sharpening
%     [m,n,nb] = size(MLow_resampled);

    tic; % a time recorder
    
    if strcmp(method,'NONE')
        solutionPan = x_MLow_resampled;%% bi-cubic expansion
    elseif strcmp(method,'GS')||strcmp(method,'GS3')||strcmp(method,'Brovey')||strcmp(method,'GSLocal')
        [solutionPan,residual_low,var_mask] = pansharpenCS(x_MLow, x_MLow_resampled, RGB_MLow,RGB_MLow_resampled,P,PLow,method,ResolutionRatio,windowsize,beta);
    else
        fprintf('you are off the route\n');
    end
    fprintf('Processing time %5.3f...minutes\n',toc/60);

    %% write the results 
    geotiffwrite(sprintf('%sPXS_%s_%s_%s.tif', dataDir,outprefix,method,resample),int16(solutionPan),R,...
        'GeoKeyDirectoryTag',Pinfo.GeoTIFFTags.GeoKeyDirectoryTag);
    
   %% evaluation method, see Zhang and Roy 2016
    if exist('residual_low','var')
        % relative
        residual_low_rel = int16(residual_low./max(1,PLow)*scaleF);
        residual_low2d = hyperConvert2d(residual_low_rel);
        fprintf('The mean relative error %7.6f and sqrt(sigma) = %7.6f\n',single(mean(residual_low2d))/scaleF,sqrt(var(single(residual_low2d)))/scaleF);
        % residual        
        residual_low2d = hyperConvert2d(residual_low);
        fprintf('The mean absolute error %7.6f and sqrt(sigma) = %7.6f\n',single(mean(residual_low2d))/scaleF,sqrt(var(single(residual_low2d)))/scaleF);
    end
    fprintf('Processing time %5.3f...minutes\n',toc/60);
% end


%     stretch = stretch_log(solutionPan(:,:,1:3));
%     imwrite(stretch, sprintf('%sPXS_%s_%s.tif',jpgDir,method,resample),'tif');
%     resample='even';
%     geotiffwrite(sprintf('%splain_%s.tif', dataDir,resample),MLow_resampled,R,...
%         'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
%     jpgDir = sprintf('%s\\jpg\\', dataDir);
%     if (~exist(jpgDir))    mkdir(jpgDir);    end
    
%     if ~exist('lambda')
%        lambda = 1;
%        lambda = 1/100;
%        lambda = 1/1000;
%     end

%         geotiffwrite(sprintf('%sPXS_residual_%s.tif', dataDir,combine_case),int16(residual_low),R,...
%             'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
%         geotiffwrite(sprintf('%sPXS_realtive_%s.tif', dataDir,combine_case),residual_low_rel,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);      
%         [residual_low_stretch,index_sta_ref] = stretch_log_pan_RGB(residual_low_rel,350);
%         imwrite(residual_low_stretch,sprintf('%s.PXS_realtive_%s.tif', dataDir,combine_case));   
%         [residual_low_stretch,index_sta] = stretch_log_pan_RGB(residual_low,35);
%         imwrite(residual_low_stretch,sprintf('%s.PXS_residual_%s.tif', dataDir,combine_case));
