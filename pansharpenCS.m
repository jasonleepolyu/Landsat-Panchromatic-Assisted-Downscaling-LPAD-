 function [sharpenedMS, residual_low,CovVar_variance] = pansharpenCS(x_MLow,x_MLow_resampled,MLow,MLow_resampled,P,PLow,model_case,FACTOR,windowsize,beta)
%%
% x_MLow_resampled is the low resolution MS images with N bands resampled to 15 m
% MLow is the low resolution MS images with N bands at 30 m
% xPAN is the high resolution PAN image with one band at 30 m
% PLow is the degraded PAN image at 30 m (only for evaluation)
% FACTOR is the spatial resolution ratio
%%---------------------------------------------------------------------
%%
    if ~exist('beta','var')
        beta = ones(bands,1)/bands;
    end

%     level = 2;% an important parameter
    [m,n,nb] = size(MLow);
    [m,n,nb_sharpen] = size(x_MLow);
    
    %% 
    sharpenedMS = x_MLow_resampled;%% 23-taps filter expansion
    if (strcmp(beta,'Regress'))
        Plow2d = (hyperConvert2d(PLow))'; %% regress(Plow2d,Mlow2d);
        Mlow2d = (hyperConvert2d(MLow(:,:,1:3)))'; %% fitlm(Mlow2d',Plow2d','y ~ x1 + x2 + x3');
        beta = mvregress(Mlow2d,Plow2d);
        beta = beta';
        beta
    end
    
    [D,MLow_resampled_intensity,MLow_intensity,P_hist] = CS_Transform(P,PLow,MLow,MLow_resampled,beta',1:size(beta,2));%% 'PanLow_combination' is y in Eq. 7-5

   
    
    %% variance calculation
    CovVar_variance = 0;
    if strcmp(model_case,'GS')||strcmp(model_case,'GS3')
        [CC_G,Coeffi_G,CovVar] = getCoeffiGlobal(x_MLow,MLow_intensity,model_case);
    end
    if strcmp(model_case,'GSLocal')
        [AABP,ECB,CovVar] = getCoeffiLocal(x_MLow,MLow_intensity,windowsize,model_case);
        CovVar = getHigh2(CovVar,FACTOR,7);
    end
    
    %% doing the injection
	Dnew = D;%%default value 
    for k=1:nb_sharpen    
        if strcmp(model_case,'GS')||strcmp(model_case,'GS3')
            Dnew(:,:)= D(:,:).*CovVar(k);
            sharpenedMS(:,:,k) = sharpenedMS(:,:,k)+Dnew;
        end
        if strcmp(model_case,'GSLocal')
            Dnew(:,:)= D(:,:).*CovVar(:,:,k);
            sharpenedMS(:,:,k) = sharpenedMS(:,:,k)+Dnew;
        end
        if strcmp(model_case,'Brovey')
%             fprintf('What is wrong');
            xtemp = P_hist./max(MLow_resampled_intensity,1);
            sharpenedMS(:,:,k) = xtemp.*sharpenedMS(:,:,k);
        end
    end% end of the bands


    %% calculate the residuals
    residual_low = abs(PLow-MLow_intensity);
 end