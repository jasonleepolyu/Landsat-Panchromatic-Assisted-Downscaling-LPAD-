function [D,Low_combination_res,Low_combination,xPan_hist] = CS_Transform(xPAN,PLow,MLow,MLow_res,beta,bands_array)
%% GS transformation  
% MLow: low resolution MS images
% MLow_res: low resolution MS images resampled by using a 23-p filter
    [m,n,bands] = size(MLow_res);
    Low_combination_res = zeros(m,n,'single');
    [mlow,nlow,bands] = size(MLow);
    Low_combination = zeros(mlow,nlow,'single');
%% I = BX+a
    beta_len = size(beta,1);
%     size(bands_array)
    if beta_len>2
        for b=1:size(bands_array,2)
            beta(b);
            Low_combination_res = Low_combination_res+MLow_res(:,:,bands_array(b))*beta(b);
            Low_combination = Low_combination+MLow(:,:,bands_array(b))*beta(b);
        end
    else
        for b=2:bands
            Low_combination_res = Low_combination_res+MLow_res(:,:,b)*beta(b-1);
            Low_combination = Low_combination+MLow(:,:,b)*beta(b-1);
        end
    end

%% histogram matching from PAN to I (Low_combination_res)
%     LowI2d = hyperConvert2d(Low_combination_res);
%     PAN2d = hyperConvert2d(xPAN);
%     b = regress(LowI2d',PAN2d');
%     fprintf('The histogram b = %f\n',b);
%     xPan_hist = b*xPAN;

    LowI2d = hyperConvert2d(Low_combination);
    Panlow2d = hyperConvert2d(PLow);
    xPan_hist = (xPAN-mean(Panlow2d))*sqrt(var(LowI2d)/var(Panlow2d))+mean(LowI2d);    
  
    %% getting details
    D = xPan_hist - Low_combination_res;
end