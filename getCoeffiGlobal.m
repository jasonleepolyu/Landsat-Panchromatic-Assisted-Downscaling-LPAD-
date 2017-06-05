function [CC,Coeffi_x,Coeffi_y,CovVar_variance] = getCoeffiGlobal(MLow,Panlow,LOCALcase,p)
    %% here we suppose x is the Pan band and y is the MS band
    if ~exist('LOCALcase')
        LOCALcase = 'CovVar';
    end
    [ni,nj,nb] = size(MLow);
    CC = zeros(nb,1,'single');
    CovVar_variance = zeros(ni,nj,nb,'single');
    Coeffi_y = CC;
    Coeffi_x = CC;
    %% input parameters
    method_variance = 0;% the methods used to calculate the variance and covariance (used by the coefficient), which can be 
    %:0 (no calculation), 1 (normal calculation)
    if ~exist('p')
        p = 0.4;
        p = 0.5;
%         p = 0.75;
%         p = 0.9;


    end
    if strcmp(LOCALcase,'GS')||strcmp(LOCALcase,'GS3')
        p = 0.5;
    end
    
    if strcmp(LOCALcase,'CovVar_G2')
        method_variance = 1;
    end
    
    if 0 %%test by using the up-sampled images as the mean
        Combine = zeros(ni,nj,nb+1,'single'); %#ok<UNRCH>
        Combine(:,:,1) = Panlow;
        for i=1:nb
            Combine(:,:,i+1) = MLow(:,:,i);
        end
        FACTOR = 4;
        ridus_or_level = 2;
        low_pass_method = 3;%%filtering and decimation
        
        covxyT = get_variance_with_resampled_mean(Combine,FACTOR,low_pass_method,ridus_or_level);
    else
        X = hyperConvert2d(Panlow); 
        Y = hyperConvert2d(MLow);
        Y = [X' Y'];
        covxyT = cov(Y);
    end
    
    vx = covxyT(1,1);
    %% starting to calculate the coefficient
    for k=1:nb               
        vy = covxyT(k+1,k+1);
        covxy = covxyT(1,k+1);
        r = covxy/sqrt(vx*vy);
        CC(k) = r;
        Coeffi_x(k) = covxy/vy;
        Coeffi_y(k) = covxy/vx;
        if p~=0.5%%&&method==1
            r2 = r*r;
            temp = p/(1-p+(2*p-1)*r2);      % Eq. 13-3
            Coeffi_y(k) = Coeffi_y(k)*temp; % Eq. 13-3
        end
        if method_variance==1
%             covxy2_vy = covxy*Coeffi_y(k);
%             vx_covxy2_vy = vy-covxy2_vy;% Cg in Eq. 7-7 when p is equaling to 0.5
            
            beta = Coeffi_x(k);% beta in Eq. 7-7 or Eq. 4-3
%             covxy2_vy = beta*vy*beta;
            vx_covxy2_vy = vx-beta*vy*beta;% Cg in Eq. 7-7
            beta_vx_covxy2_vy_beta = beta/vx_covxy2_vy*beta;% (beta)TC-1beta in Eq. 7-7
            CovVar_variance(:,:,k) = 2*(1-p)/vy+2*p*beta_vx_covxy2_vy_beta;% Eq. 7-7
        end
    end   
    fprintf('The beta value of each band:\t');
    fprintf('%5.3f\t',p,Coeffi_y,mean(Coeffi_y));
    fprintf('\n');
    if ~strcmp(LOCALcase,'GS')&&~strcmp(LOCALcase,'GS3')
        fprintf('The CC value of each band:\n');
        fprintf('%5.3f\t',p,CC,mean(CC));
        fprintf('\n');
    end
end