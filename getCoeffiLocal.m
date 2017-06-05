function [AABP,ECB,CovVar,CovVar_variance,CorreLow] = getCoeffiLocal(MLow,PanLow,windowsize,LOCALcase)
    %% here we suppose x is the Pan band and y is the MS band
    if ~exist('LOCALcase')
        LOCALcase = 'CovVar';
    end
    clippingConstant = 3;
    clippingConstantMin = 0;
    [ni,nj,nb] = size(MLow);
    AABP = zeros(ni,nj,nb,'single');
    ECB = AABP;
    CorreLow = AABP;
    CovVar = ECB;
    CovVar_variance = ECB;
    %% input parameters
    method = 0;
    method_variance = 0;% the methods used to calculate the variance and covariance (used by the coefficient), which can be 
    %:0 (no calculation), 1 (normal calculation), or 2 (calculation with up-sampled image as the mean)
    if strcmp(LOCALcase,'AABP')||strcmp(LOCALcase,'AABP2')
        method = 1;
    elseif strcmp(LOCALcase,'CovVar')||strcmp(LOCALcase,'CovVar2')||strcmp(LOCALcase,'HRCovVar')...
            ||strcmp(LOCALcase,'HRCovVar2')||strcmp(LOCALcase,'HRCovVar_Diff')||strcmp(LOCALcase,'GSLocal')
        method = 2;
    elseif strcmp(LOCALcase,'ECB')||strcmp(LOCALcase,'ECB2')
        method = 3;
    end
    if strcmp(LOCALcase,'AABP2')||strcmp(LOCALcase,'CovVar2')||strcmp(LOCALcase,'ECB2')||strcmp(LOCALcase,'HRCovVar2')
        method_variance = 1;
    end

    if strcmp(LOCALcase,'HRCovVar_Diff')%%using the up-sampled images as the mean
        method_variance = 2;
        FACTOR = 4;
        ridus_or_level = 2;
        low_pass_method = 3;%%filtering and decimation
        [PanLow,MLow] = get_difference_with_resampled_mean(PanLow,MLow,FACTOR,low_pass_method,ridus_or_level);
    else
        [CC_Global,Coeffi_Global] = getCoeffiGlobal(MLow,PanLow,LOCALcase);
        beta_cliped = 0.5*Coeffi_Global;
    end
    p = 0.5;
    %%
    critical_min = 0;
    critical_max = 1;
%     error = 0.1;
%     critical_min = critical_min+error;
%     critical_max = critical_max-error;
    averagep = zeros(1,nb,'single');
    %% starting to calculate the coefficient
    for i=1:1:ni
         for j=1:1:nj
            starti = max(1,i-windowsize);
            endi = min(ni,i+windowsize);
            startj = max(1,j-windowsize);
            endj = min(nj,j+windowsize);
            X = PanLow(starti:endi,startj:endj);
            [tempm,tempn] = size(X);
            n = tempm*tempn;
            X = reshape(X,n,1);

            for k=1:nb               
                Y = reshape(MLow(starti:endi,startj:endj,k),n,1);

              %% different ways to calculate the coefficient
                if method_variance==2
                    vx = mean(X.*X);
                    vy = mean(Y.*Y);
                    covxy = mean(X.*Y);
                else
                    covxy = cov(X,Y);
                    vx = covxy(1,1);
                    vy = covxy(2,2);
                    covxy = covxy(1,2);
                end
              %% initialize the coefficient 
                if method>=2
                    coffi = covxy/vx;
                end
              %% calculate the coefficient 
                switch method 
                    case 1 %AABP
                    coffi = min(sqrt(vy)/(sqrt(vx)+1),clippingConstant);
                    AABP(i,j,k) = coffi;
                    if method_variance==1
                        R = covxy/(sqrt(vx*vy));
                        R1_R = 1/R-R;
                        p = R1_R/((R1_R-R)+1);% the equation below Eq 11-2
                     %% reset the p values
                        p = max(critical_min,p);
                        p = min(critical_max,p);
                     %% get the variance from p
                        beta = covxy/vy;% beta in Eq. 7-7 or Eq. 4-3
                        if abs(beta)<beta_cliped(k)
                            beta = beta_cliped(k);
                        end
                        CovVar_variance(i,j,k) = 2*(1-p)/vy+2*p*(beta/(R1_R*R*vx)*beta);% Eq. 7-7
                    end
                    case 2 %M3
                    if p~=0.5&&p~=1
                        r = covxy/(sqrt(vx*vy));
                        r2 = r*r;
                        temp = p/(1-p+(2*p-1)*r2);     % Eq. 13-3
                        coffi = coffi*temp;            % Eq. 13-3
                    elseif p==1
                        coffi = vy/covxy;
                    end
                    coffi = min(coffi,clippingConstant);
                    coffi = max(coffi,clippingConstantMin);
                    CovVar(i,j,k) = coffi;
                    if method_variance==1% this variance is ONLY valid for p = 0.5
                        CovVar_variance(i,j,k) = 1/(vy-covxy*covxy/vx);
                    end
                    case 3  %ECB
                    coffiE = coffi/CC_Global(k);
                    coffiE = min(coffiE,clippingConstant);
                    coffiE = max(coffiE,clippingConstantMin);
                    ECB(i,j,k) = coffiE;
                    if method_variance==1
                        R2 = covxy*covxy/(vx*vy);
                        one_R2 = 1-R2;
                        denominator = (one_R2-R2)+CC_Global(k);
                        p = one_R2/denominator;% the equation below Eq 11-2
                     %% reset the p values
                        p = max(critical_min,p);
                        p = min(critical_max,p);
                        averagep(k) = averagep(k)+p;
                     %%
                        vx_covxy2_vy = vx*one_R2;% Cg in Eq. 7-7
                        beta = covxy/vy;% beta in Eq. 7-7 or Eq. 4-3
                        if abs(beta)<beta_cliped(k)
                            beta = beta_cliped(k);
                        end
                        beta_vx_covxy2_vy_beta = beta/vx_covxy2_vy*beta;% (beta)TC-1beta in Eq. 7-7
                        CovVar_variance(i,j,k) = 2*(1-p)/vy+2*p*beta_vx_covxy2_vy_beta;
                    end
                end
            end%%end of band
        end
    end
%     if method~=2 %M3
%         fprintf('The average p value of each band: ');
%         fprintf('%5.3f\t',averagep/(ni*nj),mean(averagep/(ni*nj)));
%         fprintf('\n');
%     else %ECB or AABP
%         fprintf('The p value: ');
%         fprintf('%5.3f\t',p);
%         fprintf('\n');  
%     end
end