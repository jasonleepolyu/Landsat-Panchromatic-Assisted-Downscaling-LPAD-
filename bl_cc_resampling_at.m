function Output = bl_cc_resampling_at(I, scale_factor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            I: the original image or bands
%%%            scale_factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pad image by duplicate the last row and column
[nrow, ncol] = size(I);
I = [I; I(nrow, :); I(nrow, :); I(nrow, :);];  
I = [I I(:, ncol) I(:, ncol) I(:, ncol)];  

%% Initialize the output image
nrow_new = nrow*scale_factor; 
ncol_new = ncol*scale_factor;
Output = zeros(nrow_new, ncol_new);

%% Ratio of output image to the original image
y_ratio = (nrow)/(nrow*scale_factor);   
x_ratio = (ncol)/(ncol*scale_factor);  

%% BL and CC Affine transformation coefficients
% landslide  20m from 15m
% a0 = -1.346995354+12.5/20;
% a1 = 1.000008878;
% a2 = 0.000033128;
% b0 = -1.218221657+12.5/20;
% b1 = -0.000042180;
% b2 = 1.000012239;

% New_york

a0 = -0.771690110 +12.5/20; % +12.5/20  col
a1 = 0.999915075;
a2 = -0.000189602;
b0 = 1.078519619 +2.5/20; % +2.5/20 row
b1 = -0.000064901;
b2 = 0.999881731;

% crop field

% a0 = -1.007142445 +12.5/20;
% a1 = 1.000007977;
% a2 = 0.000012005;
% b0 = 0.225335187 +2.5/20;
% b1 = 0.000015280;
% b2 = 0.999996748;

% burned_area

% a0 = -0.772688990 +2.5/20;
% a1 = 1.000018270;
% a2 = 0.000036507;
% b0 = 0.489724920 +2.5/20;
% b1 = -0.000030781;
% b2 = 1.000013579;

% forest

% a0 = 0.184715799 +2.5/20;
% a1 = 0.999815008;
% a2 = -0.000100652;
% b0 = 0.221510137 +2.5/20;
% b1 = -0.000149049;
% b2 = 0.999983763;

for j = 1:nrow_new
    
    for i = 1:ncol_new
        
        m = a0 + a1 * i + a2 * j;  % col    affine transformation
        n = b0 + b1 * i + b2 * j;  % row    affine transformation
        
%         m = i;  % col    affine transformation
%         n = j;  % row    affine transformation

        if (m >=1 && m < ncol_new  && n >=1 && n < nrow_new )
        
%             YY = floor((n-1)*y_ratio+1); % row
%             y = (n-1)*y_ratio+1 - YY;
% 
%             XX = floor((m-1)*x_ratio+1); % col
%             x = (m-1)*x_ratio+1 - XX;
            
            floaty = ((n-0.5)*20 -2.5)/15+1;
            if (floaty - floor(floaty))>=0.5
                YY = floor(floaty); % row
                y = floaty - YY - 0.5;
            else
                YY = floor(floaty)-1; % row
                y = floaty - YY - 0.5 ;
            end
            
            floatx = ((m-0.5)*20 -12.5)/15+1;
            if (floatx - floor(floatx))>=0.5
                XX = floor(floatx); % col
                x = floatx - XX - 0.5;
            else
                XX = floor(floatx)-1; % col
                x = floatx - XX - 0.5;
            end
            if y<0 || x<0
                fprintf('shit this is wrong ');
            end


            if (YY >1 && YY < nrow  &&  XX >1 && XX < ncol )
                
                %% this is bi-cubic
                A = [cc_kernel(y + 1) cc_kernel(y + 0) cc_kernel(y - 1) cc_kernel(y - 2)];
                
                B = [I(YY-1, XX-1) I(YY-1, XX+0) I(YY-1, XX+1) I(YY-1, XX+2);
                     I(YY+0, XX-1) I(YY+0, XX+0) I(YY+0, XX+1) I(YY+0, XX+2);
                     I(YY+1, XX-1) I(YY+1, XX+0) I(YY+1, XX+1) I(YY+1, XX+2);
                     I(YY+2, XX-1) I(YY+2, XX+0) I(YY+2, XX+1) I(YY+2, XX+2);];
                 
                C = [cc_kernel(x + 1) cc_kernel(x + 0) cc_kernel(x - 1) cc_kernel(x - 2)]';
                                
                Output(j,i) = A*B*C; 
                
                %% this is bilinear
%                 Output(j,i) = (I(YY+0, XX+0)*(1-x)+I(YY+0, XX+1)*x)*(1-y)+ (I(YY+1, XX+0)*(1-x)+I(YY+1, XX+1)*x)*y; 
            end
        end
    end
end

end

function cc = cc_kernel(x)

    % Equation (15) in Keys's paper "Cubic Convolution Interpolation for Digital Image Processing, 1981"
    
    if abs(x)>=0 && abs(x)< 1
        cc = 1.5 * abs(x)^3 - 2.5 * abs(x)^2 + 1;
    end
    
    if abs(x)>=1 && abs(x)< 2
        cc = -0.5 * abs(x)^3 + 2.5 * abs(x)^2 - 4 * abs(x) + 2;
    end
    
    if abs(x)>=2
        cc = 0;
    end
    

end 


%                 Output(j,i) = I(YY+0, XX+0); 
%                 YY = floor(floaty); % row
%                 y = floaty - YY - 0.5;
%             YY = floor(((n)*20-12.5)/15); % row
%             y = ((n)*20-12.5)/15 - YY;
% 
%             XX = floor(((m)*20-12.5)/15); % col
%             x = ((m)*20-12.5)/15 - XX;
