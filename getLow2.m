function M = getLow2(Mtrue,FACTOR,low_pass_method)

    %%
    [m,n,bands] = size(Mtrue);
    M = zeros(m/FACTOR,n/FACTOR,bands,'single');

    if(low_pass_method==3)%% a trous wavelet filter
        B3 = [1/16 1/4 3/8 1/4 1/16];
        h = B3'*B3;
    elseif(low_pass_method==4)%% 23-taps filter
        B3 = [-0.000060081482 0 0.000807762146 0 -0.005192756653 0 .... 
            0.021809577942 0 -0.072698593239 0 ....
            0.305334091185 0.5 0.305334091185 ....
            0 -0.072698593239 0 0.021809577942 ....
            0 -0.005192756653 0 0.000807762146 0 -0.000060081482];
        h = B3'*B3;
    elseif(low_pass_method==5)%% a trous wavelet filter even
        B3 = [0.0157    0.1523    0.3320    0.3320    0.1523  0.0157];
        h = B3'*B3;
    elseif(low_pass_method==6)%% ideal filter even (Table 2 in Aiazzi et al. 2013)
        B3 = [-0.0234375 -0.0703125 0.2265625 0.8671875 0.8671875 0.2265625 -0.0703125 -0.0234375]/2;
        h = B3'*B3;
    elseif(low_pass_method==7)%% ideal filter odd (Table 2 in Aiazzi et al. 2013) 
        B3 = [-1/16 0 9/16 1 9/16 0 -1/16]/2; % half-pixel SHIFT towards lower right
        h = B3'*B3;      
    end
    
    %% start to prefilter and decimation 2 by 2
    levelT = floor(log2(FACTOR));
    for b=1:bands
        A = Mtrue(:,:,b);
        if low_pass_method==1||low_pass_method==2
            [height, width] = size(A);
            [A,indexi,indexj] = expanding(A,h);
            A = conv2(A, h, 'same');
            A = A(indexi(1:FACTOR:height),indexj(1:FACTOR:width));

        else
            for level=1:levelT
                [height, width] = size(A);
                [A,indexi,indexj] = expanding(A,h);
                A = conv2(A, h, 'same');
                A = A(indexi(1:2:height),indexj(1:2:width));

            end     
        end
        M(:,:,b) = A;
    end

end