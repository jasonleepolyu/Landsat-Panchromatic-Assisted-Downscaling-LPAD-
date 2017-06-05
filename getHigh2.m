function M = getHigh2(Mtrue,FACTOR,low_pass_method)

    [m,n,bands] = size(Mtrue);
    M = zeros(m*FACTOR,n*FACTOR,bands,'single');
    if(low_pass_method==4)%% ideal filter odd 23-taps filter
        B3 = [-0.000060081482 0 0.000807762146 0 -0.005192756653 0 .... 
            0.021809577942 0 -0.072698593239 0 ....
            0.305334091185 0.5 0.305334091185 ....
            0 -0.072698593239 0 0.021809577942 ....
            0 -0.005192756653 0 0.000807762146 0 -0.000060081482];
        h = B3'*B3;
    elseif(low_pass_method==6)%% ideal filter even (Table 2 in Aiazzi et al. 2013)
        B3 = [-0.0234375 -0.0703125 0.2265625 0.8671875 0.8671875 0.2265625 -0.0703125 -0.0234375]/2; %% no SHIFT
        h = B3'*B3;
    elseif(low_pass_method==7)%% ideal filter odd (Table 2 in Aiazzi et al. 2013) 
        B3 = [-1/16 0 9/16 1 9/16 0 -1/16]/2; % half-pixel SHIFT towards lower right
        h = B3'*B3;
    end
    
   %% start to prefilter and decimation 2 by 2
    levelT = floor(log2(FACTOR));
    for b=1:bands
        A = Mtrue(:,:,b);
        for level=1:levelT
            h2 = h*4;
            [height, width] = size(A);

            [A,indexi,indexj] = expanding(A,h);
            Atemp = A;
            [height2, width2] = size(A);
            A = zeros(height2*2,width2*2);
            A(1:2:height2*2,1:2:width2*2) = Atemp;
            A = conv2(A, h2, 'same');
            A = A((indexi(1)*2-1):(indexi(height)*2), (indexj(1)*2-1):(indexj(width)*2));
        end            
        M(:,:,b) = A;
    end
%%
end