function [A,indexi,indexj] = expanding(Mtrue,g,ridus_or_level)
%This function do the expanding
%   Detailed explanation goes here
    if ~exist('ridus_or_level')
        ridus_or_level = 1;
    end
    [m,n,nb] = size(Mtrue);
    [gm,gn] = size(g);
    gm = gm*ridus_or_level;
    gn = gn*ridus_or_level;
%     A = Mtrue;
    A = zeros(m+2*gm,n+2*gn,nb,'single');
    indexi = gm+1:gm+m;
    indexj = gn+1:gn+n;
    A(indexi,indexj,:) = Mtrue;
%%
%     onesa = ones(gm,1,'single');
%     onesb = ones(1,gn,'single');
%     for k=1:nb
%         A(1:gm,indexj,k) = onesa*Mtrue(1,:,k);
%         A(m+gm+1:m+2*gm,indexj,k) = onesa*(Mtrue(m,:,k));
%         A(:,1:gn,k) = A(:,gn+1,k)*onesb;
%         A(:,n+gn+1:n+2*gn,k) = A(:,gn+n,k)*onesb;
%     end
%%
    for i=1:gm
        A(i,indexj,:) = squeeze(Mtrue(1,:,:));
    end
    for i=m+gm+1:m+2*gm
        A(i,indexj,:) = squeeze(Mtrue(m,:,:));
    end
    for i=1:gn
        A(:,i,:) = squeeze(A(:,gn+1,:));
    end
    for i=n+gn+1:n+2*gn
        A(:,i,:) = squeeze(A(:,gn+n,:));
    end
%%
%     filename = 'A.tif';
%     dataDirLandsat = 'E:\\HongKong\\tif\\';
%     LandsatFile1 = sprintf('%s\\ETM_20011120_1200_1200.tif',dataDirLandsat);
%     infoLandsat = geotiffinfo(LandsatFile1);
%     [Landsat1,RLand] = geotiffread(LandsatFile1);
%     RLand.RasterSize = [m+2*gm n+2*gn];
%     geotiffwrite(filename,A,RLand,'GeoKeyDirectoryTag',infoLandsat.GeoTIFFTags.GeoKeyDirectoryTag);

end

