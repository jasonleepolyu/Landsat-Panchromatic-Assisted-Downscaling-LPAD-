%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% cubic convolution resampling           %%%%%
%%% Zhongbin Li, SDSU, SD, USA, May-20-2017  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

%% Input original image
% info = geotiffinfo('landslide/30m_B1_B7.tif');
[I0, R] = geotiffread('30m_B1_B7.tif');

info = geotiffinfo('20m_B4.tif');
[temp, R] = geotiffread('20m_B4.tif');

b1 = I0(:,:,1);
b2 = I0(:,:,2);
b3 = I0(:,:,3);
b4 = I0(:,:,4);
b5 = I0(:,:,5);
b6 = I0(:,:,6);
b7 = I0(:,:,7);
%% Input scale factor
K = 30/20; % 30m to 20m or 15m to 20m

%% cubic convolution resampling
b1 = bl_cc_resampling_at_30(b1,K);
b2 = bl_cc_resampling_at_30(b2,K);
b3 = bl_cc_resampling_at_30(b3,K);
b4 = bl_cc_resampling_at_30(b4,K);
b5 = bl_cc_resampling_at_30(b5,K);
b6 = bl_cc_resampling_at_30(b6,K);
b7 = bl_cc_resampling_at_30(b7,K);

%% Output
Output(:,:,1) = b1;
Output(:,:,2) = b2;
Output(:,:,3) = b3;
Output(:,:,4) = b4;
Output(:,:,5) = b5;
Output(:,:,6) = b6;
Output(:,:,7) = b7;
[row, col, band] = size(b5)
R.RasterSize = [row col];
% geotiffwrite('test_cc_registered_20m_from_30m.tif',int16(Output * 10000),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite('test_bl_registered_20m_from_30m.tif',int16(Output * 10000),R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);

% Output = uint8(Output(:,:,[7,5,4]));
% figure, imshow(Output(:,:,[7,5,4]));



