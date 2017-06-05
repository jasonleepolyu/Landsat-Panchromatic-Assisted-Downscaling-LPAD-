function [M] = hyperConvert2d(M)
% HYPERCONVERT2D Converts an HSI cube to a 2D matrix
% Converts a 3D HSI cube (m x n x p) to a 2D matrix of points (p X N)
% where N = mn
%
% Usage
%   [M] = hyperConvert2d(M)
% Inputs
%   M - 3D HSI cube (m x n x p)
% Outputs
%   M - 2D data matrix (p x N)

% if (ndims(M) ~= 3)
%     error('Input image must be m x n x p.');
% end
% if(ndims(M) == 3)
%     [h, w, numBands] = size(M);
%     M = reshape(M, w*h, numBands).';
% else
%     [h, w, numBands] = size(M);
%     M = reshape(M, w*h, numBands).';
% end
    [h, w, numBands] = size(M);
    M = reshape(M, w*h, numBands).';

return;