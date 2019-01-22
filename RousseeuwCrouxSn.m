function [Sn, x_j] = RousseeuwCrouxSn(X)
% Compute the measure of scale 'Sn', from Rousseeuw & Croux (1993).
%
%   X should be an Nx1 column vector, or a an NxM matrix, where M is the
%   number of data dimensions (e.g., Nx2 for xy data). 
%
%   A robust alternative to MAD for statistical outlier identification.
%   Unlike MAD, Sn does not make an assumption of symmetry, so in
%   principle should be more robust to skewed distributions.
%
%   The outputs of this function have been validated against equivalent
%   function in MAPLE(tm).
%
% Example:          
%                   % basic example
%                   X = [1 5 2 2 7 4 1 6]';
%                   Sn = RousseeuwCrouxSn(X) % should give 3.015
%
%                   % use Sn to identify statistical outliers
%                	X = [1 5 2 2 7 50 1 5]';
%                  	[Sn, x_j] = RousseeuwCrouxSn(X);
%                	outliers = X(x_j/Sn > 3) % NB: criterion typically 2 or 3
%
%                   % multidimensional data
%                	X = [1 5 2 2 7 50 1 5; 2 6 3 1 6 48 2 7]';
%                  	[Sn, x_j] = RousseeuwCrouxSn(X);
%                	outliers = X(x_j/Sn > 3, :) % NB: criterion typically 2 or 3
%
% Requires:         none
%
% See also:         mad.m [Statistics Toolbox]
%
% Author(s):        Pete R Jones <petejonze@gmail.com>
% 
% Version History:  19/04/2016	PJ  Initial version
%                   09/08/2018	PJ  Added support for multi-dimensional inputs
%                   22/01/2019	PJ  Added error checking for row vectors and simplified 
%                                               
%   
% Copyright 2019 : P R Jones
% *********************************************************************
% 

    % (defensive) convert row vector to column
    if size(X,1)==1, X = X'; end
    
    % get number of elements
    n = size(X,1);
    
    % Set c: bias correction factor for finite sample size. NB: the values
    % used here match those used in the MAPLE implementation of Sn. For
    % more regarding the computation of the finite sample correction 
    % factors, see Pison, Aelst & Willems (2002), Metrika, 55(1), 111-123.
    if n < 10
        cc = [NaN 0.743 1.851 0.954 1.351 0.993 1.198 1.005 1.131];
        c = cc(n);
    elseif mod(n,2)==1  % n is odd
        c = n/(n-.9);
    else                % n is even
        c = 1;
    end
    
    % compute median distance of each element to all other elements
    x_j = nan(n,1);
    for i = 1:n
        X_other = X([1:i-1 i+1:end], :); % get all values except the current one
        d = sqrt(sum(bsxfun(@minus, X_other, X(i,:)).^2, 2)); % compute distance
        x_j(i) = median(d); % compute median distance
    end

    % compute median of all median differences, and apply finite sample
    % correction, c
    Sn = c * median(x_j);
end