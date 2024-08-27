% SPDX-FileCopyrightText: 2021-2024 Helmholtz-Zentrum hereon GmbH
% SPDX-FileContributor: Carsten Lemmen  <carsten.lemmen@hereon.de>
% SPDX-License-Identifier: CC0-1.0 

function s = nanstd(x,dim)

% Find NaNs and set them to zero
x(x>1E+9) = NaN;
x(x<-1E+9) = NaN;

nans = isnan(x);
x(nans) = 0;

if nargin == 1 % let sum deal with figuring out which dimension to use
    % Count up non-NaNs.
    n = sum(~nans);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x) ./ n;
    s = sqrt(sum(x.*x) ./ n - m.*m);

else
    % Count up non-NaNs.
    n = sum(~nans,dim);
    n(n==0) = NaN; % prevent divideByZero warnings
    % Sum up non-NaNs, and divide by the number of non-NaNs.
    m = sum(x,dim) ./ n;
    s = sqrt(sum(x.*x,dim) ./ n - m.*m);

end
