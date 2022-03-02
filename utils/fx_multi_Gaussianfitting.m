function [amp, parameter, yhat] = fx_multi_Gaussianfitting(numGaussians, x, y, centers, sigmas, amplitudes)
% using multiply Gaussian distributino to fit the 1D data
  
% Put all the parameters into a table for convenience in looking at, and using, the results.
tActual = table((1:numGaussians)', amplitudes(:), centers(:), sigmas(:), 'VariableNames', {'Number', 'Amplitude', 'Mean', 'Width'});
% Now sort parameters in order of increasing mean, just so it's easier to think about (though it's not required).
tActual = sortrows(tActual, 3);
tActual.Number = (1:numGaussians)'; % Unsort the first column of numbers.

%----------------------------------------------------------------------------------------------------------------------------------
% Now we have our test signal and we can begin....
% Fit Gaussian Peaks:
% Initial Gaussian Parameters
initialGuesses = [tActual.Mean(:), tActual.Width(:)];
% Add a little noise so that our first guess is not dead on accurate.
initialGuesses = initialGuesses + 2 * rand(size(initialGuesses));
startingGuesses = reshape(initialGuesses', 1, []);

global c

% 	warning off
% t and y must be row vectors.
tFit = reshape(x, 1, []);
y    = reshape(y, 1, []);

%-------------------------------------------------------------------------------------------------------------------------------------------
% Perform an iterative fit using the FMINSEARCH function to optimize the height, width and center of the multiple Gaussians.
options         = optimset('MaxFunEvals', 10^12);  % Determines how close the model must fit the data
% First, set some options for fminsearch().
options.TolFun  = 1e-5;
options.TolX    = 1e-5;
options.MaxIter = 10000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HEAVY LIFTING DONE RIGHT HERE:
% Run optimization
[parameter, ~, ~, ~] = fminsearch(@(lambda)( ...
    fitgauss(lambda, tFit, y)), startingGuesses, options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

amp = c;

% Get the means and widths.
means  = parameter(1 : 2 : end);
widths = parameter(2 : 2 : end);
% Now plot results.
yhat = zeros(1, length(tFit));
numGaussians = length(c);
for k = 1 : numGaussians
    % Get each component curve.
    thisEstimatedCurve = c(k) .* gaussian(tFit, means(k), widths(k));
    yhat = yhat + thisEstimatedCurve;
end

%=======================================================================================================================================================
function theError = fitgauss(lambda, t, y)
    % Fitting function for multiple overlapping Gaussians, with statements
    % added (lines 18 and 19) to slow the progress and plot each step along the
    % way, for educational purposes.
    % Author: T. C. O'Haver, 2006
    A = zeros(length(t), round(length(lambda) / 2));
    for j = 1 : length(lambda) / 2
        if lambda(2 * j - 1) >= 180 % ensure the centers between [0 180] 
            lambda(2 * j - 1) = 180-1;
        elseif lambda(2 * j - 1) < 0
            lambda(2 * j - 1) = 0+1;
        end
        A(:,j) = gaussian(t, lambda(2 * j - 1), lambda(2 * j))';
    end
    c = A \ y';
    z = A * c;
    theError = norm(z - y');
    theError = theError + 1e-4 * sum(c); % plus penalty terms to avoid too large c
    if sum(c < 0) > 0
        theError = theError + 1000000;
    end
end % of fitgauss()


%=======================================================================================================================================================
function g = gaussian(x, peakPosition, width)
%  gaussian(x,pos,wid) = gaussian peak centered on pos, half-width=wid
%  x may be scalar, vector, or matrix, pos and wid both scalar
%  T. C. O'Haver, 1988
% Examples: gaussian([0 1 2],1,2) gives result [0.5000    1.0000    0.5000]
% plot(gaussian([1:100],50,20)) displays gaussian band centered at 50 with width 20.
g = exp(-((x - peakPosition) ./ (0.60056120439323 .* width)) .^ 2);
end % of gaussian()

end

