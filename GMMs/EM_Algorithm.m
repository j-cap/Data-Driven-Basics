% EM_Algorithm - Expecation-Maximization Algorithm
%
% According to Deisenroth "Mathematics for Machine Learning" p. 318
% Currently brute force iteration 
%        --> potential for optimization
%
% Optional file header info (to give more details about the function than in the H1 line)
% Optional file header info (to give more details about the function than in the H1 line)
% Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [output1, output2] = functionName(input1, input2, input3)
%
% Inputs:
%    input1 - Description (Dimension, Unit, Range, ...)
%    input2 - Description (Dimension, Unit, Range, ...)
%    input3 - Description (Dimension, Unit, Range, ...)
%
% Outputs:
%    output1 - Description (Dimension, Unit, Range, ...)
%    output2 - Description (Dimension, Unit, Range, ...)
%
% Example:
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Subfunctions: none
%
% Dependencies:
%    Matlab release: 2019b
%    other m-files: none
%    MAT-files: none
%    Toolboxes: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%
% This function is part of: project / package etc.
%
% Author:  Jakob Weber
% email:   jakob.weber@ait.ac.at
% Company: Austrian Institute of Technology GmbH
%          Complex Dynamical Systems
%          Center for Vision, Automation & Control
%          http://www.ait.ac.at
%
% Version: x.y.z - 2021-07-26

% Change log:
% x.y.z - 2021-07-26 - author:
% - added important feature, s. issue #34
% - fixed bug #2

%------------- BEGIN CODE ------------------------------------------------------
function [mu_k, sigma_k, weights_k, R_NK] = EM_Algorithm(X, MU, SIGMA, WEIGHTS)

% E-step    
R_NK = calcResponsibilities(X, WEIGHTS, MU, SIGMA);
    
% M-step
N_k = sum(R_NK, 1); % calc column sum of responsibilities
mu_k = sum(R_NK .* X') ./ N_k; % update mean values

% update the covariances
sigma_k = zeros(1,3);
for k=1:length(WEIGHTS)
    iter = 0;
    for n=1:length(X)
        iter = iter + (R_NK(n,k) * (X(n) - mu_k(k)) * (X(n) - mu_k(k)));
    end
    if iter == 0
        iter = eps;
    end
    sigma_k(k) = iter / N_k(k);
end 
% update the weights
weights_k = N_k / length(X);

en
%------------- END OF CODE -----------------------------------------------------
