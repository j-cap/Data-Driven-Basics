% calcResponsibilities - Calculates responsibilities 
%
% According to Deisenroth "Mathematics for Machine Learning" p. 318
% Currently brute force iteration 
%        --> potential for optimization
%
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
function r_nk = calcResponsibilities(X, WEIGHTS, MU, SIGMA)

    r_nk = zeros(length(X), length(WEIGHTS));
    for i=1:length(X)
        for j=1:length(WEIGHTS)
            numerator = WEIGHTS(j) * norm_pdf(X(i), MU(j), SIGMA(j));
            denominator = 0;
            for k=1:length(WEIGHTS)
                denominator = denominator + WEIGHTS(k) * norm_pdf(X(i), MU(k), SIGMA(k));
            end
            r_nk(i,j) = numerator / denominator;
        end
    end
    
end
%------------- END OF CODE -----------------------------------------------------
