% sGMM - Gaussian Mixture Models 
%
% Configuration / Parameters:
%    config1 - Description (Dimension, Unit, Range, ...)
%    config2 - Description (Dimension, Unit, Range, ...)
%    config3 - Description (Dimension, Unit, Range, ...)
%
% Results / Outputs:
%    output1 - Description (Dimension, Unit, Range, ...)
%    output2 - Description (Dimension, Unit, Range, ...)
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
% email:   Jakob.Weber@ait.ac.at
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

%% Create some training data 
% p(x|theta) = 0.33 N(x|-4,1) + 0.33 N(x|0, 0.2) + 0.33 N(x|8, 3);

mu = [-4, 0, 8];
sigmas = [1, 0.45, 3];
weights = ones(3,1) * 1/3;

p1 = @(x) weights(1)*norm_pdf(x, mu(1), sigmas(1));
p2 = @(x) weights(2)*norm_pdf(x, mu(2), sigmas(2));
p3 = @(x) weights(3)*norm_pdf(x, mu(3), sigmas(3));

xplot = linspace(-10, 15, 1000);
X = [-3, -2.5, -1, 0, 2, 4, 5];
Y = p1(X) + p2(X) + p3(X) ;

fig = figure(); fig.Position = [-1600, 200, 1200, 700]; hold on;
ylim([-0.05, 0.4]); xlim([-10, 15]); grid on;
scatter(X, zeros(size(X)), 60, 'black', 'o', 'filled')
plot(xplot, p1(xplot), '-.blue');
plot(xplot, p2(xplot), '-.red');
plot(xplot, p3(xplot), '-.green');
plot(xplot, p1(xplot)+p2(xplot)+p3(xplot), 'black');


%%
mu_k = mu; 
sigma_k = sigmas;
weights_k = weights;
max_iter = 3;
MUs = zeros(max_iter, length(mu_k));
SIGMAs = zeros(max_iter, length(sigma_k));
WEIGHTs = zeros(max_iter, length(weights));
for i=1:max_iter
   [mu_k, sigma_k, weights_k, R_NK] = EM_Algorithm(X, mu_k, sigma_k, weights_k)
   MUs(i, :) = mu_k;
   SIGMAs(i,:) = sigma_k;
   WEIGHTs(i,:) = weights_k;
end
%%
fig = figure(); fig.Position = [-1600, 200, 1200, 700]; hold on;
ylim([-0.05, 0.4]); xlim([-10, 15]); grid on;

idx = 1;
p1fit = norm_pdf(xplot, MUs(idx, 1), SIGMAs(idx, 1))*WEIGHTs(idx, 1);
p2fit = norm_pdf(xplot, MUs(idx, 2), SIGMAs(idx, 2))*WEIGHTs(idx, 2);
p3fit = norm_pdf(xplot, MUs(idx, 3), SIGMAs(idx, 3))*WEIGHTs(idx, 3);

scatter(X, zeros(size(X)), 60, 'black', 'o', 'filled')
plot(xplot, p1fit, '-.blue');
plot(xplot, p2fit, '-.red');
plot(xplot, p3fit, '-.green');
plot(xplot, p1fit+p2fit+p3fit, 'black');


%%

d


%%

 

%------------- END OF CODE -----------------------------------------------------
