%%
%
%25.08.2021
%
% Chapter 2 of Bishop: Pattern Recognition and Machine Learning
%
%
%% Bernoulli Distribution

Bernoulli =@(x, mu) mu.^x .* (1 - mu).^(1 - x);

x = randi(2, 1000, 1) - 1;
y = Bernoulli(x, 1);

plot(y);


%%