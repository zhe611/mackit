function data = mixture_normal(means, sigmas, weights, n)
% This function is to generate random variables from a Gaussian mixture distribution
% Input arguments:
%   means:   vector of means of the Gaussian components, e.g., [mu1, mu2, mu3]
%   sigmas:  vector of standard deviations of the Gaussian components,
%            e.g., [sigma1, sigma2, sigma3]
%   weights: vector of mixture weights (summing to one), e.g., [w1, w2, w3]
%   n:       number of samples to generate
% Output:
%   data:    random samples drawn from the Gaussian mixture distribution

    % Validate input arguments
    if length(means) ~= length(sigmas) || length(means) ~= length(weights)
        error('The lengths of the mean, standard deviation, and weight vectors must be equal.');
    end
    
    if abs(sum(weights) - 1) > 1e-10
        error('The sum of the weights must be equal to one.');
    end
    
    % Generate component indices
    k = length(means); % number of components
    component_idx = randsample(1:k, n, true, weights);
    
    % Generate Gaussian samples for each component
    data = zeros(n, 1);
    for i = 1:k
        mask = (component_idx == i);
        num_samples = sum(mask);
        if num_samples > 0
            data(mask) = normrnd(means(i), sigmas(i), num_samples, 1);
        end
    end
end
