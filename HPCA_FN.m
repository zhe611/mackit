function rhat = HPCA_FN(X, rmax, Method, tau, scale_est, threshold, L_init, F_init, maxiter_HPCA, maxiter_HLM, eps)

% This function is to estimate factor numbers via rank minimization corresponding to Huber Principal Component Analysis (HPCA).

% X is the input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.
% rmax is a positive integer indicating the maximum number of factors.
% Method="P" indicates minimizing the Huber loss of the idiosyncratic error's norm (HPCA) 
% Method="E" indicates minimizing the elementwise Huber loss.(IHR)
% tau is the user-supplied parameter for Huber loss; default is NULL, and \tau is provided by default.
% scale est is a parameter for the elementwise Huber loss.
% threshold is threshold of rank minimization; default is NULL.
% L_int is the user-supplied inital value of factors; default is the PCA estimator.
% F_int is the user-supplied inital value of factors; default is the PCA estimator.
% maxiter_HLM is the maximum number of iterations in the HPCA. The default is 100.
% maxiter_HLM is the maximum number of iterations in the iterative Huber regression algorithm. The default is 100.
% eps is stopping critetion parameter in the HPCA. The default is 1e-3.

% Reference: He, Y., Li, L., Liu, D., and Zhou, W.-X. (2025). Huber principal component analysis for large dimensional factor models. 


% Set default parameters
    if nargin < 11
        eps = 0.001;
    end
    if nargin < 10
        maxiter_HLM = 100;
    end
    if nargin < 9
        maxiter_HPCA = 100;
    end
    if nargin < 8
        F_init = [];
    end
    if nargin < 7
        L_init = [];
    end
    if nargin < 6
        threshold = [];
    end
    if nargin < 5
        scale_est = 'MAD'; %HPCA
    end
    if nargin < 4
        tau = [];
    end
    if nargin < 3
        Method = 'E';
    end

    % Get dimensions
    [T, N] = size(X);

    % Fit HPCA
    fit = HPCA(X, rmax, Method, tau, scale_est, L_init, F_init, maxiter_HPCA, maxiter_HLM, eps);

    % Compute eigenvalues
    VK = sort(diag(fit.Fhat' * fit.Fhat) / T, 'descend');

    % Determine threshold
    if isempty(threshold)
        threshold = VK(1) * (min(N, T))^(-1/3);
    end

    % Estimate rank
    rhat = sum(VK > threshold);
end

