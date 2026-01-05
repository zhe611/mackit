function result = HPCA(X, r, Method, tau, scale_est, L_init, F_init, maxiter_HPCA, maxiter_HLM, eps)

% This function is to fit the factor models via the Huber Principal Component Analysis (HPCA) method.

% X is the input matrix, of dimension \eqn{T\times N}. Each row is an observation with \eqn{N} features at time point \eqn{t}.
% r is a positive integer indicating the factor numbers.
% Method="P" indicates minimizing the Huber loss of the idiosyncratic error's norm (HPCA) 
% Method="E" indicates minimizing the elementwise Huber loss.(IHR)
% tau is the user-supplied parameter for Huber loss; default is NULL, and \tau is provided by default.
% scale est is a parameter for the elementwise Huber loss.
% L_int is the user-supplied inital value of factors; default is the PCA estimator.
% F_int is the user-supplied inital value of factors; default is the PCA estimator.
% maxiter_HLM is the maximum number of iterations in the HPCA. The default is 100.
% maxiter_HLM is the maximum number of iterations in the iterative Huber regression algorithm. The default is 100.
% eps is stopping critetion parameter in the HPCA. The default is 1e-3.

% Reference: He, Y., Li, L., Liu, D., and Zhou, W.-X. (2025). Huber principal component analysis for large dimensional factor models. 

    % Set default parameters
    if nargin<10
        eps = 10^(-8);
    end
    if nargin<9
        maxiter_HLM = 100;
    end
    if nargin<8
        maxiter_HPCA = 100;
    end
    if nargin<7
        F_init = [];
    end
    if nargin<6
        L_init = []; 
    end
    if nargin<5
        scale_est = 'MAD';
    end
    if nargin<4
        tau = [];
    end
    if nargin<3
        Method = 'E'; 
    end

    % Get dimensions
    [T, N] = size(X);

    if strcmp(Method, 'E') %HPCA
        % Method E: Robust PCA with Huber loss
        if isempty(tau)
            tau = 1.345;
        end

        if isempty(L_init) && isempty(F_init)
            % Initialize with PCA
            [U, S, V] = svd(X, 'econ');
            F_init = U(:, 1:r) * S(1:r, 1:r);
            L_init = V(:, 1:r);
        end

        if isempty(L_init)
            L_init = zeros(N, r);
            for i = 1:N
                L_init(i, :) = myrlm(F_init, X(:, i)).coefficients;
            end
            [F_init, L_init] = FLregularize(L_init, F_init, N, T, r);
        end

        if isempty(F_init)
            F_init = zeros(T, r);
            for t = 1:T
                F_init(t, :) = myrlm(L_init, X(t, :)').coefficients;
            end
            [F_init, L_init] = FLregularize(L_init, F_init, N, T, r);
        end

        CC_init = F_init * L_init';

        F_update = F_init;
        L_update = L_init;
        m = 0;

        while true
            % Update F
            for t = 1:T
                F_update(t, :) = myrlm(L_init, X(t, :)').coefficients;
            end

            % Update L
            for i = 1:N
                L_update(i, :) = myrlm(F_update, X(:, i)).coefficients;
            end

            CC_update = F_update * L_update';

            % Check maximum iterations
            if m >= maxiter_HPCA
                warning('HPCA failed to converge in %d steps', maxiter_HPCA);
                [F_update, L_update] = FLregularize(L_update, F_update, N, T, r);
                result.Fhat = F_update;
                result.Lhat = L_update;
                result.iter = m;
                return;
            end

            % Check convergence
            if SC(CC_update, CC_init, eps)
                [F_update, L_update] = FLregularize(L_update, F_update, N, T, r);
                result.Fhat = F_update;
                result.Lhat = L_update;
                result.iter = m;
                return;
            else
                L_init = L_update;
                F_init = F_update;
                CC_init = CC_update;
            end
            m = m + 1;
        end
    elseif strcmp(Method, 'P') %IHR 
        
        % Method P: Projection-based robust PCA
        if isempty(L_init) && isempty(F_init)
            [U, S, V] = svd(X, 'econ');
            L_init = V(:, 1:r);
            F_init = U(:, 1:r) * S(1:r, 1:r);
        end

        if isempty(L_init)
            [~, ~, V] = svd(X, 'econ');
            L_init = V(:, 1:r);
        end

        if isempty(F_init)
            [U, S, ~] = svd(X, 'econ');
            F_init = U(:, 1:r) * S(1:r, 1:r);
        end

        Ehat_init = X - F_init * L_init';

        m = 1;
        Loss = [];

        while true
            Ehat_l2 = vecnorm(Ehat_init, 2, 2);

            if isempty(tau)
                [~, Ehat_order] = sort(Ehat_l2);
                omega = zeros(T, 1);
                T1 = floor(T / 2);
                tau = Ehat_l2(Ehat_order(T1));

                omega(Ehat_order(1:T1)) = 0.5;
                for t = (T1 + 1):T
                    frac2 = sum(X(Ehat_order(t), :).^2) - sum((X(Ehat_order(t), :) * L_init).^2) / N;
                    omega(Ehat_order(t)) = tau / (2 * sqrt(frac2));
                end

                Sigmahat = zeros(N, N);
                for t = 1:T
                    Sigmahat = Sigmahat + omega(t) * (X(t, :)' * X(t, :));
                end
                Sigmahat = Sigmahat / T;
            end

            [V, ~] = eig(Sigmahat);
            L_update = sqrt(N) * V(:, end - r + 1:end);
            F_update = X * L_update / N;
            Ehat_update = X - F_update * L_update';

            % Check maximum iterations
            if m >= maxiter_HPCA
                warning('HPCA failed to converge in %d steps', maxiter_HPCA);
                result.Fhat = F_update;
                result.Lhat = L_update;
                result.iter = m;
                return;
            end

            % Check convergence
            if SC(Ehat_update, Ehat_init, eps)
                result.Fhat = F_update;
                result.Lhat = L_update;
                result.iter = m;
                return;
            else
                L_init = L_update;
                F_init = F_update;
                Ehat_init = Ehat_update;
                tau = [];
            end
            m = m + 1;
        end
    else
        error('Unknown method');
        
    end
end




