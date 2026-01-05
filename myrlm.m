function result = myrlm(X, y, tau, loss, scale_est, beta_init, w_init, maxiter, tol)
    % Set default parameters
    if nargin<9
        tol = 1e-8;
    end
    if nargin<8
        maxiter = 100;
    end
    if nargin<7
        w_init = [];
    end
    if nargin<6
        beta_init = [];
    end
    if nargin<5
        scale_est = 'MAD';
    end
    if nargin<4
        loss = 'huber';
    end
    if nargin<3
        tau = 1.345;
    end

    % Initialize weights
    if isempty(w_init)
        w = ones(length(y), 1);
    else
        w = w_init;
    end

    % Initialize coefficients
    if isempty(beta_init)
        beta_init = lscov(X, y, w);  % Weighted least squares
    end

    % Initialize residuals
    epsilon_init = y - X * beta_init;

    % Iteration
    m = 0;
    Loss = [];
    while true
        % Compute robust standard deviation
        sigmahat = Robust_sigma(epsilon_init, scale_est);

        % Compute weights
        w = Psi_Type(epsilon_init / sigmahat, tau, loss);

        % Update coefficients
        beta_update = lscov(X, y, w);  % Weighted least squares

        % Compute loss
        Loss(m + 1) = sum((beta_update - beta_init).^2) / max(1e-20, sum(beta_init.^2));

        % Check convergence
        if Loss(m + 1) < tol
            result.coefficients = beta_update;
            result.m = m;
            result.Loss = Loss;
            return;
        else
            beta_init = beta_update;
            epsilon_init = y - X * beta_update;
            m = m + 1;
        end

        % Check maximum iterations
        if m > maxiter
            warning('rlm failed to converge in %d steps', maxiter);
            result.coefficients = beta_update;
            result.m = m;
            result.Loss = Loss;
            return;
        end
    end
end

