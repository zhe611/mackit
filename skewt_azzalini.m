function Y = skewt_azzalini(xi, omega, alpha, nu, n, m)
% Generate random variables from the Azzalini skew-t distribution
%
% Inputs:
%   xi    : Location parameter
%   omega : Scale parameter (omega > 0)
%   alpha : Skewness parameter
%   nu    : Degrees of freedom (nu > 0)
%   n     : Number of rows (sample size)
%   m     : Number of columns (default = 1)
%
% Outputs:
%   Y     : n-by-m matrix following ST(xi, omega^2, alpha, nu)

    if nargin < 6
        m = 1;
    end

    % 1. Convert the skewness parameter alpha to the latent truncation
    %    coefficient delta used in the hidden truncation representation
    delta = alpha / sqrt(1 + alpha^2);
    
    % 2. Generate the core standard skew-normal component: SN(0, 1, alpha)
    U0 = randn(n, m);
    U1 = randn(n, m);
    Z = delta * abs(U0) + sqrt(1 - delta^2) * U1;  % Z ~ SN(0, 1, alpha)
    
    % 3. Generate an independent chi-square random variable
    V = chi2rnd(nu, n, m);
    
    % 4. Construct the standard Azzalini skew-t random variable: ST(0, 1, alpha, nu)
    T = Z ./ sqrt(V / nu);
    
    % 5. Apply the location (xi) and scale (omega) transformation
    Y = xi + omega * T;

end
