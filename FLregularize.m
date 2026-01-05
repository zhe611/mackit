function [Fhat, Lhat] = FLregularize(L_init, F_init, N, T, r)
    % Regularize F and L
    [U, S, V] = svd(F_init * L_init', 'econ');
    Fhat = U(:, 1:r) * S(1:r, 1:r);
    Lhat = V(:, 1:r);
end

