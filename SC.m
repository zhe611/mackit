function res = SC(mu_update, mu_init, tol)
    % Stopping criterion
    mu_inl2 = norm(mu_init, 'fro');
    if mu_inl2 == 0
        res = norm(mu_update - mu_init, 'fro');
    else
        res = norm(mu_update - mu_init, 'fro') / mu_inl2;
    end
    res = (res <= tol);
end