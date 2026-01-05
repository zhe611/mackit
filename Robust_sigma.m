function sigmahat = Robust_sigma(x, type)
    % Compute robust standard deviation
    if strcmp(type, 'MAD')
        sigmahat = median(abs(x)) / 0.6745;
    elseif strcmp(type, 'const')
        sigmahat = 1;
    else
        error('Unknown scale_est type');
    end
end