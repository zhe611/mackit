function w = Psi_Type(x, tau, loss)
    % Compute weights
    if strcmp(loss, 'huber')
        w = min(1, tau ./ abs(x));
    else
        error('Unknown loss type');
    end
end