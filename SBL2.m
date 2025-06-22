function mean = SBL2(Phi, t, sigma2, M)

    % Initialization
    alpha = ones(M, 1) * 1e-3; % large initial precision
    mean = zeros(M, 1);
    max_itr = 150;

    for itr = 1:max_itr
        % Prior covariance
        A = diag(alpha);

        % E-step: compute posterior mean and covariance
        Sn = inv(sigma2^-1 * (Phi' * Phi) + A);
        mean_new = sigma2^-1 * Sn * Phi' * t; 

        % M-step: update alpha
        alpha_new = 1 ./ (mean_new.^2 + diag(Sn));

        % Convergence check
        if norm(mean_new - mean) < 1e-4
            break;
        end

        mean = mean_new;
        alpha = alpha_new;
    end
end
