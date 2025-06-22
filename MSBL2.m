% EM algorithm Function

function MEAN = MSBL2(Phi, T, sigma2, M, L)
    % SIGMA = inv(sigma2.^-1*(Phi' * Phi) + inv(GAMMA));
    % Using woodbury matrix identity, for writing SIGMA as below

    % Initilaization
    gamma = ones(M,1)*1e-3;
    MEAN = ones(M,L);
    
    max_itr = 150;
    for itr = 1:max_itr
         GAMMA = diag(gamma);
         SIGMA = GAMMA - GAMMA * Phi' * inv(sigma2 * eye(size(Phi,1)) + Phi * GAMMA * Phi') * Phi * GAMMA;
         MEAN_new = (SIGMA * Phi' * T)/sigma2;
         gamma_new = zeros(M,1);
         for i = 1:M
             row_norm = 0;
             for j = 1:L
                 row_norm = row_norm + MEAN(i,j).^2;
             end
             gamma_new(i) = row_norm./L + SIGMA(i,i); % gamma -> M x 1
         end
         if(norm(MEAN_new - MEAN)<1e-4) 
              break
         end
         MEAN = MEAN_new;
         gamma = gamma_new;
     end
end
        