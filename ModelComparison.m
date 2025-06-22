% We have compared the performance of MSBL with SBL 
clear;
close all;
clc;
format long;
%---------------------------
% variables that will remain constant in simulations
% channel setup (Monte carlo loops)
max_avg = 100;
% Select what you want to vary with NMSE
SNR = 1;

if SNR == 1
    snr_db0 = 0:5:30;
    loop = length(snr_db0);
    N0 = 40;
    M0 = 20;
    L0 = 10;
    D0 = 7;
end

% initialization of error metrics
nmse_sbl = zeros(loop,1);
nmse_msbl = zeros(loop,1);
%----------------------------------------
%========================================
% select the algorithms you want to run
sbl = 1;
msbl = 1;
%========================================
%----------------------------------------

% Progress bar
wt = waitbar(0,'Initializing...');
curr_iter = 0;

for loop_iter = 1:loop
    nmse_sum_sbl = 0;
    nmse_sum_msbl = 0;
    for mc_iter = 1:max_avg
        if SNR == 1
            snr = 10.^(snr_db0(loop_iter)/10);
            N = N0;
            M = M0;
            L = L0;
            D = D0;
        end
        % W_true
        W_true = zeros(M,L);
        idx = randperm(M,D0);
        W_true(idx,:) = randn(D0,L);
            
        % Phi (Sensing matrix)
        Phi = randn(N,M);
        
        % noise
        sigma2 = 1/snr;
        En = sqrt(sigma2)*randn(N,L);
        
        % Observed data
        T = signal_generation(Phi, W_true, En);

        if sbl == 1
            w_est = zeros(M,L);
            for l = 1 : L
                mean = SBL2(Phi, T(:, l), sigma2, M);
                w_est(:,l) = mean;
            end
            nmse_sum_sbl = nmse_sum_sbl + norm(w_est - W_true)^2 / norm(W_true)^2;
        end
        if msbl == 1
            MEAN = MSBL2(Phi, T, sigma2, M, L);
            W_est = MEAN;
            nmse_sum_msbl = nmse_sum_msbl + norm(W_est - W_true,'fro')^2 / norm(W_true,'fro')^2;
        end
        %% -----------------------------------
        % Update waitbar
        curr_iter = 1 + curr_iter;
        waitbar(curr_iter/(max_avg*loop),wt,sprintf('%0.1f%% done',curr_iter/(max_avg*loop)*100))
    end
    nmse_sbl(loop_iter) = nmse_sum_sbl / max_avg;
    nmse_msbl(loop_iter) = nmse_sum_msbl / max_avg;
    
    % Plot results
    if SNR == 1
        if sbl == 1
            semilogy(snr_db0,nmse_sbl,'-o','linewidth',2,'DisplayName','SBL');
            hold on
        end
        if msbl == 1
            semilogy(snr_db0,nmse_msbl,'-x','linewidth',2,'DisplayName','MSBL');
            hold on
        end
        hold off
        xlabel('SNR (dB)');
    end
    ylabel('NMSE');
    legend show
    grid on 
end
delete(wt)


