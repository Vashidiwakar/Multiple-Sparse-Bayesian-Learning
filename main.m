%{
Dimensions for matrices
T -> N x L
Phi -> N x M
W -> M x L
En -> N x L (noise)
N << M
%}
clc;
clear;
close all;
M = 40;
N = 20;
L = 10;
D0 = 7;

% W_true
W_true = zeros(M,L);
idx = randperm(M,D0);
W_true(idx,:) = randn(D0,L);
    
% Phi
Phi = randn(N,M);

 % noise
 SNR_dB = 0:5:30;
 SNR_lin = 10.^(SNR_dB/10);
 % NMSE initialization
 NMSE = zeros(size(SNR_lin));
    
  for i = 1:length(SNR_lin)
     nmse_sum = 0;
     for j = 1:30
         sigma2 = 1/SNR_lin(i);
         En = sqrt(sigma2)*randn(N,L);
         T = signal_generation(Phi, W_true, En);
            
         MEAN = MSBL2(Phi, T, sigma2, M, L);
         W_est = MEAN;
         nmse_temp = norm(W_est - W_true,'fro')^2 / norm(W_true,'fro')^2;
         nmse_sum = nmse_sum + nmse_temp;
     end
     NMSE(i) = nmse_sum/30;
 end
% plot of SNR_dB vs NMSE
figure;
semilogy(SNR_dB, NMSE, 'r-o',LineWidth=1.5);
xlabel('SNR(dB)');
ylabel('NMSE');
title('NMSE vs SNR for Simultaneous SBL (EM)');
grid on;
