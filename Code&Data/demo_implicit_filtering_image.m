%clear all; clc; close all;

addpath('./data')
load('volume3DFSE.mat') % 3D FSE SPGR dataset (1.5 T) used in in vivo experiment of '....'

% Implicit Filtering's parameter definition
TR = 9; % Repetition time [ms]
ini = [0.2, 500]; % K [] and T1 [ms] initial constant maps 

% Options structure
%options_imfi = struct('Direct', 1);
options_imfi.k = @(k) 1 / (2^k); % Example of defining the sequence k
options_imfi.c = 0.5; % Armijo parameter
options_imfi.rho = 0.5; % Parameter in (0, 1)
options_imfi.amax = 100; % Maximum backtracking parameter
options_imfi.Tol = 1e-6; % Convergence tolerance
options_imfi.MaxIter = 30; % Maximum number of iterations

% Call Implicit Filtering
time = tic;
[K, T1, min_snm_ynm_values_imfi] = implicit_filtering_optimization(im, alpha, TR, options_imfi);
time = toc(time);

% Display final computational time
disp(['Dataset computation time = ', num2str(time), ' s']);

% Plot the minimum values of snm-ynm against the iteration number for NOVIFAST
figure;
plot(1:length(min_snm_ynm_values_imfi), min_snm_ynm_values_imfi, 'b.-');
xlabel('Iteration');
ylabel('Minimum Value of snm-ynm');
title('Minimum Value of snm-ynm Convergence (Implicit Filtering)');
grid on;