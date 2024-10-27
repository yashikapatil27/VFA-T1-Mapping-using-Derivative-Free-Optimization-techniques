%clear all; clc; close all;

addpath('./data')
load('volume3DFSE.mat') % 3D FSE SPGR dataset (1.5 T) used in in vivo experiment of '....'

% Nelder Mead's parameter definition
TR = 9; % Repetition time [ms]
ini = [0.2, 500]; % K [] and T1 [ms] initial constant maps
options = struct('Direct', 2); % If field 'Direct' is given, it means the code is run in a blind mode, i.e., no convergence criterion. 
options.Tol = 0;
options.MaxIter = 5;

% Call Nelder Mead
time = tic;
[K, T1, min_snm_ynm_values_nm] = nelder_mead_optimization(im, alpha, TR, options);
time = toc(time);

% Display final computational time
disp(['Dataset computation time = ', num2str(time), ' s']);

% Plot the minimum values of snm-ynm against the iteration number for NOVIFAST
figure;
plot(1:length(min_snm_ynm_values_nm), min_snm_ynm_values_nm, 'b.-');
xlabel('Iteration');
ylabel('Minimum Value of snm-ynm');
title('Minimum Value of snm-ynm Convergence (Nelder Mead)');
grid on;
