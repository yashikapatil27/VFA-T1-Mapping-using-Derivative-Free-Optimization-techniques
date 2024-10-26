clear all; clc; close all; 
addpath('./data') 
load('volume3DFSE.mat') % 3D FSE SPGR dataset (1.5 T) used in in vivo experiment 

% NOVIFAST's parameter definition 
TR = 9; % Repetition time [ms] 
ini = [0.2, 500]; % K [] and T1 [ms] initial constant maps for NOVIFAST  
options_novifast = struct('Direct', 2); % If field 'Direct' is given, it means NOVIFAST is run in a blind mode, i.e., no convergence criterion. Just 5 iterations are done. 

% Call NOVIFAST and measure time 
time_novifast = tic;
[K_novifast, T1_novifast] = novifast_image(im, alpha, TR, options_novifast); 
time_novifast = toc(time_novifast); 

% Nelder-Mead's parameter definition 
options_nelder = struct('Direct', 2); % Run in a blind mode, i.e., no convergence criterion. 
options_nelder.MaxIter = 100; 

% Call Nelder-Mead and measure time 
time_nelder = tic; 
[K_nelder, T1_nelder] = nelder_mead_optimization(im, alpha, TR, options_nelder); 
time_nelder = toc(time_nelder); 

% IF's parameter definition 
%options_if = struct('Direct', 2); % Run in a blind mode, i.e., no convergence criterion. 
options_if=struct('Direct',2); %If field 'Direct' is given, it means NOVIFAST is run in a blind mode, i.e., no convergence criterion. Just 5 iter are done. 
options_if.step_size = 0.01; % Set the step size value according to your needs 
options_if.max_backtracking = 20; % Set the maximum number of backtracking steps 
options_if.rho = 0.5;
options_if.c = 0.5; % Set the Armijo condition parameter (c) value 
options_if.MaxIter = 100; 

options_imfi = struct('Direct', 1);
options_imfi.k = @(k) 1 / (2^k); % Example of defining the sequence k
options_imfi.c = 0.5; % Armijo parameter
options_imfi.rho = 0.5; % Parameter in (0, 1)
options_imfi.amax = 100; % Maximum backtracking parameter
options_imfi.Tol = 1e-6; % Convergence tolerance
options_imfi.MaxIter = 100; % Maximum number of iterations

% Call IF and measure time 
time_if = tic; 
[K_if, T1_if] = implicit_filtering_optimization(im, alpha, TR, options_imfi); 
time_if = toc(time_if); 


% Display computation times 
disp(['Computation time for NOVIFAST: ', num2str(time_novifast), ' seconds']); 
disp(['Computation time for Nelder-Mead: ', num2str(time_nelder), ' seconds']); 
disp(['Computation time for Implicit-filtering: ', num2str(time_if), ' seconds']); 
%disp(['Computation time for Implicit-filtering: ', num2str(time_if), ' seconds']); 

% Visualization 
warning off 
figure(1)
if numel(size(im)) == 4 % If im is a volume, visualize middle slice 
nslice = 18; 
imshow(squeeze(T1_novifast(:,:,nslice)), [500, 5000]) % A single slice is visualized 
title('Estimated in vivo T_1 map [ms]') 
colorbar 
elseif numel(size(im)) == 3 
nslice = 1; 
imshow(T1_novifast, [500, 5000]) 
title('Estimated in vivo T_1 map [ms]') 
colorbar 
end 
strg = ['Slice nz = ', num2str(nslice)]; 
text(104, 11, strg, 'fontsize', 14, 'Color', 'red') 
strg2 = ['NOVIFAST computation time = ', num2str(time_novifast), ' s']; 
text(68, 244, strg2, 'fontsize', 14, 'Color', 'green') 
% Add computation time for Nelder-Mead to the figure 
strg3 = ['Nelder-Mead computation time = ', num2str(time_nelder), ' s']; 
text(68, 224, strg3, 'fontsize', 14, 'Color', 'yellow') 
% Add computation time for IF to the figure 
strg4 = ['IF computation time = ', num2str(time_if), ' s']; 
text(68, 204, strg4, 'fontsize', 14, 'Color', 'cyan')  