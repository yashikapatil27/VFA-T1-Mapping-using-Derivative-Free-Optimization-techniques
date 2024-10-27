function [ K, T1 ] = if_image( Im, alpha, TR, options, varargin )
    
    % Check input
    sizeim = size(Im);
    nrows = sizeim(1);
    ncols = sizeim(2);
    
    if numel(sizeim) == 3
        nslices = 1;
        nalpha = sizeim(3);
    else
        nslices = sizeim(3);
        nalpha = sizeim(4);
    end
    
    N = numel(alpha);
    
    if nalpha ~= N
        error('Dimensions do not match');
    end
    
    if isempty(nslices)
        nslices = 1;
    end
    
    if ~isfield(options, 'Direct')
        if ~isfield(options, 'MaxIter')
            options.MaxIter = 10; % Default
        elseif options.MaxIter < 1 || options.MaxIter > 200
            error('options: MaxIter should be set to a value between 1 and 200');
        end
        
        if ~isfield(options,'Tol')
            options.Tol = 1e-6; % Default
        elseif options.Tol < 0 || options.Tol > 1e-2
            error('options: Tol should be set to a value between 0 and 1e-2');
        end
        modeDirect = false;
    elseif options.Direct < 1 || options.Direct > 200
        error('options: Directiter should be set to a value between 1 and 200');
    else
        modeDirect = true;
    end
    
    if isempty(varargin)
        fprintf('User did not provide initial values neither a mask. Using default parameters...\n')
        ini = [0.5, 500];
        th = 0.05 * max(max(max(Im(:)))); % Intensity values smaller than 5% of the maximum value of the SPGR dataset are left out
        if nslices == 1
            mask = squeeze(Im(:,:,1)) > th;
        else
            mask = squeeze(Im(:,:,:,1)) > th;
        end
    elseif ~isvector(varargin{1})
         fprintf('User did not provide initial values. Using default parameters... \n');
         ini = [0.5, 500];
         mask = varargin{1};
    else
        ini = varargin{1};
        if length(varargin) == 1
            fprintf('User did not provide a mask. Using default parameters... \n');
            th = 0.05 * max(max(max(Im(:)))); % Intensity values smaller than 5% of the maximum value of the SPGR dataset are left out
            if nslices == 1
                mask = squeeze(Im(:,:,1)) > th;
            else
                mask = squeeze(Im(:,:,:,1)) > th;
            end
        else
            fprintf('User did provide initial values and a mask \n')
            mask = varargin{2};
        end
    end

    pm = find(mask);
    M = numel(pm);

    % Implicit Filtering begins here
    
    % Pre-computing
    K = squeeze(zeros(nrows, ncols, nslices));
    T1 = squeeze(zeros(nrows, ncols, nslices));
    
    alphanm = alpha * ones(1, M);
    y = reshape(Im(:), nrows * ncols * nslices, N);
    ynm = y(pm, :)';
    pnm = sind(alphanm);
    qnm = cosd(alphanm);
    
    % Initialization
    ini_K = ini(1);
    ini_T1 = ini(2);
    c1m = ini_K * (1 * exp(-TR / ini_T1)) * ones(1, M);
    c2m = exp(-TR / ini_T1) * ones(1, M);
    k = 0;
    done = false;
    
    % Set initial point
    x = ini;
    
    while ~done
        increment_k = false;
        % Compute f(x) and gradient
        f_x = objective_function(x, Im, alpha, TR);
        grad_f_x = gradient_function(x, Im, alpha, TR);
        
        if norm(grad_f_x) <= options.step_size
            increment_k = true;
        else
            m = 0;
            % Find the smallest integer m
            while m < options.max_backtracking
                if objective_function(x - options.rho^m * grad_f_x, Im, alpha, TR) <= f_x - options.c * options.rho^m * norm(grad_f_x)^2
                    break;
                end
                m = m + 1;
            end
            
            % If no such m exists, increment k
            if m == options.max_backtracking
                increment_k = true;
            else
                % Update x
                x = x - options.rho^m * grad_f_x;
            end
        end
        
        if increment_k
            done = true;
        else
            % Update k
            k = k + 1;
            xk = x;
        end
        
        if modeDirect && k == options.MaxIter
            done = true;
        end
    end

    % Update K and T1
    K(:, :, :) = xk(1);
    T1(:, :, :) = xk(2);
end

function obj_val = objective_function(parameters, Im, alpha, TR)
    % Extract parameters
    K = parameters(1);
    T1 = parameters(2);
    
    % Compute K and T1 maps using Novifast algorithm with provided parameters
    % For simplicity, let's assume a simple calculation
    % You should replace this with the actual calculation based on your problem
    % For example, you might use the Novifast algorithm to compute K and T1 maps
    
    % Example calculation: compute the sum of squared differences between
    % the actual and predicted signal intensities
    
    % Compute predicted signal intensities using the computed K and T1 maps
    predicted_signal = K * (1 - exp(-TR ./ T1));
    
    % Compute the sum of squared differences between the predicted and actual signal intensities
    % Here, Im is the actual signal intensities obtained from the imaging data
    % alpha is the flip angle, TR is the repetition time
    % You might need to reshape Im and alpha depending on your data format
    % For simplicity, let's assume Im and alpha are in the correct format
    diff_signal = predicted_signal - Im;
    sum_squared_diff = sum(diff_signal(:).^2);
    
    % Assign the sum of squared differences as the objective value
    obj_val = sum_squared_diff;
end

function grad_f = gradient_function(parameters, Im, alpha, TR)
    % Compute gradient of the objective function with respect to parameters
    
    % Compute numerical gradient using central finite differences
    h = 1e-6; % Step size for finite differences
    f_val = objective_function(parameters, Im, alpha, TR);
    grad_f = zeros(size(parameters));
    
    for i = 1:numel(parameters)
        perturbed_params_plus = parameters;
        perturbed_params_plus(i) = perturbed_params_plus(i) + h;
        perturbed_params_minus = parameters;
        perturbed_params_minus(i) = perturbed_params_minus(i) - h;
        
        grad_f(i) = (objective_function(perturbed_params_plus, Im, alpha, TR) - objective_function(perturbed_params_minus, Im, alpha, TR)) / (2 * h);
    end
end
