function [K, T1] = implicit_filtering_image(Im, alpha, TR, options, varargin)
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
    
    disp(['Size of alpha: ', num2str(size(alpha))]); % Add this line
    disp(['Size of TR: ', num2str(size(TR))]); % Add this line
    
    % Set default options
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
        ini=[0.5,500];
        th=0.05*max(max(max(Im(:)))); %Intensity values smaller than 5% of the maximum value of the SPGR dataset are left out
        if nslices==1
            mask = squeeze(Im(:,:,1))>th;
        else
            mask = squeeze(Im(:,:,:,1))>th;
        end
    elseif ~isvector(varargin{1})
         fprintf('User did not provide initial values. Using default parameters... \n');
         ini=[0.5,500];
         mask=varargin{1};
    else
        ini=varargin{1};
        if length(varargin)==1
            fprintf('User did not provide a mask. Using default parameters... \n');
            th=0.05*max(max(max(Im(:)))); %Intensity values smaller than 5% of the maximum value of the SPGR dataset are left out
            if nslices==1
                mask = squeeze(Im(:,:,1))>th;
            else
                mask = squeeze(Im(:,:,:,1))>th;
            end
        else
            fprintf('User did provide initial values and a mask \n')
            mask=varargin{2};
        end
    end

    pm=find(mask);
    M=numel(pm);
    
    % Define the objective function
    objective_function = @(parameters) sum((compute_predicted_signal(parameters, alpha, TR, sizeim) - Im).^2);

    % Initialize K and T1
    K = zeros(nrows, ncols, nslices);
    T1 = zeros(nrows, ncols, nslices);
    
    % Initialize simplex with initial parameter values
    simplex = repmat(ini', 1, N + 1);
    
    % Initialize objective values
    obj_values = zeros(1, N + 1);
    for i = 1:N+1
        obj_values(i) = objective_function(simplex(:, i));
    end
    
    k = 1;
    
    % Algorithm parameters
    c = options.C; % Armijo parameter c
    rho = options.Rho; % Armijo parameter ρ
    amax = options.Amax; % Maximum backtracking parameter amax
    
    % Main loop
    while k <= options.MaxIter
        x = simplex(:, 1); % Initial point
        increment_k = false;
    
        % Backtracking line search
        while true
            % Compute gradient
            grad_f = compute_gradient(x, Im, alpha, TR);
    
            % Check termination condition
            if norm(grad_f) < options.Tol
                increment_k = true;
                break;
            end
    
            m = 0;
            while true
                % Backtracking
                x_new = x - rho^m * grad_f;
                if objective_function(x_new) <= objective_function(x) - c * rho^m * norm(grad_f)^2
                    x = x_new;
                    break;
                elseif m >= amax
                    increment_k = true;
                    break;
                else
                    m = m + 1;
                end
            end
    
            if increment_k
                break;
            end
        end
    
        % Update simplex
        simplex(:, 1) = x;
    
        % Check termination
        if increment_k
            break;
        else
            k = k + 1;
        end
    
        % Print sizes for debugging
        disp('Size of simplex(:, i):');
        disp(size(simplex(:, i))); % Add this line
        disp('Size of Im:');
        disp(size(Im)); % Add this line
    
        % Compute objective values
        for i = 1:N+1
            obj_values(i) = objective_function(simplex(:, i));
        end
    end


    % Set final K and T1
    K(:, :, :) = reshape(simplex(1, :), [nrows, ncols, nslices]);
    T1(:, :, :) = reshape(simplex(2, :), [nrows, ncols, nslices]);
end

function predicted_signal = compute_predicted_signal(parameters, alpha, TR, sizeIm)
    K = parameters(1);
    T1 = parameters(2);
    
    % Compute the predicted signal
    predicted_signal_slice = K * (1 - exp(-TR / T1)) .* alpha;
    
    % Reshape predicted_signal_slice to match the size of Im
    predicted_signal = repmat(predicted_signal_slice, [1, sizeIm(2), sizeIm(3)]);
end


function grad_f = compute_gradient(parameters, Im, alpha, TR)
    h = 1e-6; % Step size for numerical differentiation
    grad_f = zeros(2, 1);
    for i = 1:2
        perturbed_params = parameters;
        perturbed_params(i) = perturbed_params(i) + h;
        grad_f(i) = (objective_function(perturbed_params, Im, alpha, TR) - objective_function(parameters, Im, alpha, TR)) / h;
    end
end

function obj_val = objective_function(parameters, Im, alpha, TR)
    predicted_signal = compute_predicted_signal(parameters, alpha, TR, size(Im));
    obj_val = sum((predicted_signal - Im).^2);
end
