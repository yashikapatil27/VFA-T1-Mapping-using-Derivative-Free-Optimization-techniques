function [K, T1] = pattern_search_optimization(Im, alpha, TR, options, varargin)
    
    %% Check input
    sizeim = size(Im);
    disp(['Size of Im: ', num2str(sizeim)]);
    nrows = sizeim(1);
    ncols = sizeim(2);
    
    % Check if the input image is 3D or 4D
    if numel(sizeim) == 3
        nslices = 1;
        nalpha = sizeim(3);
    else
        nslices = sizeim(3);
        nalpha = sizeim(4);
    end
    
    % Check FA
    N = numel(alpha);
    
    if nalpha ~= N
        error('Dimensions do not match');
    end
    
    if isempty(nslices)
        nslices = 1;
    end
    
    % Checks the options structure passed to the function and sets the convergence parameters accordingly.
    if ~isfield(options, 'MaxIter')
        options.MaxIter = 10; % Default
    elseif options.MaxIter < 1 || options.MaxIter > 200
        error('options: MaxIter should be set to a value between 1 and 200');
    end
    
    if ~isfield(options, 'Tol')
        options.Tol = 1e-6; % Default
    elseif options.Tol < 0 || options.Tol > 1e-2
        error('options: Tol should be set to a value between 0 and 1e-2');
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
    
    %% Pattern search optimization begins here
    
    % pre-computing
    K = zeros(nrows, ncols, nslices);
    T1 = zeros(nrows, ncols, nslices);
    
    alphanm = alpha * ones(1, M);
    y = reshape(Im(:), nrows*ncols*nslices, N);
    ynm = y(pm,:)';
    pnm = sind(alphanm);
    qnm = cosd(alphanm);

    % initialization
    Kini = ini(1);
    T1ini = ini(2);
    c1m = Kini * (1 * exp(-TR/T1ini)) * ones(1, M);
    c2m = exp(-TR/T1ini) * ones(1, M);
    k = 0;
    gamma = 0.1; % Initial step length
    gamma_tol = options.Tol; % Convergence tolerance
    theta_max = 0.5; % Contraction parameter
    
    % Initial direction set D0
    D = eye(length(ini)); % Identity matrix for initial direction set

    % iterative process
    while gamma > gamma_tol && k < options.MaxIter
        if k == 0
            % Evaluate the objective function at the initial point
            f_xk = objective_function(ini); % Implement this function
            
            % Initialize the search direction set
            Dk = D;
        else
            % Evaluate the objective function at xk + gamma * pk
            f_xk_gamma_pk = objective_function(xk + gamma * pk); % Implement this function
            
            % Check if f(xk + gamma * pk) < f(xk) - rho(gamma)
            if f_xk_gamma_pk < f_xk - rho(gamma)
                xk = xk + gamma * pk;
                gamma = phi * gamma;
            else
                gamma = theta * gamma;
            end
        end
        
        k = k + 1;
    end

    % Function to evaluate the objective function
    function f_val = objective_function(x)
        % Extract parameters from x
        Kini = x(1);
        T1ini = x(2);
        
        % Call novifast_image with the current parameters
        [K_est, T1_est] = novifast_image(Im, alpha, TR, options, [Kini, T1ini]);
        
        % Define ground truth (you need to provide this)
        % For demonstration, assuming K_gt and T1_gt are the ground truth maps
        % Compute error between estimated maps and ground truth
        error_K = sum(abs(K_est(:) - K_gt(:)));
        error_T1 = sum(abs(T1_est(:) - T1_gt(:)));
        
        % Total error as the sum of errors in K and T1
        f_val = error_K + error_T1;
    end
    
    % Function to calculate the sufficient decrease function rho
    function rho_val = rho(gamma)
        % Define a simple sufficient decrease function
        rho_val = gamma^2; % Adjust as needed
    end

end
