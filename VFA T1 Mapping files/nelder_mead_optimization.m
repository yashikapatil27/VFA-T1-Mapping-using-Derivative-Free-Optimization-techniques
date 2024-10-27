function [ K, T1, min_snm_ynm_values ] = nelder_mead_optimization( Im, alpha, TR, options, varargin )
    
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
    if ~isfield(options, 'Direct') % Checks if the field `'Direct'` is not present in the options structure.
        if ~isfield(options, 'MaxIter')
            options.MaxIter = 10; % Default
        elseif options.MaxIter < 1 || options.MaxIter > 200
            error('options: Maxiter should be set to a value between 1 and 200');
        end
        
        if ~isfield(options, 'Tol')
            options.Tol = 1e-6; % Default
        elseif options.Tol < 0 || options.Tol > 1e-2
            error('options: Tol should be set to a value between 0 and 1e-2');
        end
        modeDirect = false; % Indicates that the convergence mode is not set to direct.
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
    
    %% Nelder Mead begins here

    % Pre-computing
    K = zeros(nrows, ncols, nslices); % Initialize K with the specified size
    T1 = zeros(nrows, ncols, nslices); % Initialize T1 with the specified size

    alphanm = alpha * ones(1, M);
    y = reshape(Im(:), nrows * ncols * nslices, N);
    ynm = y(pm, :)';
    pnm = sind(alphanm);
    qnm = cosd(alphanm);

    % Initialization
    Kini = ini(1);
    T1ini = ini(2);
    c1m = Kini * (1 * exp(-TR/T1ini)) * ones(1, M);
    c2m = exp(-TR/T1ini) * ones(1, M);

    % Initialize min_snm_ynm_values array to store minimum values of snm-ynm
    min_snm_ynm_values = [];

    % Nelder mead step instead of iterative process
    [K, T1, min_snm_ynm_values] = nelder_mead_step(K, T1, c1m, c2m, ynm, pnm, qnm, alpha, N, TR, options, modeDirect, mask, pm);

end

function [ K, T1, min_snm_ynm_values ] = nelder_mead_step(K, T1, c1m, c2m, ynm, pnm, qnm, alpha, N, TR, options, modeDirect, mask, pm)

    pm = find(mask);
    M = numel(pm);
    % Initialize the simplex with the provided c1m and c2m
    simplex = [c1m(:), c2m(:)];
    % Number of variables (dimension)
    n = size(simplex, 2);
    % Number of vertices (including the best and worst points)
    m = n + 1;
    
    % Initialize objective function values for the simplex
    f_values = zeros(m, 1);
    for i = 1:m
        % Evaluate the objective function for each vertex of the simplex
        f_values(i) = objective_function(simplex(i, :), c1m, c2m, ynm, pnm, qnm, alpha, TR, M);
    end
    
    % Main loop of the Nelder-Mead algorithm
    iter = 0;
    done = false;
    min_snm_ynm_values = [];

    while iter < options.MaxIter && ~done % Termination based on maximum number of iterations
        % Update c2m_old at the beginning of each iteration
        c2m_old = c2m;
        
        % Sort the vertices of the simplex based on the objective function values
        [f_values, order] = sort(f_values);
        simplex = simplex(order, :);
        
        % Compute the centroid of all vertices except the worst one
        centroid = mean(simplex(1:end-1, :));
        
        % Reflection: Compute the reflected point and evaluate the objective function
        x_reflected = centroid + (centroid - simplex(end, :));
        f_reflected = objective_function(x_reflected, c1m, c2m, ynm, pnm, qnm, alpha, TR, M);
        % Reflection
        disp('Reflection');
        disp(['f_reflected: ', num2str(f_reflected)]);

        if f_reflected < f_values(end-1) && f_reflected >= f_values(1)
            % Replace the worst point with the reflected point
            simplex(end, :) = x_reflected;
            f_values(end) = f_reflected;
            % Update c2m with the second element of the reflected point
            c2m = x_reflected(:, 2);
        elseif f_reflected < f_values(1)
            % Expansion
            disp('Expansion');
            x_expanded = centroid + 2 * (x_reflected - centroid);
            f_expanded = objective_function(x_expanded, c1m, c2m, ynm, pnm, qnm, alpha, TR, M);
            if f_expanded < f_reflected
                % Replace the worst point with the expanded point
                simplex(end, :) = x_expanded;
                f_values(end) = f_expanded;
                % Update c2m with the second element of the expanded point
                c2m = x_expanded(:, 2);
            else
                % Replace the worst point with the reflected point
                simplex(end, :) = x_reflected;
                f_values(end) = f_reflected;
                % Update c2m with the second element of the reflected point
                c2m = x_reflected(:, 2);
            end
        else
            % Contraction
            disp('Contraction');
            x_contracted = centroid + 0.5 * (simplex(end, :) - centroid);
            f_contracted = objective_function(x_contracted, c1m, c2m, ynm, pnm, qnm, alpha, TR, M);
        
            if f_contracted < f_values(end)
                % Replace the worst point with the contracted point
                simplex(end, :) = x_contracted;
                f_values(end) = f_contracted;
                % Update c2m with the second element of the contracted point
                c2m = x_contracted(:, 2);
            else
                % Shrink the simplex towards the best point
                disp('Shrink');
                for i = 2:m
                    % Shrink the simplex
                    simplex(i, :) = 0.5 * (simplex(i, :) + simplex(1, :));
                end
            end
            % Update min_snm_ynm_values with the minimum value of snm-ynm
            min_snm_ynm = min(abs(c1m .* pnm ./ (1 - c2m_old .* qnm) - ynm), [], 'all');
            min_snm_ynm_values = [min_snm_ynm_values, min_snm_ynm];
        end


        K(pm) = simplex(1, 1);
        T1(pm) = -TR / log(simplex(1, 2));

        

        % Check convergence based on stopping criteria
        if modeDirect
            % Mode with no-convergence criterion
            if iter == options.Direct
                done = true;
            end
        else
            % Mode with convergence criterion
            rel_err = (norm(c2m(isfinite(c2m)) - c2m_old(isfinite(c2m_old)), 1)) / norm(c2m(isfinite(c2m)), 1); %Relative l1 norm for c2 (Convergence is controlled by c2 only)  
            fprintf('rel_err:%f\n', rel_err);
    
            if rel_err < options.Tol || iter >= options.MaxIter
                fprintf('Convergence criteria triggered. ');
                done = true; 
            end
        end
        

        iter = iter + 1; % Increment iteration count
        disp(['Minimum value of snm-ynm = ', num2str(min_snm_ynm_values(end))]);
    end

    % Return the best solution found
    K = K;
    T1 = T1;
end


function f_value = objective_function(x, c1m, c2m, ynm, pnm, qnm, alpha, TR, M)

    % Evaluate the objective function (sum of squared differences)

    % Extract parameters from x
    N = numel(alpha); 

    % Compute snm based on the current parameters
    denm = 1 - c2m .* qnm;
    snm = c1m .* pnm ./ denm; 

    % Compute the sum of squared differences between ynm and snm
    squared_diff = sum((ynm - snm).^2, 1);

    % Compute the objective function value
    f_value = sum(squared_diff);
end
