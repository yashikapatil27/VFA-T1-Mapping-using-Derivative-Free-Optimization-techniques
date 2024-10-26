function [K, T1] = conjugate_direction(Im, alpha, TR, options, varargin)
    
    %% check input
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

    % Initialize variables
    x = zeros(1, N); % Initial point
    %p = eye(N); % Initial directions
    p = zeros(N, 1);
    
    k = 1;
    done = false;

    while ~done
        % Compute x_k+1 as the minimizer of f along the line x_k + alpha * p_k
        alpha = compute_alpha(x, p, Im, alpha, TR);
        x_new = x + alpha * p(:, end);

        % Update z_j and p_j
        z = x_new;
        %fprintf("Z: %e", z);
        p(:, 1:end-1) = p(:, 2:end);
        p(:, end) = z - x;

        % Compute alpha_n
        alpha_n = compute_alpha(z, p(:, end), Im, alpha, TR);

        % Update x_k+1
        x = z + alpha_n * p(:, end);

        % Update iteration count
        k = k + 1;

        % Convergence test
        if k >= options.MaxIter
            done = true;
        end
    end

    % Compute K and T1 using the final x
    T1 = compute_T1(x, Im, alpha, TR);
    K = compute_K(x, Im, alpha, TR);
end

% Compute alpha for minimizing f along the line
function alpha = compute_alpha(x, p, Im, alpha, TR)
    % Here you need to compute alpha using a line search method
    alpha = 0.1;
end

% Compute T1 using the final x
function T1 = compute_T1(x, Im, alpha, TR)
    % Here you need to compute T1 using the final x
    T1 = zeros(size(Im, 1), size(Im, 2), size(Im, 3));
end

% Compute K using the final x
function K = compute_K(x, Im, alpha, TR)
    % Here you need to compute K using the final x
    K = zeros(size(Im, 1), size(Im, 2), size(Im, 3));
end
