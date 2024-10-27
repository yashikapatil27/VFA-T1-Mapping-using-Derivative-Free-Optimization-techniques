function [K, T1, min_snm_ynm_values] = implicit_filtering_optimization(Im, alpha, TR, options_imfi, varargin)     
    %% check input
    sizeim=size(Im);
    disp(['Size of Im: ', num2str(sizeim)]);
    nrows=sizeim(1);
    ncols=sizeim(2);

    % Check if the input image is 3D or 4D
    if numel(sizeim)==3
        nslices=1;
        nalpha=sizeim(3);
    else
        nslices=sizeim(3);
        nalpha=sizeim(4);
    end

    % Check FA
    N=numel(alpha);         

    if nalpha~=N
        error('Dimensions do not match');
    end

    if isempty(nslices)
        nslices=1;
    end

    disp(options_imfi);  % Print the contents of options_imfi to inspect its structure

    % Checks the options structure passed to the function and sets the convergence parameters accordingly.
    if ~isfield(options_imfi, 'Direct')          % Checks if the field `'Direct'` is not present in the options structure.
        if ~isfield(options_imfi, 'MaxIter')
            options_imfi.MaxIter = 10; %default
            modeDirect=false
        elseif options_imfi.MaxIter<1 || options_imfi.MaxIter>200
            error('options: Maxiter should be set to a value between 1 and 200');
        end
        
        if ~isfield(options_imfi,'Tol')
            options_imfi.Tol = 1e-6; %default
        elseif options_imfi.Tol<0 || options_imfi.Tol>1e-2
            error('options: Tol should be set to a value between 0 and 1e-2');
        end
        modeDirect=false;       % modeDirect to false. This indicates that the convergence mode is not set to direct.
    elseif options_imfi.Direct<1 || options_imfi.Direct>200
        error('options: Directiter should be set to a value between 1 and 200');
    else
        modeDirect=true;
    end

    if isempty(varargin)
        fprintf('User did not provide initial values neither a mask. Using default parameters...\n')
        ini=[0.2,500];
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

    % Pre-computing 
    K = zeros(nrows, ncols, nslices); 
    T1 = zeros(nrows, ncols, nslices);    
    alphanm = alpha * ones(1, M);     
    y = reshape(Im(:), nrows * ncols * nslices, N); 
    ynm = y(pm, :)'; 
    pnm = sind(alphanm); 
    qnm = cosd(alphanm); 

    % Initialization 
    
    x = ini; 
    k = 1; 
    %increment_k = false; 
    % initialization
    Kini=ini(1);
    T1ini=ini(2);
    c1m_old=Kini*(1*exp(-TR/T1ini))*ones(1,M);
    c2m_old=exp(-TR/T1ini)*ones(1,M);
    done = false; 
    iter = 0; 
    
    % Initialize min_snm_ynm_values array to store minimum values of snm-ynm
    min_snm_ynm_values = [];

    % Implicit Filtering loop 
    
    while iter < options_imfi.MaxIter && ~done 
    
    % Compute f(x) and ∇f(x)
    %f_x = objective_function(x, alpha, M, ynm);
    [f_x, min_snm_ynm_values] = objective_function(x, alpha, M, TR, ynm, min_snm_ynm_values, iter);

    %grad_f_x = finite_difference_gradient(x, alpha, M, ynm); 
    grad_f_x = finite_difference_gradient(x, alpha, M, TR, ynm, min_snm_ynm_values, iter);

    % Check Armijo condition 
    if norm(grad_f_x) <= options_imfi.k(k) 
        done = true; 
    else 
    
        % Find smallest integer m between 0 and amax 
        m = 0; 
        x_temp = x; 
        f_x_temp = f_x; 
        
        while m < options_imfi.amax 
            % Compute new x and f(x) using backtracking line search 

            x_temp = x - options_imfi.rho^m * grad_f_x; 
            %[f_x_temp, min_snm_ynm_values, c1, c2] = objective_function(x_temp, alpha, M, ynm, min_snm_ynm_values, iter);
            [f_x_temp, min_snm_ynm_values, c2m, c2m_old] = objective_function(x_temp, alpha, M, TR, ynm, min_snm_ynm_values, iter);
            if f_x_temp <= f_x - options_imfi.c * options_imfi.rho^m * (grad_f_x' * grad_f_x) 
                break; 
            end 
            m = m + 1; 
        
        end 
        if m == options_imfi.amax 
            done = true; 
        else    
            x = x_temp; 
        end    
    end 
    
    % Update iteration count 
    iter = iter + 1; 
    
    % Print minimum value of snm-ynm for each iteration
    disp(['Minimum value of snm-ynm = ', num2str(min_snm_ynm_values(end))]);
    %stopping 
    %min_snm_ynm = min(abs(c1m .* pnm ./ (1 - c2m .* qnm) - ynm), [], 'all'); 
    %min_snm_ynm_values = [min_snm_ynm_values, min_snm_ynm]; 
    
    if modeDirect %mode with no-convergence criterion 
        if isfield(options_imfi, 'Direct') && k == options_imfi.Direct 
            done = true; 
        end 
    else 
        rel_err = (norm(c2m(isfinite(c2m)) - c2m_old(isfinite(c2m_old)), 1)) / norm(c2m(isfinite(c2m)), 1); %Relative l1 norm for c2 (Convergence is controlled by c2 only)  
        fprintf('rel_err: %f, k: %d\n', rel_err, k);

        if rel_err < options_imfi.Tol || k >= options_imfi.MaxIter
            fprintf('Convergence criteria triggered. rel_err: %f, k: %d\n', rel_err, k);
            done = true; 
        end

    end  

    
    % Update K and T1 
    K(pm) = x(1); 
    T1(pm) = x(2); 

    
    % Update min_snm_ynm_values with the minimum value of snm-ynm
    %min_snm_ynm = min(abs(c1m .* pnm ./ (1 - c2m .* qnm) - ynm), [], 'all');
    %min_snm_ynm_values = [min_snm_ynm_values, min_snm_ynm];
    
    
    end 

    % Check the size of T1 
    disp(['Size of T1: ', num2str(size(T1))]); 
    % Print statements moved outside the loop
    %disp(['Objective Function - Iteration ', num2str(iter), ': Minimum value of snm-ynm = ', num2str(min_snm_ynm)]);
    %disp(['Finite Difference Gradient - Iteration ', num2str(iter), ': Calculating gradient for parameter ', num2str(1i)]);
end 

function [f_value, min_snm_ynm_values, c2m, c2m_old] = objective_function(x, alpha, M, TR, ynm, min_snm_ynm_values, iter) 
    % Evaluate the objective function (sum of squared differences)
    % Extract parameters from x 
    
    %c1 = x(1); 
    %c2 = x(2); 
    N=numel(alpha);         

    Kini = x(1);
    T1ini = x(2);

    c1m=Kini*(1*exp(-TR/T1ini))*ones(1,M);
    c2m=exp(-TR/T1ini)*ones(1,M);
    c2m_old=c2m;

    c1m=repmat(c1m,[N,1]);
    c2m=repmat(c2m,[N,1]);
    
    alphanm=alpha*ones(1,M);
    pnm=sind(alphanm); 
    qnm=cosd(alphanm); 
    
    % Compute snm based on the current parameters 
    denm = 1 - c2m .* qnm;
    snm = c1m.* pnm ./ denm; 
    
    % Compute the sum of squared differences between ynm and snm 
    squared_diff = sum((ynm - snm).^2, 1); 

    % Compute minimum value of snm-ynm
    min_snm_ynm = min(abs(c1m .* pnm ./ (1 - c2m .* qnm) - ynm), [], 'all'); 
    
    % Store minimum value in the list
    min_snm_ynm_values = [min_snm_ynm_values, min_snm_ynm]; 

    %disp(['Iteration ', num2str(iter), ': Minimum value of snm-ynm = ', num2str(min_snm_ynm)]);
    %disp(['Objective Function - Iteration ', num2str(iter), ': Minimum value of snm-ynm = ', num2str(min_snm_ynm)]);
    
    % Compute the objective function value 
    f_value = sum(squared_diff); 
end 

function grad = finite_difference_gradient(x, alpha, M, TR, ynm, min_snm_ynm_values, iter) 
    % Compute gradient using finite differences 
    grad = zeros(size(x)); 
    h = 1e-6; % Step size for finite differences 

    %for i = 1:numel(x)
    %    x_plus_h = x; 
    %    x_plus_h(i) = x_plus_h(i) + h; 
    %    grad(i) = (objective_function(x_plus_h, alpha, M, ynm, min_snm_ynm_values, iter) - objective_function(x, alpha, M, ynm, min_snm_ynm_values, iter)) / h; 
    %end 

    for i = 1:numel(x)
        x_plus_h = x; 
        x_plus_h(i) = x_plus_h(i) + h; 
        
        obj_func_before = objective_function(x, alpha, M, TR, ynm, min_snm_ynm_values, iter);
        obj_func_after = objective_function(x_plus_h, alpha, M, TR, ynm, min_snm_ynm_values, iter);
        
        grad(i) = (obj_func_after - obj_func_before) / h; 
        
        %disp(['Finite Difference Gradient - Iteration ', num2str(iter), ': Calculating gradient for parameter ', num2str(i)]);
    end 
end 

 

 