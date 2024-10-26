function [ K, T1 ] = nelder_mead_image(Im, alpha, TR, options, varargin)

    % Initialization
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
            options.MaxIter = 10; %default
        elseif options.MaxIter < 1 || options.MaxIter > 200
            error('options: MaxIter should be set to a value between 1 and 200');
        end

        if ~isfield(options,'Tol')
            options.Tol = 1e-6; %default
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
        ini = [0.5,500];
        th = 0.05 * max(max(max(Im(:)))); %Intensity values smaller than 5% of the maximum value of the SPGR dataset are left out
        if nslices == 1
            mask = squeeze(Im(:,:,1)) > th;
        else
            mask = squeeze(Im(:,:,:,1)) > th;
        end
    elseif ~isvector(varargin{1})
        fprintf('User did not provide initial values. Using default parameters... \n');
        ini = [0.5,500];
        mask = varargin{1};
    else
        ini = varargin{1};
        if length(varargin) == 1
            fprintf('User did not provide a mask. Using default parameters... \n');
            th = 0.05 * max(max(max(Im(:)))); %Intensity values smaller than 5% of the maximum value of the SPGR dataset are left out
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

    % Nelder Mead begins here

    % Pre-computing
    K = squeeze(zeros(nrows, ncols, nslices));
    T1 = squeeze(zeros(nrows, ncols, nslices));

    % Row vector alphanm with the same size as alpha but with repeated values of alpha
    alphanm = alpha * ones(1, M);
    y = reshape(Im(:), nrows*ncols*nslices, N);
    ynm = y(pm,:)';
    pnm = sind(alphanm);
    qnm = cosd(alphanm);

    % initialization
    Kini = ini(1);
    T1ini = ini(2);
    %c1m = repmat(Kini * (1 - exp(-TR / T1ini)), [N, 1]);  % Initialize c1m as a matrix with dimensions 11x522467
    c1m = Kini * (1 - exp(-TR / T1ini)) * ones(1, M);
    %c2m = exp(-TR / T1ini) * ones(1, M);
    c2m = exp(-TR / T1ini) * ones(1, M);


    % Initialize simplex with initial parameter values
    simplex = repmat(ini', 1, N + 1);
    
    k = 0;
    done = false;
    
    while ~done
        c1m_old = c1m;
        c2m_old = c2m;
        % Compute predicted signal intensities for all points in the simplex
        predicted_signals = zeros(M, N + 1);
        for i = 1:N+1
            K = simplex(1, i);
            T1 = simplex(2, i);
            
            % Compute c1m and c2m
            c1m_i = (K * (1 - exp(-TR / T1)) * ones(1, M))';    % Initialize c1m as a matrix with dimensions 11x522467
            c2m_i = (exp(-TR / T1) * ones(1, M))';  
            %denm_i = (1 - c2m_i .* qnm)';  
            denm_i = (1 - (c2m_i .* qnm)')';

            
            % Ensure pnm, c1m, and denm have compatible sizes
            if size(pnm, 2) == M
                pnm_i = pnm;
            else
                pnm_i = pnm(:, 1);  % Use the first column of pnm if necessary
            end

            % Print sizes for debugging
            disp(['Size of pnm_i: ', num2str(size(pnm_i))]);
            disp(['Size of c1m_i: ', num2str(size(c1m_i))]);
            disp(['Size of denm_i: ', num2str(size(denm_i))]);

                 
            % Replicate c1m_i_replicated along the second dimension to match the size of pnm_i
            %c1m_i_replicated = repmat(c1m_i, 1, size(pnm_i, 2));
            %disp(['Size of c1m_i_replicated: ', num2str(size(c1m_i_replicated))]);
            disp(['Size of pnm_i: ', num2str(size(pnm_i))]);
            disp(['Size of denm_i: ', num2str(size(denm_i))]);

            % Compute predicted signals
            predicted_signals(:, i) = (c1m_i .* pnm_i) ./ denm_i;
        end
        
        % Compute objective values
        obj_values = sum((predicted_signals - ynm).^2, 1);
        
        [best_val, best_idx] = min(obj_values);
        [worst_val, worst_idx] = max(obj_values);
        
        % Print best and worst point values
        disp(['Iteration ', num2str(k)]);
        disp(['Best point: ', num2str(simplex(:, best_idx)'), ', Objective value: ', num2str(best_val)]);
        disp(['Worst point: ', num2str(simplex(:, worst_idx)'), ', Objective value: ', num2str(worst_val)]);
        
        % Compute centroid
        centroid = mean(simplex(:, [1:worst_idx-1, worst_idx+1:end]), 2);
        
        % Reflect worst point through centroid
        reflected_point = centroid + (centroid - simplex(:, worst_idx));
        c1_reflected = reflected_point(1) * ones(1, M);
        c2_reflected = reflected_point(2) * ones(1, M);
        denm_reflected = 1 - c2_reflected .* qnm;
        predicted_signal_reflected = c1_reflected .* pnm ./ denm_reflected;
        reflected_val = sum((predicted_signal_reflected - ynm).^2);
        
        % Print reflected point and its value
        disp(['Reflected point: ', num2str(reflected_point'), ', Objective value: ', num2str(reflected_val)]);
        
        % If reflected point is better than the second worst, try expanding
        if best_val <= reflected_val && reflected_val < obj_values(end)
            expanded_point = centroid + 2 * (reflected_point - centroid);
            c1_expanded = expanded_point(1) * ones(1, M);
            c2_expanded = expanded_point(2) * ones(1, M);
            denm_expanded = 1 - c2_expanded .* qnm;
            predicted_signal_expanded = c1_expanded .* pnm ./ denm_expanded;
            expanded_val = sum((predicted_signal_expanded - ynm).^2);
            
            % Print expanded point and its value
            disp(['Expanded point: ', num2str(expanded_point'), ', Objective value: ', num2str(expanded_val)]);
            
            if expanded_val < reflected_val
                simplex(:, end) = expanded_point;
            else
                simplex(:, end) = reflected_point;
            end
        % If reflected point is worse than worst, try contracting
        elseif reflected_val >= obj_values(end)
            contracted_point = centroid + 0.5 * (simplex(:, worst_idx) - centroid);
            c1_contracted = contracted_point(1) * ones(1, M);
            c2_contracted = contracted_point(2) * ones(1, M);
            denm_contracted = 1 - c2_contracted .* qnm;
            predicted_signal_contracted = c1_contracted .* pnm ./ denm_contracted;
            contracted_val = sum((predicted_signal_contracted - ynm).^2);
            
            % Print contracted point and its value
            disp(['Contracted point: ', num2str(contracted_point'), ', Objective value: ', num2str(contracted_val)]);
            
            if contracted_val < obj_values(end)
                simplex(:, end) = contracted_point;
            % If contraction is not better, perform shrinkage
            else
                for i = 2:size(simplex, 2)
                    simplex(:, i) = 0.5 * (simplex(:, 1) + simplex(:, i));
                end
            end
        end
        
        k = k + 1;
        
        % Check termination criteria
        if modeDirect
            if k == options.Direct
                done = true;
            end
        else
            rel_err = norm(simplex(:, worst_idx) - simplex(:, best_idx), 1) / norm(simplex(:, best_idx), 1);
            if rel_err < options.Tol || k >= options.MaxIter
                done = true;
            end
        end
    end
    
    % Extract optimized parameter values from best point in the simplex
    K(:, :, :) = simplex(1, best_idx);
    T1(:, :, :) = simplex(2, best_idx);

end
