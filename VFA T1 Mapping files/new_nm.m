function [ K, T1 ] = new_nm( Im, alpha, TR, options, varargin)
    
    
    % Check input
    sizeim = size(Im);
    disp(['Size of Im: ', num2str(sizeim)]);
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

    % NOVIFAST begins here
    
    % pre-computing
    % K and T1 initialized with zeros. Matrices with dim: [nrows, ncols,
    % nslices]
    K=squeeze(zeros(nrows,ncols,nslices));
    T1=squeeze(zeros(nrows,ncols,nslices));
    
    % row vector alphanm with the same size as alpha but with repeated values of alpha
    alphanm=alpha*ones(1,M);
    
    % Initialize simplex with initial parameter values
    %simplex = repmat(ini', 1, N + 1);
    % Initialize simplex with random samples
    %num_samples = N; % Number of random samples equals the number of parameters
    %simplex = zeros(length(ini), N + 1); % Initialize simplex matrix
    
    %for i = 1:N
    %    % Generate random samples within the range of parameter values
    %    random_sample = rand(size(ini)) .* (max(ini) - min(ini)) + min(ini);
    %    simplex(:, i+1) = random_sample; % Add random sample to simplex matrix
    %end
    
    %simplex(:, 1) = ini; % Set the first column of the simplex matrix to the initial parameters

    % Initialize simplex with Latin Hypercube Sampling (LHS)
    num_samples = N; % Number of LHS samples equals the number of parameters
    lhs_samples = lhsdesign(num_samples, length(ini)); % Generate LHS samples
    
    % Scale the LHS samples to the range of parameter values
    scaled_lhs_samples = lhs_samples .* (max(ini) - min(ini)) + min(ini);
    
    % Ensure that the number of elements in ini matches the number of samples
    if length(ini) ~= size(scaled_lhs_samples, 2)
        error('Number of elements in ini must match the number of samples');
    end
    
    % Initialize simplex with LHS samples
    simplex = [ini; scaled_lhs_samples'];



    %perturbation_factor = 0.1; % Adjust the perturbation factor as needed
    %perturbed_ini = ini .* (1 + perturbation_factor * rand(size(ini)) - perturbation_factor / 2);
    %simplex = repmat(perturbed_ini', 1, N + 1);

    
    
    % Initialize objective values
    obj_values = zeros(1, N + 1);
    for i = 1:N+1
        obj_values(i) = objective_function(simplex(:, i), Im, alphanm, TR);
    end
    
    k = 0;
    done = false;
    
    while ~done
        [best_val, best_idx] = min(obj_values);
        [worst_val, worst_idx] = max(obj_values);
        
        % Initialize centroid
        centroid = mean(simplex(:, [1:worst_idx-1, worst_idx+1:end]), 2);
        
        % Check the size of alpha
        size_alpha = size(alpha);
        
        % Check the size of centroid
        size_centroid = size(centroid);
        
        % Check the size of simplex(:, worst_idx)
        size_simplex_worst = size(simplex(:, worst_idx));
        
        % Display the sizes
        %disp(['Size of alpha: ', num2str(size_alpha)]);
        %disp(['Size of centroid: ', num2str(size_centroid)]);
        %disp(['Size of simplex(:, worst_idx): ', num2str(size_simplex_worst)]);

        
        %reflected_point = centroid + alpha * (centroid - simplex(:, worst_idx));
        %reflected_point = centroid + alpha .* (centroid - simplex(:, worst_idx));
        %reflected_point = centroid + alpha .* (centroid' - simplex(:, worst_idx));
        reflected_point = centroid + alpha' .* (centroid - simplex(:, worst_idx));

        reflected_val = objective_function(reflected_point, Im, alpha, TR);

        
        if best_val <= reflected_val && reflected_val < obj_values(end)
            simplex(:, end) = reflected_point;
            obj_values(end) = reflected_val;
        elseif reflected_val < best_val
            expanded_point = centroid + 2 * (reflected_point - centroid);
            expanded_val = objective_function(expanded_point, Im, alpha, TR);
            if expanded_val < reflected_val
                simplex(:, end) = expanded_point;
                obj_values(end) = expanded_val;
            else
                simplex(:, end) = reflected_point;
                obj_values(end) = reflected_val;
            end
        else
            if obj_values(end) <= reflected_val && reflected_val < obj_values(worst_idx)
                contracted_point = centroid + 0.5 * (reflected_point - centroid);
                contracted_val = objective_function(contracted_point, Im, alpha, TR);
                if contracted_val <= reflected_val
                    simplex(:, end) = contracted_point;
                    obj_values(end) = contracted_val;
                else
                    %contracted_point = centroid + 0.5 * (simplex(:, worst_idx) - centroid);
                    %contracted_val = objective_function(contracted_point, Im, alpha, TR);
                   
                    %if contracted_val < obj_values(end)
                    %    simplex(:, end) = contracted_point;
                    %    obj_values(end) = contracted_val;
                    %else
                    %    for i = 2:size(simplex, 2)
                    %        simplex(:, i) = 0.5 * (simplex(:, 1) + simplex(:, i));
                    %        obj_values(i) = objective_function(simplex(:, i), Im, alpha, TR);
                    %    end
                    %end
                    for i = 2:size(simplex, 2)
                        simplex(:, i) = 0.5 * (simplex(:, 1) + simplex(:, i));
                        obj_values(i) = objective_function(simplex(:, i), Im, alpha, TR);
                    end
                end
            end
        end
        
        k = k + 1;
        
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
    
    % Compute T1 and K based on best point in the simplex
    K(:, :, :) = simplex(1, best_idx);
    T1(:, :, :) = simplex(2, best_idx);
    % Compute T1 and K based on best point in the simplex
    %K(:, :, :) = reshape(simplex(1, best_idx), size(K));
    %T1(:, :, :) = reshape(simplex(2, best_idx), size(T1));


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
