function [ K, T1, min_snm_ynm_values ] = novifast_image( Im, alpha, TR, options , varargin )
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
    
    disp(options);
    % Checks the options structure passed to the function and sets the convergence parameters accordingly.
    if ~isfield(options, 'Direct')          % Checks if the field `'Direct'` is not present in the options structure.
        if ~isfield(options, 'MaxIter')
            options.MaxIter = 10; %default
        elseif options.MaxIter<1 || options.MaxIter>200
            error('options: Maxiter should be set to a value between 1 and 200');
        end

        if ~isfield(options,'Tol')
            options.Tol = 1e-6; %default
        elseif options.Tol<0 || options.Tol>1e-2
            error('options: Tol should be set to a value between 0 and 1e-2');
        end
        modeDirect=false;       % modeDirect to false. This indicates that the convergence mode is not set to direct.
    elseif options.Direct<1 || options.Direct>200
        error('options: Directiter should be set to a value between 1 and 200');
    else
        modeDirect=true;
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
    %% NOVIFAST begins here

    % pre-computing
    K=squeeze(zeros(nrows,ncols,nslices));
    T1=squeeze(zeros(nrows,ncols,nslices));

    % Initialize min_snm_ynm_values array to store minimum values of snm-ynm
    min_snm_ynm_values = [];

    alphanm=alpha*ones(1,M);
    y=reshape(Im(:),nrows*ncols*nslices,N);
    ynm=y(pm,:)';
    pnm=sind(alphanm);
    qnm=cosd(alphanm);

    % initialization
    Kini=ini(1);
    T1ini=ini(2);
    c1m=Kini*(1*exp(-TR/T1ini))*ones(1,M);
    c2m=exp(-TR/T1ini)*ones(1,M);
    k=0;
    done=false;


    % iterative process
    while ~done

        c2m_old=c2m;
        c1m=repmat(c1m,[N,1]);
        c2m=repmat(c2m,[N,1]);
        denm=1-c2m.*qnm;
        snm=c1m.*pnm./denm;

        %definition of vectors
        A=ynm.*qnm./denm;
        Ahat=snm.*qnm./denm;
        B=pnm./denm;
        Z=ynm./denm;


        %definition of inner products
        BB=sum(B.^2,1);
        AAhat=sum(A.*Ahat,1);
        BAhat=sum(B.*Ahat,1);
        BA=sum(B.*A,1);
        BZ=sum(B.*Z,1);
        ZAhat=sum(Z.*Ahat,1);

        %calculation of c1m and c2m
        detm=BB.*AAhat- BAhat.*BA;
        c1m=(BZ.*AAhat - ZAhat.*BA)./detm;
        c2m=(BB.*ZAhat - BAhat.*BZ)./detm;
        k=k+1;

        min_snm_ynm = min(abs(snm - ynm), [], 'all');
        min_snm_ynm_values = [min_snm_ynm_values, min_snm_ynm];

        disp(['Iteration ', num2str(k), ': Minimum value of snm-ynm = ', num2str(min_snm_ynm)]);
        
        %stopping
        if modeDirect %mode with no-convergence criterion
            if k==options.Direct
                done=true;
            end
        else
            rel_err=( norm(c2m(isfinite(c2m))-c2m_old(isfinite(c2m_old)),1) )/ norm(c2m(isfinite(c2m)),1);  %Relative l1 norm for c2 (Convergence is controlled by c2 only)
            if rel_err< options.Tol || k>=options.MaxIter %mode with convergence criterion
                fprintf("Convergence Criteria Triggered.")
                done=true;
            end
        end  
    end

    Km=c1m./(1-c2m);
    T1m=-TR./log(c2m);

    %K and T1 maps
    T1(pm)=T1m;
    K(pm)=Km;

end