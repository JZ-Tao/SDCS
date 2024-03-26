function [Opts, arg_str] = initSDCoSOptByTestCase(varargin)
Opts.normalize_type = 0;
Opts.solver = 'Fista';
Opts.NL_sigma = 0.02;
Opts.NL_weight = 0;%[0.92];
Opts.Sylvester = false;
Opts.TV_lambda = 0.0005;
Opts.TV.MAXITER = 3;
Opts.TV.type = 'iso';
Opts.n_iter = 150;
Opts.cluster_on = false;
Opts.cluster_optim = false;
Opts.cluster_method = 'block';
Opts.cluster_p = 1;
Opts.BM.method = 'BM_org';
Opts.outer_diff = 0;
Opts.SS.add_rest_sub = 0;
Opts.SS.prePCA = 0;
Opts.SS.n_sub = 10;
Opts.SS.basis_type = 'SVD'; 
Opts.cascade = 0;
Opts.v_field_switch = 0;
Opts.tol_switch = 0;
Opts.tol = [0.01 0.99];
Opts.BP = 0;
Opts.NL_Phase = BM3DProfile.ALL_STAGES;
arg_str = [];
nn = nargin;
Opts.iter_tol = 1e-5;
for k = 1:nargin
    for j = 1:length(varargin{k})
        arg = varargin{k}{j};
        if iscell(arg) arg = arg{1}; end
        switch(arg)
            case 'iter_tol_off'
                Opts.iter_tol = 1e-12;
            case 'iter_tol_on'
                Opts.iter_tol = 1e-5;
            case 'BP'
                Opts.BP = 1;
            case 'SSBP'
                Opts.BP = 2;
            case 'ADMM'
                Opts.solver = 'ADMM';
            case 'Fista'
                Opts.solver = 'Fista';
            case 'ADMM_sy'
                Opts.solver = 'ADMM_sy';
            case 'dTV'
                Opts.v_field_switch = 1;
            case 'cascade' % cascading TV and NL
                Opts.cascade = 1;
            case 'NL_ALL'
                Opts.NL_Phase = BM3DProfile.ALL_STAGES;
            case 'NL_TH'
                Opts.NL_Phase = BM3DProfile.HARD_THRESHOLDING;      
            case 'NL_only'
                Opts.NL_weight = 1;
            case 'TV_only'
                Opts.NL_weight = 0;
            case 'X_on'
                Opts.outer_diff = 1;
            case 'rest_sub'
                Opts.SS.add_rest_sub = 1;
            case 'prePCA'
                Opts.SS.prePCA = 1;
            case 'Hysime'
                Opts.SS.basis_type = 'Hysime';
            case 'VCA'
                Opts.SS.basis_type = 'VCA';
            case 'SVD'
                Opts.SS.basis_type = 'SVD';        
            case 'NLPCA'
                Opts.SS.basis_type = 'NLPCA';   
            case 'PPA'
                Opts.SS.basis_type = 'PPA';
            case 'DRR'
                Opts.SS.basis_type = 'DRR';
            case 'Sylvester'
                Opts.Sylvester = true;
            case 'block'
                Opts.cluster_method = 'block';
            case 'kmeans'
                Opts.cluster_method = 'kmeans';
            case 'super'
                Opts.cluster_method = 'super';
            case 'cluster_on'    
                Opts.cluster_on = true;
            case 'cluster_optim'
                Opts.cluster_optim = true;
            case 'N'
                Opts.normalize_type = 1;
            otherwise
                str_TV_NL = regexp(arg, 'TV_NL_(\d*)','tokens');
                if ~isempty(str_TV_NL)
                    Opts.NL_weight = str2double(str_TV_NL{1}{1})/100;
                    arg_str = [arg_str '_' arg];
                    continue;
                end
                str_MAXITER = regexp(arg, 'MAXITER_(\d*)','tokens');
                if ~isempty(str_MAXITER)
                    Opts.MAXITER = str2double(str_MAXITER{1}{1});
                    arg_str = [arg_str '_' arg];
                    continue;
                end
                str_lambda = regexp(arg, 'lambda_(.*)$','tokens');
                if ~isempty(str_lambda)
                    Opts.TV_lambda = str2double(str_lambda{1}{1});
                    arg_str = [arg_str '_' arg];
                    continue;
                end
                str_n_sub = regexp(arg, 'n_sub(\d*)','tokens');
                if ~isempty(str_n_sub)
                    Opts.SS.n_sub = str2double(str_n_sub{1}{1});
                    arg_str = [arg_str '_' arg];
                    continue;
                end
                str_n_iter = regexp(arg, 'n_iter(\d*)','tokens');
                if ~isempty(str_n_iter)
                    Opts.n_iter = str2double(str_n_iter{1}{1});
                    arg_str = [arg_str '_' arg];
                    continue;
                end
                str_BM = regexp(arg, 'BM_(.*)','tokens');
                if ~isempty(str_BM)
                    Opts.BM.method = arg;
                    arg_str = [arg_str '_' arg];
                    continue;
                end
                str_sigma = regexp(arg, 'sigma_(.*)','tokens');
                if ~isempty(str_sigma)
                    Opts.NL_sigma = str2double(str_sigma{1}{1});
                    arg_str = [arg_str '_' arg];
                    continue;
                end
                str_cluster_p = regexp(arg, 'cluster_p(\d*)','tokens');
                if ~isempty(str_cluster_p)
                    Opts.cluster_p = str2double(str_cluster_p{1}{1});
                    arg_str = [arg_str '_' arg];
                    continue;
                end
                str_iter_tol = regexp(arg, 'iter_tol_(.*)$','tokens');
                if ~isempty(str_iter_tol)
                    Opts.iter_tol = str2double(str_iter_tol{1}{1});
                    arg_str = [arg_str '_' arg];
                    continue;
                end
                error(['invalid arg:' arg]);
        end        
        arg_str = [arg_str '_' arg];
    end
end
  
if Opts.cluster_p > 1
    Opts.cluster_on = true;
end
Opts.TV_switch = 0; 
Opts.NL_switch = 0; 
if Opts.NL_weight == 0 || Opts.cascade == 1

    Opts.method = 'TV';
elseif Opts.NL_weight == 1

    Opts.method = 'BM3D';
else
    Opts.method = 'TV_NL';
end


