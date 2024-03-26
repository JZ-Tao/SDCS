function X_out = SDCS_FCSA(X_LR_obs, ratio, lambda, Opts)
% The code is modified from Chen's SIRF and Beck's tv_fista.
NL_weight = Opts.NL_weight;

n_iter = Opts.n_iter;
MTF_info = Opts.MTF_info;
TV_lambda = lambda;

DP = Opts.DP;

X = interpWrapper(X_LR_obs, ratio, MTF_info);
method = Opts.method;
Y = X;


[n_r, n_c, n_b] = size(X);
t_new = 1;
L = 1;
P1 = [];
P2 = [];
P = zeros([n_r n_c n_b 2]); P0 = P;
N_pix = n_r*n_c*n_b;

itr = 1;
residual = inf;
tol = Opts.iter_tol;
block_init = 0;
n_NL_runs = 0;

while((residual>tol) && (itr<=n_iter))
    % store the old value of the iterate and the t-constant    
    X_old = X;
    t_old = t_new;
    % gradient step
    D_LR = conv_downsample_MTF(X,ratio,MTF_info) - X_LR_obs;
    D = interpWrapper(D_LR, ratio, MTF_info);
    Y = Y - D/L;
    ztilde = Y-DP;
    
    [M_Z, PS] = mapminmax(im2mat(ztilde),0,1);
    ztilde = mat2im(M_Z, size(ztilde,1));

    if Opts.cascade 
        if residual>tol*1.5 && itr < n_iter*0.9 && n_NL_runs == 0
            [X, P1, P2]=denoise_TV_MT(ztilde, 2*TV_lambda./L, P1, P2,Opts.TV);
        else
            if block_init == 0
                block_init = 1;
                [X, BlockMatches] = BM3D(ztilde, TV_lambda,'np',Opts.NL_Phase, {1,Opts.NL_Phase == BM3DProfile.ALL_STAGES});
            else
                X = BM3D(ztilde, TV_lambda, 'np', Opts.NL_Phase, BlockMatches);
            end
            if n_NL_runs == 0
                disp('Final denoiser');
            end
            n_NL_runs = n_NL_runs+1;
            if n_NL_runs > 0.1*n_iter
                break;
            end
        end
    else
        switch method
            case 'TV'
                [X, P1, P2]=denoise_TV_MT(ztilde, 2*TV_lambda./L, P1, P2,Opts.TV);
            case 'TV_NL'
                [z1, P1, P2]=denoise_TV_MT(ztilde, 2*TV_lambda./L, P1, P2,Opts.TV);
                X = NL_weight*wrapper_BM3D(ztilde, 2*TV_lambda./L) + (1-NL_weight)*z1;
            case 'BM3D'
                X=wrapper_BM3D(ztilde, 2*TV_lambda./L);
            otherwise
                X      = denoise(ztilde,2*TV_lambda./L);
        end
    end
    X = mat2im(mapminmax('reverse', im2mat(X), PS), size(X,1));
    X = X+DP;
    % updating t and Y
    t_new = (1+sqrt(1+4*t_old^2))/2;
    Y = X+((t_old-1)/t_new)*(X-X_old);
    
    %calculate residual
    residual = (1/sqrt(N_pix))*(sqrt(sum(sum(sum((X-X_old).^2))))); %fun_val_old-fun_val;
    
    if Opts.print==true
        fprintf('%3g \t %3.5e \n', itr, residual);
    end
    itr = itr+1;
end
X_out = X;

end