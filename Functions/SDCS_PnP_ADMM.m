function out = SDCS_PnP_ADMM(y,ratio,lambda,Opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A modified code by Jingzhe Tao from:
% out = PlugPlayADMM_super(y,h,K,lambda,method,opts)
%
% inversion step: x=argmin_x(||Ax-y||^2+rho/2||x-(v-u)||^2)
% denoising step: v=Denoise(x+u)
%       update u: u=u+(x-v)
%
%Input:           y    -  the observed gray scale image
%                 h    -  blur kernel
%                 K    -  downsampling factor
%              lambda  -  regularization parameter
%              method  -  denoiser, e.g., 'BM3D'
%       opts.rho       -  internal parameter of ADMM {1}
%       opts.gamma     -  parameter for updating rho {1}
%       opts.maxitr    -  maximum number of iterations for ADMM {20}
%       opts.tol       -  tolerance level for residual {1e-4}   
%       ** default values of opts are given in {}. 
%
%Output:          out  -  recovered gray scale image 
%         
%
%Xiran Wang and Stanley Chan
%Copyright 2016
%Purdue University, West Lafayette, In, USA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check defaults
if ~isfield(Opts,'rho')
    Opts.rho = 1;
end
if ~isfield(Opts,'n_iter')
    Opts.n_iter = 20;
end
if ~isfield(Opts,'iter_tol')
    Opts.iter_tol = 1e-4;
end
if ~isfield(Opts,'gamma')
    Opts.gamma=1;
end
if ~isfield(Opts,'print')
    Opts.print = false;
end

% set parameters
n_iter   = Opts.n_iter;
tol       = Opts.iter_tol;
gamma     = Opts.gamma;
rho       = Opts.rho;

%initialize variables
[rows_in,cols_in, band_in] = size(y);
rows      = rows_in*ratio;
cols      = cols_in*ratio;
N         = rows*cols*band_in;

% method = Opts.method;
MTF_info = Opts.MTF_info;
h = MTF_info.BluKer;

[G,Gt]    = defGGt_mod_HS(MTF_info,ratio);
GGt       = constructGGt_mod_HS(h,ratio,rows,cols, band_in);
Gty       = Gt(y);
%v         = imresize(y,K);
v = interpWrapper(y, ratio, MTF_info);
x         = v;
u         = zeros(size(v));
residual  = inf;
method = Opts.method;
NL_weight = Opts.NL_weight;
DP = Opts.DP;

%set function handle for denoiser
switch method
    case 'TV_NL'
        denoise=@(in,sigma) NL_weight*wrapper_BM3D(in,sigma) + (1-NL_weight)*wrapper_TV(in,sigma);
    case 'BM3D'
        denoise=@wrapper_BM3D;
    case 'TV'
        denoise=@wrapper_TV;
    case 'NLM'
        denoise=@wrapper_NLM;
    case 'RF'
        denoise=@wrapper_RF;
    case 'FFDNet'
        denoise=@wrapper_FFDNet;
    otherwise
        error('unknown denoiser \n');
end

BlockMatches = [];

% main loop
if Opts.print==true
    fprintf('Plug-and-Play ADMM --- SDCS \n');
    fprintf('Denoiser = %s \n\n', method);
    fprintf('itr \t ||x-xold|| \t ||v-vold|| \t ||u-uold|| \n');
end

itr = 1;
P1 = [];
P2 = [];

n_NL_runs = 0;
block_init = 0;
while(residual>tol && (itr<=n_iter))
    %store x, v, u from previous iteration for psnr residual calculation
    x_old=x;
    v_old=v;
    u_old=u;

    %inversion step
    xtilde = v-u;
    rhs = Gty + rho*xtilde;
    x = real((rhs - Gt(ifft2(fft2(G(rhs))./(GGt + rho))))/rho);

    %denoising step
    vtilde = x+u;
    ztilde = vtilde - DP;
    [M_Z, PS] = mapminmax(im2mat(ztilde),0,1);
    ztilde = mat2im(M_Z, size(ztilde,1));    

    sigma  = sqrt(lambda/rho);
    if Opts.cascade 
        if residual>tol*1.5 && itr < n_iter*0.9 && n_NL_runs == 0
            [z, P1, P2]=denoise_TV_MT(ztilde, 2*sigma, P1, P2,Opts.TV);
        else
            if block_init == 0
                block_init = 1;
                [z, BlockMatches] = BM3D(ztilde, sigma,'np', Opts.NL_Phase, {1,Opts.NL_Phase == BM3DProfile.ALL_STAGES});
            else
                z = BM3D(ztilde, sigma, 'np', Opts.NL_Phase, BlockMatches);
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
                [z, P1, P2] = denoise_TV_MT(ztilde, 2*sigma, P1, P2,Opts.TV);
            case 'TV_NL'
                [z1, P1, P2] = denoise_TV_MT(ztilde, 2*sigma, P1, P2,Opts.TV);
                z = NL_weight*wrapper_BM3D(ztilde, 2*sigma) + (1-NL_weight)*z1;
            case 'BM3D'
                z = wrapper_BM3D(ztilde, 2*sigma);
            otherwise
                z = denoise(ztilde,sigma);
        end
    end
    z = mat2im(mapminmax('reverse', im2mat(z), PS), size(z,1));
    v = z+DP;
    %update langrangian multiplier
    u      = u + (x-v);
    
    %update rho
    rho=rho*gamma;
    
    %calculate residual
    residualx = (1/sqrt(N))*(sqrt(sum(sum(sum((x-x_old).^2)))));
    residualv = (1/sqrt(N))*(sqrt(sum(sum(sum((v-v_old).^2)))));
    residualu = (1/sqrt(N))*(sqrt(sum(sum(sum((u-u_old).^2)))));
    
    residual = residualx + residualv + residualu;

    if Opts.print==true % && mod(itr,10) == 0
        fprintf('%3g \t %3.5e \n', itr, residual);
    end
    
    itr=itr+1;
end

out=real(x);
end