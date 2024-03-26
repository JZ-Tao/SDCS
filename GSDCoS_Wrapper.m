function I_HS_fused = GSDCoS_Wrapper(I_HS_LR, I_PAN, ratio, Opts)

basis_type = Opts.SS.basis_type; 
add_rest_sub = Opts.SS.add_rest_sub;
n_sub = Opts.SS.n_sub;
MTF_info = Opts.MTF_info;
L = Opts.L;
max_v = max(2^L-1, 1);
I_HS_LR = I_HS_LR./max_v;
I_PAN = I_PAN./max_v;

DP = bandMatching(I_HS_LR, I_PAN, ratio, Opts);

I_HS = interpWrapper(I_HS_LR, ratio, MTF_info);
DP_LR = conv_downsample_MTF(DP, ratio, MTF_info);
[n_row, n_col, n_band] = size(I_HS);

if Opts.outer_diff
    X_LR = I_HS_LR - DP_LR;
    X = interpWrapper(X_LR, ratio, MTF_info);
else
    X = I_HS;
    X_LR = I_HS_LR;
end
    
if Opts.normalize_type == 1
    [M_X_LR, PS] = mapminmax(im2mat(X_LR),0,1);
    X_LR = mat2im(M_X_LR, size(X_LR,1));
end
if Opts.cluster_on
    cluster_p = Opts.cluster_p;
    C_map = ClusteringHS(X_LR, Opts.cluster_method, cluster_p);
    C_map_U = imresize(C_map, ratio, 'nearest');    
    n_cluster = max(C_map(:));
    color_map = GetColorMap(n_cluster);
    class_map = GetClassMap(C_map-1, color_map);
    figure, imshow(class_map);

    [E, iE, E_r, iE_r, ~] = ClusterwiseSubspace(C_map, X_LR, basis_type, n_cluster, Opts.SS);
    [iCETrans, CETrans] = ClusterizeTransSelector(C_map, C_map_U, n_cluster, basis_type, n_sub, size(I_HS,3));
else
    [E, iE, E_r, iE_r, ~] = genSubspace(X_LR, basis_type, Opts.SS, Opts.tag);
    [iCETrans, CETrans]= subspaceTransSelector(basis_type, n_sub);
end

DP_sub = iCETrans(iE, DP);

X_LR_sub = iCETrans(iE, X_LR);

if add_rest_sub
    % TODO: Only valid for PCA.
    if ~isempty(iE_r) && ~isempty(E_r)%strcmp(basis_type, 'SVD') || strcmp(basis_type, 'Hysime')
        I_X_r_sub = iCETrans(iE_r, X);
    else
        I_X_r_sub = [];
    end
end
if Opts.outer_diff
    Opts.DP = 0;
else
    Opts.DP = DP_sub;
    DP = 0;
end

lambda = Opts.TV_lambda;
switch(Opts.solver)
    case 'ADMM'
        funcWrapper = @PnP_ADMM;
    case 'Fista'
        funcWrapper = @DCoS_FCSA;
end


if Opts.normalize_type ~= 0
    [M_X_LR_sub, PS2] = mapminmax(im2mat(X_LR_sub),0,1);
    X_LR_sub = mat2im(M_X_LR_sub, size(X_LR_sub,1));
end
if ~Opts.cluster_on || ~Opts.cluster_optim
    X_sub_hat = funcWrapper(X_LR_sub, ratio, lambda, Opts);
elseif strcmp(Opts.cluster_method, 'block')
    block_func = @(block_struct) block_optimize(block_struct.data, ratio, lambda, Opts);
    %
    X_sub_hat = blockproc(X_LR_sub,[cluster_p cluster_p],block_func,'UseParallel',true);

else
    X_sub_hat_all = zeros(n_row, n_col, n_sub, n_cluster);
    parfor ii = 1 : n_cluster
        t = X_LR_sub.*(C_map == ii);
        X_sub_hat_all(:,:,:,ii) = funcWrapper(t, ratio, lambda, Opts);
    end
    X_sub_hat = sum(X_sub_hat_all,4);
end
if Opts.normalize_type ~= 0
    X_sub_hat = mat2im(mapminmax('reverse', im2mat(X_sub_hat), PS2), size(X_sub_hat,1));
end

if add_rest_sub
    for i=1:size(I_X_r_sub,3)
      [thr,sorh,keepapp] = ddencmp('den','wv',I_X_r_sub(:,:,i));
      I_X_r_sub(:,:,i) = wdencmp('gbl',I_X_r_sub(:,:,i),'sym4',4,thr,sorh,keepapp);
    end
    % 
    X_f_sub = cat(3, X_sub_hat, I_X_r_sub);
    E_f = [E E_r];
    imX_sub = CETrans(E_f, X_f_sub);
else
    E_f = E;
    imX_sub = CETrans(E, X_sub_hat);
end
    
if Opts.normalize_type == 1
    imX_sub = mat2im(mapminmax('reverse', im2mat(imX_sub), PS), size(imX_sub,1));
elseif Opts.normalize_type == 2
    imX_sub = imX_sub.*min_max_r(3) + min_max_r(1);
end
I_HS_fused = (imX_sub + DP).*max_v;


    function X_sub_hat_i = block_optimize(X_LR_sub, ratio, lambda, Opts)
        X_sub_hat_i = funcWrapper(X_LR_sub, ratio, lambda, Opts);
    end    

end