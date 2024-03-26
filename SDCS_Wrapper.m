function I_HS_fused = SDCS_Wrapper(I_HS_LR, I_PAN, ratio, Opts)

basis_type = Opts.SS.basis_type; 
n_sub = Opts.SS.n_sub;


L = Opts.L;
max_v = max(2^L-1, 1);
I_HS_LR = I_HS_LR./max_v;
I_PAN = I_PAN./max_v;
% Band matching
disp('Band matching calculation in progress...');
DP = bandMatching(I_HS_LR, I_PAN, ratio, Opts);
disp('Band matching is complete.');
X_LR = I_HS_LR;

% Subspace processing
[E, iE, ~, ~, ~] = genSubspace(X_LR, basis_type, Opts.SS, Opts.tag);
[iCETrans, CETrans]= subspaceTransSelector(basis_type, n_sub);

DP_sub = iCETrans(iE, DP); 
X_LR_sub = iCETrans(iE, X_LR);

Opts.DP = DP_sub;
DP = 0;


lambda = Opts.TV_lambda;
switch(Opts.solver)
    case 'ADMM'
        funcWrapper = @SDCS_PnP_ADMM;
    case 'Fista'
        funcWrapper = @SDCS_FCSA;
end

X_sub_hat = funcWrapper(X_LR_sub, ratio, lambda, Opts);
imX_sub = CETrans(E, X_sub_hat);
I_HS_fused = (imX_sub + DP).*max_v;

end