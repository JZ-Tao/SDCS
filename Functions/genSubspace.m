function [E, iE, E_r, iE_r, eig] = genSubspace(X, basis_type, Opts, tag)
if isstruct(Opts)
    if ~isfield(Opts, 'n_sub'), Opts.n_sub = 10; end
    n_sub = Opts.n_sub;
else
    n_sub = Opts;
end
if ~exist('tag', 'var')
    tag = [];
end
if size(X,3) ~= 1
    m_X = im2mat(X);
else
    m_X = X; % Already in 2D matrix form.
end
eig = 0;
E_r = [];
iE_r = [];
switch  basis_type 
%     case 'SPCA_PSP'
%         [U1,V1,U2,V2] = Sparse_PCA(m_X, lambda, K, iter1, iter2, tol1, tol2)
%     case 'SPCA_NIPM'
%         [F, adj_var, cum_var] = sparsePCA(m_X, card, num_comp, num_runs, verbosity);
    case 'VCA'
%     Find endmembers with VCA (pick the one with smallest volume from 20 
%     runs of the algorithm)
        max_vol = 0;
        vol = zeros(1, 20);
        for idx_VCA = 1:20
            E_aux = VCA(m_X,'Endmembers',n_sub,'SNR',0,'verbose','off');
            vol(idx_VCA) = abs(det(E_aux'*E_aux));
            if vol(idx_VCA) > max_vol
                E = E_aux;
                max_vol = vol(idx_VCA);
            end   
        end
        iE = pinv(E);	%((E'*E)^-1)*E'

    case 'Hysime'
    	[w,Rn] = estNoise(m_X);
        [n_sub_opt, E_opt] = hysime(m_X,w,Rn,'off');%diag(sigma2y_real)
        
        E = E_opt(:, 1:n_sub);
        E_r = E_opt(:, n_sub+1:end);
        iE = E'; %pinv(E);	%((E'*E)^-1)*E'
        iE_r = E_r';
        disp('n_sub by hysime:');

    case 'SVD'

        [E_full,eig_full,~] = svd(m_X*m_X'/size(m_X,2));
        E = E_full(:, 1:n_sub);
        E_r = E_full(:, n_sub+1:end);

        eig_full = diag(eig_full);
        eig = eig_full(1:n_sub);
        iE = E';
        iE_r = E_r';
%     case 'NLPCA'
%         prePCA = 1;
%         file_name = ['NLPCA_net_' tag '_sub_' num2str(n_sub) '_prePCA_' num2str(prePCA) '.mat'];
%         if 1 %~exist(file_name,'file')
%             n_band = size(m_X,1);
%             if prePCA
%                 lf = floor(n_band*0.5);
%                 lh = min(2+n_sub*2, floor(lf*0.5));
%                 [~, iE, ~] = nlpca(m_X, n_sub,'pre_pca','yes','units_per_layer',[lf, lh, n_sub, lh, lf]);
%             else
%                 lf = n_band;
%                 lh = min(2+n_sub*2, floor(lf*0.5));
%                 [~, iE, ~] = nlpca(m_X, n_sub,'units_per_layer',[lf, lh, n_sub, lh, lf]);
%             end
%             save(file_name, 'iE');
%         else
%             iE = importdata(file_name);
%         end
%         E = iE;
%     case 'PPA'
%         order = 2;
%         [~, iE] = PPA(m_X,order,n_sub);
%         %Morden = 10;
%         %[~, iE] = PPA_auto(m_X,[],Morden,n_sub);
%         
%         E = iE;
%     case 'DRR'
%         file_name = ['DRR_' tag '_sub_' num2str(n_sub) '.mat'];
%         if 1%~exist(file_name,'file')
%             [~, iE] = DRR_method(m_X,n_sub);
%             save(file_name, 'iE');
%         else
%             iE = importdata(file_name);
%         end
%         E = iE;
end


