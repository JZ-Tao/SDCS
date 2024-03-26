function X_sub = iETransfrom(iE, X, basis_type, seg_indices)
    if ~exist('seg_indices','var')
        seg_indices = [];
    end
    if ~strcmp(basis_type, 'SegPCA')
        X_sub = mat2im(iE*im2mat(X), size(X,1));
    else
        X_sub = [];
        for i=1:length(seg_indices)
            X_sub_seg = mat2im(iE{i}*im2mat(X(:,:,seg_indices{i})), size(X,1));
            X_sub = cat(3,X_sub, X_sub_seg);
        end
    end