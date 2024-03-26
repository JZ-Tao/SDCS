function X = ETransfrom(E, X_sub, basis_type, size_x, seg_indices)
    if ~exist('seg_indices','var')
        seg_indices = [];
    end
    if ~strcmp(basis_type, 'SegPCA')
        X = mat2im(E*im2mat(X_sub), size_x(1));
    else
        X = zeros(size_x);
        n_seg_sub = cell(1, size(X_sub,3)+1);
        n_seg_sub{1} = 0;
        for i=1:length(seg_indices)
            n_seg_sub{i+1} = size(E{i},2);
            X_seg = mat2im(E{i}*im2mat(X_sub(:,:,n_seg_sub{i}+1:n_seg_sub{i}+n_seg_sub{i+1})), size_x(1));
            X(:,:,seg_indices{i}) = X_seg;
        end
    end
end

