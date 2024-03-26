function [iETrans, ETrans] = subspaceTransSelector(type, n_sub)
if ~exist('type','var')
    type = 'linear';
end
switch type
    case {'linear','PCA','SVD','VCA','Hysime'}
        iETrans = @ (iE, x) mat2im(iE*im2mat(x), size(x,1));
        ETrans = iETrans;
%     case 'NLPCA'
%         iETrans = @ (iE, x) mat2im(nlpca_get_components(iE, im2mat(x)), size(x,1));
%         ETrans = @ (iE, x) mat2im(nlpca_get_data(iE, im2mat(x)), size(x,1));
%     case 'PPA'
%         iETrans = @ (iE, x) mat2im(apply_PPA(im2mat(x), iE), size(x,1));
%         ETrans = @ (iE, x) mat2im(inv_PPA(im2mat(x), iE), size(x,1));
%     case 'DRR'
% %         iETrans = @ (iE, x) mat2im(apply_DRR_method(im2mat(x), iE), size(x,1));
%          ETrans = @ (iE, x) mat2im(inv_DRR_method(im2mat(x), iE), size(x,1));
%         iETrans = @ (iE, x) TransDRR(iE, x, n_sub);
%         ETrans = @ (iE, x) iTransDRR(iE, x, n_band);
end
% 
% function x_s = TransDRR(iE, x, n_sub)
%     t = apply_DRR_method(im2mat(x), iE);
%     x_s = mat2im(t(1:n_sub,:), size(x,1));
% end
% 
% function x_s = iTransDRR(iE, x)
%     x_m = im2mat(x);
%     x_s = mat2im(inv_DRR_method([x_m; zeros(n_band-size(x_m,1),size(x_m,2))], iE), size(x,1));
% end

end