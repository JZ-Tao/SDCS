function [G,Gt] = defGGt_mod_HS(MTF_info, ratio)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A modified code by Jingzhe Tao from:
% Operators for super-resolution
% Stanley Chan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G  = @(x) fdown(x, MTF_info, ratio);
Gt = @(x) upf(x, MTF_info, ratio);
end

function y = fdown(x,MTF_info,ratio)
% tmp = imfilter(x,h,'circular');
% y = downsample2(tmp,ratio);
y = conv_downsample_MTF(x,ratio,MTF_info);

end

function y = upf(x,MTF_info,ratio)

y = interpWrapper(x, ratio, MTF_info);

% [r, c, b] = size(x);
% y = zeros(ratio.*r, ratio.*c, b);
% y(MTF_info.start_pos(1):ratio:end,MTF_info.start_pos(2):ratio:end,:) = x;
% y = imfilter(y, MTF_info.BluKer, 'circular');
end
