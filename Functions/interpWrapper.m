function Y = interpWrapper(X, ratio, MTF_info)
if ~exist('MTF_info','var')
    MTF_info.start_pos = [1,1];
end
% if ~exist('tap','var')
%     tap = -1;
% end

Y = interp23tapGeneral(X, ratio, MTF_info.start_pos);
