function Out = QualityIndices(I_HS,I_REF,ratio, flag_cut, dim_cut)
%--------------------------------------------------------------------------
% Quality Indices
%
% USAGE
%   Out = QualityIndices(I_HS,I_REF,ratio)
%
% INPUT
%   I_HS  : target HS data (rows,cols,bands)
%   I_REF : reference HS data (rows,cols,bands)
%   ratio : GSD ratio between HS and MS imagers
%	is_per_band : 0 for the default index calculation / 1 for index per band.
% OUTPUT
%   Out.cc   : CC
%   Out.sam  : SAM
%   Out.rmse : RMSE
%   Out.ergas: ERGAS
%	Out.psnr : PSNR
%--------------------------------------------------------------------------
if ~exist('flag_cut','var')
    flag_cut = 0;
end
if ~exist('dim_cut','var')
    dim_cut = 0;
end

%Remove border from the analysis
if flag_cut
    I_HS  = I_HS(1+dim_cut:end-dim_cut,1+dim_cut:end-dim_cut,:);
    I_REF = I_REF(1+dim_cut:end-dim_cut,1+dim_cut:end-dim_cut,:);
end

cc = CC(I_HS,I_REF);
Out.cc = mean(cc);
[angle_SAM,map] = SAM(I_HS,I_REF);
Out.sam = angle_SAM;
Out.sam_map = map;
Out.rmse = RMSE(I_HS,I_REF);
Out.ergas = ERGAS(I_HS,I_REF,ratio);
psnr = PSNR(I_REF,I_HS);
Out.psnrall = psnr.all;
Out.psnr = psnr.ave;


disp(['CC   : ' num2str(Out.cc)]);
disp(['SAM  : ' num2str(Out.sam)]);
disp(['RMSE : ' num2str(Out.rmse)]);
disp(['ERGAS: ' num2str(Out.ergas)]);
disp(['PSNR: ' num2str(Out.psnr)]);
