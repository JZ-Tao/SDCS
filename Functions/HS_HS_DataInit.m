function Dataset = HS_HS_DataInit(data_root_path, Dataset)

%% import HS data

file_full_path = [data_root_path '/HS/' Dataset.filename];
if ~exist(file_full_path,'file')
    error(['file: ' file_full_path ' not exist.']);
else
    data = double(importdata(file_full_path));
    max_d = max(data(:));
    if ~isfield(Dataset, 'L_org')
        nn = 0;
        while(max_d > 2^nn-1)
            nn = nn+1;
        end
        Dataset.L_org = nn;
    end
    if Dataset.cropped
        pos = Dataset.cropped_pos;
        sz = Dataset.cropped_size;
        data = data(pos(1):pos(1)+sz(1)-1, ...
            pos(2):pos(2)+sz(2)-1, :);
    end
    if ~isfield(Dataset, 'n_row'), Dataset.n_row = size(data,1); end
    if ~isfield(Dataset, 'n_column'), Dataset.n_column = size(data,2); end
    if ~isfield(Dataset, 'n_band'), Dataset.n_band = size(data,3); end
    if ~isfield(Dataset, 'band_set')
        band_idx_1 = floor(Dataset.n_band/2); % simply half-half
        Dataset.band_set = {(1:band_idx_1)'; (band_idx_1+1:Dataset.n_band)'}; 
    end

    data = data(1:Dataset.n_row, 1:Dataset.n_column, :);
    Dataset.min_v = min(data(:));
    Dataset.max_v = max(data(:));


    switch(Dataset.scaling_type)
        case 0
            Dataset.scaling = 1;
            Dataset.REF = data;
        case 1
            Dataset.scaling = 2^(Dataset.L_org)-1;
            Dataset.REF = data/Dataset.scaling;
            Dataset.REF(Dataset.REF<0) = 0;
            Dataset.REF(Dataset.REF>Dataset.scaling) = Dataset.scaling;
    end
end

if size(Dataset.band_set,1) == 2
    % HR VNIR
    Dataset.PAN = Dataset.REF(:,:,Dataset.band_set{1});
    % HR SWIR
    Dataset.REF = Dataset.REF(:,:,Dataset.band_set{2});
else
    error('Incorrect band set setting. Must be 2 cells');
    % Generate band sets (VNIR | SWIR) according to corr. coef, or simply
    % divide into two halfs.
%     HS_ccvals = zeros(Dataset.n_band, Dataset.n_band);
%     for i = 1:Dataset.n_band
%         for j = 1:Dataset.n_band
%             a = Dataset.REF(:,:,i);
%             b = Dataset.REF(:,:,j);
%             cc = corrcoef(a(:), b(:));
%             HS_ccvals(i,j) = cc(1,2);
%         end
%     end
%     figure, imshow(HS_ccvals, []);
end
%% Generating the spatial degraded data from the reference
% LR SWIR
[Dataset.HS]=conv_downsample_MTF(Dataset.REF,Dataset.ratio,Dataset.MTF_info);
%% Note that there may be no spectral degradation part.
Dataset.sp_overlap = intersect(Dataset.band_set{1}, Dataset.band_set{2});
% Using Yokoya's code for R estimating. Not directly used.
R_est = estR(Dataset.HS,Dataset.PAN);
Dataset.R = R_est;

