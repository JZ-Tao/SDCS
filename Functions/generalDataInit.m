function Dataset = generalDataInit(data_root_path, Dataset)

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
    Dataset.n_row = size(data,1);
    Dataset.n_column = size(data,2);
    if ~isfield(Dataset, 'n_band'), Dataset.n_band = size(data,3); end
    if ~isfield(Dataset, 'band_set'), Dataset.band_set = 1:Dataset.n_band; end
   
    data = data(1:Dataset.n_row, 1:Dataset.n_column, Dataset.band_set);
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
        otherwise
            error('Unknown scaling type');
    end
end

%% Generating the spatial degraded data from the reference
[Dataset.HS]=conv_downsample_MTF(Dataset.REF,Dataset.ratio,Dataset.MTF_info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generating the spectral degraded data from the reference

if isempty(Dataset.srf_filename)
    R = zeros(1, Dataset.n_band);
    n_ov_band = length(Dataset.sp_overlap);
    R(Dataset.sp_overlap) = ones(1,n_ov_band)/n_ov_band;
else
    sp_data = double(importdata([data_root_path '/SRF/' Dataset.srf_filename]));
    [~, valid_bands] = intersect(sp_data(:,1), Dataset.sp_range);
    no_wa = length(valid_bands);
    % Spline interpolation
    xx  = linspace(1, no_wa, Dataset.n_band);
    x = 1:no_wa;

    if strcmp(Dataset.case_type, 'HS-PAN')
        n_pan_band = 1;
        R = spline(x, sp_data(valid_bands, 2), xx);
    elseif strcmp(Dataset.case_type, 'HS-MS')
        n_pan_band = size(sp_data,2)-2;
        if n_pan_band <= 0
            error('invalid sp data. Not enough data columns!');
        end
        R = zeros(n_pan_band, Dataset.n_band);
        for i = 1:n_pan_band
            R(i,:) = spline(x, sp_data(valid_bands,i+2), xx);
        end
    end
    if isempty(valid_bands)
        warning('no overlap bands.');
        Dataset.sp_overlap = 1:Dataset.n_band; 
    else
        Dataset.sp_overlap = valid_bands;
    end
    Dataset.sp_overlap = repmat(Dataset.sp_overlap,[n_pan_band 1]);

    R = bsxfun(@rdivide, R, sum(R,2));
end
% SRF used for spectral degradation, i.e. true SRF.
Dataset.R = R;

% Spectral degradation
Dataset.PAN = mat2im(R*im2mat(Dataset.REF), size(Dataset.REF,1));

