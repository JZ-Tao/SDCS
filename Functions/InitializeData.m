function [Datasets] = InitializeData(data_root_path, on_test_datasets)
%--------------------------------------------------------------------------
% INPUTS
%   data_root_path: path of dataset.
%   on_test_datasets: dataset tag list.
% OUTPUTS
%   dataset: the sturcture array containing the information of initiali-
%   zed datasets. A record of dataset includes the following fields:
%       -- tag: A keyword for corrent dataset record, eg. 'indian'.
%       -- filename: The name of the data file, eg. 'Indian_pines.mat'.
%       -- n_row: Rows (ie. width) of the selected data, eg. 256.
%       -- n_column: Column (ie. height) of the selected data, eg. 256.
%       -- n_band: Number of selected bands, eg. 192.
%       -- band_set: The indices of the selected bands, eg. 1:192.
%       -- ratio: The GSD ratio between LR and HR data, eg. 4.
%       -- REF: The reference data (ie. the original HR HS data) with a 
%       dynamic range scaled to [0,1].
%       -- max_v / min_v: For storing the original dynamic range of data.
%       -- HS: the simulated LR data by spatial degradation from REF.
%       -- KerBlu: The blur kernel of spatial downsampling.
%       -- PAN: The simulated HR Pan data by spectral degradation from REF.
% By Jingzhe Tao, Apl. 2018

%--------------------------------------------------------------------------
global_sensor = 'none'; % 'none' | 'GaussKernel'
n_data = 0;

%% Moffett Field

if numel(find(ismember(on_test_datasets, 'Moffett_field')))
    n_data = n_data+1;
    Datasets{n_data}.tag = 'Moffett Field HS-PAN';
    Datasets{n_data}.filename = 'Moffett_field.mat';
    Datasets{n_data}.band_set = (1:176)';
    Datasets{n_data}.sp_range = 400:2500;
    Datasets{n_data}.srf_filename = '';
    Datasets{n_data}.sp_overlap = 1:41;
    Datasets{n_data}.ratio = 6;
    Datasets{n_data}.scaling_type = 1;
    Datasets{n_data}.sensor = global_sensor;
    init_func{n_data} = @(x) generalDataInit(data_root_path, x);
    Datasets{n_data}.case_type = 'HS-PAN';
    Datasets{n_data}.L_org = 13;
    Datasets{n_data}.cropped = 1;
    Datasets{n_data}.cropped_size = [390, 180];
    Datasets{n_data}.cropped_pos = [1, 1];
    Datasets{n_data}.color_band_idx = [29, 19, 9];
    Datasets{n_data}.PAN_color_band_idx = 1;
end

%% Hydice

if numel(find(ismember(on_test_datasets, 'Washington_MS'))) % 1280x307
    n_data = n_data+1;
    Datasets{n_data}.tag = 'Washington DC HS-MS';
    Datasets{n_data}.filename = 'dc.mat';
    Datasets{n_data}.n_band = 191;
    Datasets{n_data}.band_set = (1:191)';
    Datasets{n_data}.sp_range = 400:2500;
    Datasets{n_data}.srf_filename = 'quickbird.mat';
    Datasets{n_data}.ratio = 5;
    Datasets{n_data}.scaling_type = 1;
    Datasets{n_data}.sensor = global_sensor;
    init_func{n_data} = @(x) generalDataInit(data_root_path, x);
    Datasets{n_data}.case_type = 'HS-MS';   
    Datasets{n_data}.cropped = 1;
    Datasets{n_data}.cropped_size = [480, 300];
    Datasets{n_data}.cropped_pos = [780, 5];
    Datasets{n_data}.color_band_idx = [50, 36, 15]; %[60, 27, 17];
    Datasets{n_data}.PAN_color_band_idx = [3 2 1];
end

if numel(find(ismember(on_test_datasets, 'Washington_HS'))) % 1280x307
    n_data = n_data+1;
    Datasets{n_data}.tag = 'Washington DC HS-HS';
    Datasets{n_data}.filename = 'dc.mat'; 
    Datasets{n_data}.n_band = 191;
    Datasets{n_data}.band_set = {(1:78)'; (79:191)'};
    Datasets{n_data}.sp_range = 400:2500;
    Datasets{n_data}.srf_filename = 'quickbird.mat';
    Datasets{n_data}.ratio = 4;
    Datasets{n_data}.scaling_type = 1;
    Datasets{n_data}.sensor = global_sensor;
    init_func{n_data} = @(x) HS_HS_DataInit(data_root_path, x);
    Datasets{n_data}.case_type = 'HS-HS'; 
    Datasets{n_data}.cropped = 1;
    Datasets{n_data}.cropped_size = [480, 300];
    Datasets{n_data}.cropped_pos = [780, 5];
    Datasets{n_data}.color_band_idx = [50, 36, 15]; %[60, 27, 17];
    Datasets{n_data}.PAN_color_band_idx = [50, 36, 15];
end

%% Rosis Pavia

if numel(find(ismember(on_test_datasets, 'Pavia'))) % 610x340
    n_data = n_data+1;
    Datasets{n_data}.tag = 'Pavia University HS-PAN';
    Datasets{n_data}.filename = 'PaviaU.mat';
    Datasets{n_data}.n_band = 103;
    Datasets{n_data}.band_set = (1:103)';
    Datasets{n_data}.sp_range = 430:840;
    Datasets{n_data}.srf_filename = 'ikonos_spec_resp.mat';
    Datasets{n_data}.ratio = 6;
    Datasets{n_data}.scaling_type = 1;
    Datasets{n_data}.sensor = 'GaussKernel';
    init_func{n_data} = @(x) generalDataInit(data_root_path, x);
    Datasets{n_data}.case_type = 'HS-PAN';   
    Datasets{n_data}.cropped = 1;
    Datasets{n_data}.cropped_size = [600, 330]; % 320, 320
    Datasets{n_data}.cropped_pos = [5, 10]; % 50, 10
    Datasets{n_data}.color_band_idx = [60, 35, 8];
    Datasets{n_data}.PAN_color_band_idx = 1;
end


if numel(find(ismember(on_test_datasets, 'Pavia_MS'))) % 610x340
    n_data = n_data+1;
    Datasets{n_data}.tag = 'Pavia University HS-MS';
    Datasets{n_data}.filename = 'PaviaU.mat';
    Datasets{n_data}.n_band = 103;
    Datasets{n_data}.sp_range = 430:840;
    Datasets{n_data}.srf_filename = 'ikonos_spec_resp.mat';
    Datasets{n_data}.ratio = 5;
    Datasets{n_data}.scaling_type = 1;
    Datasets{n_data}.sensor = 'GaussKernel';
    init_func{n_data} = @(x) generalDataInit(data_root_path, x);
    Datasets{n_data}.case_type = 'HS-MS';
    Datasets{n_data}.cropped = 1;
    Datasets{n_data}.cropped_size = [600, 330];
    Datasets{n_data}.cropped_pos = [5, 10];
    Datasets{n_data}.color_band_idx = [60, 35, 8];
    Datasets{n_data}.PAN_color_band_idx = [3 2 1];
end

%% Chikusei
if numel(find(ismember(on_test_datasets, 'Chikusei_HS')))
    n_data = n_data+1;
    Datasets{n_data}.tag = 'Chikusei HS-HS';
    Datasets{n_data}.filename = 'Chikusei_729_512.mat';
    Datasets{n_data}.srf_filename = '';
    Datasets{n_data}.n_band = 128;
    Datasets{n_data}.band_set = {(1:68)'; (69:128)'};
    Datasets{n_data}.sp_range = 360:1020;
    Datasets{n_data}.ratio = 4;
    Datasets{n_data}.scaling_type = 1;
    Datasets{n_data}.sensor = global_sensor;
    init_func{n_data} = @(x) HS_HS_DataInit(data_root_path, x);
    Datasets{n_data}.case_type = 'HS-HS';
    Datasets{n_data}.cropped = 0;
    Datasets{n_data}.cropped_size = [320, 320];
    Datasets{n_data}.cropped_pos = [512, 729]; 
    Datasets{n_data}.color_band_idx = [56 35 17];
    Datasets{n_data}.PAN_color_band_idx = [56 35 17];
end

if n_data == 0
    error('Name of dataset with no matches.')
end
%% import HS data

for i = 1:n_data
    Datasets{i}.MTF_info = SensorMTFInfo(Datasets{i}.ratio, Datasets{i}.sensor);
    Datasets{i} = init_func{i}(Datasets{i});
end



