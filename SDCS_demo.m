% A demo for SDCS.
% Subspace Dynamic Combined Sparsity-Based Hyper-Sharpening for Diverse
% Auxiliary Images, TGRS, 2024. DOI: 10.1109/TGRS.2024.3382402
% By Jingzhe Tao, Mar. 2024
clc;
clear;
close all;
file_path = matlab.desktop.editor.getActive;
cd(fileparts(file_path.Filename));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_run_times = 1;
L = 1;
flag_cut = 1;
dim_cut = 1;
data_root_path = 'Datasets';
disp('Data initialization is in progress...');
on_test_datasets = {'Chikusei_HS'}; % 'Moffett_field' | 'Chikusei_HS'
Datasets = InitializeData(data_root_path, on_test_datasets);
Dataset = Datasets{1};
disp('Data initialization is complete.');
I_REF = Dataset.REF;
I_HS_LR = Dataset.HS;
I_PAN = Dataset.PAN; % Auxiliary image
ratio = Dataset.ratio;
MTF_info = Dataset.MTF_info;
tag = Dataset.tag;
clear Dataset
k = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
test_cases = {
    % Each line corresponds to a set of SDCS parameter combinations.
    % From top to bottom are SDCS (the proposed), SDCS-TV, SDCS-NL, SDCS-TVNL.
    {{'ADMM'},{'cascade'},{'NL_TH'},{'BM_MIX'},{'SVD'},{'lambda_1e-3'},{'n_sub10'},{'iter_tol_5e-3'}},...  
%     {{'Fista'},{'TV_only'},{'BM_FG'},{'SVD'},{'lambda_1e-3'},{'n_sub10'},{'iter_tol_5e-3'}},...
%     {{'ADMM'},{'NL_only'},{'BM_MIX'},{'SVD'},{'lambda_1e-3'},{'n_sub10'},{'iter_tol_5e-3'}},...
%     {{'ADMM'},{'TV_NL_92'},{'BM_MIX'},{'SVD'},{'lambda_1e-3'},{'n_sub10'},{'iter_tol_5e-3'}},...
};
n_test = length(test_cases);
for t = 1:n_test
    [Opts, arg_str] = initSDCSOptByTestCase(test_cases{t});
    Opts.tag = tag;
    Opts.MTF_info = MTF_info;
    Opts.L = L;
    Opts.arg_str = arg_str;
    Opts.print = true;
    k = k+1;

    Method_{k} = ['SDCS_' arg_str];
    disp(['Computing ' Method_{k}]);
    t2 = tic;
    for i = 1:n_run_times
        I_{k} = SDCS_Wrapper(I_HS_LR,I_PAN, ratio, Opts);
    end
    Time_{k} = toc(t2)/n_run_times;
    fprintf('Elaboration time %s: %.4f [sec]\n',Method_{k},Time_{k});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QI_ = cell(1, k);
for i = 1:k
    disp(['Quality indices of method ' Method_{i} ' :']);
    QI_{i} = QualityIndices(I_{i}, I_REF, ratio, flag_cut, dim_cut);
end
