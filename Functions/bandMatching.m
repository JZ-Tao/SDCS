function D_MS = bandMatching(HS, I_MS, ratio, Opts)
% I_MS forms a band match with respect to HS in any of the following ways: 
%  -- BM_org: band duplication (HS-PAN) or band selection (other).
%  -- BM_FG: Gain factor estimated using full-scale regression (from TIP
%  2018 GLP-REG method).
%  -- BM_MIX: Consensus enhancement (proposed).
%  -- Other eligible sharpening methods.
method = Opts.BM.method;
MTF_info = Opts.MTF_info;
if size(I_MS, 3) == 1
    fusion_type = 'HS-PAN';
else
    fusion_type = 'HS-MS';  % HS-MS or HS-HS
end
n_hs_band = size(HS, 3);
if size(HS,1) == size(I_MS,1)
    I_HS = HS;
    start_pos = MTF_info.start_pos;
  
    I_HS_LR = I_HS(start_pos(1):ratio:end, start_pos(2):ratio:end,:);
else
    I_HS_LR = HS;
    I_HS = interpWrapper(I_HS_LR, ratio,MTF_info);
end

basic_method = {'BM_org', 'BM_HM', 'BM_G', 'BM_FG', 'BM_SFG'};
if ismember(method, basic_method) % strcmp(method, 'BM_org') || strcmp(method, 'BM_HM')  || strcmp(method, 'BM_G')  || strcmp(method, 'BM_FG')
    if strcmp(fusion_type, 'HS-PAN')
        D_MS = repmat(I_MS, [1,1,n_hs_band]);
    else
        
        I_MS_LR = conv_downsample_MTF(I_MS, ratio, MTF_info);
        n_ms_band = size(I_MS_LR,3);
        A = zeros(size(I_MS, 3),n_hs_band);
        for i = 1:n_ms_band
            ms = I_MS_LR(:,:,i);
            for j = 1:n_hs_band
                hs = I_HS_LR(:,:,j);
                cc = corrcoef(ms(:),hs(:));
                A(i,j) = cc(1,2);
            end
        end
        [~,indices] = max(A,[],1);
        D_MS = zeros(size(I_HS));
        for j = 1:n_ms_band
            idx_tmp = find(indices == j);
            if ~isempty(idx_tmp) 
                D_MS(:,:,idx_tmp) = repmat(I_MS(:,:,j), [1 1 length(idx_tmp)]);
            end
        end
    end
    if strcmp(method, 'BM_org')
        return
    else
        D_MS_LR = conv_downsample_MTF(D_MS, ratio, MTF_info);
        D_MS_LP = interpWrapper(D_MS_LR, ratio, MTF_info);
        switch(method)
            case 'BM_HM'
                hm_mode = 1;
                D_MS = PanHistEqualization(I_HS, D_MS, D_MS_LP, hm_mode);
            case 'BM_G'
               for i = 1:n_hs_band
                   M_i = I_HS(:,:,i);
                   PL_i = D_MS_LP(:,:,i);
                   covMP = cov(M_i(:), PL_i(:));
                   G_i = covMP(1,2)/var(PL_i(:));
                   D_MS(:,:,i) = G_i.*D_MS(:,:,i);
               end
            case 'BM_FG'
               for i = 1:n_hs_band
                    M_i = I_HS(:,:,i);
                    PL_i = D_MS_LP(:,:,i);
                    P_i = D_MS(:,:,i);
                    CMSPAN = cov(M_i(:), P_i(:));    
                    CPANPANLR = cov(PL_i(:), P_i(:));
                    G_i = CMSPAN(1,2)./CPANPANLR(1,2);
                    D_MS(:,:,i) = G_i.*P_i;
               end
            case 'BM_SFG'
                n_cluster = 20;
                Im2seg = D_MS;
                if size(I_HS,3) ~= 1
                    [E,~] = svds(im2mat(Im2seg),1); Im = mat2im(E'*im2mat(Im2seg), size(Im2seg,1));
                else
                    Im = Im2seg;
                end
                Im = mapminmax(Im,0,1);
                labels = mex_ers(Im,n_cluster)+1;
                if size(Im,1) == size(I_HS,1)
                    C_map = reshape(labels,[size(Im,1) size(Im,2)]);
                    C_map_U = imresize(C_map, ratio, 'nearest');
                else
                    C_map_U = reshape(labels,[size(Im,1) size(Im,2)]);
                    C_map = imresize(C_map_U, 1/ratio, 'nearest');
                end
                Class = unique(C_map(:));
                n_C = size(Class,1);
                for i = 1:n_hs_band
                    M_i = I_HS(:,:,i);
                    PL_i = D_MS_LP(:,:,i);
                    P_i = D_MS(:,:,i);
                    cG_i = zeros(size(D_MS,1), size(D_MS,2));

                    for ii = 1 : n_C
                   %x_s = x_s + Trans(T{ii}, x.*(C_map == ii));
                        idx = C_map == Class(ii);
                        CM_i = M_i(idx);
                        CPL_i = PL_i(idx);
                        CP_i = P_i(idx);
                        CCMSPAN = cov(CM_i(:), CP_i(:));    
                        CCPANPANLR = cov(CPL_i(:), CP_i(:));
                        cG_i = cG_i + idx.*(CCMSPAN(1,2)./CCPANPANLR(1,2));
                    end
                    
                    CMSPAN = cov(M_i(:), P_i(:));    
                    CPANPANLR = cov(PL_i(:), P_i(:));
                    gG_i = CMSPAN(1,2)./CPANPANLR(1,2);
                    
                    
                    h = fspecial('gaussian',[5,5], 0.5);
                    cG_i = imfilter(cG_i,h, 'replicate');
                    tau = 0.5;
                    G_i = (1-tau)*gG_i+tau*cG_i;
                    D_MS(:,:,i) = G_i.*P_i;
                end       
        end
        return;
    end
end
if strcmp(method, 'BM_REG')
    I_MS_LR = conv_downsample_MTF(I_MS, ratio, MTF_info);
    v_hs = im2mat(I_HS_LR)';
    v_ms = im2mat(I_MS_LR)';
    v_MS = im2mat(I_MS)';
    A = [v_ms ones(size(v_hs,1),1)];    
    X = (A'*A)\A'*v_hs;
    A_HR = [v_MS ones(size(v_MS,1),1)];
    D_MS = mat2im((A_HR*X)', size(I_MS,1));
end
if strcmp(method, 'BM_GSA')
    D_MS = GSA_wrapper(I_HS_LR,I_MS,ratio, MTF_info);
    return;
end
if strcmp(method, 'BM_GLP_HS')
    mode = 2;
    D_MS = MTF_GLP_wrapper(I_HS_LR,I_MS,ratio,mode,MTF_info,Opts.L);
    return;
end
if strcmp(method, 'BM_CNMF_n')
    % Add noise
    % Adding Gaussian noise seems to help improve the quality, possibly
    % because the process allows the data to better conform to the method's
    % imaging model assumption that Y=MX+N, with N being the Gaussian noise.
    SNRh = 30; % SNR (in dB) for the hyperspectral image
    SNRm = 40; % SNR (in dB) for the multispectral/panchromatic image
    sigmah = sqrt(sum(I_HS_LR(:).^2)/(10^(SNRh/10))/numel(I_HS_LR));
    I_HS_LR_n = I_HS_LR + sigmah*randn(size(I_HS_LR));
    sigmam = sqrt(sum(I_MS(:).^2)/(10^(SNRm/10))/numel(I_MS));
    I_MS_n = I_MS + sigmam*randn(size(I_MS));
    % The number of end-members is limited to 10 to improve computational
    % efficiency.
    D_MS = CNMF_fusion_mod(I_HS_LR_n,I_MS_n,MTF_info,10);
    return;
end
if strcmp(method, 'BM_CNMF')
    D_MS = CNMF_fusion_mod(I_HS_LR,I_MS,MTF_info);
    return;
end
if strcmp(method, 'BM_CNMF_Sub')
    if size(I_MS, 3)>10
        D_MS = CNMF_fusion_mod(I_HS_LR,reduceBandsByAverage(I_MS, 10),MTF_info);
    else
        D_MS = CNMF_fusion_mod(I_HS_LR,I_MS,MTF_info);
    end
    return;
end
if strcmp(method, 'BM_FUSE')
    [R,~] = estR_new(I_HS_LR,conv_downsample_MTF(I_MS, ratio, MTF_info));
    R = R(:,1:end-1);
    D_MS = BayesianFusion(I_HS_LR,I_MS,R, MTF_info.BluKer,ratio,'Gaussian',MTF_info.start_pos);
end

if strcmp(method, 'BM_MIX') || strcmp(method, 'BM_MIX_REG') || strcmp(method, 'BM_MIX_CE')
    % BM_MIX: Consensus enhancement
    methods = {'BM_CNMF_n','BM_GLP_HS', 'BM_FG'};
    n_method = length(methods);
    Opts_sub = Opts;

    D_MSi = zeros(size(I_MS,1),size(I_MS,2),size(I_HS_LR,3),n_method);
    for i = 1:n_method
        Opts_sub.BM.method = methods{i};
        D_MSi(:,:,:,i) = bandMatching(I_HS_LR, I_MS, ratio, Opts_sub);
    end
    D_MS = zeros(size(I_HS));
    if strcmp(method, 'BM_MIX_REG')
        Opts_sub.BM.method = 'BM_REG';
        for i = 1:size(I_HS_LR,3)
            Bi = squeeze(D_MSi(:,:,i,:));     
            D_MS(:,:,i) = bandMatching(I_HS_LR(:,:,i), Bi, ratio, Opts_sub);
        end
    else
        for i = 1:size(I_HS_LR,3)
            Bi = squeeze(D_MSi(:,:,i,:));
            [E,~] = svds(im2mat(Bi),1); E = E./sum(E);
            Bi_sub = (E'*im2mat(Bi));
            D_MS(:,:,i) = mat2im(Bi_sub, size(Bi,1));
        end
        if strcmp(method, 'BM_MIX_CE')
            [rows, cols, bands] = size(I_HS);
            max_itr = 5;
            n_method = length(methods);
            gamma = 0.5;
            rho = 1;
            
            w = D_MSi;%(:,:,:,1:end);
            clear D_MSi;
            zhat = D_MS;
            Fv    = zeros(rows,cols,bands, n_method);
            for itr = 1:max_itr
                %==== Update v ====%
                Gw = zhat;
                Gw = repmat(Gw,[1 1 1 n_method]);
                v  = 2*Gw - w;

                %==== Update w ====%
                Fv(:,:,:,1) = (D_MS + rho*v(:,:,:,1))/(1+rho);
                for i=1:n_method-1
                    Opts_sub.BM.method = methods{i};
                    Fv(:,:,:,i+1) = bandMatching(v(:,:,:,i+1), I_MS, ratio, Opts_sub);
                end
                w = (1-gamma)*w + gamma*(2*Fv - v);

                %==== Compute zhat ====%
                zhat = sum(repmat(reshape(E',[1 1 1 n_method]),[rows,cols,bands,1]).*w,4);
            end
            
        end
    end
    return;
end
