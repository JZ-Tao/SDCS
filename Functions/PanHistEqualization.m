function P = PanHistEqualization(M, P, PL, mode)
if mode == 0
    return;
end
n_band = size(M, 3);
% global, mixed Res (MTF-GLP etc. default)
if mode == 1
    for i = 1:n_band
        P(:,:,i) = (P(:,:,i) - mean2(P(:,:,i))).*(std2(M(:,:,i))./std2(PL(:,:,i))) + mean2(M(:,:,i));
    end
    return;
end
% local HM
if mode == 2
    w = 3;
    for i = 1:n_band
        [U_M, SD_M] = LocalStatics(M(:,:,i), w); 
        [U_P, SD_P] = LocalStatics(P(:,:,i), w); 
        [U_PL, SD_PL] = LocalStatics(PL(:,:,i), w); 
        P(:,:,i) = (P(:,:,i) - U_P).*(SD_M./SD_PL) + U_M;
        %P(:,:,i) = (P(:,:,i) - mean2(P(:,:,i))).*(std2(M(:,:,i))./std2(PL(:,:,i))) + mean2(M(:,:,i));
    end
    return;
end
% Local Regress
if mode == 3
    w = 3;
    for i = 1:n_band
        [G, O] = regionRegress(M(:,:,i), PL(:,:,i), w); 
        P(:,:,i) = G.*P(:,:,i)+O;
        %P(:,:,i) = (P(:,:,i) - mean2(P(:,:,i))).*(std2(M(:,:,i))./std2(PL(:,:,i))) + mean2(M(:,:,i));
    end
    
    return;
end


