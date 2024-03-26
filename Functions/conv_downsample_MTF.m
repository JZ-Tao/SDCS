function Y = conv_downsample_MTF(X, ratio, MTF_info)

    start_pos = MTF_info.start_pos;
    BluKer = MTF_info.BluKer;
    
    
    X_LP=imfilter(X, BluKer, 'circular');
    
     Y=X_LP(start_pos(1):ratio:end, start_pos(2):ratio:end,:);
 
