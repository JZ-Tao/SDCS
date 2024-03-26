function MTF_info = SensorMTFInfo(ratio, sensor)

MTF_info.start_pos = [1,1]; % The starting point of downsampling

switch sensor

    case 'none'
        GNyq = 0.3;
        N = 31;
        fcut = 1/ratio;

        alpha = sqrt((N*(fcut/2))^2/(-2*log(GNyq)));
        H = fspecial('gaussian', N, alpha);
        Hd = H./max(H(:));
        h = fwind1(Hd,kaiser(N));
        MTF_info.BluKer = real(h);
        MTF_info.start_pos = [floor(ratio/2)+1, floor(ratio/2)+1]; % The starting point of downsampling     

    case 'GaussKernel'
        size_kernel=[9 9];
        sig = (1/(2*(2.7725887)/ratio^2))^0.5;
        MTF_info.BluKer = fspecial('gaussian',[size_kernel(1) size_kernel(2)],sig);
end