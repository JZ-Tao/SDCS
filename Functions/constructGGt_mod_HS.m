function GGt = constructGGt_mod_HS(h, K, rows, cols, band)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A modified code by Jingzhe Tao from:
% Eigen-decomposition for super-resolution
% Stanley Chan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hth = conv2(h,rot90(h,2));
%%%%%%%%%%%%%%%%%%
% Tap-23
% Other interpolating convolution kernels can be used if necessary.
L = 44;
g = K.*fir1(L,1./K);
g2 = conv2(g',g);
hth = conv2(g2,h);
% hth = conv2(h,h);
%%%%%%%%%%%%%%%%%%
yc = ceil(size(hth,1)/2);  % mark the center coordinate
xc = ceil(size(hth,2)/2);

L = floor(size(hth,1)/K);  % width of the new filter 
                           %  = (1/k) with of the original filter
                           
g = zeros(L,L);            % initialize new filter
for i=-floor(L/2):floor(L/2)
    for j=-floor(L/2):floor(L/2)
        g(i+floor(L/2)+1,j+floor(L/2)+1) = hth(yc+K*i, xc+K*j);
    end
end

GGt = abs(fft2(g,rows/K,cols/K));
GGt = repmat(GGt, [1 1 band]);