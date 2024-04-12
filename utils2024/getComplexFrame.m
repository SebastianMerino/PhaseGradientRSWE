function [uFrame,uFilt] = getComplexFrame(u,fv,cMin,dinf)
% ====================================================================== %
% Inputs:  
%          u           : 3D matriz (x,z,t)
%          fv         : vibration frequency
%          cMin         : minimum SWS to filter
%          dinf        : structure with spatial and temporal resolutions
%
% Outputs: 
%          uFrame    : complex PV frame
%          uFilt     : filtered PV 
% ====================================================================== %
kmax = 2*pi*fv/cMin;
uFilt = zeros(size(u));
for it = 1:size(u,3)
    uFilt(:,:,it) = filt2D(u(:,:,it),dinf.dx,dinf.dz,kmax);
end
% figure,
% for it = 1:20
%     imagesc(x,z,uFilt(:,:,it), std(uFilt(:))*[-1 1])
%     axis image
%     pause(0.1)
% end
fs = dinf.PRFe;
t = (0:size(u,3)-1)/fs;
t = reshape(t,[1,1,length(t)]);
uFrame = sum(uFilt.*exp(-1j*2*pi*fv*t),3);
end

% ====================================================================== %
% ====================================================================== %

function uFilt  = filt2D(u,dx,dz,kmax)
% ====================================================================== %
% Functions that low-pass filters spatial frequencies until kmax
% ====================================================================== %
[M,N] = size(u);
[fx,fz] = freqspace([M,N],'meshgrid'); % -1 to 1
kx = fx*2*pi/dx/2;
kz = fz*2*pi/dz/2;

% Freq filter
n = 4;
kmod = sqrt(kx.^2 + kz.^2);
H = 1./(1 + (kmod/kmax).^(2*n));
uFT = fftshift(fft2(u));
uFilt = ifft2(ifftshift(uFT.*H));
% figure, imagesc(kx(1,:),kz(:,1),H)
% figure, imagesc(kx(1,:),kz(:,1),db(uFT))
end
