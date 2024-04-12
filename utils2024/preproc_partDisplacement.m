function [uFilt] = preproc_partDisplacement(u,fv,cMin,dinf,opts)
% ====================================================================== %
% Inputs:  
%          u           : 3D matriz (x,z,t)
%          fv         : vibration frequency
%          cMin         : minimum SWS to filter
%          dinf        : structure with spatial and temporal resolutions
%          opts      :
%               opts.lpf : 2D LOW PASS FILTER 
%               opts.lpf : 2D PASS BAND FILTER (especify opts.cMaX)
% Outputs: 
%          uFilt     : filtered PV 
% ====================================================================== %
    
    % PEAK IDENTIFICATIOn
    Ufft = (fft(u, [], 3)); % do fft pixel-wise 3rd dimension time
    powerSpectrum = abs(Ufft).^2; % power spectrum 

    [~, maxFreqIndex] = max(powerSpectrum, [], 3); % identify peaks

    peak_Complex = zeros(size(u,1), size(u,2));
    peak_phase = peak_Complex;
    peak_magnitude = peak_Complex;

    % Extract the magnitude and phase for the frequency with the maximum amplitude
    for i = 1:size(Ufft,1)
        for j = 1:size(Ufft,2)
            maxIndex = maxFreqIndex(i,j);
            peak_magnitude(i,j) = abs(  Ufft(i,j,maxIndex) );
            peak_phase(i,j) = angle(  Ufft(i,j,maxIndex) );
            peak_Complex(i,j) = peak_magnitude(i,j) * exp(1i * peak_phase(i,j));
        end
    end
        
    % 2D filtering
    kmax = 2*pi*fv/cMin;
    if strcmp(opts.typeFilter, 'lpf')
        uFilt  = filt2D_lpf(peak_Complex,dinf.dx,dinf.dz,kmax);
    end
    if strcmp(opts.typeFilter, 'bpf')
        kmin = 2*pi*fv/opts.cMax;
        uFilt  = filt2D_bpf(peak_Complex,dinf.dx,dinf.dz,kmin,kmax);
    end
end

% ====================================================================== %
% ====================================================================== %

function uFilt  = filt2D_lpf(u,dx,dz,kmax)
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

function uFilt = filt2D_bpf(u, dx, dz, k_min, k_max)
% ====================================================================== %
% Function that band-pass filters spatial frequencies between k_min and k_max
% ====================================================================== %
[M, N] = size(u);
[fx, fz] = freqspace([M, N], 'meshgrid');  % -1 to 1
kx = fx * 2 * pi / dx / 2;
kz = fz * 2 * pi / dz / 2;

% Frequency filter
n = 4;
kmod = sqrt(kx.^2 + kz.^2);
% Band-pass filter: 1 within the band, 0 outside
H = 1 ./ (1 + (kmod / k_max).^(2 * n)) - 1 ./ (1 + (kmod / k_min).^(2 * n));
uFT = fftshift(fft2(u));
uFilt = ifft2(ifftshift(uFT .* H));

% Optional: Uncomment to visualize the filter response and the magnitude spectrum
% figure, imagesc(kx(1, :), kz(:, 1), H), title('Filter Response'), xlabel('kx'), ylabel('kz');
% figure, imagesc(kx(1, :), kz(:, 1), db(abs(uFT))), title('Magnitude Spectrum'), xlabel('kx'), ylabel('kz');
end