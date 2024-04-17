% SCRIP TO test in vivio data
setup,
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Elastrography\reverberant\CIRS_phantom\L7-4'];
resultsDir = fullfile(baseDir, 'results');
[~,~,~] = mkdir(resultsDir);

%% Setting parameters

% Kernel size and step
% w_kernel = [15 15];
% kernel_lsq = [11,11];
% kernel_tv = [11,11];
stride = 1;

% Plotting const
cm = 1e2;
sws_range = [0,4];

% Reg constants
tol = 1e-5;

subfolders = dir(fullfile(baseDir,'data*'));
%% Loading file
%iAcq = 3;
for iAcq = 1:length(subfolders) 
% loading data
name = [subfolders(iAcq).name,'.mat'];
data = load(fullfile(baseDir,subfolders(iAcq).name,name));
u = data.u;
dinf = data.dinf;
f_v = data.f_vib(1);
x = (0:size(u,2)-1)*dinf.dx - dinf.offset_x;
z = (0:size(u,1)-1)*dinf.dz;
Bmode = db(data.IQBmodeData);
Bmode = Bmode - max(Bmode(:));
xBm = (0:size(Bmode,2)-1)*dinf.dx - dinf.offset_x;
zBm = (0:size(Bmode,1)-1)*dinf.dz;
   
% %% Taking peak
% uFrame = zeros(size(u, [1 2]));
% mask = zeros(size(uFrame));
% for ii = 1:size(u,1)
%     for jj = 1:size(u,2)
%         signal = squeeze(u(ii,jj,:));
%         U = fft(signal);
%         [~,iMax] = max(abs(U));
%         uFrame(ii,jj) = U(iMax);
%         % mask(ii,jj) = snr(signal)>=10; 
%     end
% end
% 
% figure('Units','centimeters', 'Position',[5 5 20 10]), 
% tiledlayout(1,2, 'TileSpacing','compact', 'Padding','compact')
% t1 = nexttile;
% imagesc(xBm*cm,zBm*cm, real(uFrame), std(uFrame(:))*[-1 1])
% colormap(t1,"parula")
% axis image
% title('PV')
% 
% t1 = nexttile;
% imagesc(x*cm,z*cm, mask, [0 1])
% colormap(t1,"parula")
% axis image
% title('SNR mask')
% 
% %% 2D bandpass filter
% cMin = 0.3; cMax = 5.1; 
% kMin = 2*pi*f_v/cMax; kMax = 2*pi*f_v/cMin;
% uFilt = filt2D_bpf(uFrame, dinf.dx, dinf.dz, kMin, kMax);
% 
% load phasemap.mat
% cm = 1e2;
% figure('Units','centimeters', 'Position',[5 5 20 10]), 
% tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact')
% t1 = nexttile;
% imagesc(xBm*cm,zBm*cm, Bmode, [-50 0])
% colormap(t1,"gray")
% axis image
% title('Bmode')
% 
% t1 = nexttile;
% imagesc(x*cm,z*cm, real(uFilt), std(uFilt(:))*[-1 1])
% colormap(t1,"parula")
% axis image
% title('PV')
% t2 = nexttile;
% imagesc(x*cm,z*cm, (angle(uFilt)))
% colormap(t2,phasemap)
% axis image
% title('Phase')

%% Other function
[uFiltv2] = fun_JO_v1(u, f_v, dinf);
figure('Units','centimeters', 'Position',[5 5 20 10]), 
tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact')
t1 = nexttile;
imagesc(xBm*cm,zBm*cm, Bmode, [-50 0])
colormap(t1,"gray")
axis image
title('Bmode')

t1 = nexttile;
imagesc(x*cm,z*cm, real(uFiltv2), std(uFiltv2(:))*[-1 1])
colormap(t1,"parula")
axis image
title('PV')
t2 = nexttile;
imagesc(x*cm,z*cm, (angle(uFiltv2)))
colormap(t2,phasemap)
axis image
title('Phase')

%% Fitting
roi_size = [1.34,1.34]*1e-3; %% ROI size, depth x width
overlap = 0.8;

nz = round(roi_size(1)/dinf.dz);
nx = round(roi_size(2)/dinf.dx);
wz = round(roi_size(1)/dinf.dz*(1-overlap));
wx = round(roi_size(2)/dinf.dx*(1-overlap));

z_axis = (0:nz-1)*dinf.dz; 
x_axis = (0:nx-1)*dinf.dx; 
[X,Z] = meshgrid(x_axis,z_axis);
A = [X(:) Z(:) ones(length(x_axis)*length(z_axis),1)];

ix0 = 1:wx:size(u,2)-nx;
iz0 = 1:wz:size(u,1)-nz;
N = length(ix0);
M = length(iz0);
phase_grad = zeros(M,N);
for ii = 1:M
    for jj = 1:N
        uRoi = uFiltv2(iz0(ii):iz0(ii)+nz-1, ix0(jj):ix0(jj)+nx-1);
        % maskRoi = mask(iz0(ii):iz0(ii)+nz-1, ix0(jj):ix0(jj)+nx-1);
        % W = diag(maskRoi(:));
        phase = unwrap(angle(uRoi),[],2);
        b = phase(:);
        % results = (W*A)\(W*b);
        results = A\b;
        phase_grad(ii,jj) = sqrt(results(1).^2 + results(2).^2);
    end
end

%%
sws = 2*pi*f_v./phase_grad;
% sws_filt = imgaussfilt(sws,9);
h = fspecial('gaussian',100,7);
sws_filt = imfilter(sws,h);
figure('Units','centimeters', 'Position',[5 5 20 10]), 
tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact')
t1 = nexttile;
imagesc(x*cm,z*cm, phase_grad)
colormap(t1,"parula")
axis image
title('|\Delta\phi|')

t1 = nexttile;
imagesc(x*cm,z*cm, sws/2, sws_range)
colorbar
colormap(t1,"turbo")
axis image
title('SWS')

t1 = nexttile;
imagesc(x*cm,z*cm, sws_filt/2, sws_range)
colorbar
colormap(t1,"turbo")
axis image
title('SWS')

%% PG-LSQ + MED
constant_lsq = 0.33;
kernel_lsq = [11 11];

extended_u = padarray(uFiltv2,(kernel_lsq-1)/2,'symmetric');
og_size = size(u, [1 2]);
% Constructing matrices
tic
[Ax_large, Az_large, bx_large, bz_large, size_out] = getmat_pg_v3(...
    extended_u, kernel_lsq, dinf, og_size, stride);
toc

% Direct inversion
tic
[results_x,~] = minres(Ax_large'*Ax_large,Ax_large'*bx_large);
[results_z,~] = minres(Az_large'*Az_large,Az_large'*bz_large);
results_x = reshape(results_x,[size_out,2]);
results_z = reshape(results_z,[size_out,2]);
toc

% Post-processing
grad_x = results_x(:,:,1);
grad_z = results_z(:,:,1);
phase_grad_2 = (grad_x.^2 + grad_z.^2)/constant_lsq;
wl_samples = floor (2.5/f_v/dinf.dx/stride(1))*2+1;
% wl_samples = 1;
k2_med = medfilt2(phase_grad_2,[wl_samples wl_samples],'symmetric');
k = sqrt(k2_med);
sws_lsq = (2*pi*f_v)./k;
% 
% h = fspecial('gaussian',100,7);
% sws_filt = imfilter(sws_lsq,h);

%% PG-TV
constant_tv = 0.7;
kernel_tv = [15 15]; stride = 1;
extended_u = padarray(uFiltv2,(kernel_tv-1)/2,'symmetric');

[grad_abs, size_out] = pg_norm(extended_u, kernel_tv, dinf, og_size, stride);
grad_abs = grad_abs/sqrt(constant_tv);
M = size_out(1); N = size_out(2);
mu = 10.^4; 

tol = 1e-5;
mask = ones(size(M*N));
minimask = ones(M*N,1);

tic
denoised_grad = IRLS_TV(grad_abs(:),speye(M*N),mu,M,N,tol,mask,minimask);
denoised_grad = reshape(denoised_grad, size_out);
toc
k = denoised_grad;
sws_tv = (2*pi*f_v)./k;   


    %%
    xSWS = x(1:stride:end);
    zSWS = z(1:stride:end);
    figure('Units','centimeters', 'Position',[5 5 30 10]),
    tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact')
    nexttile,
    imagesc(x*cm, z*cm, Bmode);
    colormap(gca, gray); 
    axis image;
    title(['Bmode, freq=',num2str(f_v)])
    xlabel('Lateral [cm]'), ylabel('Axial [cm]'),

%     nexttile,
%     imagesc(xSWS*cm, zSWS*cm, sws_cf, sws_range);
%     colorbar
%     colormap(gca, turbo); 
%     axis image;
%     title('SWS PG-CF')
%     xlabel('Lateral [cm]')
%     
    nexttile,
    imagesc(xSWS*cm, zSWS*cm, sws_lsq, sws_range);
    colorbar
    colormap(gca, turbo); 
    axis image;
    title('SWS PG-LSQ')
    xlabel('Lateral [cm]')

    nexttile,
    imagesc(xSWS*cm, zSWS*cm, sws_tv, sws_range);
    colorbar
    colormap(gca, turbo); 
    axis image;
    title('SWS PG-TV')
    xlabel('Lateral [cm]')
    %%
    save_all_figures_to_directory(resultsDir,[name(1:end-4),'_fig']);
    close all,
    
    % save(fullfile(resultsDir,name),'sws_cf',"sws_lsq","sws_lsq","uFilt")
    save(fullfile(resultsDir,name),"sws_lsq","sws_lsq","uFilt")

end
%%
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

% Optional: Uncomment to visualize the filter response and the 
% magnitude spectrum
% figure, imagesc(kx(1, :), kz(:, 1), H), 
% title('Filter Response'), xlabel('kx'), ylabel('kz');
% figure, imagesc(kx(1, :), kz(:, 1), db(abs(uFT))), 
% title('Magnitude Spectrum'), xlabel('kx'), ylabel('kz');
end