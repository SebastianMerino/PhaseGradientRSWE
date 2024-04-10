clear, clc
load phasemap.mat

% file_dir = 'P:\rswe\dataold\Data800Hz-10000ondas\R-FIELD_inc_1.mat';
file_dir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Elastrography\reverberant\800.mat'];
w_kernel = [15 15];
constant = 0.33;
stride = round(w_kernel)/5;

% plotting const
cm = 1e2;
sws_range = [2,5];

data = load(file_dir);
u = data.pv_complexZ;
dinf.dx = data.x(2) - data.x(1);
dinf.dz = data.z(2) - data.z(1);
f_v = data.freq;
x = data.x; z = data.z;
og_size = size(u);

%% Observing data
figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, real(u(:,:,1)));
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('PV')
xlabel('Lateral'), ylabel('Axial'),

nexttile,
imagesc(x*cm, z*cm, angle(u));
colormap(gca, phasemap); % Apply jet colormap to the current axes
colorbar;
axis image;
title('Phase')
xlabel('Lateral'), ylabel('Axial'),

extended_u = padarray(u,[(w_kernel(1)-1)/2 (w_kernel(2)-1)/2],'symmetric');

% %% 2D filtering -> NOT EFECTIVE
% dx = x(2)-x(1);
% cmin = 0.5; cmax = 5;
% kmin = 2*pi*f_v/cmax;
% kmax = 2*pi*f_v/cmin;
% knyq = 2*pi/dx /2;
% 
% NFFT = 21;
% [f1,f2] = freqspace(NFFT,'meshgrid');
% r = sqrt(f1.^2 + f2.^2);
% Hd = ones(NFFT); 
% Hd((r<kmin/knyq)|(r>kmax/knyq)) = 0;
% 
% win = fspecial('gaussian',NFFT,3);
% win = win ./ max(win(:));
% 
% h = fwind2(Hd,win);
% % freqz2(h)
% 
% filt_u = filter2(h,extended_u);
% filt_u = filt_u(8:end-7,8:end-7);
% 
% figure('Units','centimeters', 'Position',[5 5 20 10]),
% tiledlayout(1,2)
% nexttile,
% imagesc(x*cm, z*cm, real(filt_u(:,:,1)));
% colormap(gca, parula); % Apply jet colormap to the current axes
% colorbar;
% axis image;
% title('PV')
% xlabel('Lateral'), ylabel('Axial'),
% 
% nexttile,
% imagesc(x*cm, z*cm, angle(filt_u));
% colormap(gca, phasemap); % Apply jet colormap to the current axes
% colorbar;
% axis image;
% title('Phase')
% xlabel('Lateral'), ylabel('Axial'),

%% ================ SYSTEM WITH MEAN FILTER ============
% testing modyfied system
clear, clc
load phasemap.mat

% file_dir = 'P:\rswe\dataold\Data800Hz-10000ondas\R-FIELD_inc_1.mat';
file_dir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Elastrography\reverberant\800.mat'];
w_kernel = [5 5];
constant = 0.3;
stride = round(w_kernel/5);
%the median window contains at least a wavelenght

% plotting const
cm = 1e2;
sws_range = [2,5];

data = load(file_dir);
u = data.pv_complexZ;
dinf.dx = data.x(2) - data.x(1);
dinf.dz = data.z(2) - data.z(1);
f_v = data.freq;
x = data.x; z = data.z;
og_size = size(u);
med_wind = floor (2.5/f_v/dinf.dx/stride(1))*2+1; 

%% Observing data
figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, real(u(:,:,1)));
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('PV')
xlabel('Lateral'), ylabel('Axial'),

nexttile,
imagesc(x*cm, z*cm, angle(u));
colormap(gca, phasemap); % Apply jet colormap to the current axes
colorbar;
axis image;
title('Phase')
xlabel('Lateral'), ylabel('Axial'),


%% Getting system
extended_u = padarray(u,[(w_kernel(1)-1)/2 (w_kernel(2)-1)/2],'symmetric');
tic
[Ax_large, Az_large, bx_large, bz_large, size_out] = getmat_pg_v3(extended_u,...
    w_kernel, dinf, og_size, stride);
toc

%% Direct inversion
tic
results_x = minres(Ax_large'*Ax_large,Ax_large'*bx_large);
results_z = minres(Az_large'*Az_large,Az_large'*bz_large);
toc
results_x = reshape(results_x,[size_out, 2]);
results_z = reshape(results_z,[size_out, 2]);


%% Unwrap x system
figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, (results_x(:,:,1)), [-2000 2000]);
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('k_x')
xlabel('Lateral'), ylabel('Axial'),

nexttile,
imagesc(x*cm, z*cm, (results_x(:,:,2)));
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('c')
xlabel('Lateral'), ylabel('Axial'),


%% Unwrap z system
figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, (results_z(:,:,1)), [-2000 2000]);
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('k_z')
xlabel('Lateral'), ylabel('Axial'),

nexttile,
imagesc(x*cm, z*cm, (results_z(:,:,2)));
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('c')
xlabel('Lateral'), ylabel('Axial'),

%% SWS
grad_x = results_x(:,:,1);
grad_z = results_z(:,:,1);
phase_grad_2 = (grad_x.^2 + grad_z.^2)/constant;
% ----- MedFilt  ----
% k2_med = medfilt2(phase_grad_2,[med_wind med_wind],'symmetric');
% h = fspecial('average', [60 60]);
% k2_med = filter2(h,phase_grad_2);
k2_med = imgaussfilt(phase_grad_2,20);
k = sqrt(k2_med);
% --------------------
sws_matrix = (2*pi*f_v)./k;   

figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, real(u(:,:,1)));
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('k_x')
xlabel('Lateral'), ylabel('Axial'),

nexttile,
imagesc(x*cm, z*cm, sws_matrix);
colormap(gca, turbo); % Apply jet colormap to the current axes
colorbar;
axis image;
title('SWS')
xlabel('Lateral'), ylabel('Axial'),