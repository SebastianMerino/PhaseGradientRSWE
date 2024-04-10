% testing modyfied system
clear, clc
load phasemap.mat

% file_dir = 'P:\rswe\dataold\Data800Hz-10000ondas\R-FIELD_inc_1.mat';
file_dir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Elastrography\reverberant\800.mat'];
w_kernel = [15 15];
constant = 0.3;
stride = round(w_kernel)/3;
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
k2_med = medfilt2(phase_grad_2,[med_wind med_wind],'symmetric');
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
imagesc(x*cm, z*cm, sws_matrix, sws_range);
colormap(gca, turbo); % Apply jet colormap to the current axes
colorbar;
axis image;
title('SWS')
xlabel('Lateral'), ylabel('Axial'),

%% Denoiding gradient squared
M = size_out(1);
N = size_out(2);
mu = 10.^(8); % mu=0.5
tol = 1e-5;
mask = ones(size(bx_large));
minimask = ones(M*N,1);

tic
denoised_grad = IRLS_TV(phase_grad_2(:),speye(M*N),mu,M,N,tol,mask,minimask);
denoised_grad = reshape(denoised_grad, size_out);
toc
k = sqrt(denoised_grad);
% k = denoised_grad;
% --------------------
sws_denoised = (2*pi*f_v)./k;   



figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, sws_matrix, sws_range);
colormap(gca, turbo); % Apply jet colormap to the current axes
colorbar;
axis image;
title('SWS')
xlabel('Lateral'), ylabel('Axial'),
nexttile,
imagesc(x*cm, z*cm, sws_denoised, sws_range);
colormap(gca, turbo); % Apply jet colormap to the current axes
colorbar;
axis image;
title('SWS')
xlabel('Lateral'), ylabel('Axial'),


%% Denoising gradient
mu = 10.^(4.5); % mu=0.5
tic
denoised_grad = IRLS_TV(sqrt(phase_grad_2(:)),speye(M*N),mu,M,N,tol,mask,minimask);
denoised_grad = reshape(denoised_grad, size_out);
toc
k = denoised_grad;
% --------------------
sws_denoised = (2*pi*f_v)./k;   



figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, sws_matrix, sws_range);
colormap(gca, turbo); % Apply jet colormap to the current axes
colorbar;
axis image;
title('SWS')
xlabel('Lateral'), ylabel('Axial'),
nexttile,
imagesc(x*cm, z*cm, sws_denoised, sws_range);
colormap(gca, turbo); % Apply jet colormap to the current axes
colorbar;
axis image;
title('SWS')
xlabel('Lateral'), ylabel('Axial'),

%% Denoising sws
sws_noisy = (2*pi*f_v)./sqrt(phase_grad_2);   
mu = 10.^(1); % mu=0.5
tic
sws_denoised = IRLS_TV(sws_noisy(:),speye(M*N),mu,M,N,tol,mask,minimask);
sws_denoised = reshape(sws_denoised, size_out);
toc

figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, sws_matrix, sws_range);
colormap(gca, turbo); % Apply jet colormap to the current axes
colorbar;
axis image;
title('SWS')
xlabel('Lateral'), ylabel('Axial'),
nexttile,
imagesc(x*cm, z*cm, sws_denoised, sws_range);
colormap(gca, turbo); % Apply jet colormap to the current axes
colorbar;
axis image;
title('SWS')
xlabel('Lateral'), ylabel('Axial'),