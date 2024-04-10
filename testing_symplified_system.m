% testing modyfied system
clear, clc
load phasemap.mat

% file_dir = 'P:\rswe\dataold\Data800Hz-10000ondas\R-FIELD_inc_1.mat';
file_dir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Elastrography\reverberant\700.mat'];
w_kernel = [15 15];
constant = 0.3;
stride = round(w_kernel/5);
%the median window contains at least a wavelenght

% plotting const
cm = 1e2;
sws_range = [2,5];

data = load(file_dir);
u = data.pv_complexZ(:,:,1);
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
[Ax_large, Az_large, bx_large, bz_large, size_out] = getmat_pg_vSimp_nophi0(extended_u,...
    w_kernel, dinf, og_size, stride);
toc

%% Direct inversion
tic
results_x = minres(Ax_large'*Ax_large,Ax_large'*bx_large);
results_z = minres(Az_large'*Az_large,Az_large'*bz_large);
toc
results_x = reshape(results_x,size_out);
results_z = reshape(results_z,size_out);


%% both systems
figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, (results_x), [-2000 2000]);
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('k_x')
xlabel('Lateral'), ylabel('Axial'),

nexttile,
imagesc(x*cm, z*cm, (results_z), [-2000 2000]);
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('k_z')
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

%% Filtering each axis

grad_x = medfilt2(abs(results_x(:,:,1)),[med_wind med_wind],'symmetric');
grad_z = medfilt2(abs(results_z(:,:,1)),[med_wind med_wind],'symmetric');
phase_grad_2 = (grad_x.^2 + grad_z.^2)/constant;
k2_med = phase_grad_2;
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

%% NEW SYSTEM WITH ABSOLUTE VALUES
extended_u = padarray(u,[(w_kernel(1)-1)/2 (w_kernel(2)-1)/2],'symmetric');
tic
[Ax_large, Az_large, bx_large, bz_large, size_out] = getmat_pg_v4(extended_u,...
    w_kernel, dinf, og_size, stride);
toc

tic
results_x = minres(Ax_large'*Ax_large,Ax_large'*bx_large);
results_z = minres(Az_large'*Az_large,Az_large'*bz_large);
toc

results_x = reshape(results_x,size_out);
results_z = reshape(results_z,size_out);

figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, (results_x), [0 2000]);
colormap(gca, parula); 
colorbar;
axis image;
title('k_x')
xlabel('Lateral'), ylabel('Axial'),

nexttile,
imagesc(x*cm, z*cm, (results_z), [0 2000]);
colormap(gca, parula);
colorbar;
axis image;
title('k_z')
xlabel('Lateral'), ylabel('Axial'),

%% Regularized solutions for x and z
M = size_out(1);
N = size_out(2);
mu = 2;% good for 15 kernel
% mu = 10^(-2);
tol = 1e-5;
mask = ones(size(bx_large));
minimask = ones(M*N,1);

tic
results_x = IRLS_TV(bx_large,Ax_large,mu,M,N,tol,mask,minimask);
results_z = IRLS_TV(bz_large,Az_large,mu,M,N,tol,mask,minimask);
toc
results_x = reshape(results_x,size_out);
results_z = reshape(results_z,size_out);


figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, (results_x), [0 2000]);
colormap(gca, parula); 
colorbar;
axis image;
title('k_x')
xlabel('Lateral'), ylabel('Axial'),

nexttile,
imagesc(x*cm, z*cm, (results_z), [0 2000]);
colormap(gca, parula);
colorbar;
axis image;
title('k_z')
xlabel('Lateral'), ylabel('Axial'),

%% Filtering each axis
med_wind = 1;
grad_x = medfilt2(abs(results_x(:,:,1)),[med_wind med_wind],'symmetric');
grad_z = medfilt2(abs(results_z(:,:,1)),[med_wind med_wind],'symmetric');
phase_grad_2 = (grad_x.^2 + grad_z.^2)/constant;
k2_med = phase_grad_2;
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
imagesc(x*cm, z*cm, sws_matrix, [2 5]);
colormap(gca, turbo); % Apply jet colormap to the current axes
colorbar;
axis image;
title('SWS')
xlabel('Lateral'), ylabel('Axial'),