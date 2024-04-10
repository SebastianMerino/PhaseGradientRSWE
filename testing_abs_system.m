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

%% FORMING NEW SYSTEM
tic
[Ax_large, Az_large, bx_large, bz_large, size_out] = getmat_pg_v3(extended_u,...
    w_kernel, dinf, og_size, stride);
toc

% Direct inversion
tic
results_x = minres(Ax_large'*Ax_large,Ax_large'*bx_large);
results_z = minres(Az_large'*Az_large,Az_large'*bz_large);
toc
results_x = reshape(results_x,[size_out, 2]);
results_z = reshape(results_z,[size_out, 2]);

% Modifying system for x
kx = results_x(:,:,1);
bx_mask = ( sign(kx(:)).*ones(1,w_kernel(1)*w_kernel(2)) )';
bx_new = bx_large.*bx_mask(:);

% Modifying system for z
kz = results_z(:,:,1);
bz_mask = ( sign(kz(:)).*ones(1,w_kernel(1)*w_kernel(2)) )';
bz_new = bz_large.*bz_mask(:);

% Inverting again
tic
results_x = minres(Ax_large'*Ax_large,Ax_large'*bx_new);
results_z = minres(Az_large'*Az_large,Az_large'*bz_new);
toc
results_x = reshape(results_x,[size_out, 2]);
results_z = reshape(results_z,[size_out, 2]);

% SWS
grad_x = results_x(:,:,1);
grad_z = results_z(:,:,1);
phase_grad_2 = (grad_x.^2 + grad_z.^2)/constant;
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

%% ========================= REGULARIZATION ========================= 
%% Solving Inverse x problem
mu1 = 10.^(0);
mu2 = 10.^(0);
tol = 1e-5;
mask = ones(size(bx_large));

M = size_out(1);
N = size_out(2);
A1 = Ax_large(:,1:M*N);
A2 = Ax_large(:,M*N+1:end);
tic
[kxx,cx] = AlterOpti_ADMM(A1,A2,bx_new,mu1,mu2,M,N,tol,mask);
toc
kxx = reshape(kxx,size_out);
cx = reshape(cx,size_out);
%%
figure('Units','centimeters', 'Position',[5 5 20 10]),
tl = tiledlayout(1,2);
nexttile,
imagesc(x*cm, z*cm, (kxx), [0 2000]);
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('k_x')
xlabel('Lateral'), ylabel('Axial'),

nexttile,
imagesc(x*cm, z*cm, abs(cx));
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('c')
xlabel('Lateral'), ylabel('Axial'),

%% Solving Inverse z problem
A1 = Az_large(:,1:M*N);
A2 = Az_large(:,M*N+1:end);
mask = ones(size(bz_large));

%
tic
[kzz,cz] = AlterOpti_ADMM(A1,A2,bz_new,mu1,mu2,M,N,tol,mask);
toc
kzz = reshape(kzz,size_out);
cz = reshape(cz,size_out);


figure('Units','centimeters', 'Position',[5 5 20 10]),
tl = tiledlayout(1,2);
nexttile,
imagesc(x*cm, z*cm, (kzz), [0 2000]);
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('k_z')
xlabel('Lateral'), ylabel('Axial'),

nexttile,
imagesc(x*cm, z*cm, abs(cz));
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('c')
xlabel('Lateral'), ylabel('Axial'),

%%
phase_grad_2 = (kxx.^2 + kzz.^2)/constant;
k = sqrt(phase_grad_2);
% --------------------
sws_tv = (2*pi*f_v)./k;   

figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, real(u(:,:,1)));
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('u')
xlabel('Lateral'), ylabel('Axial'),

nexttile,
imagesc(x*cm, z*cm, sws_tv, sws_range);
colormap(gca, turbo); % Apply jet colormap to the current axes
colorbar;
axis image;
title('SWS')
xlabel('Lateral'), ylabel('Axial'),
