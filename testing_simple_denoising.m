% testing modyfied system
clear, clc
load phasemap.mat

file_dir = 'P:\rswe\dataold\Data800Hz-10000ondas\R-FIELD_inc_1.mat';
% file_dir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
%     'Elastrography\reverberant\new\sim\6_1.mat'];
w_kernel = [15 15];
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


%% =================================================================== %%
%% NEW APPROACH
extended_u = padarray(u,[(w_kernel(1)-1)/2 (w_kernel(2)-1)/2],'symmetric');
stride = 3;
[grad_abs, size_out] = pg_norm(extended_u, w_kernel, dinf, og_size, stride);
sws_denoised = (2*pi*f_v)./grad_abs; 

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
imagesc(x*cm, z*cm, sws_denoised);
colormap(gca, turbo); % Apply jet colormap to the current axes
colorbar;
axis image;
title('SWS')
xlabel('Lateral'), ylabel('Axial'),

%%
med_wind = floor (2.5/f_v/dinf.dx/stride(1))*2+1;
k2_med = medfilt2(grad_abs,[med_wind med_wind],'symmetric');
k = sqrt(k2_med);
sws_direct = (2*pi*f_v)./k;

%%
M = size_out(1); N = size_out(2);
mu = 10.^5; % mu=0.5

tol = 1e-5;
mask = ones(size(M*N));
minimask = ones(M*N,1);

tic
denoised_grad = IRLS_TV(grad_abs(:),speye(M*N),mu,M,N,tol,mask,minimask);
denoised_grad = reshape(denoised_grad, size_out);
toc
k = denoised_grad;
% --------------------
sws_denoised = (2*pi*f_v)./k;   

%%

figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, sws_direct);
colormap(gca, turbo); % Apply jet colormap to the current axes
colorbar;
axis image;
title('SWS')
xlabel('Lateral'), ylabel('Axial'),
nexttile,
imagesc(x*cm, z*cm, sws_denoised);
colormap(gca, turbo); % Apply jet colormap to the current axes
colorbar;
axis image;
title('SWS')
xlabel('Lateral'), ylabel('Axial'),


