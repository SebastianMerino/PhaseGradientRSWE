clear, clc

baseDir = 'P:\rswe\dataold\Data800Hz-10000ondas';
w_kernel = [15,15];
constant = 0.33;
stride = 3;

% plotting const
cm = 1e2;
sws_range = [2,5];

data = load(fullfile(baseDir,'R-FIELD_inc_1.mat'));
u = data.pv_complexZ;
dinf.dx = data.x(2) - data.x(1);
dinf.dz = data.z(2) - data.z(1);
f_v = data.freq;
x = data.x; z = data.z;
og_size = size(u);
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
imagesc(x*cm, z*cm, abs(results_x(:,:,1)), [0 4000]);
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('k_x')
xlabel('Lateral'), ylabel('Axial'),

nexttile,
imagesc(x*cm, z*cm, abs(results_x(:,:,2)));
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('c')
xlabel('Lateral'), ylabel('Axial'),


%% Unwrap z system
figure('Units','centimeters', 'Position',[5 5 20 10]),
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, abs(results_z(:,:,1)));
colormap(gca, parula); % Apply jet colormap to the current axes
colorbar;
axis image;
title('k_z')
xlabel('Lateral'), ylabel('Axial'),

nexttile,
imagesc(x*cm, z*cm, abs(results_z(:,:,2)));
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
med_wind = floor (2.5/f_v/dinf.dx )*2+1; %the median window contains at least a wavelenght
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


%% Solving Inverse x problem
A1 = Ax_large(:,1:167*167);
A2 = Ax_large(:,167*167+1:end);
mu1 = 10.^(-1);
mu2 = 10.^(0);
M = size_out(1);
N = size_out(2);
tol = 1e-2;
mask = ones(size(bx_large));
%
tic
[kxx,cx] = AlterOpti_ADMM(A1,A2,bx_large,mu1,mu2,M,N,tol,mask);
toc
kxx = reshape(kxx,size_out);
cx = reshape(cx,size_out);
%%
figure('Units','centimeters', 'Position',[5 5 20 10]),
tl = tiledlayout(1,2);
nexttile,
imagesc(x*cm, z*cm, abs(kxx), [0 4000]);
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

%%
sws_x = (2*pi*f_v)./kxx;
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
imagesc(x*cm, z*cm, sws_x, sws_range);
colormap(gca, turbo); % Apply jet colormap to the current axes
colorbar;
axis image;
title('SWS')
xlabel('Lateral'), ylabel('Axial'),
