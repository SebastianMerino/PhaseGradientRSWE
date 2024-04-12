% SCRIP TO GENERATE SWS MAPS OF HOMOGENEOUS REVERBERANT FIELDS
% Creation: 26/03/2024 (EMZ)
setup,
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography' ...
    '\reverberant\new'];
dataDir = fullfile(baseDir,'sim');
resultsDir = fullfile(baseDir,'results');
[~,~,~]= mkdir(dataDir);
[~,~,~]= mkdir(resultsDir);

%% Setting parameters
% Kernel size and step
w_kernel = [15 15];
constant = 0.3;
stride = round(w_kernel/5);

% Plotting const
cm = 1e2;
sws_range = [2,5];

% Reg constants
mu = 2;
tol = 1e-5;

% file management
files = dir(fullfile(dataDir,'*.mat'));
%% Looping each file

iSim = 5;

% loading data
data = load(fullfile(dataDir,files(iSim).name));
u = data.pv_complexZ(:,:,1);
dinf.dx = data.x(2) - data.x(1);
dinf.dz = data.z(2) - data.z(1);
f_v = data.freq;
x = data.x; z = data.z;

% Constructing matrices
tic
og_size = size(u);
extended_u = padarray(u,[(w_kernel(1)-1)/2 (w_kernel(2)-1)/2],'symmetric');
[A_complex, b_complex, size_out] = getmat_pg_v5(extended_u, ...
    w_kernel, dinf, og_size, stride);
toc

% Direct inversion
tic
results_xz = minres(A_complex'*A_complex,A_complex'*b_complex);
results_xz = reshape(results_xz,size_out);
toc

% Post-processing
med_wind = floor (2.5/f_v/dinf.dx/stride(1))*2+1;
grad_x = real(results_xz(:,:,1));
grad_z = imag(results_xz(:,:,1));
phase_grad_2 = (grad_x.^2 + grad_z.^2)/constant;
% med_wind = 1;
k2_med = medfilt2(phase_grad_2,[med_wind med_wind],'symmetric');
k = sqrt(k2_med);
sws_direct = (2*pi*f_v)./k;   
%%
xSWS = x(1:stride:end);
zSWS = z(1:stride:end);

figure('Units','centimeters', 'Position',[5 5 30 10]),
tiledlayout(1,3)
nexttile,
imagesc(x*cm, z*cm, real(u(:,:,1)));
colormap(gca, parula); 
axis image;
title('PV')
xlabel('Lateral [cm]'), ylabel('Axial [cm]'),
% hold on
% contour(xSWS*cm, zSWS*cm, roi_inc, 'k--')
% contour(xSWS*cm, zSWS*cm, roi_back, 'k--')
% hold off

nexttile,
imagesc(xSWS*cm, zSWS*cm, sws_direct, sws_range);
colormap(gca, turbo); 
axis image;
title('SWS PG')
xlabel('Lateral [cm]')
% nexttile,
% imagesc(xSWS*cm, zSWS*cm, sws_reg, sws_range);
% colormap(gca, turbo); 
% c = colorbar;
% c.Label.String = 'SWS [m/s]';
% axis image;
% title('SWS PG-TV')
% xlabel('Lateral [cm]'),

