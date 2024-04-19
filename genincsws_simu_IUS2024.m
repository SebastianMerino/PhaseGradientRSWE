% SCRIP TO GENERATE SWS MAPS OF HOMOGENEOUS REVERBERANT FIELDS
% Creation: 26/03/2024 (EMZ)

clc, clear all, close all;
addpath(genpath(pwd));

%% GENERATE DATA SWS INCLUSION ONE FRAME WITH PHASE ESTIMATOR QR SOLVER v1.0

addpath(genpath(pwd));
fprintf('-------QR solver 1.0-------\n')

% 500Hz 2.5m/s 1l = 5mm = 51pix 4.5m/s 1l = 9mm = 91pix
% 600Hz 2.5m/s 1l = 4.17mm = 41pix 4.5m/s 1l = 7.5mm = 75pix
% 700Hz 2.5m/s 1l = 3.57mm = 35pix 4.5m/s 1l = 6.42mm = 65pix

% 800Hz 2.5m/s 1l = 3.125mm = 33pix 4.5m/s 1l = 5.625mm = 57pix
% 900Hz 2.5m/s 1l = 2.77mm = 29pix 4.5m/s 1l = 5mm = 51pix
% 1000Hz 2.5m/s 1l = 2.5mmm = 27pix 4.5m/s 1l = 4.5mm = 47pix

nWaves = 10e3; % number of waves
v_freq = [500, 600, 700, 800, 900, 1000];
% v_freq = [500];
nFields = 1;

window = 15; %11 pixels as described in paper
stride = 3;
nFields = 1;

w_kernel = [window, window];
tic;


pathdata = './dataold/';
pathout = './out/old/';

if ~exist("pathout","dir"); mkdir(pathout); end

for ii = 1:length(v_freq)
    freq = v_freq(ii)
    pathfreq_in = [pathdata,'Data', num2str(freq),'Hz-',num2str(nWaves),'ondas/'];
    pathfreq_out = [pathout, 'Out', num2str(freq),'Hz/'];

    if ~exist(pathfreq_out,"dir"); mkdir(pathfreq_out); end

    for field = 1:nFields
        name = ['R-FIELD_inc_',num2str(field),'.mat'];
        R_Field = load([pathfreq_in, name]);
        dinf.dx = min(diff(R_Field.x));
        dinf.dz = min(diff(R_Field.z));
        frame = R_Field.pv_complexZ(:,:,1); % number of frame
        
        %frame = (frame'); %transpose for Z (vertical-axial) X(horizontal-lateral)
        
        og_size = size(frame);
        mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
    
        % FUNCION SWS MAP (Con linearizacion)
        [grad_z,grad_x,k,sws_matrix] = phase_estimator_QR_kernel(mirror_frame, w_kernel, freq, dinf, og_size, 0.33, stride);
    
        % EMPAQUETAR RESULTADOS
        
        pg_QRv1.grad_z = grad_z;
        pg_QRv1.grad_x = grad_x;
        pg_QRv1.grad_k = k;
        pg_QRv1.sws_matrix = sws_matrix;
        pg_QRv1.sws_abs_3D(:,:,ii) = sws_matrix;
       
        
    % Save
%     pathQR = [pathfreq_out, 'PhaseGradientQRv1/'];
%     if ~exist(pathQR,"dir"); mkdir(pathQR); end
%     save([pathQR, 'SWS_PG_QRv1_homo_',num2str(field),'.mat'],'pg_QRv1');

    end
end
toc
fprintf('---------------------------\n')

%%

%% PLOT QR v1 KERNEL

caxis_sws = [2 5];
figure, % SWS
sgtitle(['SWS PGvKernel w=', num2str(window(1)),', str=', num2str(stride)])
set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
% set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.55 0.25])
for ii = 1 : length(v_freq)
    freq = v_freq(ii);
%     subplot (2, 3, ii)
    subplot (2, 3, ii) % IUS ABSTRACT
    imagesc(pg_QRv1.sws_abs_3D(:,:,ii), caxis_sws)
%         colormap ("jet");
    colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar, axis("tight")
    title (['SWS f_v = ', num2str(freq) ])

end

%%
tol = 1e-3;
mask = ones(M*N,1);
isotropic = true;
colMajor = false;


%% DENOISING 

% GENERATE DATA SWS INCLUSION ONE FRAME WITH PHASE ESTIMATOR QR SOLVER v1.0

addpath(genpath(pwd));
fprintf('-------QR solver 1.0-------\n')
% resol
% 1 pix = 0.154mm
% 11 pixe 1.69mm

% 500Hz 2.5m/s 1l = 5mm = 51pix 4.5m/s 1l = 9mm = 91pix
% 600Hz 2.5m/s 1l = 4.17mm = 41pix 4.5m/s 1l = 7.5mm = 75pix
% 700Hz 2.5m/s 1l = 3.57mm = 35pix 4.5m/s 1l = 6.42mm = 65pix

% 800Hz 2.5m/s 1l = 3.125mm = 33pix 4.5m/s 1l = 5.625mm = 57pix
% 900Hz 2.5m/s 1l = 2.77mm = 29pix 4.5m/s 1l = 5mm = 51pix
% 1000Hz 2.5m/s 1l = 2.5mmm = 27pix 4.5m/s 1l = 4.5mm = 47pix

nWaves = 10e3; % number of waves
v_freq = [500, 600, 700, 800, 900, 1000];
v_freq = [1000];
nFields = 1;

window = 15; %11 pixels as described in paper
stride = 3;
nFields = 1;

w_kernel = [window, window];
tic;


pathdata = './dataold/';
pathout = './out/old/';

if ~exist("pathout","dir"); mkdir(pathout); end

for freq = v_freq
   
    pathfreq_in = [pathdata,'Data', num2str(freq),'Hz-',num2str(nWaves),'ondas/'];
    pathfreq_out = [pathout, 'Out', num2str(freq),'Hz/'];

    if ~exist(pathfreq_out,"dir"); mkdir(pathfreq_out); end

    name = ['R-FIELD_inc_',num2str(field),'.mat'];
    R_Field = load([pathfreq_in, name]);
    dinf.dx = min(diff(R_Field.x));
    dinf.dz = min(diff(R_Field.z));
    frame = R_Field.pv_complexZ(:,:,1); % number of frame
    
    %frame = (frame'); %transpose for Z (vertical-axial) X(horizontal-lateral)
    
    og_size = size(frame);
    mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');

    % FUNCION SWS MAP (Con linearizacion)
    [grad_abs, size_out] = pg_norm(mirror_frame, w_kernel, dinf, og_size, stride);
     M = size_out(1); N = size_out(2);

     pathMat = [pathfreq_out, 'Matrix/'];
    if ~exist(pathMat,"dir"); mkdir(pathMat); end
    save([pathMat, 'PG_abs_str',num2str(stride),'.mat'],'grad_abs','size_out', ...
        'window', 'stride');



end
toc
fprintf('---------------------------\n')


%% ORIGINAL 3D CREATION

v_freq = [500, 600, 700, 800, 900, 1000];
numChannels = length(v_freq);
stride = 3;

M = length(1:stride:501); % size_out(1)
N = length(1:stride:501); % size_out(2)

grad_abs_3D = zeros(M, N, numChannels); 
sws_abs_3D = grad_abs_3D;
for ii = 1 : numChannels

    freq = v_freq(ii);
    my_obj = load(['./out/old/Out', num2str(freq), 'Hz/Matrix/PG_abs_str' num2str(stride), '.mat']);

    og.grad_abs_3D (:, :, ii) = my_obj.grad_abs;     
    og.sws_abs_3D(:,:, ii) = 2*pi*freq ./ og.grad_abs_3D (:, :, ii);

end

%% PLOT ORIGINAL
caxis_sws = [0 5];
figure, % SWS
sgtitle('SWS')
set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
for ii = 1 : numChannels
    freq = v_freq(ii);
    subplot (2, 3, ii)
    imagesc(og.sws_abs_3D(:,:,ii), caxis_sws)
        colormap ("jet");
%     colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS f_v = ', num2str(freq) ])

end

figure, % grad phi
sgtitle('|\nabla\phi|')
for ii = 1 : numChannels
    freq = v_freq(ii);
    subplot (2, 3, ii)
    imagesc(og.grad_abs_3D(:,:,ii))
    colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['\nabla\phi f_v = ', num2str(freq) ])

end
%% MEDIAN FILTER

for ii = 1 : numChannels
    freq = v_freq(ii); 
    med_wind = [15];
    medf.grad_abs_3D(:,:,ii) = medfilt2(og.grad_abs_3D(:,:,ii),[med_wind med_wind],'symmetric');
    medf.sws_abs_3D(:,:,ii) = 2*pi*freq ./ medf.grad_abs_3D (:, :, ii);

end

%% PLOT MED FILT 
caxis_sws = [0 5];
figure, % SWS
sgtitle('SWS Median filter')
set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
for ii = 1 : numChannels
    freq = v_freq(ii);
    subplot (2, 3, ii)
    imagesc(medf.sws_abs_3D(:,:,ii), caxis_sws)
        colormap ("jet");
%     colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS f_v = ', num2str(freq) ])

end

%% AVERAGE FILTER

for ii = 1 : numChannels
    freq = v_freq(ii);
    avg_kernel = ones(7, 7) / 49;  % Create a 7x7 averaging filter kernel
    avef.grad_abs_3D(:,:,ii) = filter2(avg_kernel, og.grad_abs_3D(:,:,ii), 'same');
    avef.sws_abs_3D(:,:,ii) = 2*pi*freq ./ avef.grad_abs_3D(:, :, ii);
end

% PLOT AVERAGE FILT 
caxis_sws = [2 5];
figure, % SWS
sgtitle('SWS Average filter')
set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
for ii = 1 : numChannels
    freq = v_freq(ii);
    subplot (2, 3, ii)
    imagesc(avef.sws_abs_3D(:,:,ii), caxis_sws)
%         colormap ("jet");
    colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS f_v = ', num2str(freq) ])

end
%% TOTAL VARIATION
M = my_obj.size_out(1);
N = my_obj.size_out(2);

mu = 1e4;
for ii = 1 : numChannels
    freq = v_freq(ii);

    my_grad = og.grad_abs_3D (:, :, ii);

    [grad_tv] = IRLS_TV_simple(my_grad(:),speye(M*N),mu,M,N,1e-4,ones(size(M*N)),ones(M*N,1));

    tv.grad_abs_3D(:,:,ii) = reshape(grad_tv, [ M N ] );

    tv.sws_abs_3D(:,:,ii) = (2*pi*freq)./tv.grad_abs_3D(:,:,ii);

end

%% PLOT  TV
figure, % SWS
caxis_sws = [2 5];

% sgtitle(['SWS TV, \mu=', num2str(mu)])
sgtitle(['\bfSWS TV'])
set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
for ii = 1 : numChannels
    freq = v_freq(ii);
    subplot (2, 3, ii)
    imagesc(x*1e2, z*1e2, tv.sws_abs_3D(:,:,ii), caxis_sws)
%     colormap ("jet");
    colormap("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS f_v = ', num2str(freq) ])

end

% figure, % grad phi
% sgtitle('|\nabla\phi| TV')
% for ii = 1 : numChannels
%     freq = v_freq(ii);
%     subplot (2, 3, ii)
%     imagesc(tv.grad_abs_3D(:,:,ii))
%     colormap ("turbo");
%     xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
%     title (['\nabla\phi f_v = ', num2str(freq) ])
% 
% end

%% TOTAL NUCLEAR VARIATION
M = my_obj.size_out(1);
N = my_obj.size_out(2);

bestmu = 10^4.1;

besttau = 0.001;
maxIter = 1000;
stableIter = 50;
tol = 5e-4; % tolerance error

weightEstimators = ones(1, length(v_freq));

clear tnv
[tnv.grad_abs_3D, cost, error, fid, reg] = pdo_den_wtnv(og.grad_abs_3D, bestmu, besttau, maxIter, tol, stableIter, weightEstimators); 


tnv.sws_abs_3D =  2*pi* reshape(v_freq, [1, 1, numChannels]) ./ tnv.grad_abs_3D; % more elegant

%% PLOT TNV

caxis_sws = [2 5];
figure, % SWS
sgtitle('\bfSWS TNV')
set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
for ii = 1 : numChannels
    freq = v_freq(ii);
    subplot (2, 3, ii)
    imagesc(x*1e2, z*1e2,tnv.sws_abs_3D(:,:,ii), caxis_sws), axis('tight')
%     colormap ("jet");
    colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS f_v = ', num2str(freq) ])

end

% figure, % grad phi
% sgtitle('|\nabla\phi| TNV')
% for ii = 1 : numChannels
%     freq = v_freq(ii);
%     subplot (2, 3, ii)
%     imagesc(tnv.grad_abs_3D(:,:,ii))
%     colormap ("turbo");
%     xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
%     title (['\nabla\phi f_v = ', num2str(freq) ])
% 
% end

%%
%% METRICS


mm = 1e3;

% Parameters to select ROI
L = 7.5;
cx = 0; cz = 0;
d = 5;

%% REGION METRICS
%%%%%%%%%%% OLD FORM %%%%%%%%%%%%%%%%
% [back,inc] = getRegionMasks(x*mm,z*mm,cx,cz,L,d,L);

%%%%%%% METRICS MASK %%%%%%%

%%%%%%%%%%%%%%%%%% INCLUSION %%%%%%%%%%%%%%%%%%
cx1 = -5e-3; % [m] 
cz1 = -5e-3; % [m] 

cx2 = 5e-3; % [m]
cz2 = 5e-3; % [m]

%%%%%%%%%%%%%%%%%% BACKGROUND %%%%%%%%%%%%%%%%%%
dist_x = 10e-3; %[m]
dist_z = 20e-3; %[m]

bx1_l = -20e-3; %[m]
bx2_l = bx1_l + dist_x;
bz1 = -10e-3; %[m]
bz2 = bz1 + dist_z;

bx1_r = 10e-3; %[m]
bx2_r = bx1_r + dist_x;

[~, inc] = mask_rect(x,z, cx1, cx2, cz1, cz2,  pv_complexZ(:,:,1)); %inclusion

[~, bl_mask] = mask_rect(x,z, bx1_l, bx2_l, bz1, bz2,  pv_complexZ(:,:,1)); %background l
[~, br_mask] = mask_rect(x,z, bx1_r, bx2_r, bz1, bz2,  pv_complexZ(:,:,1));  %background r
back = or(bl_mask, br_mask);

figure, 
subplot(1,2,1), imagesc(x*1e3, z*1e3, inc_mask), title('Inclusion'), colorbar;
subplot(1,2,2), imagesc(x*1e3, z*1e3, back_mask), title('Background'), colorbar;
%%

%%%%%%% METRICS MASK %%%%%%%
gt_inc = 4.5; % [m/s]
gt_back = 2.5; % [m/s]

clear MetricsAveF MetricsTV MetricsTNV
for ii = 1 : numChannels
    freq = v_freq(ii);

    sws = bigImg( avef.sws_abs_3D(:,:,ii) , pv_complexZ(:,:,1) );
    MetricsAveF(ii) = get_metrics_v2(sws,inc,back,'avef',freq, gt_inc, gt_back);

    sws =  bigImg( tv.sws_abs_3D(:,:,ii) , pv_complexZ(:,:,1) );
    MetricsTV(ii) = get_metrics_v2(sws,inc,back,'tv',freq, gt_inc, gt_back);

    sws = bigImg( tnv.sws_abs_3D(:,:,ii), pv_complexZ(:,:,1) );
    MetricsTNV(ii) = get_metrics_v2(sws,inc,back,'tnv',freq, gt_inc, gt_back);
end

T = [struct2table(MetricsAveF);
    struct2table(MetricsTV);
    struct2table(MetricsTNV)];
writetable(T,fullfile(pathout,'results.xlsx'),'WriteRowNames',true);
close all

%% Plots

figure,
hold on
errorbar(v_freq,T.mean_inc(1:6),T.std_inc(1:6), 'r')
errorbar(v_freq,T.mean_inc(7:12),T.std_inc(7:12), 'g')
errorbar(v_freq,T.mean_inc(13:18),T.std_inc(13:18), 'b')
yline(4.5, 'k--', 'Ground Truth'); % 'k--' denotes a black dashed line

hold off
legend('AveF','TV','TNV', 'Location','northwest')
grid on
xlim([450 1050])
title('SWS dispersion inclusion')
ylabel('SWS [m/s]'), xlabel('Frequency [Hz]')

figure,
hold on

errorbar(v_freq,T.mean_back(1:6),T.std_back(1:6), 'r')
errorbar(v_freq,T.mean_back(7:12),T.std_back(7:12), 'g')
errorbar(v_freq,T.mean_back(13:18),T.std_back(13:18), 'b')
yline(2.5, 'k--', 'Ground Truth'); % 'k--' denotes a black dashed line
hold off
legend('AveF','TV','TNV', 'Location','northwest')
grid on
xlim([450 1050])
title('SWS dispersion background')
ylabel('SWS [m/s]'), xlabel('Frequency [Hz]')

figure,
hold on
plot(v_freq,T.cnr(1:6), 'ro-')
plot(v_freq,T.cnr(7:12), 'go-')
plot(v_freq,T.cnr(13:18), 'bo-')

hold off
legend('AveF','TV','TNV', 'Location','northwest')
grid on
xlim([450 1050])
ylim([1 7])
title('CNR')
xlabel('Frequency [Hz]')

%% Figures
% PLotting constants
zlim_mm = [-25 25];
caxis_sws = [2 5];
fontSize = 20;
fontText = 20;
lw = 2;
mm = 1e3;

% Upper left corner of each background rectangle
x0 = cx - L/2; z0 = cz-L/2;
xb1 = x0 - d - L/2;
xb2 = x0 + L + d;
figure('Position',[0 100 2600 900]), % SWS
tiledlayout(3,numChannels, 'TileSpacing','loose', 'Padding','compact')
% sgtitle('SWS TNV')
% set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.75 0.35])

for ii = 1 : numChannels
    freq = v_freq(ii);
%     subplot (2, 3, ii)
    % subplot (1, numChannels+1, ii+1) % abstract IUS203
    nexttile;
    imagesc(x*mm,z*mm,medf.sws_abs_3D(:,:,ii), caxis_sws), axis('tight')
    axis image
%     rectangle('Position',[x0, z0,L,L], 'LineStyle','-', 'LineWidth',2)
%     rectangle('Position',[xb1, z0,L/2,L], 'LineStyle','-', 'LineWidth',2)
%     rectangle('Position',[xb2, z0,L/2,L], 'LineStyle','-', 'LineWidth',2)

    % INCLUSION REGION METRICS
    rectangle('Position',mm*[cx1, cz1, cx2-cx1, cz2-cz1],'LineWidth', lw), hold on;
    % BACKGROUND REGION METRICS
    % LEFT
    rectangle('Position', mm*[bx1_l bz1, bx2_l-bx1_l, bz2-bz1],'EdgeColor','k','LineWidth',lw,'LineStyle','-');
    % RIGHT
    rectangle('Position', mm*[bx1_r bz1, bx2_r-bx1_r, bz2-bz1],'EdgeColor','k','LineWidth',lw,'LineStyle','-');


    colormap ("turbo");
    title (['PG, f = ', num2str(freq),'Hz'])
    ylim(zlim_mm)
    xlim(zlim_mm)
    if ii==1; ylabel('Axial [mm]'); end
    text(0,12,['\bfCNR = ',num2str(T.cnr(ii),2)], ...
        'HorizontalAlignment', 'center', 'FontSize', fontText)
    text(0,17,['\bfBias_{in} = ',num2str(T.bias_inc(ii),2), '%'], ...
        'HorizontalAlignment', 'center', 'FontSize', fontText)
    text(0,22,['\bfBias_{bg} = ',num2str(T.bias_back(ii),2), '%'], ...
        'HorizontalAlignment', 'center', 'FontSize', fontText)
    set(gca, 'FontSize',fontSize)

end
c = colorbar;
c.Label.String = 'SWS [m/s]';

%nexttile;
%axis off

for ii = 1 : numChannels
    freq = v_freq(ii);
%     subplot (2, 3, ii)
    % subplot (1, numChannels+1, ii+1) % abstract IUS203
    nexttile;
    imagesc(x*mm,z*mm,tv.sws_abs_3D(:,:,ii), caxis_sws), axis('tight')
    axis image
%     rectangle('Position',[x0, z0,L,L], 'LineStyle','-', 'LineWidth',2)
%     rectangle('Position',[xb1, z0,L/2,L], 'LineStyle','-', 'LineWidth',2)
%     rectangle('Position',[xb2, z0,L/2,L], 'LineStyle','-', 'LineWidth',2)

    % INCLUSION REGION METRICS
    rectangle('Position',mm*[cx1, cz1, cx2-cx1, cz2-cz1],'LineWidth', lw), hold on;
    % BACKGROUND REGION METRICS
    % LEFT
    rectangle('Position', mm*[bx1_l bz1, bx2_l-bx1_l, bz2-bz1],'EdgeColor','k','LineWidth',lw,'LineStyle','-');
    % RIGHT
    rectangle('Position', mm*[bx1_r bz1, bx2_r-bx1_r, bz2-bz1],'EdgeColor','k','LineWidth',lw,'LineStyle','-');

    colormap ("turbo");
    xlabel('Lateral [mm]'),
    if ii==1; ylabel('Axial [mm]'); end
    title (['PG-TV, f = ', num2str(freq),'Hz'])
    ylim(zlim_mm)
    text(0,12,['\bfCNR = ',num2str(T.cnr(ii+6),2)], ...
        'HorizontalAlignment', 'center', 'FontSize', fontText)
    text(0,17,['\bfBias_{in} = ',num2str(T.bias_inc(ii+6),2), '%'], ...
        'HorizontalAlignment', 'center', 'FontSize', fontText)
    text(0,22,['\bfBias_{bg} = ',num2str(T.bias_back(ii+6),2), '%'], ...
        'HorizontalAlignment', 'center', 'FontSize', fontText)
    set(gca, 'FontSize',fontSize)
end
c = colorbar;
c.Label.String = 'SWS [m/s]';


for ii = 1 : numChannels
    freq = v_freq(ii);
%     subplot (2, 3, ii)
    % subplot (1, numChannels+1, ii+1) % abstract IUS203
    nexttile;
    imagesc(x*mm,z*mm,tnv.sws_abs_3D(:,:,ii), caxis_sws), axis('tight')
    axis image
%     rectangle('Position',[x0, z0,L,L], 'LineStyle','-', 'LineWidth',2)
%     rectangle('Position',[xb1, z0,L/2,L], 'LineStyle','-', 'LineWidth',2)
%     rectangle('Position',[xb2, z0,L/2,L], 'LineStyle','-', 'LineWidth',2)

    % INCLUSION REGION METRICS
    rectangle('Position',mm*[cx1, cz1, cx2-cx1, cz2-cz1],'LineWidth', lw), hold on;
    % BACKGROUND REGION METRICS
    % LEFT
    rectangle('Position', mm*[bx1_l bz1, bx2_l-bx1_l, bz2-bz1],'EdgeColor','k','LineWidth',lw,'LineStyle','-');
    % RIGHT
    rectangle('Position', mm*[bx1_r bz1, bx2_r-bx1_r, bz2-bz1],'EdgeColor','k','LineWidth',lw,'LineStyle','-');

    colormap ("turbo");
    xlabel('Lateral [mm]'),
    if ii==1; ylabel('Axial [mm]'); end
    title (['PG-TNV, f = ', num2str(freq),'Hz'])
    ylim(zlim_mm)
    text(0,12,['\bfCNR = ',num2str(T.cnr(ii+12),2)], ...
        'HorizontalAlignment', 'center', 'FontSize', fontText)
    text(0,17,['\bfBias_{in} = ',num2str(T.bias_inc(ii+12),2), '%'], ...
        'HorizontalAlignment', 'center', 'FontSize', fontText)
    text(0,22,['\bfBias_{bg} = ',num2str(T.bias_back(ii+12),2), '%'], ...
        'HorizontalAlignment', 'center', 'FontSize', fontText)

    set(gca, 'FontSize',fontSize)
end
c = colorbar;
c.Label.String = 'SWS [m/s]';



%% TOTAL NUCLEAR VARIATION grid search
M = my_obj.size_out(1);
N = my_obj.size_out(2);

tau_vector = [0.001];
mu_vector = logspace(log10(0.1),log10(10^5),10); % logarithmic
maxIter = 1000;
stableIter = 50;
tol = 5e-4; % tolerance error

weightEstimators = ones(1, length(v_freq));


% van quedando mu = 10^3 10^3.67

% mu_vector = 10^3.5;
% tau_vector = [0.1 0.05 0.01 0.005 0.001 0.0005 0.0001];


v_freq = 500:100:1000;
numChannels = length(v_freq);

for u = 1:length(mu_vector)
    bestmu = mu_vector(u);
    for t = 1:length(tau_vector)
        besttau = tau_vector(t);

        clear tnv
        [tnv.grad_abs_3D, cost, error, fid, reg] = pdo_den_wtnv(og.grad_abs_3D, bestmu, besttau, maxIter, tol, stableIter, weightEstimators); 

%         for ii = 1 : numChannels
%         
%             freq = v_freq_best(ii);
%         
%             tnv.sws_abs_3D(:,:, ii) = (2*pi*freq)./tnv.grad_abs_3D(:,:,ii);
%         
%         end

        tnv.sws_abs_3D =  2*pi* reshape(v_freq, [1, 1, numChannels]) ./ tnv.grad_abs_3D; % more elegant

        caxis_sws = [2 5];
        figure, % SWS
        sgtitle(['SWS TNV, \mu=10^{', num2str(log10(bestmu)), '} \tau=', num2str(besttau)]);
        set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
        for ii = 1 : numChannels
            freq = v_freq(ii);
            subplot (2, 3, ii)
            imagesc(tnv.sws_abs_3D(:,:,ii), caxis_sws), axis('tight')
%             colormap ("jet");
            colormap ("turbo");
            xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
            title (['SWS f_v = ', num2str(freq) ])
        
        end

    end
end
