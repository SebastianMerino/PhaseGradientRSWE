%% GENERATE NORM PHASE GRADIENT STORE MATRIX

addpath(genpath(pwd));
fprintf('-------Norm Phase gradient -------\n')


list_data = [10 11 13];
v_freq = [400 600 900];

window = 15; %11 pixels as described in paper
stride = 1;
w_kernel = [window, window];
constant = 1; % 1.16 gives good results
phaseExtrac = 'JO';
tic;

% pathdata = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\reverberant\CIRS_phantom\L7-4';
pathdata = 'C:\Users\smerino.C084288\Documents\MATLAB\Datasets\RSWE-PG\Data_3_PUCP\CIRS_phantom\L7-4';
pathout = fullfile(pathdata,'AromCode');

if ~exist("pathout","dir"); mkdir(pathout); end

for ii = 1: length(v_freq)
   
    freq = v_freq(ii);
    sampleCase = list_data(ii);

    pathfreq_in = fullfile(pathdata,['data_', num2str(sampleCase)]);
    pathfreq_out = fullfile(pathout, ['data_', num2str(sampleCase)]);

    if ~exist(pathfreq_out,"dir"); mkdir(pathfreq_out); end

    name = ['data_',num2str(sampleCase),'.mat'];
    structdata = load(fullfile(pathfreq_in, name));

    dinf.dx = structdata.dinf.dx;
    dinf.dz = structdata.dinf.dz; % MISTAKE
    dinf.PRFe = structdata.dinf.PRFe;


    if strcmp(phaseExtrac, 'JO')
        folderFigout = 'fig_uPhaseExtrJO';
        [u_out] = fun_JO_v1(structdata.u, freq, dinf);
    end
    if strcmp(phaseExtrac, 'LIM')
        folderFigout = 'fig_uPhaseExtrLIM';
        u2 = signal_period(freq, dinf.PRFe, structdata.u);   

        % Temperal filtering process, a bandpass FIR filter is used around 
        % +- 20 Hz the vibration frequency
        cs_min = 1; % [m/s]
        cs_max = 5;   % [m/s]
        f_tol = 10;   % [Hz] tolerance for +/-2*f_tol
        % This function uses spatial_fil_phase_extrac inside
        [u_new, u_out, Frames1] = u_filt(u2, freq, f_tol, dinf, cs_min, cs_max);           

    end

    % figure, tiledlayout(1,2)
    % nexttile,
    % imagesc(x,z,structdata.u(:,:,end), std(structdata.u(:))*[-1 1])
    % axis image
    % nexttile,
    % imagesc(x,z,real(u_out), std(u_out(:))*[-1 1])
    % axis image


    frame = u_out;
    og_size = size(frame);
    mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
   
    [grad_abs, size_out] = pg_norm(mirror_frame, w_kernel, dinf, og_size, stride);

    % Save
    pathMat = fullfile(pathfreq_out, ['MatrixW',num2str(w_kernel(1))] );
    if ~exist(pathMat,"dir"); mkdir(pathMat); end
    save(fullfile(pathMat, ['PG_abs_str',num2str(stride),'data', num2str(sampleCase), '.mat']) ...
        ,'grad_abs','size_out', ...
            'window', 'stride', 'freq');
end
toc
fprintf('---------------------------\n')

dinf.dz = structdata.dinf.dz;
x = (0:size(grad_abs,2)-1)*dinf.dx*stride;
z = (0:size(grad_abs,1)-1)*dinf.dz*stride;
x = x - mean(x);
%% ORIGINAL 3D CREATION

list_data = [10 11 13];
v_freq = [400 600 900];

v_freq_best = v_freq;
list_data_best = list_data;


numChannels = length(v_freq_best);
stride = 1;
window = 15;

M = length(1:stride:290); % size_out(1)
N = length(1:stride:248); % size_out(2)

grad_abs_3D = zeros(M, N, numChannels); 
sws_abs_3D = grad_abs_3D;
clear og;
for ii = 1 : numChannels

    freq = v_freq_best(ii);
    sample = list_data_best(ii);
%     my_obj = load(['./out/JO/', num2str(freq), 'Hz/Matrix', '/PG_abs_str' num2str(stride),'data',num2str(sample) ,'.mat']);
    my_obj = load(fullfile(pathout,['data_', num2str(sample)], ...
        ['MatrixW', num2str(window)], ...
        ['PG_abs_str' num2str(stride),'data',num2str(sample) ,'.mat']));
    og.grad_abs_3D (:, :, ii) = my_obj.grad_abs;     
    og.sws_abs_3D(:,:, ii) = 2*pi*freq ./ og.grad_abs_3D (:, :, ii);

end

%% PLOT ORIGINAL
caxis_sws = [1 6];
figure, % SWS
sgtitle(['SWS w=', num2str(window),', str=', num2str(stride)])
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.75 0.35])
for ii = 1 : numChannels
    freq = v_freq_best(ii);
    subplot (1, numChannels, ii)
    imagesc(x,z, og.sws_abs_3D(:,:,ii), caxis_sws)
    axis image
%         colormap ("jet");
    colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS f_v = ', num2str(freq) ])

end

% figure, % grad phi
% sgtitle('|\nabla\phi|')
% for ii = 1 : numChannels
%     freq = v_freq_best(ii);
%     subplot (1, 4, ii)
%     imagesc(og.grad_abs_3D(:,:,ii))
%     colormap ("turbo");
%     xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
%     title (['\nabla\phi f_v = ', num2str(freq) ])
% 
% end
%% MEDIAN FILTER
for ii = 1 : numChannels
    freq = v_freq_best(ii); 
    med_wind = [my_obj.window];
    medf.grad_abs_3D(:,:,ii) = medfilt2(og.grad_abs_3D(:,:,ii),[med_wind med_wind],'symmetric');
    medf.sws_abs_3D(:,:,ii) = 2*pi*freq ./ medf.grad_abs_3D (:, :, ii) * constant;

end

% PLOT MED FILT 
caxis_sws = [1 4];
figure, % SWS
sgtitle('SWS Median filter')
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.75 0.35])
for ii = 1 : numChannels
    freq = v_freq_best(ii);
    subplot (1, numChannels, ii)
    imagesc(x,z,medf.sws_abs_3D(:,:,ii), caxis_sws)
%         colormap ("jet");
    colormap ("turbo");
    axis image
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS f_v = ', num2str(freq) ])

end

%% AVERAGE FILTER

for ii = 1 : numChannels
    freq = v_freq(ii);
    avg_kernel = ones(7, 7) / 49;  % Create a 7x7 averaging filter kernel
    avef.grad_abs_3D(:,:,ii) = filter2(avg_kernel, og.grad_abs_3D(:,:,ii), 'same');
    avef.sws_abs_3D(:,:,ii) = 2*pi*freq ./ avef.grad_abs_3D(:, :, ii) * constant;
end

% PLOT AVERAGE FILT 
caxis_sws = [0 4];
figure, % SWS
sgtitle('SWS Average filter')
set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.3])
for ii = 1 : numChannels
    freq = v_freq(ii);
    subplot (1, 3, ii)
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
clear tv;
for ii = 1 : numChannels
    freq = v_freq_best(ii);

    my_grad = og.grad_abs_3D (:, :, ii);

    [grad_tv] = IRLS_TV(my_grad(:),speye(M*N),mu,M,N,1e-4,ones(size(M*N)),ones(M*N,1));

    tv.grad_abs_3D(:,:,ii) = reshape(grad_tv, [ M N ] );

    tv.sws_abs_3D(:,:,ii) = (2*pi*freq)./tv.grad_abs_3D(:,:,ii)* constant;

end

% PLOT  TV
figure, % SWS
caxis_sws = [1 4];

sgtitle('SWS TV')
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.75 0.35])
for ii = 1 : numChannels
    freq = v_freq_best(ii);
    subplot (1, numChannels, ii)
    imagesc(x,z,tv.sws_abs_3D(:,:,ii), caxis_sws)
%     colormap ("jet");
    colormap("turbo");
    axis image
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS f_v = ', num2str(freq) ])

end

% figure, % grad phi
% sgtitle('|\nabla\phi| TV')
% for ii = 1 : numChannels
%     freq = v_freq_best(ii);
%     subplot (2, 3, ii)
%     imagesc(x,z,tv.grad_abs_3D(:,:,ii))
%     colormap ("turbo");
%     axis image
%     xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
%     title (['\nabla\phi f_v = ', num2str(freq) ])
% 
% end

%% TOTAL NUCLEAR VARIATION
M = my_obj.size_out(1);
N = my_obj.size_out(2);

bestmu = 10^3.5;

besttau = 10.^-2.5;
maxIter = 1000;
stableIter = 50;
tol = 10e-4; % tolerance error

weightEstimators = ones(1, length(v_freq));
clear tnv
[tnv.grad_abs_3D, cost, error, fid, reg] = pdo_den_wtnv(og.grad_abs_3D, bestmu, besttau, maxIter, tol, stableIter, weightEstimators); 

% 
% for ii = 1 : numChannels
% 
%     freq = v_freq_best(ii);
% 
%     tnv.sws_abs_3D(:,:, ii) = (2*pi*freq)./tnv.grad_abs_3D(:,:,ii);
% 
% end

tnv.sws_abs_3D =  2*pi* reshape(v_freq_best, [1, 1, numChannels]) ./ tnv.grad_abs_3D * constant; % more elegant


% tnv.sws_abs_3D_big = bigImg(tnv.sws_abs_3D, pg_QRv1.sws_abs_3D);
%% PLOT TNV

caxis_sws = [1.01 4];
figure, % SWS
sgtitle('SWS TNV')
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.75 0.35])

for ii = 1 : numChannels
    freq = v_freq_best(ii);
%     subplot (2, 3, ii)
    subplot (1, numChannels, ii) % abstract IUS203

    imagesc(x,z,tnv.sws_abs_3D(:,:,ii), caxis_sws), axis('tight')
    axis image
%     colormap ("jet");
    colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS f_v = ', num2str(freq) ])

end


%% METRICS
mm = 1e3;
Bmode = db(structdata.IQBmodeData);
Bmode = Bmode - max(Bmode(:));
xBm = (0:size(Bmode,2)-1)*dinf.dx*stride;
xBm = xBm - mean(xBm);
zBm = (0:size(Bmode,1)-1)*dinf.dz*stride;

% REGION METRICS
%%%%%%%%%%% OLD FORM %%%%%%%%%%%%%%%%
% [back,inc] = getRegionMasks(x*mm,z*mm,cx,cz,L,d,L);

%%%%%%% METRICS MASK %%%%%%%

%%%%%%%%%%%%%%%%%% INCLUSION %%%%%%%%%%%%%%%%%%
cx = 0.5e-3; % [m] 
cz = 22.5e-3; % [m] 
Lx = 7.5e-3;
Lz = 7.5e-3;

%%%%%%%%%%%%%%%%%% BACKGROUND %%%%%%%%%%%%%%%%%%
dist_x = 12e-3; %[m]
Lzb = 15e-3;

cx1b = cx + dist_x;
cx2b = cx - dist_x;

inc = mask_rect_v2(x, z, cx, cz, Lx, Lz);
bl_mask = mask_rect_v2(x, z, cx1b, cz, Lx, Lzb);
br_mask = mask_rect_v2(x, z, cx2b, cz, Lx, Lzb);
back = or(bl_mask, br_mask);

figure('Position',[100 100 400 250]), 
tiledlayout(1,2)
nexttile,
imagesc(xBm*mm,zBm*mm, Bmode, [-50 0])
% c =colorbar(t1,'westoutside');
% colorbar
colormap gray
hold on
contour(x*mm,z*mm,back,1,'k--', 'LineWidth',1);
contour(x*mm,z*mm,inc,1,'--', 'LineWidth',1);
hold off
title('B-mode')
xlabel('Lateral [mm]'), ylabel('Axial [mm]'),
axis image
ylim([4 46])

t2= nexttile;
imagesc(xBm*mm,zBm*mm, real(frame), [-1 1]*1e-5)
% c =colorbar(t1,'westoutside');
% colorbar
colormap(t2,parula)
% hold on
% contour(x*mm,z*mm,back,1,'k--', 'LineWidth',1);
% contour(x*mm,z*mm,inc,1,'--', 'LineWidth',1);
% hold off
title('Particle velocity')
xlabel('Lateral [mm]'), ylabel('Axial [mm]'),
axis image
ylim([4 46])

saveas(gcf,'./figs/Bm.png')
%%
% [back,inc] = getRegionMasks(x*mm,z*mm,cx,cz,L,d,L);
gtInc = 3.65;
gtBack = 2.6;

clear MetricsAveF MetricsTV MetricsTNV
for ii = 1 : numChannels
    freq = v_freq(ii);

    sws = avef.sws_abs_3D(:,:,ii);
    MetricsAveF(ii) = get_metrics_v2(sws,inc,back,'avef',freq, gtInc, gtBack);

    sws =  tv.sws_abs_3D(:,:,ii);
    MetricsTV(ii) = get_metrics_v2(sws,inc,back,'tv',freq, gtInc, gtBack);

    sws = tnv.sws_abs_3D(:,:,ii);
    MetricsTNV(ii) = get_metrics_v2(sws,inc,back,'tnv',freq, gtInc, gtBack);
end

T = [struct2table(MetricsAveF);
    struct2table(MetricsTV);
    struct2table(MetricsTNV)];
writetable(T,fullfile(pathout,'results.xlsx'),'WriteRowNames',true);
close all



%% Figures
% PLotting constants
zlim_mm = [2 43];
caxis_sws = [1.5 4.5];
fontSize = 10;
fontText = 10;
lw = 1.5;
mm = 1e3;

% Upper left corner of each background rectangle
x0 = cx - Lx/2; z0 = cz-Lz/2;
z0b = cz-Lzb/2;
xb1 = x0 - dist_x;
xb2 = x0 + dist_x;

figure('Position',[100 100 650 600]), % SWS
tiledlayout(3,numChannels, 'TileSpacing','compact', 'Padding','loose')
% sgtitle('SWS TNV')
% set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.75 0.35])

for ii = 1 : numChannels
    freq = v_freq(ii);
%     subplot (2, 3, ii)
    % subplot (1, numChannels+1, ii+1) % abstract IUS203
    nexttile;
    imagesc(x*mm,z*mm,medf.sws_abs_3D(:,:,ii), caxis_sws), axis('tight')
    axis image

    % % INCLUSION REGION METRICS
    % rectangle('Position',mm*[x0, z0,Lx,Lz],'LineWidth', lw), hold on;
    % % BACKGROUND REGION METRICS
    % % LEFT
    % rectangle('Position', mm*[xb1, z0b,Lx,Lzb],'EdgeColor','k','LineWidth',lw,'LineStyle','-');
    % % RIGHT
    % rectangle('Position', mm*[xb2, z0b,Lx,Lzb],'EdgeColor','k','LineWidth',lw,'LineStyle','-');
    % 

    colormap ("turbo");
    title (['f = ', num2str(freq),'Hz'])
    ylim(zlim_mm)
    % xlim(zlim_mm)
    if ii==1
        ylabel('Axial [mm]'); 
        text(-35,22.5,'\bfPG', 'Rotation',90, ...
        'HorizontalAlignment', 'center', 'FontSize', fontText)
    end    
    % text(0,10,sprintf('\\bfCNR = %.2f', T.cnr(ii)), ...
    %     'HorizontalAlignment', 'center', 'FontSize', fontText)
    % text(0,34, sprintf('\\bfSWS_{in} = %.2f \\pm %.2f', ...
    %     T.mean_inc(ii), T.std_inc(ii)),...
    %     'HorizontalAlignment', 'center', 'FontSize', fontText)
    % text(0,38, sprintf('\\bfSWS_{bg} = %.2f \\pm %.2f', ...
    %     T.mean_back(ii), T.std_back(ii)),...
    %     'HorizontalAlignment', 'center', 'FontSize', fontText)
    set(gca, 'FontSize',fontSize)
    % axis off

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

    % % INCLUSION REGION METRICS
    % rectangle('Position',mm*[x0, z0,Lx,Lz],'LineWidth', lw), hold on;
    % % BACKGROUND REGION METRICS
    % % LEFT
    % rectangle('Position', mm*[xb1, z0b,Lx,Lzb],'EdgeColor','k','LineWidth',lw,'LineStyle','-');
    % % RIGHT
    % rectangle('Position', mm*[xb2, z0b,Lx,Lzb],'EdgeColor','k','LineWidth',lw,'LineStyle','-');

    colormap ("turbo");
    % xlabel('Lateral [mm]'),
    if ii==1
        ylabel('Axial [mm]'); 
        text(-35,22.5,'\bfPG-TV', 'Rotation',90, ...
        'HorizontalAlignment', 'center', 'FontSize', fontText)
    end    % title (['PG-TV, f = ', num2str(freq),'Hz'])
    ylim(zlim_mm)
    % text(0,10,sprintf('\\bfCNR = %.2f', T.cnr(ii+3)), ...
    %     'HorizontalAlignment', 'center', 'FontSize', fontText)
    % text(0,34, sprintf('\\bfSWS_{in} = %.2f \\pm %.2f', ...
    %     T.mean_inc(ii+3), T.std_inc(ii+3)),...
    %     'HorizontalAlignment', 'center', 'FontSize', fontText)
    % text(0,38, sprintf('\\bfSWS_{bg} = %.2f \\pm %.2f', ...
    %     T.mean_back(ii+3), T.std_back(ii+3)),...
    %     'HorizontalAlignment', 'center', 'FontSize', fontText)
    set(gca, 'FontSize',fontSize)
    % axis off
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

    % % INCLUSION REGION METRICS
    % rectangle('Position',mm*[x0, z0,Lx,Lz],'LineWidth', lw), hold on;
    % % BACKGROUND REGION METRICS
    % % LEFT
    % rectangle('Position', mm*[xb1, z0b,Lx,Lzb],'EdgeColor','k','LineWidth',lw,'LineStyle','-');
    % % RIGHT
    % rectangle('Position', mm*[xb2, z0b,Lx,Lzb],'EdgeColor','k','LineWidth',lw,'LineStyle','-');

    colormap ("turbo");
    xlabel('Lateral [mm]'),
    if ii==1
        ylabel('Axial [mm]'); 
        text(-35,22.5,'\bfPG-TNV', 'Rotation',90, ...
        'HorizontalAlignment', 'center', 'FontSize', fontText)
    end
    % title (['PG-TNV, f = ', num2str(freq),'Hz'])

    ylim(zlim_mm)
    % text(0,10,sprintf('\\bfCNR = %.2f', T.cnr(ii+6)), ...
    %     'HorizontalAlignment', 'center', 'FontSize', fontText)
    % text(0,34, sprintf('\\bfSWS_{in} = %.2f \\pm %.2f', ...
    %     T.mean_inc(ii+6), T.std_inc(ii+6)),...
    %     'HorizontalAlignment', 'center', 'FontSize', fontText)
    % text(0,38, sprintf('\\bfSWS_{bg} = %.2f \\pm %.2f', ...
    %     T.mean_back(ii+6), T.std_back(ii+6)),...
    %     'HorizontalAlignment', 'center', 'FontSize', fontText)

    set(gca, 'FontSize',fontSize)
    % axis off
end
c = colorbar;
c.Label.String = 'SWS [m/s]';

% saveas(gcf,'./figs/SWSphantoms.png')


%% Plots
lw = 1.5;

figure,
hold on
errorbar(v_freq_best,T.mean_inc(1:3),T.std_inc(1:3), 'ro-', ...
    "MarkerFaceColor",'r', 'LineWidth', lw)
errorbar(v_freq_best,T.mean_inc(4:6),T.std_inc(4:6), 'bs-', ...
    "MarkerFaceColor",'b', 'LineWidth', lw)
errorbar(v_freq_best,T.mean_inc(7:9),T.std_inc(7:9), 'k*-', ...
    "MarkerFaceColor",'k', 'LineWidth', lw)
yline(gtInc, 'k--', 'GT_{in}'); % 'k--' denotes a black dashed line
% hold off
% legend('PG','PG-TV','PG-TNV', 'Location','northwest')
% grid on
% xlim([300 1000])
% ylim([2.3,4.8])
% title('SWS dispersion inclusion')
% ylabel('SWS [m/s]'), xlabel('Frequency [Hz]')
% 
% figure,
% hold on
errorbar(v_freq_best,T.mean_back(1:3),T.std_back(1:3), 'ro-', ...
    "MarkerFaceColor",'r', 'LineWidth', lw)
errorbar(v_freq_best,T.mean_back(4:6),T.std_back(4:6),  'bs-', ...
    "MarkerFaceColor",'b', 'LineWidth', lw)
errorbar(v_freq_best,T.mean_back(7:9),T.std_back(7:9), 'k*-', ...
    "MarkerFaceColor",'k', 'LineWidth', lw)
yline(gtBack, 'k--', 'GT_{bg}'); % 'k--' denotes a black dashed line
hold off
legend('PG','PG-TV','PG-TNV', 'Location','northwest')
grid on
xlim([300 1000])
ylim([1.8,4.8])
title('SWS dispersion')
ylabel('SWS [m/s]'), xlabel('Frequency [Hz]')

figure,
hold on
plot(v_freq_best,T.cnr(1:3), 'ro-', ...
    "MarkerFaceColor",'r', 'LineWidth', lw)
plot(v_freq_best,T.cnr(4:6),'bs-', ...
    "MarkerFaceColor",'b', 'LineWidth', lw)
plot(v_freq_best,T.cnr(7:9), 'k*-', ...
    "MarkerFaceColor",'k', 'LineWidth', lw)

hold off
legend('PG','PG-TV','PG-TNV', 'Location','northwest')
grid on
xlim([300 1000])
ylim([1 7])
title('CNR')
xlabel('Frequency [Hz]')

save_all_figures_to_directory('figs','fig','.svg')

%%
writetable(T,'inVivo.xlsx')
%% =====================================================================
% %% TOTAL NUCLEAR VARIATION grid search
% M = my_obj.size_out(1);
% N = my_obj.size_out(2);
% 
% % tau_vector = [0.001];
% % mu_vector = logspace(log10(0.1),log10(10^5),10); % logarithmic
% mu_vector = 10.^(3:0.5:4);
% tau_vector = 10.^(-3.5:0.5:-2.5);
% maxIter = 1000;
% stableIter = 50;
% tol = 10e-4; % tolerance error
% 
% weightEstimators = ones(1, length(v_freq));
% 
% 
% % van quedando mu = 10^3 10^3.67
% 
% % mu_vector = 10^3.5;
% % tau_vector = [0.1 0.05 0.01 0.005 0.001 0.0005 0.0001];
% 
% 
% list_data = [10 11 13];
% v_freq = [400 600 900];
% 
% v_freq_best = v_freq;
% list_data_best = list_data;
% numChannels = length(v_freq_best);
% 
% for u = 1:length(mu_vector)
%     bestmu = mu_vector(u);
%     for t = 1:length(tau_vector)
%         besttau = tau_vector(t);
% 
%         clear tnv
%         [tnv.grad_abs_3D, cost, error, fid, reg] = pdo_den_wtnv(og.grad_abs_3D, bestmu, besttau, maxIter, tol, stableIter, weightEstimators); 
% 
% %         for ii = 1 : numChannels
% %         
% %             freq = v_freq_best(ii);
% %         
% %             tnv.sws_abs_3D(:,:, ii) = (2*pi*freq)./tnv.grad_abs_3D(:,:,ii);
% %         
% %         end
% 
%         tnv.sws_abs_3D =  2*pi* reshape(v_freq_best, [1, 1, numChannels]) ./ tnv.grad_abs_3D; % more elegant
% 
%         caxis_sws = [1 4];
%         figure, % SWS
%         sgtitle(['SWS TNV, \mu=10^{', num2str(log10(bestmu)), '} \tau=10^{', num2str(log10(besttau)),'}']);
%         set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
%         for ii = 1 : numChannels
%             freq = v_freq_best(ii);
%             subplot (2, 3, ii)
%             imagesc(tnv.sws_abs_3D(:,:,ii), caxis_sws), axis('tight')
% %             colormap ("jet");
%             colormap ("turbo");
%             xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
%             title (['SWS f_v = ', num2str(freq) ])
% 
%         end
% 
%     end
% end



% %% PLOT TNV grid search
% 
% caxis_sws = [0 5];
% figure, % SWS
% sgtitle('SWS TNV')
% set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
% for ii = 1 : numChannels
%     freq = v_freq_best(ii);
%     subplot (2, 3, ii)
%     imagesc(tnv.sws_abs_3D(:,:,ii), caxis_sws), axis('tight')
%     colormap ("jet");
% %     colormap ("turbo");
%     xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
%     title (['SWS f_v = ', num2str(freq) ])
% 
% end
% 
% figure, % grad phi
% sgtitle('|\nabla\phi| TNV')
% for ii = 1 : numChannels
%     freq = v_freq_best(ii);
%     subplot (2, 3, ii)
%     imagesc(tnv.grad_abs_3D(:,:,ii))
%     colormap ("turbo");
%     xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
%     title (['\nabla\phi f_v = ', num2str(freq) ])
% 
% end
