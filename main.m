% SCRIP TO GENERATE SWS MAPS OF HOMOGENEOUS REVERBERANT FIELDS
% Creation: 26/03/2024 (EMZ)
setup,
% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography' ...
%     '\reverberant\new'];
baseDir = 'P:\smerino\reverberant\simu_inc';
dataDir = fullfile(baseDir,'raw');
resultsDir = fullfile(baseDir,'results');
[~,~,~]= mkdir(dataDir);
[~,~,~]= mkdir(resultsDir);

%% Setting parameters
% Kernel size and step
kernel_lsq = [15,15];
kernel_tv = [15,15];

% constant = 0.3;
constant_LSQ = 0.3;
constant_TV = 0.7; %0.7;

stride = 3;

% Plotting const
cm = 1e2;
sws_range = [2,5];

% file management
files = dir(fullfile(dataDir,'*.mat'));
grad_abs_all = zeros(167,167,length(files));
% grad_abs_all = zeros(501,501,length(files));
%% Looping each file

for iSim = 1:length(files) 
    % loading data
    data = load(fullfile(dataDir,files(iSim).name));
    u = data.pv_complexZ(:,:,1);
    dinf.dx = data.x(2) - data.x(1);
    dinf.dz = data.z(2) - data.z(1);
    f_v = data.freq;
    x = data.x; z = data.z;
    og_size = size(u);
    %% PG-LSQ + MED
    extended_u = padarray(u,(kernel_lsq-1)/2,'symmetric');

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
    phase_grad_2 = (grad_x.^2 + grad_z.^2)/constant_LSQ;
    wl_samples = floor (2.5/f_v/dinf.dx/stride(1))*2+1;
    k2_med = medfilt2(phase_grad_2,[wl_samples wl_samples],'symmetric');
    k = sqrt(k2_med);
    sws_lsq = (2*pi*f_v)./k;
    %% PG-TV
    extended_u = padarray(u,(kernel_tv-1)/2,'symmetric');
%     overlap = 0.8;
%     wl_samples = floor (2.5/f_v/dinf.dx/stride(1))*2+1;
%     stride_tv = round((1-overlap)*wl_samples);
    stride_tv = stride;
    [grad_abs, size_out] = pg_norm(extended_u, kernel_tv, dinf, og_size, stride_tv);
    grad_abs = grad_abs/sqrt(constant_TV);
    grad_abs_all(:,:,iSim) = grad_abs;
    M = size_out(1); N = size_out(2);
    mu = 10.^4.5; 
    
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
    [X,Z] = meshgrid(xSWS,zSWS);
    a = 10e-3/sqrt(2);
    roi_z = abs(Z)<a;
    roi_inc = abs(X)<a & roi_z;
    roi_back = (abs(X)>2*a & abs(X)<3*a) & roi_z;
    inc_ideal = sqrt(X.^2+Z.^2)<10e-3;

    figure('Units','centimeters', 'Position',[5 5 30 10]),
    tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact')
    nexttile,
    imagesc(x*cm, z*cm, real(u), 2*std(u(:))*[-1 1]);
    colormap(gca, parula); 
    axis image;
    title(['PV, freq=',num2str(f_v)])
    xlabel('Lateral [cm]'), ylabel('Axial [cm]'),
        hold on
    contour(xSWS*cm, zSWS*cm, roi_inc, 'k--')
    contour(xSWS*cm, zSWS*cm, roi_back, 'k--')
    hold off
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
    xSWS = x(1:stride:end);
    zSWS = z(1:stride:end);
    [X,Z] = meshgrid(xSWS,zSWS);
    a = 10e-3/sqrt(2);
    roi_z = abs(Z)<a;
    roi_inc = abs(X)<a & roi_z;
    roi_back = (abs(X)>2*a & abs(X)<3*a) & roi_z;
    resultsPG(iSim) = get_metrics(sws_lsq,roi_inc,roi_back,'PG',f_v);

    xSWS = x(1:stride_tv:end);
    zSWS = z(1:stride_tv:end);
    [X,Z] = meshgrid(xSWS,zSWS);
    a = 10e-3/sqrt(2);
    roi_z = abs(Z)<a;
    roi_inc = abs(X)<a & roi_z;
    roi_back = (abs(X)>2*a & abs(X)<3*a) & roi_z;
    resultsPGTV(iSim) = get_metrics(sws_tv,roi_inc,roi_back,'PGTV',f_v);
end
%%
results1 = struct2table(resultsPG);
results2 = struct2table(resultsPGTV);

T = [results1;results2];
writetable(T,fullfile(resultsDir,'results.xlsx'),'WriteRowNames',true);
save_all_figures_to_directory(resultsDir,'fig')
close all

%%
function metrics = get_metrics(img,roi_inc,roi_back,method,freq)
    metrics.mean_inc = mean(img(roi_inc));
    metrics.mean_back =  mean(img(roi_back));
    metrics.std_inc = std(img(roi_inc));
    metrics.std_back = std(img(roi_back));
    metrics.cnr = abs(metrics.mean_inc - metrics.mean_back)/...
        sqrt(metrics.std_inc.^2 + metrics.std_back.^2);
    metrics.method = method;
    metrics.freq = freq;
end