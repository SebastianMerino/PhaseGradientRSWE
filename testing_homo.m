% SCRIP TO GENERATE SWS MAPS OF HOMOGENEOUS REVERBERANT FIELDS
% Creation: 26/03/2024 (EMZ)
setup,
% baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography' ...
%     '\reverberant\new_homo'];
baseDir = 'P:\smerino\reverberant\simu_homo';
dataDir = fullfile(baseDir,'raw');
resultsDir = fullfile(baseDir,'results');
[~,~,~]= mkdir(dataDir);
[~,~,~]= mkdir(resultsDir);

c_back = 3.5; % background SWS [m/s]
% %% Generating simulations
% nWaves = 1e3; % number of waves
% v_freq = 300:100:1000; % [Hz]
% nFields = 1; % number of fields to generate at each freq
% % vid_duration = 25e-3; % time of field generation [s]
% iFile = 1;
% for freq = v_freq
%     for field = 1:nFields
%         tic;
%         [pv_complexZ,x,z] = simulate_reverberant_complexZ_vhomo(nWaves, ...
%             freq, c_back);
% 
%         name = [num2str(iFile),'_',num2str(field),'.mat'];
%         save(fullfile(dataDir,name),...
%             'pv_complexZ','x','z','freq','c_back');
%         toc
%         iFile = iFile+1;
%     end
% end

%% Setting parameters
% Kernel size and step
kernel_lsq = [15,15];
kernel_tv = [15,15];

constant = 1;
stride = 3;

% Plotting const
cm = 1e2;
sws_range = [2,5];


% file management
files = dir(fullfile(dataDir,'*.mat'));
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
    fprintf('\nFrecuencia = %d Hz\n',f_v)
    %% PG-LSQ + MED
    extended_u = padarray(u,(kernel_lsq-1)/2,'symmetric');

    % Constructing matrices
    % tic
    [Ax_large, Az_large, bx_large, bz_large, size_out] = getmat_pg_v3(...
        extended_u, kernel_lsq, dinf, og_size, stride);
    % toc
    
    % Direct inversion
    %tic
    [results_x,~] = minres(Ax_large'*Ax_large,Ax_large'*bx_large);
    [results_z,~] = minres(Az_large'*Az_large,Az_large'*bz_large);
    results_x = reshape(results_x,[size_out,2]);
    results_z = reshape(results_z,[size_out,2]);
    %toc
    
    % Post-processing
    grad_x = results_x(:,:,1);
    grad_z = results_z(:,:,1);
    phase_grad_2 = (grad_x.^2 + grad_z.^2)/constant;
    wl_samples = floor (2.5/f_v/dinf.dx/stride(1))*2+1;
    k2_med = medfilt2(phase_grad_2,[wl_samples wl_samples],'symmetric');
    k = sqrt(k2_med);
    sws_lsq = (2*pi*f_v)./k;
    c_est = mean(sws_lsq(:));
    fprintf('Constante LSQ: \t%.4f\n',(c_back/c_est)^2)
    %% PG-TV
    overlap = 0.8;
    wl_samples = floor (2.5/f_v/dinf.dx/stride(1))*2+1;
    stride_tv = round((1-overlap)*wl_samples);
    kernel_tv = [wl_samples wl_samples];
    extended_u = padarray(u,(kernel_tv-1)/2,'symmetric');

    [grad_abs, size_out] = pg_norm(extended_u, kernel_tv, dinf, og_size, stride_tv);
    % grad_abs = grad_abs/sqrt(constant);
    M = size_out(1); N = size_out(2);
    mu = 10.^5; 
    
    tol = 1e-5;
    mask = ones(size(M*N));
    minimask = ones(M*N,1);
    
    
    denoised_grad = IRLS_TV(grad_abs(:),speye(M*N),mu,M,N,tol,mask,minimask);
    denoised_grad = reshape(denoised_grad, size_out);
    
    k = denoised_grad;
    sws_tv = (2*pi*f_v)./k;  
    c_est = mean(sws_tv(:));
    fprintf('Constante TV: \t%.4f\n\n',(c_back/c_est)^2)

    %% PLOT
    xSWS = x(1:stride:end);
    zSWS = z(1:stride:end);

    figure('Units','centimeters', 'Position',[5 5 30 10]),
    tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact')
    nexttile,
    imagesc(x*cm, z*cm, real(u), 2*std(u(:))*[-1 1]);
    colormap(gca, parula); 
    axis image;
    title(['PV, freq=',num2str(f_v)])
    xlabel('Lateral [cm]'), ylabel('Axial [cm]'),

%     
    nexttile,
    imagesc(xSWS*cm, zSWS*cm, sws_lsq);
    colorbar
    colormap(gca, turbo); 
    axis image;
    title('SWS PG-LSQ')
    xlabel('Lateral [cm]')

    nexttile,
    imagesc(xSWS*cm, zSWS*cm, sws_tv);
    colorbar
    colormap(gca, turbo); 
    axis image;
    title('SWS PG-TV')
    xlabel('Lateral [cm]')


end
figDir =fullfile(baseDir,'fig');
[~,~,~] = mkdir(figDir);
save_all_figures_to_directory(figDir,'fig')
