% SCRIP TO test in vivio data
setup,
% baseDir = 'P:\rswe\DATA_2_PUCP';
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Elastrography\reverberant\in_vivo'];
dataDir = baseDir;
resultsDir = 'P:\rswe\DATA_2_PUCP\smerino\best_v2';
[~,~,~] = mkdir(resultsDir);

%% Setting parameters
list_data = [1 2 17 20 23];

% Kernel size and step
% w_kernel = [15 15];
kernel_lsq = [11,11];
kernel_tv = [11,11];
constant = 1;
stride = 1;

% Plotting const
cm = 1e2;
sws_range = [0,6];

% Reg constants
tol = 1e-5;


%% Looping each file
% iAcq = length(files) -1;
for iAcq = 1:length(list_data) 
    % loading data
    name = ['data',num2str(list_data(iAcq)),'.mat'];
    data = load(fullfile(dataDir,name));
    u = data.u;
    dinf.dx = data.dinf.dx;
    dinf.dz = data.dinf.dz;
    f_v = data.freq;
    clear data
    x = (0:size(u,2)-1)*dinf.dx; z = (0:size(u,1)-1)*dinf.dz;
    %% Pre-processing
    tic
    cMin = 1; opts.cMax = 5; 
    opts.typeFilter = 'bpf';
    uFilt = preproc_partDisplacement(u,f_v,cMin,dinf,opts);
    og_size = size(uFilt);
    toc

    %% Curve fitting
    wl_samples = floor (2.5/f_v/dinf.dx/2)*2+1;
    kernel_cf = [wl_samples,wl_samples];
    extended_u = padarray(uFilt,(kernel_cf-1)/2,'symmetric');
    tic
    [~,~,~,~,~,sws_cf] = theoretical_fitting(extended_u,kernel_cf,...
        f_v,dinf,og_size);
%     toc
    %% PG-LSQ + MED
    extended_u = padarray(uFilt,(kernel_lsq-1)/2,'symmetric');

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
    phase_grad_2 = (grad_x.^2 + grad_z.^2)/constant;
    wl_samples = floor (2.5/f_v/dinf.dx/stride(1))*2+1;
    k2_med = medfilt2(phase_grad_2,[wl_samples wl_samples],'symmetric');
    k = sqrt(k2_med);
    sws_lsq = (2*pi*f_v)./k;

    %% PG-TV
    extended_u = padarray(uFilt,(kernel_tv-1)/2,'symmetric');

    [grad_abs, size_out] = pg_norm(extended_u, kernel_tv, dinf, og_size, stride);
    grad_abs = grad_abs/sqrt(constant);
    M = size_out(1); N = size_out(2);
    mu = 10.^4; 
    
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
    figure('Units','centimeters', 'Position',[5 5 30 10]),
    tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact')
    nexttile,
    imagesc(x*cm, z*cm, real(uFilt), 2*std(uFilt(:))*[-1 1]);
    colormap(gca, parula); 
    axis image;
    title(['PV, freq=',num2str(f_v)])
    xlabel('Lateral [cm]'), ylabel('Axial [cm]'),

%     nexttile,
%     imagesc(xSWS*cm, zSWS*cm, sws_cf, sws_range);
%     colorbar
%     colormap(gca, turbo); 
%     axis image;
%     title('SWS PG-CF')
%     xlabel('Lateral [cm]')
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
    save_all_figures_to_directory(resultsDir,[num2str(f_v),'fig']);
    close all,
    
    % save(fullfile(resultsDir,name),'sws_cf',"sws_lsq","sws_lsq","uFilt")
    save(fullfile(resultsDir,name),"sws_lsq","sws_lsq","uFilt")

end
%%
