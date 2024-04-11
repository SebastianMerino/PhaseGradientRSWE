% SCRIP TO GENERATE SWS MAPS OF HOMOGENEOUS REVERBERANT FIELDS
% Creation: 26/03/2024 (EMZ)
setup,
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography' ...
    '\reverberant\new'];
dataDir = fullfile(baseDir,'sim');
resultsDir = fullfile(baseDir,'results');
[~,~,~]= mkdir(dataDir);
[~,~,~]= mkdir(resultsDir);
%% Generating simulations
nWaves = 1e3; % number of waves
v_freq = 300:100:1000; % [Hz]
c_back = 2.5; % background SWS [m/s] 
c_inc = 4.5; % inclusion SWS [m/s]
nFields = 1; % number of fields to generate at each freq
% vid_duration = 25e-3; % time of field generation [s]
iFile = 1;
for freq = v_freq
    for field = 1:nFields
        tic;        
        [pv_complexZ,x,z] = simulate_reverberant_complexZ_vFast(nWaves, ...
            freq, c_back, c_inc);
        
        name = [num2str(iFile),'_',num2str(field),'.mat'];
        save(fullfile(dataDir,name),...
            'pv_complexZ','x','z','freq','c_back','c_inc');
        toc
        iFile = iFile+1;
    end
end

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

for iSim = 1:length(files) 
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
    [Ax_large, Az_large, bx_large, bz_large, size_out] = getmat_pg_v4(...
        extended_u, w_kernel, dinf, og_size, stride);
    toc
    
    % Direct inversion
    tic
    results_x = minres(Ax_large'*Ax_large,Ax_large'*bx_large);
    results_z = minres(Az_large'*Az_large,Az_large'*bz_large);
    results_x = reshape(results_x,size_out);
    results_z = reshape(results_z,size_out);
    toc
    
    % Post-processing
    med_wind = floor (2.5/f_v/dinf.dx/stride(1))*2+1;
    grad_x = results_x(:,:,1);
    grad_z = results_z(:,:,1);
    phase_grad_2 = (grad_x.^2 + grad_z.^2)/constant;
    med_wind = 1;
    k2_med = medfilt2(phase_grad_2,[med_wind med_wind],'symmetric');
    k = sqrt(k2_med);
    sws_direct = (2*pi*f_v)./k;   

    %% Regularized systems
    tic
    M = size_out(1);
    N = size_out(2);
    mask = ones(size(bx_large));
    minimask = ones(M*N,1);
    results_x = IRLS_TV(bx_large,Ax_large,mu,M,N,tol,mask,minimask);
    results_z = IRLS_TV(bz_large,Az_large,mu,M,N,tol,mask,minimask);
    results_x = reshape(results_x,size_out);
    results_z = reshape(results_z,size_out);
    toc

    % Post-processing
    grad_x = results_x(:,:,1);
    grad_z = results_z(:,:,1);
    phase_grad_2 = (grad_x.^2 + grad_z.^2)/constant;
    k2_med = phase_grad_2;
    k = sqrt(k2_med);
    sws_reg = (2*pi*f_v)./k;   
    
    % Plot
    xSWS = x(1:stride:end);
    zSWS = z(1:stride:end);
    [X,Z] = meshgrid(xSWS,zSWS);
    a = 10e-3/sqrt(2);
    roi_z = abs(Z)<a;
    roi_inc = abs(X)<a & roi_z;
    roi_back = (abs(X)>2*a & abs(X)<3*a) & roi_z;
    inc_ideal = sqrt(X.^2+Z.^2)<10e-3;

    figure('Units','centimeters', 'Position',[5 5 30 10]),
    tiledlayout(1,3)
    nexttile,
    imagesc(x*cm, z*cm, real(u(:,:,1)));
    colormap(gca, parula); 
    axis image;
    title('PV')
    xlabel('Lateral [cm]'), ylabel('Axial [cm]'),
    hold on
    contour(xSWS*cm, zSWS*cm, roi_inc, 'k--')
    contour(xSWS*cm, zSWS*cm, roi_back, 'k--')
    hold off
    
    nexttile,
    imagesc(xSWS*cm, zSWS*cm, sws_direct, sws_range);
    colormap(gca, turbo); 
    axis image;
    title('SWS PG')
    xlabel('Lateral [cm]')
    nexttile,
    imagesc(xSWS*cm, zSWS*cm, sws_reg, sws_range);
    colormap(gca, turbo); 
    c = colorbar;
    c.Label.String = 'SWS [m/s]';
    axis image;
    title('SWS PG-TV')
    xlabel('Lateral [cm]'),

    %%
    resultsPG(iSim) = get_metrics(sws_direct,roi_inc,roi_back,'PG',f_v);
    resultsPGTV(iSim) = get_metrics(sws_reg,roi_inc,roi_back,'PGTV',f_v);
end
%%
results1 = struct2table(resultsPG);
results2 = struct2table(resultsPGTV);

T = [results1;results2];
writetable(T,fullfile(resultsDir,'results.xlsx'),'WriteRowNames',true);
save_all_figures_to_directory(resultsDir,'fig')


%%
function metrics = get_metrics(img,roi_inc,roi_back,method,freq)
    metrics.mean_inc = mean(img(roi_inc), 'omitmissing');
    metrics.mean_back =  mean(img(roi_back), 'omitmissing');
    metrics.std_inc = std(img(roi_inc), 'omitmissing');
    metrics.std_back = std(img(roi_back), 'omitmissing');
    metrics.cnr = abs(metrics.mean_inc - metrics.mean_back)/...
        sqrt(metrics.std_inc.^2 + metrics.std_back.^2);
    metrics.method = method;
    metrics.freq = freq;
end