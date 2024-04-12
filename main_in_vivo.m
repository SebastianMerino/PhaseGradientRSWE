% SCRIP TO test in vivio data
setup,
% baseDir = 'P:\rswe\DATA_2_PUCP';
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\' ...
    'Elastrography\reverberant\in_vivo'];
dataDir = baseDir;
matrixDir = fullfile(baseDir,'matrices');
figDir = fullfile(baseDir,'figures');

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
% iAcq = length(files) -1;
for iAcq = 1:length(files) 
    % loading data
    data = load(fullfile(dataDir,files(iAcq).name));
    u = data.u;
    dinf.dx = data.dinf.dx;
    dinf.dz = data.dinf.dz;
    f_v = data.freq;
    x = (0:size(u,2)-1)*dinf.dx; z = (0:size(u,1)-1)*dinf.dz;

    %% Constructing matrices
    tic
    og_size = size(u, [1,2]);
    [Ax_large, Az_large, bx_large, bz_large, size_out] = getmat_pg_invivo(...
        u, w_kernel, dinf, og_size, stride);
    toc
    %%
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
    k2_med = medfilt2(phase_grad_2,[med_wind med_wind],'symmetric');
    k = sqrt(k2_med);
    sws_direct = (2*pi*f_v)./k;   

    %%
    xSWS = x(1:stride:end);
    zSWS = z(1:stride:end);
    figure('Units','centimeters', 'Position',[5 5 20 10]),
    tiledlayout(1,2)
    nexttile,
    imagesc(x*cm, z*cm, real(u(:,:,1)), std(u(:))*[-1 1]);
    colormap(gca, parula); 
    axis image;
    title(['PV, freq=',num2str(f_v)])
    xlabel('Lateral [cm]'), ylabel('Axial [cm]'),
    
    nexttile,
    imagesc(xSWS*cm, zSWS*cm, sws_direct, [0 6]);
    colorbar
    colormap(gca, turbo); 
    axis image;
    title('SWS PG')
    xlabel('Lateral [cm]')

    %%
            cMin = 1; opts.cMax = 5; 
        opts.typeFilter = 'bpf';
        uFilt = preproc_partDisplacement(structdata.u,freq,cMin,dinf,opts);
 
        frame = uFilt;
        og_size = size(frame);
        mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
       
        [grad_abs, size_out] = pg_norm(mirror_frame, w_kernel, dinf, og_size, stride);

%%
    save_all_figures_to_directory(figDir,[files(iAcq).name(1:end-4),'fig']);
    close all,
end
%%
%     nexttile,
%     imagesc(xSWS*cm, zSWS*cm, sws_reg, sws_range);
%     colormap(gca, turbo); 
%     c = colorbar;
%     c.Label.String = 'SWS [m/s]';
%     axis image;
%     title('SWS PG-TV')
%     xlabel('Lateral [cm]'),
% 
%     %%
%     resultsPG(iSim) = get_metrics(sws_direct,roi_inc,roi_back,'PG',f_v);
%     resultsPGTV(iSim) = get_metrics(sws_reg,roi_inc,roi_back,'PGTV',f_v);
% end
% %%
% results1 = struct2table(resultsPG);
% results2 = struct2table(resultsPGTV);
% 
% T = [results1;results2];
% writetable(T,fullfile(resultsDir,'results.xlsx'),'WriteRowNames',true);
% save_all_figures_to_directory(resultsDir,'fig')
% 
% 
% %%
% function metrics = get_metrics(img,roi_inc,roi_back,method,freq)
%     metrics.mean_inc = mean(img(roi_inc), 'omitmissing');
%     metrics.mean_back =  mean(img(roi_back), 'omitmissing');
%     metrics.std_inc = std(img(roi_inc), 'omitmissing');
%     metrics.std_back = std(img(roi_back), 'omitmissing');
%     metrics.cnr = abs(metrics.mean_inc - metrics.mean_back)/...
%         sqrt(metrics.std_inc.^2 + metrics.std_back.^2);
%     metrics.method = method;
%     metrics.freq = freq;
% end