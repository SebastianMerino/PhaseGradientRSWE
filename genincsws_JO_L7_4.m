% NEW DATA 3 PUCP L7-4 
% DATA 10 400-500
% DATA 11 600-700
% DATA 13 900-1000
% DATA 14 900-1000



%% GENERATE DATA SWS INCLUSION ONE FRAME WITH PHASE ESTIMATOR QR SOLVER kernel V1.0

addpath(genpath(pwd));
fprintf('-------QR solver 1.0 (kernel) -------\n')

list_data = [10 11 13 14];
v_freq = [400 600 900 900];

window = 15; %11 pixels as described in paper
stride = 1;


w_kernel = [window, window];
tic;

pathdata = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\reverberant\CIRS_phantom\L7-4';
pathout = fullfile(pathdata,'AromCode');

if ~exist("pathout","dir"); mkdir(pathout); end

preproc = 'new';

for ii = 1: length(v_freq)
   
    freq = v_freq(ii);
    sampleCase = list_data(ii);
    pathfreq_in = fullfile(pathdata,['data_', num2str(sampleCase)]);
    pathfreq_out = fullfile(pathout,['data_', num2str(sampleCase)]);

    if ~exist(pathfreq_out,"dir"); mkdir(pathfreq_out); end

    name = ['data_',num2str(sampleCase),'.mat'];
    structdata = load(fullfile(pathfreq_in, name));

    dinf.dx = structdata.dinf.dx;
    dinf.dz = structdata.dinf.dx;

    cMin = 1; opts.cMax = 5; 
    opts.typeFilter = 'bpf';

    if strcmp(preproc, 'new')
        uFilt = preproc_partDisplacement(structdata.u, freq, cMin, dinf, opts);
        frame = uFilt;
    end
    if strcmp(preproc, 'old')
        [uFrame,uFilt] = getComplexFrame(structdata.u, freq, cMin, dinf);
        frame = uFrame;
    end
    
    og_size = size(frame);
    mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');

    % FUNCION SWS MAP (Con linearizacion)
    [grad_z,grad_x,k,sws_matrix] = phase_estimator_QR_kernel(mirror_frame, w_kernel, freq, dinf, og_size, 1, stride);

    % EMPAQUETAR RESULTADOS
    if strcmp(preproc, 'new')
        pg_QRv1.sws_abs_3D(:,:,ii) = sws_matrix;
        pg_QRv1.grad_abs_3D(:,:,ii) = sqrt(grad_z.^2 + grad_x.^2);    
        pg_QRv1.k(:,:,ii) = k;
    end

    if strcmp(preproc, 'old')
        pg_QRv1.sws_abs_3D(:,:,ii) = sws_matrix;
        pg_QRv1.grad_abs_3D(:,:,ii) = sqrt(grad_z.^2 + grad_x.^2);    
        pg_QRv1.k(:,:,ii) = k;
    end

%         og_size = size(frame);
%         mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
%     
%         % FUNCION SWS MAP (Con linearizacion)
%         [grad_z,grad_x,k,sws_matrix] = phase_estimator_QR_kernel(mirror_frame, w_kernel, freq, dinf, og_size, 0.33, stride);
%     
%         % EMPAQUETAR RESULTADOS
%         
%         pg_QRv1.grad_z = grad_z;
%         pg_QRv1.grad_x = grad_x;
%         pg_QRv1.grad_k = k;
%         pg_QRv1.sws_matrix = sws_matrix;
%        
        
    % Save
%         pathQR = [pathfreq_out, 'PhaseGradientQRv1/'];
%         if ~exist(pathQR,"dir"); mkdir(pathQR); end
%         save([pathQR, 'SWS_PG_QRv1_jo_',num2str(sampleCase),'.mat'],'pg_QRv1');


end
toc
fprintf('---------------------------\n')

%%
%% PLOT QR v1 KERNEL
numChannels = length(v_freq);
caxis_sws = [1 4];
figure, % SWS
sgtitle(['SWS PGvKernel w=', num2str(window(1)),', str=', num2str(stride)])
% set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.55 0.25])
for ii = 1 : numChannels
    freq = v_freq(ii);
%     subplot (2, 3, ii)
    subplot (1, 4, ii) % IUS ABSTRACT
    imagesc(pg_QRv1.sws_abs_3D(:,:,ii), caxis_sws)
        colormap ("jet");
%     colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar, axis("tight")
    title (['SWS f_v = ', num2str(freq) ])

end

%% THERORICAL FITTING 

addpath(genpath(pwd));
fprintf('-------Curve Fitting Classical -------\n')
% GENERATE DATA SWS INCLUSION ONE FRAME CURVE FITTING

% DATA 10 400-500
% DATA 11 600-700
% DATA 13 900-1000
% DATA 14 900-1000

list_data = [10 11 13 14];
v_freq = [400 600 900 900];

% lambda = cs/freq


% 3.65m/s - 400Hz -> 9.125mm
% 3.65m/s - 600Hz -> 6.08mm
% 3.65m/s - 900Hz -> 4.05mm

% Resoluc dx 0.29mm/pixel
% Resoluc dz 0.15mm/pixel

% dx = 0.29 mm/px
% 400 Hz | 2.6m/s 1l = 6.5mm = 23px | 3.65m/s 1l = 9.125mm = 33px
% 600 Hz | 2.6m/s 1l = 4.3mm = 15px | 3.65m/s 1l = 6.08mm = 21px
% 900 Hz | 2.6m/s 1l = 2.8mm = 11px | 3.65m/s 1l = 4.05mm = 15px

% dz = 0.15 mm/px
% 400 Hz | 2.6m/s 1l = 6.5mm = 43px | 3.65m/s 1l = 9.125mm = 61px
% 600 Hz | 2.6m/s 1l = 4.3mm = 29px | 3.65m/s 1l = 6.08mm = 41px
% 900 Hz | 2.6m/s 1l = 2.8mm = 19px | 3.65m/s 1l = 4.05mm = 27px


% 800Hz 2.5m/s 1l = 3.125mm = 33pix 4.5m/s 1l = 5.625mm = 57pix
% 900Hz 2.5m/s 1l = 2.77mm = 29pix 4.5m/s 1l = 5mm = 51pix
% 1000Hz 2.5m/s 1l = 2.5mmm = 27pix 4.5m/s 1l = 4.5mm = 47pix


window_v = 1*[61 41 27] + 0; %1 lamda pixels
preproc = 'old';

pathdata = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\reverberant\CIRS_phantom\L7-4';
pathout = fullfile(pathdata,'AromCode');

if ~exist("pathout","dir"); mkdir(pathout); end

cont = 1;
for ii = 1: length(v_freq)
    tic;

    window = window_v(cont);
    freq = v_freq(ii);
    sampleCase = list_data(ii);

    pathfreq_in = fullfile(pathdata,['data_', num2str(sampleCase)]);
    pathfreq_out = fullfile(pathout,['data_', num2str(sampleCase)]);

    if ~exist(pathfreq_out,"dir"); mkdir(pathfreq_out); end

    name = ['data_',num2str(sampleCase),'.mat'];
    structdata = load(fullfile(pathfreq_in, name));

    dinf.dx = structdata.dinf.dx;
    dinf.dz = structdata.dinf.dx;

    cMin = 1; opts.cMax = 5; 
    opts.typeFilter = 'bpf';

    if strcmp(preproc, 'new')
        uFilt = preproc_partDisplacement(structdata.u, freq, cMin, dinf, opts);
        frame = uFilt;
    end
    if strcmp(preproc, 'old')
        [uFrame,uFilt] = getComplexFrame(structdata.u, freq, cMin, dinf);
        frame = uFrame;
    end
    
    og_size = size(frame);
    mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
        
    % FUNCTION SWS MAP (CURVE FITTING)
    [k_z,R_ax,k_x,R_lat,k,sws_matrix] = theoretical_fitting(mirror_frame,[window window],freq,dinf,og_size);
%     k_z = 1; R_ax = 1; k_x=1; R_lat =1; k = 1; sws_matrix=1;

    % Pack results
    theoricalFitting.k_z = k_z;
    theoricalFitting.R_ax = R_ax;
    theoricalFitting.k_x = k_x;
    theoricalFitting.R_lat = R_lat;
    theoricalFitting.k = k;
    theoricalFitting.sws_matrix = sws_matrix;
    theoricalFitting.window = window;

    tf.sws_abs_3D(:,:,ii) = sws_matrix;
    tf.k(:,:,ii) = sqrt(k_z.^2 + k_x.^2);    
    
    % Save
    pathCurveFitting = [pathfreq_out, '/CurveFitting/'];
    if ~exist(pathCurveFitting,"dir"); mkdir(pathCurveFitting); end

    save([pathCurveFitting, 'SWS_TF_',num2str(sampleCase),'.mat'],'theoricalFitting');
   
    cont = cont + 1;
    t = toc;
    fprintf('Time f_v = %d Hzis %f sec\n', freq, t);
   
end

%%

%%

%% WAVE APPROXIMATIOm

addpath(genpath(pwd));
fprintf('-------WAVE APPROXIMATION -------\n')
% GENERATE DATA SWS INCLUSION ONE FRAME CURVE FITTING

% DATA 10 400-500
% DATA 11 600-700
% DATA 13 900-1000
% DATA 14 900-1000

list_data = [10 11 13 14];
v_freq = [400 600 900 900];

% lambda = cs/freq


% 3.65m/s - 400Hz -> 9.125mm
% 3.65m/s - 600Hz -> 6.08mm
% 3.65m/s - 900Hz -> 4.05mm

% Resoluc dx 0.29mm/pixel
% Resoluc dz 0.15mm/pixel

% dx = 0.29 mm/px
% 400 Hz | 2.6m/s 1l = 6.5mm = 23px | 3.65m/s 1l = 9.125mm = 33px
% 600 Hz | 2.6m/s 1l = 4.3mm = 15px | 3.65m/s 1l = 6.08mm = 21px
% 900 Hz | 2.6m/s 1l = 2.8mm = 11px | 3.65m/s 1l = 4.05mm = 15px

% dz = 0.15 mm/px
% 400 Hz | 2.6m/s 1l = 6.5mm = 43px | 3.65m/s 1l = 9.125mm = 61px
% 600 Hz | 2.6m/s 1l = 4.3mm = 29px | 3.65m/s 1l = 6.08mm = 41px
% 900 Hz | 2.6m/s 1l = 2.8mm = 19px | 3.65m/s 1l = 4.05mm = 27px


% 800Hz 2.5m/s 1l = 3.125mm = 33pix 4.5m/s 1l = 5.625mm = 57pix
% 900Hz 2.5m/s 1l = 2.77mm = 29pix 4.5m/s 1l = 5mm = 51pix
% 1000Hz 2.5m/s 1l = 2.5mmm = 27pix 4.5m/s 1l = 4.5mm = 47pix


window_v = 1*[61 41 27 27] + 0; %1 lamda pixels

shift = 7; %2pixels
est_z = 10;     % - Axial estimator (usually 10)
est_x = 5;  % - Lateral estimator (usually 5)


pathdata = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\reverberant\CIRS_phantom\L7-4';
pathout = fullfile(pathdata,'AromCode');

if ~exist("pathout","dir"); mkdir(pathout); end

cont = 4;
for ii = 4: length(v_freq)
    tic;

    window = window_v(cont);
    freq = v_freq(ii);
    sampleCase = list_data(ii);

    pathfreq_in = fullfile(pathdata,['data_', num2str(sampleCase)]);
    pathfreq_out = fullfile(pathout,['data_', num2str(sampleCase)]);

    if ~exist(pathfreq_out,"dir"); mkdir(pathfreq_out); end

    name = ['data_',num2str(sampleCase),'.mat'];
    structdata = load(fullfile(pathfreq_in, name));

    dinf.dx = structdata.dinf.dx;
    dinf.dz = structdata.dinf.dx;

    cMin = 1; opts.cMax = 5; 
    opts.typeFilter = 'bpf';
    uFilt = preproc_partDisplacement(structdata.u, freq, cMin, dinf, opts);

    frame = real(uFilt);
    og_size = size(frame);
    mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
        
    sws_matrix = sws_generator(mirror_frame,[window window],freq,shift,dinf,og_size,est_z,est_x);

    % Pack results
    wa.sws_abs_3D(:,:,ii) = sws_matrix;
    wa.shift = shift;

    % Save
    pathWaveApprox = [pathfreq_out, '/WaveApprox/'];
    if ~exist(pathWaveApprox,"dir"); mkdir(pathWaveApprox); end

    save([pathWaveApprox, 'SWS_WA_sh',num2str(shift),'.mat'],'wa');
   
    cont = cont + 1;
    t = toc;
    fprintf('Time f_v = %d Hzis %f sec\n', freq, t);
   
end

%% PLOT WAVE APPROXIMATION
numChannels = length(v_freq);
caxis_sws = [1 4];
figure, % SWS
sgtitle(['SWS WAprox shift=', num2str(wa.shift)])
% set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.55 0.25])
for ii = 1 : numChannels
    freq = v_freq(ii);
%     subplot (2, 3, ii)
    subplot (1, 4, ii) % IUS ABSTRACT
    imagesc(wa.sws_abs_3D(:,:,ii), caxis_sws)
        colormap ("jet");
%     colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar, axis("tight")
    title (['SWS f_v = ', num2str(freq) ])

end

%% DENOISING NORM PHASE GRADIENT

% GENERATE DATA SWS INCLUSION ONE FRAME WITH PHASE ESTIMATOR QR SOLVER v1.0

addpath(genpath(pwd));
fprintf('-------Validation PreProcv1.0 -------\n')

% 500Hz 2.5m/s 1l = 5mm = 51pix 4.5m/s 1l = 9mm = 91pix
% 600Hz 2.5m/s 1l = 4.17mm = 41pix 4.5m/s 1l = 7.5mm = 75pix
% 700Hz 2.5m/s 1l = 3.57mm = 35pix 4.5m/s 1l = 6.42mm = 65pix

% 800Hz 2.5m/s 1l = 3.125mm = 33pix 4.5m/s 1l = 5.625mm = 57pix
% 900Hz 2.5m/s 1l = 2.77mm = 29pix 4.5m/s 1l = 5mm = 51pix
% 1000Hz 2.5m/s 1l = 2.5mmm = 27pix 4.5m/s 1l = 4.5mm = 47pix


list_data = [10 11 13 14];
v_freq = [400 600 900 900];

window = 15; %11 pixels as described in paper
stride = 1;


w_kernel = [window, window];
tic;


pathdata = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\reverberant\CIRS_phantom\L7-4';
pathout = fullfile(pathdata,'AromCode');

if ~exist("pathout","dir"); mkdir(pathout); end
contFig = 1;


for ii = 1: length(v_freq)
   
    freq = v_freq(ii);
    sampleCase = list_data(ii);
    pathfreq_in = fullfile(pathdata,['data_', num2str(sampleCase)]);
    pathfreq_out = fullfile(pathout,['data_', num2str(sampleCase)]);

    if ~exist(pathfreq_out,"dir"); mkdir(pathfreq_out); end

        name = ['data_',num2str(sampleCase),'.mat'];
        structdata = load(fullfile(pathfreq_in, name));

        dinf.dx = structdata.dinf.dx;
        dinf.dz = structdata.dinf.dx;
        dinf.PRFe = structdata.dinf.PRFe;
  
        cMin = 1; opts.cMax = 5; 
        opts.typeFilter = 'bpf';
        uFilt = preproc_partDisplacement(structdata.u,freq,cMin,dinf,opts);
        magnitude = abs(uFilt);
        phase = angle(uFilt);
   
        figure(contFig), 
        set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
        sgtitle(['Data', num2str(sampleCase), ' f_v=', num2str(freq)])
        subplot(121), imagesc(phase), title('Phase'), axis("tight")
        subplot(122), imagesc(magnitude), title('Magn'), axis("tight")
        
  
        frame = uFilt;
        og_size = size(frame);
        mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
    
        % FUNCION SWS MAP (Con linearizacion)
        [grad_abs, size_out] = pg_norm(mirror_frame, w_kernel, dinf, og_size, stride);
       
        sws_abs = (2*pi*freq)./grad_abs;
        
        figure(contFig+1), 
        set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
        sgtitle(['Data', num2str(sampleCase), ' Original, f_v=', num2str(freq)])
        subplot(121), imagesc(grad_abs), colorbar, title('|\nabla\phi|'), axis("tight")
        subplot(122), imagesc(sws_abs, sws_range), colormap ("turbo"), colorbar, title('sws'), axis("tight")

        mu = 10^4;
        M = size_out(1); N = size_out(2);
        [b_opt] = IRLS_TV(grad_abs(:),speye(M*N),mu,M,N,1e-4,ones(size(M*N)),ones(M*N,1));

        grad_abs_opt = reshape(b_opt, size(grad_abs));
        sws_abs_opt = (2*pi*freq)./grad_abs_opt;

        figure(contFig+2),
        set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
        sgtitle(['Data', num2str(sampleCase), ' TV, f_v=', num2str(freq)])
        subplot(121), imagesc(grad_abs_opt), colorbar, title('|\nabla\phi_{opt}|'), axis("tight")
        subplot(122), imagesc(sws_abs_opt, sws_range), colormap ("turbo"), colorbar, title('sws_{opt}'), axis("tight")

        contFig = contFig + 3;
        
        folderFigout = 'fig_uPreprocv1';
        pathFigout = ['./DATA_3_PUCP/CIRS_PHANTOM/L7-4/', folderFigout];
        if ~exist(pathFigout,"dir"); mkdir(pathFigout); end
        save_all_figures_to_directory(pathFigout);

%         close all;
end



toc
fprintf('---------------------------\n')

%%

%% DENOISING NORM PHASE GRADIENT

% GENERATE DATA SWS INCLUSION ONE FRAME WITH PHASE ESTIMATOR QR SOLVER v1.0

addpath(genpath(pwd));
fprintf('-------Validation PhaseExtraction vJ.O. -------\n')

list_data = [10 11 13 14];
v_freq = [400 600 900 900];

window = 15; %11 pixels as described in paper
stride = 1;


w_kernel = [window, window];
sws_range = [0, 4];
tic;


pathdata = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\reverberant\CIRS_phantom\L7-4';
pathout = fullfile(pathdata,'AromCode');

if ~exist("pathout","dir"); mkdir(pathout); end
contFig = 1;

phaseExtrac = 'JO';

for ii = 1: length(v_freq)
   
    freq = v_freq(ii);
    sampleCase = list_data(ii);
    pathfreq_in = fullfile(pathdata,['data_', num2str(sampleCase)]);
    pathfreq_out = fullfile(pathout,['data_', num2str(sampleCase)]);

    if ~exist(pathfreq_out,"dir"); mkdir(pathfreq_out); end

        name = ['data_',num2str(sampleCase),'.mat'];
        structdata = load(fullfile(pathfreq_in, name));

        dinf.dx = structdata.dinf.dx;
        dinf.dz = structdata.dinf.dx;
        dinf.PRFe = structdata.dinf.PRFe;
  
        cMin = 1; opts.cMax = 5; 
        opts.typeFilter = 'bpf';

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

        

        magnitude = abs(u_out);
        phase = angle(u_out);
   
        figure(contFig), 
        set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
        sgtitle(['Data', num2str(sampleCase), ' f_v=', num2str(freq)])
        subplot(121), imagesc(phase), title('Phase'), axis("tight")
        subplot(122), imagesc(magnitude), title('Magn'), axis("tight"), colormap('default')
        
  
        frame = u_out;
        og_size = size(frame);
        mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
    
        % FUNCION SWS MAP (Con linearizacion)
        [grad_abs, size_out] = pg_norm(mirror_frame, w_kernel, dinf, og_size, stride);
       
        sws_abs = (2*pi*freq)./grad_abs;
        
        figure(contFig+1), 
        set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
        sgtitle(['Data', num2str(sampleCase), ' Original, f_v=', num2str(freq)])
        subplot(121), imagesc(grad_abs), colorbar, title('|\nabla\phi|'), axis("tight")
        subplot(122), imagesc(sws_abs, sws_range), colormap ("turbo"), colorbar, title('sws'), axis("tight")

        mu = 10^4;
        M = size_out(1); N = size_out(2);
        [b_opt] = IRLS_TV(grad_abs(:),speye(M*N),mu,M,N,1e-4,ones(size(M*N)),ones(M*N,1));

        grad_abs_opt = reshape(b_opt, size(grad_abs));
        sws_abs_opt = (2*pi*freq)./grad_abs_opt;

        figure(contFig+2),
        set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.45 0.4])
        sgtitle(['Data', num2str(sampleCase), ' TV, f_v=', num2str(freq)])
        subplot(121), imagesc(grad_abs_opt), colorbar, title('|\nabla\phi_{opt}|'), axis("tight")
        subplot(122), imagesc(sws_abs_opt, sws_range), colormap ("turbo"), colorbar, title('sws_{opt}'), axis("tight")

        contFig = contFig + 3;
        

        pathFigout = ['./DATA_3_PUCP/CIRS_PHANTOM/L7-4/', folderFigout];
        if ~exist(pathFigout,"dir"); mkdir(pathFigout); end
        save_all_figures_to_directory(pathFigout);

%         close all;
end



toc
fprintf('---------------------------\n')

%%
%% GENERATE DATA SWS INCLUSION ONE FRAME WITH PHASE ESTIMATOR QR SOLVER v1.0

addpath(genpath(pwd));
fprintf('-------QR solver Kernel 1.0-------\n')

list_data = [10 11 13 14];
v_freq = [400 600 900 900];

window = 11;
w_kernel = [window, window];
stride = 1;
tic;


pathdata = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\reverberant\CIRS_phantom\L7-4';
pathout = fullfile(pathdata,'AromCode');

if ~exist("pathout","dir"); mkdir(pathout); end


% all channels
v_freq_best = v_freq;
list_data_best = list_data;

phaseExtrac = 'JO';
for ii = 1: length(v_freq_best)
   
    freq = v_freq_best(ii);
    sampleCase = list_data_best(ii);
    pathfreq_in = '';
    pathfreq_out = [pathout, num2str(freq),'Hz/'];

    if ~exist(pathfreq_out,"dir"); mkdir(pathfreq_out); end

        name = ['data',num2str(sampleCase),'.mat'];
        structdata = load(fullfile(pathfreq_in, name));

        dinf.dx = structdata.dinf.dx;
        dinf.dz = structdata.dinf.dx;
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
 
        frame = u_out;
        og_size = size(frame);
        mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
    
        % FUNCION SWS MAP (Con linearizacion)
        [grad_z,grad_x,k,sws_matrix] = phase_estimator_QR_kernel(mirror_frame, w_kernel, freq, dinf, og_size, 1, stride);
    
        % EMPAQUETAR RESULTADOS
        
        pg_QRv1.sws_abs_3D(:,:,ii) = sws_matrix;
        pg_QRv1.grad_abs_3D(:,:,ii) = sqrt(grad_z.^2 + grad_x.^2);    
        pg_QRv1.k(:,:,ii) = k;
    
end
toc
fprintf('---------------------------\n')

%% PLOT QR v1 KERNEL
numChannels = length(v_freq_best);
caxis_sws = [1 6];
figure, % SWS
sgtitle(['SWS PGvKernel w=', num2str(window(1)),', str=', num2str(stride)])
% set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.55 0.25])
for ii = 1 : numChannels
    freq = v_freq_best(ii);
%     subplot (2, 3, ii)
    subplot (1, 4, ii) % IUS ABSTRACT
    imagesc(pg_QRv1.sws_abs_3D(:,:,ii), caxis_sws)
%         colormap ("jet");
    colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar, axis("tight")
    title (['SWS f_v = ', num2str(freq) ])

end
%%
figure, % grad phi
sgtitle('|\nabla\phi|')
for ii = 1 : numChannels
    freq = v_freq_best(ii);
    subplot (2, 3, ii)
    imagesc(pg_QRv1.k(:,:,ii))
    colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['\nabla\phi f_v = ', num2str(freq) ])

end


%% ======================================================================
%% ======================================================================
%% ======================================================================
%% GENERATE NORM PHASE GRADIENT STORE MATRIX

addpath(genpath(pwd));
fprintf('-------Norm Phase gradient -------\n')


list_data = [10 11 13];
v_freq = [400 600 900];

window = 15; %11 pixels as described in paper
stride = 1;
w_kernel = [window, window];

phaseExtrac = 'JO';
tic;

pathdata = 'C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography\reverberant\CIRS_phantom\L7-4';
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
    medf.sws_abs_3D(:,:,ii) = 2*pi*freq ./ medf.grad_abs_3D (:, :, ii);

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

    tv.sws_abs_3D(:,:,ii) = (2*pi*freq)./tv.grad_abs_3D(:,:,ii);

end

%% PLOT  TV
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

tnv.sws_abs_3D =  2*pi* reshape(v_freq_best, [1, 1, numChannels]) ./ tnv.grad_abs_3D; % more elegant


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

% Parameters to select ROI
L = 7.5;
cx = 0.5; cz = 22.5;
d = 5;

%%
[back,inc] = getRegionMasks(x*mm,z*mm,cx,cz,L,d,L);

for ii = 1 : numChannels
    freq = v_freq_best(ii);

    sws = medf.sws_abs_3D(:,:,ii);
    MetricsMedF(ii) = get_metrics(sws,inc,back,'medf',freq);

    sws = tv.sws_abs_3D(:,:,ii);
    MetricsTV(ii) = get_metrics(sws,inc,back,'tv',freq);

    sws = tnv.sws_abs_3D(:,:,ii);
    MetricsTNV(ii) = get_metrics(sws,inc,back,'tnv',freq);
end

T = [struct2table(MetricsMedF);
    struct2table(MetricsTV);
    struct2table(MetricsTNV)];
writetable(T,fullfile(pathout,'results.xlsx'),'WriteRowNames',true);
close all

%% Plots

figure,
hold on
errorbar(v_freq_best,T.mean_inc(1:3),T.std_inc(1:3), 'r')
errorbar(v_freq_best,T.mean_inc(4:6),T.std_inc(4:6), 'g')
errorbar(v_freq_best,T.mean_inc(7:9),T.std_inc(7:9), 'b')

errorbar(v_freq_best,T.mean_back(1:3),T.std_back(1:3), 'r')
errorbar(v_freq_best,T.mean_back(4:6),T.std_back(4:6), 'g')
errorbar(v_freq_best,T.mean_back(7:9),T.std_back(7:9), 'b')
hold off
legend('MedF','TV','TNV', 'Location','northwest')
grid on
xlim([300 1000])
title('SWS dispersion')
ylabel('SWS [m/s]'), xlabel('Frequency [Hz]')

figure,
hold on
plot(v_freq_best,T.cnr(1:3), 'ro-')
plot(v_freq_best,T.cnr(4:6), 'go-')
plot(v_freq_best,T.cnr(7:9), 'bo-')

hold off
legend('MedF','TV','TNV', 'Location','northwest')
grid on
xlim([300 1000])
ylim([1 7])
title('CNR')
xlabel('Frequency [Hz]')

%% Figures
% PLotting constants
zlim_mm = [5 40];
caxis_sws = [1.01 4];
fontSize = 10;
gtInc = 2.6;
gtBack = 3.65;

% Upper left corner of each background rectangle
x0 = cx - L/2; z0 = cz-L/2;
xb1 = x0 - d - L/2;
xb2 = x0 + L + d;
figure('Position',[100 100 500 450]), % SWS
tiledlayout(3,numChannels, 'TileSpacing','tight', 'Padding','compact')
% sgtitle('SWS TNV')
% set(gcf, 'units', 'Normalized', 'Position', [0 0.1 0.75 0.35])

% t1 = nexttile([3,1]);
% imagesc(xBm*mm,zBm*mm, Bmode, [-50 0])
% axis image
% % c =colorbar(t1,'westoutside');
% % colorbar
% rectangle('Position',[x0, z0,L,L], 'LineStyle','--')
% rectangle('Position',[xb1, z0,L/2,L], 'LineStyle','--')
% rectangle('Position',[xb2, z0,L/2,L], 'LineStyle','--')
% title('Bmode')
% xlabel('Lateral [mm]'), ylabel('Axial [mm]'),
% ylim(zlim_mm)
% xlim([x(1) x(end)]*mm)
% set(gca, 'FontSize',fontSize)

for ii = 1 : numChannels
    freq = v_freq_best(ii);
%     subplot (2, 3, ii)
    % subplot (1, numChannels+1, ii+1) % abstract IUS203
    nexttile;
    imagesc(x*mm,z*mm,medf.sws_abs_3D(:,:,ii), caxis_sws), axis('tight')
    axis image
    rectangle('Position',[x0, z0,L,L], 'LineStyle','--')
    rectangle('Position',[xb1, z0,L/2,L], 'LineStyle','--')
    rectangle('Position',[xb2, z0,L/2,L], 'LineStyle','--')

    colormap ("turbo");
    title (['f = ', num2str(freq),'Hz'])
    ylim(zlim_mm)
    % if ii==1; ylabel('Axial [mm]'); end
    text(0,10,['CNR = ',num2str(T.cnr(ii),2)], ...
        'HorizontalAlignment', 'center')
    biasInc = abs(T.mean_inc(ii) - gtInc)/gtInc*100;
    biasBack = abs(T.mean_back(ii) - gtBack)/gtBack*100;
    text(0,31,['Bias_{inc} = ',num2str(biasInc,2),'%'], ...
        'HorizontalAlignment', 'center')
    text(0,35,['Bias_{back} = ',num2str(biasBack,2),'%'], ...
        'HorizontalAlignment', 'center')
    if ii==1
        text(-22,22.5,'\bf PG', ...
            'HorizontalAlignment', 'center', 'Rotation',90);
    end
    axis off
    set(gca, 'FontSize',fontSize)

end
c = colorbar;
c.Label.String = 'SWS [m/s]';


for ii = 1 : numChannels
    freq = v_freq_best(ii);
%     subplot (2, 3, ii)
    % subplot (1, numChannels+1, ii+1) % abstract IUS203
    nexttile;
    imagesc(x*mm,z*mm,tv.sws_abs_3D(:,:,ii), caxis_sws), axis('tight')
    axis image
    rectangle('Position',[x0, z0,L,L], 'LineStyle','--')
    rectangle('Position',[xb1, z0,L/2,L], 'LineStyle','--')
    rectangle('Position',[xb2, z0,L/2,L], 'LineStyle','--')

    colormap ("turbo");
    % title (['PG-TV, f = ', num2str(freq),'Hz'])
    ylim(zlim_mm)
    % if ii==1; ylabel('Axial [mm]'); end
    text(0,10,['CNR = ',num2str(T.cnr(ii+3),2)], ...
        'HorizontalAlignment', 'center')
    biasInc = abs(T.mean_inc(ii+3) - gtInc)/gtInc*100;
    biasBack = abs(T.mean_back(ii+3) - gtBack)/gtBack*100;
    text(0,31,['Bias_{inc} = ',num2str(biasInc,2),'%'], ...
        'HorizontalAlignment', 'center')
    text(0,35,['Bias_{back} = ',num2str(biasBack,2),'%'], ...
        'HorizontalAlignment', 'center')
    if ii==1
        text(-22,22.5,'\bf PG-TV', ...
            'HorizontalAlignment', 'center', 'Rotation',90);
    end
    axis off
    set(gca, 'FontSize',fontSize)

end
c = colorbar;
c.Label.String = 'SWS [m/s]';
%nexttile;
%axis off

for ii = 1 : numChannels
    freq = v_freq_best(ii);
%     subplot (2, 3, ii)
    % subplot (1, numChannels+1, ii+1) % abstract IUS203
    nexttile;
    imagesc(x*mm,z*mm,tnv.sws_abs_3D(:,:,ii), caxis_sws), axis('tight')
    axis image
    rectangle('Position',[x0, z0,L,L], 'LineStyle','--')
    rectangle('Position',[xb1, z0,L/2,L], 'LineStyle','--')
    rectangle('Position',[xb2, z0,L/2,L], 'LineStyle','--')

    colormap ("turbo");
    xlabel('Lateral [mm]'),
    % if ii==1; ylabel('Axial [mm]'); end
    % title (['PG-TNV, f = ', num2str(freq),'Hz'])
    ylim(zlim_mm)
    text(0,10,['CNR = ',num2str(T.cnr(ii+6),2)], ...
        'HorizontalAlignment', 'center')
    
    biasInc = abs(T.mean_inc(ii+6) - gtInc)/gtInc*100;
    biasBack = abs(T.mean_back(ii+6) - gtBack)/gtBack*100;
    text(0,31,['Bias_{inc} = ',num2str(biasInc,2),'%'], ...
        'HorizontalAlignment', 'center')
    text(0,35,['Bias_{back} = ',num2str(biasBack,2),'%'], ...
        'HorizontalAlignment', 'center')
    if ii==1
        text(-22,22.5,'\bf PG-TNV', ...
            'HorizontalAlignment', 'center', 'Rotation',90);
    end

    axis off
    set(gca, 'FontSize',fontSize)
end
c = colorbar;
c.Label.String = 'SWS [m/s]';
% colormap(t1,gray)

%% =====================================================================
%% TOTAL NUCLEAR VARIATION grid search
M = my_obj.size_out(1);
N = my_obj.size_out(2);

% tau_vector = [0.001];
% mu_vector = logspace(log10(0.1),log10(10^5),10); % logarithmic
mu_vector = 10.^(3:0.5:4);
tau_vector = 10.^(-3.5:0.5:-2.5);
maxIter = 1000;
stableIter = 50;
tol = 10e-4; % tolerance error

weightEstimators = ones(1, length(v_freq));


% van quedando mu = 10^3 10^3.67

% mu_vector = 10^3.5;
% tau_vector = [0.1 0.05 0.01 0.005 0.001 0.0005 0.0001];


list_data = [10 11 13];
v_freq = [400 600 900];

v_freq_best = v_freq;
list_data_best = list_data;
numChannels = length(v_freq_best);

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

        tnv.sws_abs_3D =  2*pi* reshape(v_freq_best, [1, 1, numChannels]) ./ tnv.grad_abs_3D; % more elegant

        caxis_sws = [1 4];
        figure, % SWS
        sgtitle(['SWS TNV, \mu=10^{', num2str(log10(bestmu)), '} \tau=10^{', num2str(log10(besttau)),'}']);
        set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
        for ii = 1 : numChannels
            freq = v_freq_best(ii);
            subplot (2, 3, ii)
            imagesc(tnv.sws_abs_3D(:,:,ii), caxis_sws), axis('tight')
%             colormap ("jet");
            colormap ("turbo");
            xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
            title (['SWS f_v = ', num2str(freq) ])
        
        end

    end
end



%% PLOT TNV grid search

caxis_sws = [0 5];
figure, % SWS
sgtitle('SWS TNV')
set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
for ii = 1 : numChannels
    freq = v_freq_best(ii);
    subplot (2, 3, ii)
    imagesc(tnv.sws_abs_3D(:,:,ii), caxis_sws), axis('tight')
    colormap ("jet");
%     colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS f_v = ', num2str(freq) ])

end

figure, % grad phi
sgtitle('|\nabla\phi| TNV')
for ii = 1 : numChannels
    freq = v_freq_best(ii);
    subplot (2, 3, ii)
    imagesc(tnv.grad_abs_3D(:,:,ii))
    colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['\nabla\phi f_v = ', num2str(freq) ])

end
