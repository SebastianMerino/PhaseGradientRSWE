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
%% Generating simulations
nWaves = 1e3; % number of waves
v_freq = 300:100:1000; % [Hz]
nFields = 1; % number of fields to generate at each freq
% vid_duration = 25e-3; % time of field generation [s]
iFile = 1;
for freq = v_freq
    for field = 1:nFields
        tic;
        [pv_complexZ,x,z] = simulate_reverberant_complexZ_vhomo(nWaves, ...
            freq, c_back);

        name = [num2str(iFile),'_',num2str(field),'.mat'];
        save(fullfile(dataDir,name),...
            'pv_complexZ','x','z','freq','c_back');
        toc
        iFile = iFile+1;
    end
end

%% Setting parameters
% Kernel size and step
w_kernel = [15 15];
constant = 1; % 0.28 is measured
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
% for iSim = 1:length(files)
iSim = 5;

% loading data
data = load(fullfile(dataDir,files(iSim).name));
u = data.pv_complexZ(:,:,1);
dinf.dx = data.x(2) - data.x(1);
dinf.dz = data.z(2) - data.z(1);
f_v = data.freq;
x = data.x; z = data.z;
xSWS = x(1:stride:end); zSWS = z(1:stride:end);

% Constructing matrices
tic
og_size = size(u);
extended_u = padarray(u,[(w_kernel(1)-1)/2 (w_kernel(2)-1)/2],'symmetric');
[Ax_large, Az_large, bx_large, bz_large, size_out] = getmat_pg_v4(...
    extended_u, w_kernel, dinf, og_size, stride);
toc

% Direct inversion
tic
[results_x,~] = minres(Ax_large'*Ax_large,Ax_large'*bx_large);
[results_z,~] = minres(Az_large'*Az_large,Az_large'*bz_large);
results_x = reshape(results_x,size_out);
results_z = reshape(results_z,size_out);
toc

grad_x = results_x(:,:,1);
grad_z = results_z(:,:,1);
phase_grad_2 = (grad_x.^2 + grad_z.^2)/constant;

figure, tiledlayout(3,1)
nexttile, 
histogram(results_x, 100, "BinLimits",[-5000 5000], 'Normalization','pdf')
title('Hist k_x')
ylim([0 1e-3])
nexttile, 
histogram(results_z, 100, "BinLimits",[-5000 5000], 'Normalization','pdf')
title('Hist k_z')
ylim([0 1e-3])

nexttile, 
histogram(sqrt(phase_grad_2), 100, "BinLimits",[0 7500], 'Normalization','pdf')
title('Hist |k|')
ylim([0 1e-3])

% end
%% Post-processing
% med_wind = floor (2.5/f_v/dinf.dx/stride(1))*2+1;
med_wind = 1;
k2_med = medfilt2(phase_grad_2,[med_wind med_wind],'symmetric');
k = sqrt(k2_med);
sws_direct = (2*pi*f_v)./k;   

% Plotting
figure('Units','centimeters', 'Position',[5 5 30 10]),
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, real(u(:,:,1)));
colormap(gca, parula);
axis image;
title('PV')
xlabel('Lateral [cm]'), ylabel('Axial [cm]'),
nexttile,
imagesc(xSWS*cm, zSWS*cm, sws_direct, sws_range);
colormap(gca, turbo);
axis image;
title('SWS PG')
xlabel('Lateral [cm]')
%%
c_est = mean(sws_direct(:));
disp('SWS')
disp(c_est)
disp('Constant')
disp((c_back/c_est)^2)

