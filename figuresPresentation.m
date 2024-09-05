% SCRIP TO GENERATE SWS MAPS OF HOMOGENEOUS REVERBERANT FIELDS
% Creation: 26/03/2024 (EMZ)
setup,
baseDir = ['C:\Users\sebas\Documents\MATLAB\DataProCiencia\Elastrography' ...
    '\reverberant\new'];
% baseDir = 'P:\smerino\reverberant\simu_inc';
dataDir = fullfile(baseDir,'sim');
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

iSim = 3;
% loading data
data = load(fullfile(dataDir,files(iSim).name));
u = data.pv_complexZ(:,:,1);
dinf.dx = data.x(2) - data.x(1);
dinf.dz = data.z(2) - data.z(1);
f_v = data.freq;
x = data.x; z = data.z;
og_size = size(u);

[X,Z] = meshgrid(x,z);
inc = sqrt(X.^2 + Z.^2) < 10e-3;
idealSws = data.c_back*ones(og_size);
idealSws(inc) = data.c_inc;

figure('Position', [100 100 500 200]),
tiledlayout(1,2)
nexttile,
imagesc(x*1e3,z*1e3,idealSws, [2 5])
axis image
xlabel('Lateral [mm]')
ylabel('Axial [mm]')
title('SWS [m/s]')
colormap turbo
colorbar

t2 = nexttile;
imagesc(x*1e3,z*1e3,real(u))
axis image
xlabel('Lateral [mm]')
ylabel('Axial [mm]')
title('Particle velocity')
colormap(t2,parula)

saveas(gcf,'simulated.png')
