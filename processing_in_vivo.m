setup,

fprintf('-------Norm Phase gradient -------\n')

list_data = 1:25;
v_freq = [400,500,600*ones(1, 17), 700*ones(1, 3), 800*ones(1,3)];
window = 11; %11 pixels as described in paper
stride = 1;
w_kernel = [window, window];

tic;

pathdata = 'P:\rswe\DATA_2_PUCP';
pathout = 'P:\rswe\DATA_2_PUCP\smerino';

if ~exist("pathout","dir"); mkdir(pathout); end

for ii = 1: length(v_freq)
   
    freq = v_freq(ii);
    sampleCase = list_data(ii);
    pathfreq_in = pathdata;
    pathfreq_out = fullfile(pathout,[num2str(freq),'Hz']);

    if ~exist(pathfreq_out,"dir"); mkdir(pathfreq_out); end
        name = ['data',num2str(sampleCase),'.mat'];
        structdata = load(fullfile(pathfreq_in, name));

        dinf.dx = structdata.dinf.dx;
        dinf.dz = structdata.dinf.dx;
        dinf.PRFe = structdata.dinf.PRFe;
  
        cMin = 1; opts.cMax = 5; 
        opts.typeFilter = 'bpf';
        uFilt = preproc_partDisplacement(structdata.u,freq,cMin,dinf,opts);
 
        frame = uFilt;
        og_size = size(frame);
        mirror_frame = padarray(frame,[(window-1)/2 (window-1)/2],'symmetric');
       
        [grad_abs, size_out] = pg_norm(mirror_frame, w_kernel, dinf, og_size, stride);

        % Save
        pathMat = fullfile(pathfreq_out, 'MatrixW',num2str(w_kernel(1)));
        if ~exist(pathMat,"dir"); mkdir(pathMat); end
        save(fullfile(pathMat, ...
            ['PG_abs_str',num2str(stride),'data', num2str(sampleCase), '.mat']) ...
            ,'grad_abs','size_out', 'window', 'stride', 'freq');
end
toc
fprintf('---------------------------\n')


%% ORIGINAL 3D CREATION

list_data = 1:25;
v_freq = [400,500,600*ones(1, 17), 700*ones(1, 3), 800*ones(1,3)];

v_freq_best = 400:100:800;
list_data_best = [1 2 17 20 23];

numChannels = length(v_freq_best);
stride = 1;
window = 11;

M = length(1:stride:290); % size_out(1)
N = length(1:stride:248); % size_out(2)

grad_abs_3D = zeros(M, N, numChannels); 
sws_abs_3D = grad_abs_3D;
for ii = 1 : numChannels

    freq = v_freq_best(ii);
    sample = list_data_best(ii);
    my_obj = load(['./out/JO/', num2str(freq), 'Hz/MatrixW', num2str(window), '/PG_abs_str' num2str(stride),'data',num2str(sample) ,'.mat']);
    og.grad_abs_3D (:, :, ii) = my_obj.grad_abs;     
    og.sws_abs_3D(:,:, ii) = 2*pi*freq ./ og.grad_abs_3D (:, :, ii);

end

%% TOTAL VARIATION
M = my_obj.size_out(1);
N = my_obj.size_out(2);

mu = 1e4;
for ii = 1 : numChannels
    freq = v_freq_best(ii);

    my_grad = og.grad_abs_3D (:, :, ii);

    [grad_tv] = IRLS_TV_simple(my_grad(:),speye(M*N),mu,M,N,1e-4,ones(size(M*N)),ones(M*N,1));

    tv.grad_abs_3D(:,:,ii) = reshape(grad_tv, [ M N ] );

    tv.sws_abs_3D(:,:,ii) = (2*pi*freq)./tv.grad_abs_3D(:,:,ii);

end

%% PLOT  TV
figure, % SWS
caxis_sws = [0 5];

sgtitle('SWS TV')
set(gcf, 'units', 'Normalized', 'Position', [0 0 0.55 0.55])
for ii = 1 : numChannels
    freq = v_freq_best(ii);
    subplot (2, 3, ii)
    imagesc(tv.sws_abs_3D(:,:,ii), caxis_sws)
%    colormap ("jet");
    colormap("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['SWS f_v = ', num2str(freq) ])

end

figure, % grad phi
sgtitle('|\nabla\phi| TV')
for ii = 1 : numChannels
    freq = v_freq_best(ii);
    subplot (2, 3, ii)
    imagesc(tv.grad_abs_3D(:,:,ii))
    colormap ("turbo");
    xlabel('Lateral [cm]'), ylabel('Axial[cm]'), colorbar
    title (['\nabla\phi f_v = ', num2str(freq) ])

end