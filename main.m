clear,
baseDir = 'P:\rswe\dataold\Data800Hz-10000ondas';
w_kernel = [15,15];
constant = 0.33;
stride = 3;

data = load(fullfile(baseDir,'R-FIELD_inc_1.mat'));
u = data.pv_complexZ;
dinf.dx = data.x(2) - data.x(1);
dinf.dz = data.z(2) - data.z(1);
f_v = data.freq;
x = data.x; z = data.z;
og_size = size(u);
extended_u = padarray(u,[(w_kernel(1)-1)/2 (w_kernel(2)-1)/2],'symmetric');
%%
tic
[Ax_large, Az_large, bx_large, bz_large, size_out] = getmat_pg_v2(extended_u,...
    w_kernel, dinf, og_size, stride);
toc

tic
results_x = minres(Ax_large'*Ax_large,Ax_large'*bx_large);
results_z = minres(Az_large'*Az_large,Az_large'*bz_large);
toc

grad_x = reshape(results_x(1:prod(size_out)),size_out);
grad_z = reshape(results_z(prod(size_out)+1:2*prod(size_out)),size_out);


cm = 1e2;
figure,
tiledlayout(1,2)
nexttile,
imagesc(x*cm, z*cm, abs(grad_x));
colormap(gca, jet); % Apply jet colormap to the current axes
colorbar;
axis image;
title('k_x')
xlabel('Lateral'), ylabel('Axial'),

nexttile,
imagesc(x*cm, z*cm, abs(grad_z));
colormap(gca, jet); % Apply jet colormap to the current axes
colorbar;
axis image;
title('k_z')
xlabel('Lateral'), ylabel('Axial'),

