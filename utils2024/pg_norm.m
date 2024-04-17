function [grad_abs, size_out] = pg_norm(u, w_kernel, dinf, og_size, stride)
% function [Ax_large, Az_large, bx_large, bz_large]  = getmat_pg_Simp_nophi0(u, w_kernel, freq, dinf, og_size, stride)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that yields the matrixes Ax Az and vectors bx & bz of the 
% PHASE GRADIENT FRAMEWORK BY J. ORMACHEA
% 
% Inputs:  
%          u           : 2D region of interest to evaluate (previously mirror padding)
%          w_kernel    : vector with the size of the window kernel
%          f_v         : vibration frequency
%          dinf        : structure that contains the spatial resolutions
%          og_size     : vector containing the original size of the data
%          stride      : stride for window

% Outputs: 
%          Ax_large    : 
%          Az_large    : 
%          bx_large    : 
%          bz_large    : 
%          size_out    :
% Author: EMZ 4 S.Merino
%
% LATER:
%
% ------To generate a big matrix of both axial and lateral and big vector:
%
%  AA = blkdiag(Ax_large, Az_large); % faster kron(speye(2*numSubMatrices), A_small);
%  bb = [bx_large; bz_large];
%
% ------To form kz_plane kx_plane cte_plane, each one 2D matrix
%     Remember: A*results_vect = b
%     results_vect = A\b;  % could be results_x = my_TV(A, b)
%     
%     results_3D  = reshape(results_vect, [3, size_out(2), size_out(1)]); 
%     results_3D = permute(res3D_x, [3 2 1]); 
%     kx_plane = results_3D(:,:,1), kz_plane = results_3D(:,:,2), cte_plane = results_3D(:,:,3),
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    st = stride;
    res_z = dinf.dz; % Axial resolution
    res_x = dinf.dx; % Lateral resolution
    
     % Axis of the kernels
    z_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(1))*res_z; 
    x_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(2))*res_x; 
    
    %% Initializing vectors and matrixes    
    angle_u = angle(u);
    
    [X, Z] = meshgrid(x_axis,z_axis);
    Az_small = [Z(:) ones(length(x_axis)*length(z_axis),1)]; 
    Ax_small = [X(:) ones(length(x_axis)*length(z_axis),1)]; 
    [numRows, ~] = size(Az_small); % [ww, 3]

    %% HELP FOR DIMENSIONS %% % THEORY
    size_mirror = size(u); % ogsize + w_kernel - 1; %  = size (u)
    numkernels = floor( (size_mirror - w_kernel)./st + 1 ); % numKernels
    size_out = floor( (og_size - 1)./st + 1 );

    numSubMatrices = prod(numkernels);

    Az_small_nophi0 = [Z(:)];
    Ax_small_nophi0 = [X(:)];

    Az_large = kron(speye(numSubMatrices), Az_small_nophi0);
    bz_large = zeros(numSubMatrices*numRows, 1); 
    Ax_large = kron(speye(numSubMatrices), Ax_small_nophi0);
    bx_large = bz_large;
    	     
    angle_z = unwrap(angle_u,[],1);
    angle_x = unwrap(angle_u,[],2);
    cont_kernel = 1; 

    grad_abs = zeros(size_out);

    for ii = 1:st:og_size(1)
        for jj = 1:st:og_size(2) %% for faster computing pararell toolbox
            area_z = angle_z(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel
            area_x = angle_x(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel

            diff_area_z = diff(area_z, 1, 1)/res_x;
            diff_area_x = diff(area_x, 1, 2)/res_x;

            diff_area_z = [diff_area_z; diff_area_z(end, :)];  % replicate the last row
            diff_area_x = [diff_area_x, diff_area_x(:, end)];  % replicate the last column
        
            val_abs = sqrt(diff_area_z.^2 + diff_area_x.^2);
            grad_abs(ii,jj) = mean(val_abs, 'all');
%             val_abs = (diff_area_z.^2 + diff_area_x.^2);
%             grad_abs(ii,jj) = sqrt(mean(val_abs, 'all'));

        end
    end
    grad_abs = grad_abs(1:st:og_size(1),1:st:og_size(2));
    
end