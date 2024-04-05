function [grad_z,grad_x,k,sws_matrix] = phase_estimator_TV(u, w_kernel,f_v,dinf,og_size,constant, stride)
% function [grad_z,grad_x,k,sws_matrix] = phase_estimator_QR_bigmat(u, w_kernel,f_v,dinf,og_size,constant)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that yields the shear wave speed of a region with the 
% phase gradient method with QR solver MATLAB first version. 
% 
% Inputs:  
%          u           : 2D region of interest to evaluate (previously mirror padding)
%          w_kernel    : vector with the size of the window kernel
%          f_v         : vibration frequency
%          dinf        : structure that contains the spatial resolutions
%          og_size     : vector containing the original size of the data
%          lambda      : regularization coefficient

% Outputs: 
%          grad_z       : Gradient matrix for the axial direction
%          grad_x       : Gradient matrix for the lateral direction
%          k            : Total Wavenumber matrix 
%          sws_matrix   : Shear wave speed matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    

    st = stride;
    res_z = dinf.dz; % Axial resolution
    res_x = dinf.dx; % Lateral resolution
    
     % Axis of the kernels
    z_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(1))*res_z; 
    x_axis = linspace(-(w_kernel(1)-1)/2,(w_kernel(1)-1)/2,w_kernel(2))*res_x; 
    
    %% Initializing vectors and matrixes
    grad_z = zeros(og_size); 
    grad_x = grad_z; % matrixes containing the estimated wavenumber along each direction
    
    angle_u = angle(u);
    [M, N] = size(angle_u);
    %angle_u_unwrp = unwrap(angle_u, [], direction);
    %angle_u_unwrp = unwrap(angle_u, [], 2);   
    
    [X, Z] = meshgrid(x_axis,z_axis);
    A_small = [X(:) Z(:) ones(length(x_axis)*length(z_axis),1)]; 
    [numRows, numCols] = size(A_small); % [ww, 3]
%     b_small = zeros( size(w_kernel(1), w_kernel(2)) );

%     numCondA_small = cond(A_small);
%     fprintf('cond(A_small) = %f\n', numCondA_small);

    % Better pre-allocation v2.0
%    501−15+1 = 487.^2, if not mirror padding is applieds
    %% HELP FOR DIMENSIONS %% % THEORY
    size_mirror = size(u); % ogsize + w_kernel - 1; %  = size (u)
    numkernels = floor( (size_mirror - w_kernel)./st + 1 ); % numKernels
    overlap_ax1 = 1 - st/w_kernel(1); 
    overlap_la1 = 1 - st/w_kernel(2);
    size_out = floor( (og_size - 1)./st + 1 );

    numSubMatrices = prod(numkernels);

    Az_large = kron(speye(numSubMatrices), A_small);
    bz_large = zeros(numSubMatrices*numRows, 1); 
    Ax_large = Az_large;
    bx_large = bz_large;
    
%     sigma_max = svds(Az_large, 1, 'largest');
%     sigma_min = svds(Az_large, 1, 'smallest');
%     numCondA_large = sigma_max / sigma_min;
%     fprintf('cond(A_large) = %f\n', numCondA_large);

    % For concatenation v1.0
%     A_large = [];
%     b_large = [];
	     
    angle_z = unwrap(angle_u,[],1);
    angle_x = unwrap(angle_u,[],2);
    cont_kernel = 1; 
    for ii = 1:st:og_size(1)

        for jj = 1:st:og_size(2) %% for faster computing pararell toolbox
            
            area_z = angle_z(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel
            bz_small = area_z(:);
            area_x = angle_x(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel
            bx_small = area_x(:);

            %%%%%%%%%%%% BETTER EFFICIENCY v2.0 %%%%%%%%%%%%
            rowStart = (cont_kernel-1)*numRows + 1; % size ww
            rowEnd = rowStart + numRows - 1;
%             colStart = (cont_kernel-1)*numCols + 1; % size 3

%             fprintf('Index iter # %d idx [%d %d]\n', cont_kernel, rowStart, rowEnd);

            bz_large(rowStart:rowEnd) = bz_small;
            bx_large(rowStart:rowEnd) = bx_small;

%             Az_large(rowStart:rowStart+numRows-1, colStart:colStart+numCols-1) = A_small;
%             bz_large(rowStart:rowStart+numRows-1) = bz_small;
%             bx_large(rowStart:rowStart+numRows-1) = bx_small;
            
            %%%%%%%%%%%% BETTER EFFICIENT v2.0 %%%%%%%%%%%%

%             %%%%%%%%%%%% NOT SO EFFICIENT v1.0 %%%%%%%%%%%%
%             if (ii==1 && jj==1) % first window
%                 A_large = A_small;
%                 b_large = b_small; % concatenate b
%                 % concatenate A_small
%             else
%                 A_large = [A_large               sparse(size(A_large, 1), size(A_small, 2)); 
%                      sparse(size(A_small, 1), size(A_large, 2)), A_small];
%                 b_large = [b_large; b_small];
%             end
%             %%%%%%%%%%%% NOT SO EFFICIENT v1.0 %%%%%%%%%%%%

            cont_kernel = cont_kernel + 1;
        end
    end
    


%     disp(cont_kernel);

    %%%%% FOR x %%%%%
    
    % results_x = Ax_large\bx_large;  
    results_x = minres(Ax_large'*Ax_large,Ax_large'*bx_large);
    % results_z = Az_large\bz_large;
    results_z = minres(Az_large'*Az_large,Az_large'*bz_large);

    
    res3D_x  = reshape(results_x, [3, size_out(2), size_out(1)]); 
    res3D_x = permute(res3D_x, [3 2 1]); 
    
    %%%%% FOR z %%%%%
    
    res3D_z  = reshape(results_z, [3, size_out(2), size_out(1)]); 
    res3D_z = permute(res3D_z, [3 2 1]);

    
    grad_x = res3D_x(:,:,1); grad_z = res3D_z(:,:,2);
    phase_grad_2 = (grad_x.^2 + grad_z.^2)/constant;
    
    % ----- MedFilt  ----
    med_wind = floor (2.5/f_v/dinf.dx)*2+1; %the median window contains at least a wavelenght
    k2_med = medfilt2(phase_grad_2,[med_wind med_wind],'symmetric');
    k = sqrt(k2_med);
    % --------------------
    % k = sqrt(phase_grad_2);
    sws_matrix = (2*pi*f_v)./k;   
end