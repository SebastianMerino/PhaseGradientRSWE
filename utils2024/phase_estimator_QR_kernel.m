function [grad_z,grad_x,k,sws_map] = phase_estimator_QR_kernel(u, w_kernel,f_v,dinf,og_size,constant,stride)
% function [grad_z,grad_x,k,sws_map] = phase_estimator_QR_kernel(u, w_kernel,f_v,dinf,og_size,constant,stride)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function that yields the shear wave speed of a region with the 
% phase gradient method with QR solver MATLAB optimized version doing by KERNEL. 
% 
% Inputs:  
%          u           : 2D region of interest to evaluate (previously mirror padding)
%          w_kernel    : vector with the size of the window kernel
%          f_v         : vibration frequency
%          dinf        : structure that contains the spatial resolutions
%          og_size     : vector containing the original size of the data
%          constant    : constant from equations (0.33 gives good results)
%                           (kx + kz)/constant
%          stride      : stride for window movement
% Outputs: 
%          grad_z       : Gradient matrix for the axial direction
%          grad_x       : Gradient matrix for the lateral direction
%          k            : Total Wavenumber matrix 
%          sws_map      : Shear wave speed map 
% Author: EMZ optimized
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
    %angle_u_unwrp = unwrap(angle_u, [], direction);
    %angle_u_unwrp = unwrap(angle_u, [], 2);   
    

    [X,Z] = meshgrid(x_axis,z_axis);
    A = [X(:) Z(:) ones(length(x_axis)*length(z_axis),1)];
%     At = A'; 
%     AtA = A'*A;
%     numCond = cond(A);

%     fprintf('cond(A) = %f\n', numCond);

    %% HELP FOR DIMENSIONS %% % THEORY
    size_mirror = size(u); % ogsize + w_kernel - 1; %  = size (u)
    disp(size_mirror)
    disp(w_kernel)
    disp(st)
    numkernels = floor( (size_mirror - w_kernel)./st + 1 ); % numKernels
    overlap_ax1 = 1 - st/w_kernel(1);
    overlap_la1 = 1 - st/w_kernel(2);
    size_out = floor( (og_size - 1)./st + 1 );

    angle_z = unwrap(angle_u,[],1);
    angle_x = unwrap(angle_u,[],2);
    for ii = 1:st:og_size(1)
        for jj = 1:st:og_size(2) %% for faster computing pararell toolbox
            area_z = angle_z(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel
            area_x = angle_x(ii: ii+w_kernel(1)-1,jj:jj+w_kernel(2)-1); % Window kernel
        
            %% QR solver z
            b_z = area_z(:);
            results_z = A\b_z;
%             results_z = AtA\At*b_z; % same effect
            grad_z(ii,jj) = results_z(2);
            %% QR solver x
            b_x = area_x(:);
            results_x = A\b_x;
%             results = AtA\At*b_x; % same effect
            grad_x(ii,jj) = results_x(1);

        end
    end
    
    phase_grad_2 = (grad_x.^2 + grad_z.^2)/constant;

    med_wind = floor (2.5/f_v/dinf.dx)*2+1; %the median window contains at least a wavelenght
%     k2_med = medfilt2(phase_grad_2,[med_wind med_wind],'symmetric')/constant;
    k2_med = medfilt2(phase_grad_2,[med_wind med_wind],'symmetric');
    k = sqrt(k2_med);

    sws_map = (2*pi*f_v)./k;
end