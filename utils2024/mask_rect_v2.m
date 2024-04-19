function [rect_mask] = mask_rect_v2(x, z, cx, cz, Lx, Lz)
% function [im_out, rect_mask] = mask_rect(x, z, cx1, cx2, cz1, cz2, im_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DESCRIPTION : Generic function for metric calculation, extract rectangular ROI (V1.0)
% rectangular region
% INPUTS: 
%         - x   : Dim2 MATLAB (left-right) x_axis vector (only matters inital and final value)
%         - z   : Dim1 MATLAB (up-down) z_axis vector (only matters inital and final value) 
%         - cx : left position  
%         - cz : up position 
%         - Lx : 
%         - Lz: input image
% Example Region of Interest (ROI)
%
%           (cx1, cz1)  ---------  (cx2, cz1)  
%               |                       |
%               |                       |
%           (cx1, cz2)  ---------  (cx2, cz2)  
%
% OUTPUTS: 
%         - im_out     : Output image (Only ROI values, others are NaN) 
%         - rect_mask  : Mask (ROI are '1', others are '0')
% AUTHORs: Edmundo Arom Miranda 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [X_m,Z_m]= meshgrid(x,z);
    x0 = cx-Lx/2;
    z0 = cz-Lz/2;
    
    rect_mask = and( X_m >= x0 , X_m <= x0+Lx ) & ...
        and(Z_m >= z0, Z_m <= z0+Lz);
    

end

