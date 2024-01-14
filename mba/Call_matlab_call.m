% call c++ function


% input parameters
load voxel.mat;   % voxel
load bvec.mat;    % bvec_7T_b2000_2
highres_points = load('100_shell.txt')';


vo = matlab_call(bvec_7T_b2000_2', voxel, highres_points');

dwi_plot_function(vo,highres_points(1,:)',highres_points(2,:)',highres_points(3,:)');

% dwi_plot_function(voxel,bvec_7T_b2000_2(1,:)',bvec_7T_b2000_2(2,:)',bvec_7T_b2000_2(3,:)');