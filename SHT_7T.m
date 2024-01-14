%{
tic
mask = load_nii('brain_mask.nii'); 
brain_mask = mask.img;
order = [0 2 4 6 8];


nii_7T = load_nii('126426_7T_DWI_dir72_final.nii'); 
bval_7T = load('126426_7T_DWI_dir72_AP.bval');
bvec_7T = load('126426_7T_DWI_dir72_AP.bvec');

index_7T_b2000 = find(bval_7T<80);  % bval_7T<2100 & bval_7T>1900  % bval_7T<1020 & bval_7T>900
bvec_7T_b2000 = bvec_7T(:,index_7T_b2000); 
[azimuth,elevation,r] = cart2sph( bvec_7T_b2000(1,:),bvec_7T_b2000(2,:),bvec_7T_b2000(3,:) );
bvec_7T_b2000_angle = [azimuth',pi/2 - elevation'];

dwi_7T_b2000 = nii_7T.img(:,:,:,index_7T_b2000);

RISH_7T = zeros(200,200,132,length(order));

for k = 1:132
    for j = 1:200
        for i = 1:200
            if brain_mask(i,j,k) == 1
                voxel = double(squeeze(dwi_7T_b2000(i,j,k,:)));
                for r = 1:5
                    [F_N, Y_N] = leastSquaresSHT(order(r), voxel, bvec_7T_b2000_angle, 'real');
                    RISH_7T(i,j,k,r) = sum( abs(F_N(end-2*order(r):end)).*abs(F_N(end-2*order(r):end)), 1 );
                end
            end
        end
    end
end

save RISH_7T_126426_B0.mat RISH_7T;
toc
%}


%%
Subjects = {'126426','130114','130518','134627','135124',...
            '146735','165436','167440','177140','180533',...
            '193845','239136','360030','385046','401422',...
            '463040','550439','644246','654552','757764',...
            '765864','878877','905147','943862','971160',...
            '995174'};
        
highres_points = load('100_shell.txt')';
highres_points_2 = [highres_points -highres_points];

Num_of_Subs = length(Subjects);
for ii = 21:21
    tic
    mask_name = [Subjects{ii} '_mask.nii'];
    mask = load_untouch_nii(mask_name);
    brain_mask = mask.img;
    order = [0 2 4 6 8];
    
    nii_7T_name = [Subjects{ii} '_7T_DWI_dir72_final.nii'];
    nii_7T = load_untouch_nii(nii_7T_name);
    bval_7T = load([Subjects{ii} '_7T_DWI_dir72_AP.bval']);
    bvec_7T = load([Subjects{ii} '_7T_DWI_dir72_AP.bvec']);
 
    % b2000
    index_7T_b2000 = find(bval_7T<2100 & bval_7T>1900);  % bval_7T<2100 & bval_7T>1900  % bval_7T<1020 & bval_7T>900 % bval_7T<80
    bvec_7T_b2000 = bvec_7T(:,index_7T_b2000);
    bvec_7T_b2000_2 = [bvec_7T_b2000 -bvec_7T_b2000];
%     [azimuth,elevation,~] = cart2sph( bvec_7T_b2000_2(1,:),bvec_7T_b2000_2(2,:),bvec_7T_b2000_2(3,:) );
%     bvec_7T_b2000_angle = [azimuth',pi/2 - elevation'];
    dwi_7T_b2000 = nii_7T.img(:,:,:,index_7T_b2000);
    dwi_7T_b2000_2 = repmat(dwi_7T_b2000, [1 1 1 2]);
    
    SHC_7T_b2000 = zeros(200,200,132,45);
    RISH_7T_b2000 = zeros(200,200,132,length(order));
    for k = 1:132
        for j = 1:200
            for i = 1:200
                if brain_mask(i,j,k) == 1
                    voxel = double(squeeze(dwi_7T_b2000_2(i,j,k,:)));
                    
%                     F = scatteredInterpolant(bvec_7T_b2000_2(1,:)',bvec_7T_b2000_2(2,:)',bvec_7T_b2000_2(3,:)', double(voxel));
%                     F.Method = 'natural';
%                     voxel_2 = F(highres_points(1,:)',highres_points(2,:)',highres_points(3,:)');

                    voxel_2 = matlab_call_66(bvec_7T_b2000_2', voxel, highres_points');
                    
                    [azimuth, elevation,~] = cart2sph( highres_points_2(1,:),highres_points_2(2,:),highres_points_2(3,:) );
                    bvec_7T_b2000_angle = [azimuth', pi/2 - elevation'];
                    voxel_3 = [voxel_2; voxel_2];
                    
                    F_N = directSHT(8, voxel_3, bvec_7T_b2000_angle, 'real');
                    
                    temp = F_N(1);
                    RISH_7T_b2000(i,j,k,1) = sum(abs(temp).*abs(temp));
                    temp = F_N(5:9);
                    RISH_7T_b2000(i,j,k,2) = sum(abs(temp).*abs(temp));
                    temp = F_N(17:25);
                    RISH_7T_b2000(i,j,k,3) = sum(abs(temp).*abs(temp));
                    temp = F_N(37:49);
                    RISH_7T_b2000(i,j,k,4) = sum(abs(temp).*abs(temp));
                    temp = F_N(65:81);
                    RISH_7T_b2000(i,j,k,5) = sum(abs(temp).*abs(temp));
                    
                    SHC_7T_b2000(i,j,k,1) = F_N(1);
                    SHC_7T_b2000(i,j,k,2:6) = F_N(5:9);
                    SHC_7T_b2000(i,j,k,7:15) = F_N(17:25);
                    SHC_7T_b2000(i,j,k,16:28) = F_N(37:49);
                    SHC_7T_b2000(i,j,k,29:45) = F_N(65:81);
                end
            end
        end
    end
    
    save_name = [Subjects{ii} '_RISH_7T_b2000.mat'];
    save(save_name, 'RISH_7T_b2000');
    save_name2 = [Subjects{ii} '_SHC_7T_b2000.mat'];
    save(save_name2, 'SHC_7T_b2000');
    toc
    
    tic
    % b1000
    index_7T_b1000 = find(bval_7T<1080 & bval_7T>900);  % bval_7T<2100 & bval_7T>1900  % bval_7T<1020 & bval_7T>900 % bval_7T<80
    bvec_7T_b1000 = bvec_7T(:,index_7T_b1000);
    bvec_7T_b1000_2 = [bvec_7T_b1000 -bvec_7T_b1000];
%     [azimuth,elevation,~] = cart2sph( bvec_7T_b1000_2(1,:),bvec_7T_b1000_2(2,:),bvec_7T_b1000_2(3,:) );
%     bvec_7T_b1000_angle = [azimuth',pi/2 - elevation'];
    dwi_7T_b1000 = nii_7T.img(:,:,:,index_7T_b1000);
    dwi_7T_b1000_2 = repmat(dwi_7T_b1000, [1 1 1 2]);
    
    SHC_7T_b1000 = zeros(200,200,132,45);
    RISH_7T_b1000 = zeros(200,200,132,length(order));
    for k = 1:132
        for j = 1:200
            for i = 1:200
                if brain_mask(i,j,k) == 1
                    voxel = double(squeeze(dwi_7T_b1000_2(i,j,k,:)));
                         
%                     F = scatteredInterpolant(bvec_7T_b1000_2(1,:)',bvec_7T_b1000_2(2,:)',bvec_7T_b1000_2(3,:)', double(voxel));
%                     F.Method = 'natural';
%                     voxel_2 = F(highres_points(1,:)',highres_points(2,:)',highres_points(3,:)');

                    voxel_2 = matlab_call_64(bvec_7T_b1000_2', voxel, highres_points');
                    
                    [azimuth, elevation,~] = cart2sph( highres_points_2(1,:),highres_points_2(2,:),highres_points_2(3,:) );
                    bvec_7T_b1000_angle = [azimuth', pi/2 - elevation'];
                    voxel_3 = [voxel_2; voxel_2]; 
                    
                    F_N = directSHT(8, voxel_3, bvec_7T_b1000_angle, 'real');
                    
                    temp = F_N(1);
                    RISH_7T_b1000(i,j,k,1) = sum(abs(temp).*abs(temp));
                    temp = F_N(5:9);
                    RISH_7T_b1000(i,j,k,2) = sum(abs(temp).*abs(temp));
                    temp = F_N(17:25);
                    RISH_7T_b1000(i,j,k,3) = sum(abs(temp).*abs(temp));
                    temp = F_N(37:49);
                    RISH_7T_b1000(i,j,k,4) = sum(abs(temp).*abs(temp));
                    temp = F_N(65:81);
                    RISH_7T_b1000(i,j,k,5) = sum(abs(temp).*abs(temp));
                    
                    SHC_7T_b1000(i,j,k,1) = F_N(1);
                    SHC_7T_b1000(i,j,k,2:6) = F_N(5:9);
                    SHC_7T_b1000(i,j,k,7:15) = F_N(17:25);
                    SHC_7T_b1000(i,j,k,16:28) = F_N(37:49);
                    SHC_7T_b1000(i,j,k,29:45) = F_N(65:81);
                end
            end
        end
    end
    
    save_name = [Subjects{ii} '_RISH_7T_b1000.mat'];
    save(save_name, 'RISH_7T_b1000');  
    save_name2 = [Subjects{ii} '_SHC_7T_b1000.mat'];
    save(save_name2, 'SHC_7T_b1000');
    
    toc
end






