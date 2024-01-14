% 1. 7T -> 3T; 
% 2. 3T B1000/B2000 sphereical interpolation along 7T diffusion directions
% 3. Compare 3T/7T signals
% 4. optimize diffusion signals


% 7T->3T harmonization learning
e = 0.001;
order = [0 2 4 6 8];

Subjects = {'126426','130114','130518','134627','135124',...
            '146735','165436','167440','177140','180533',...
            '193845','239136','360030','385046','401422',...
            '463040','550439','644246','654552','757764'};

Num_of_Subs = length(Subjects);
        
RISH_3T_b2000_All = zeros(200,200,132,5,Num_of_Subs);
RISH_7T_b2000_All = zeros(200,200,132,5,Num_of_Subs);
RISH_3T_b1000_All = zeros(200,200,132,5,Num_of_Subs);
RISH_7T_b1000_All = zeros(200,200,132,5,Num_of_Subs);

for ii = 1:Num_of_Subs
    load([Subjects{ii} '_RISH_3T_b2000.mat']);
    RISH_3T_b2000_All(:,:,:,:,ii) = RISH_3T_b2000;
    load([Subjects{ii} '_RISH_7T_b2000.mat']);
    RISH_7T_b2000_All(:,:,:,:,ii) = RISH_7T_b2000;
    load([Subjects{ii} '_RISH_3T_b1000.mat']);
    RISH_3T_b1000_All(:,:,:,:,ii) = RISH_3T_b1000;
    load([Subjects{ii} '_RISH_7T_b1000.mat']);
    RISH_7T_b1000_All(:,:,:,:,ii) = RISH_7T_b1000;
end

% b2000
E0_3T_b2000 = zeros(200,200,132);
E2_3T_b2000 = zeros(200,200,132);
E4_3T_b2000 = zeros(200,200,132);
E6_3T_b2000 = zeros(200,200,132);
E8_3T_b2000 = zeros(200,200,132);

E0_7T_b2000 = zeros(200,200,132);
E2_7T_b2000 = zeros(200,200,132);
E4_7T_b2000 = zeros(200,200,132);
E6_7T_b2000 = zeros(200,200,132);
E8_7T_b2000 = zeros(200,200,132);

for ii = 1:Num_of_Subs
    T = squeeze( RISH_3T_b2000_All(:,:,:,:,ii) );
    E0_3T_b2000 = E0_3T_b2000 + T(:,:,:,1);
    E2_3T_b2000 = E2_3T_b2000 + T(:,:,:,2);
    E4_3T_b2000 = E4_3T_b2000 + T(:,:,:,3);
    E6_3T_b2000 = E6_3T_b2000 + T(:,:,:,4);
    E8_3T_b2000 = E8_3T_b2000 + T(:,:,:,5);
    
    T = squeeze( RISH_7T_b2000_All(:,:,:,:,ii) );
    E0_7T_b2000 = E0_7T_b2000 + T(:,:,:,1);
    E2_7T_b2000 = E2_7T_b2000 + T(:,:,:,2);
    E4_7T_b2000 = E4_7T_b2000 + T(:,:,:,3);
    E6_7T_b2000 = E6_7T_b2000 + T(:,:,:,4);
    E8_7T_b2000 = E8_7T_b2000 + T(:,:,:,5);
end

E0_3T_b2000 = E0_3T_b2000./Num_of_Subs;
E2_3T_b2000 = E2_3T_b2000./Num_of_Subs;
E4_3T_b2000 = E4_3T_b2000./Num_of_Subs;
E6_3T_b2000 = E6_3T_b2000./Num_of_Subs;
E8_3T_b2000 = E8_3T_b2000./Num_of_Subs;

E0_7T_b2000 = E0_7T_b2000./Num_of_Subs;
E2_7T_b2000 = E2_7T_b2000./Num_of_Subs;
E4_7T_b2000 = E4_7T_b2000./Num_of_Subs;
E6_7T_b2000 = E6_7T_b2000./Num_of_Subs;
E8_7T_b2000 = E8_7T_b2000./Num_of_Subs;

T_b2000_R0 = sqrt(E0_3T_b2000./(E0_7T_b2000+e));
T_b2000_R2 = sqrt(E2_3T_b2000./(E2_7T_b2000+e));
T_b2000_R4 = sqrt(E4_3T_b2000./(E4_7T_b2000+e));
T_b2000_R6 = sqrt(E6_3T_b2000./(E6_7T_b2000+e));
T_b2000_R8 = sqrt(E8_3T_b2000./(E8_7T_b2000+e));


% b1000
E0_3T_b1000 = zeros(200,200,132);
E2_3T_b1000 = zeros(200,200,132);
E4_3T_b1000 = zeros(200,200,132);
E6_3T_b1000 = zeros(200,200,132);
E8_3T_b1000 = zeros(200,200,132);

E0_7T_b1000 = zeros(200,200,132);
E2_7T_b1000 = zeros(200,200,132);
E4_7T_b1000 = zeros(200,200,132);
E6_7T_b1000 = zeros(200,200,132);
E8_7T_b1000 = zeros(200,200,132);

for ii = 1:Num_of_Subs
    T = squeeze( RISH_3T_b1000_All(:,:,:,:,ii) );
    E0_3T_b1000 = E0_3T_b1000 + T(:,:,:,1);
    E2_3T_b1000 = E2_3T_b1000 + T(:,:,:,2);
    E4_3T_b1000 = E4_3T_b1000 + T(:,:,:,3);
    E6_3T_b1000 = E6_3T_b1000 + T(:,:,:,4);
    E8_3T_b1000 = E8_3T_b1000 + T(:,:,:,5);
    
    T = squeeze( RISH_7T_b1000_All(:,:,:,:,ii) );
    E0_7T_b1000 = E0_7T_b1000 + T(:,:,:,1);
    E2_7T_b1000 = E2_7T_b1000 + T(:,:,:,2);
    E4_7T_b1000 = E4_7T_b1000 + T(:,:,:,3);
    E6_7T_b1000 = E6_7T_b1000 + T(:,:,:,4);
    E8_7T_b1000 = E8_7T_b1000 + T(:,:,:,5);
end

E0_3T_b1000 = E0_3T_b1000./Num_of_Subs;
E2_3T_b1000 = E2_3T_b1000./Num_of_Subs;
E4_3T_b1000 = E4_3T_b1000./Num_of_Subs;
E6_3T_b1000 = E6_3T_b1000./Num_of_Subs;
E8_3T_b1000 = E8_3T_b1000./Num_of_Subs;

E0_7T_b1000 = E0_7T_b1000./Num_of_Subs;
E2_7T_b1000 = E2_7T_b1000./Num_of_Subs;
E4_7T_b1000 = E4_7T_b1000./Num_of_Subs;
E6_7T_b1000 = E6_7T_b1000./Num_of_Subs;
E8_7T_b1000 = E8_7T_b1000./Num_of_Subs;

T_b1000_R0 = sqrt(E0_3T_b1000./(E0_7T_b1000+e));
T_b1000_R2 = sqrt(E2_3T_b1000./(E2_7T_b1000+e));
T_b1000_R4 = sqrt(E4_3T_b1000./(E4_7T_b1000+e));
T_b1000_R6 = sqrt(E6_3T_b1000./(E6_7T_b1000+e));
T_b1000_R8 = sqrt(E8_3T_b1000./(E8_7T_b1000+e));



% testing
Subjects = {'765864','878877','905147','943862','971160','995174'};
Num_of_Subs = length(Subjects);
for ii = 1:1
    mask_name = [Subjects{ii} '_mask.nii'];
    mask = load_untouch_nii(mask_name);
    brain_mask = mask.img;
    order = [0 2 4 6 8];   
   
    % b2000     
    bval_7T = load([Subjects{ii} '_7T_DWI_dir72_AP.bval']);
    bvec_7T = load([Subjects{ii} '_7T_DWI_dir72_AP.bvec']);
    index_7T_b2000 = find(bval_7T<2100 & bval_7T>1900);  % bval_7T<2100 & bval_7T>1900  % bval_7T<1020 & bval_7T>900 % bval_7T<80
    bvec_7T_b2000 = bvec_7T(:,index_7T_b2000);
    [azimuth,elevation,~] = cart2sph( bvec_7T_b2000(1,:),bvec_7T_b2000(2,:),bvec_7T_b2000(3,:) );
    bvec_7T_b2000_angle = [azimuth',pi/2 - elevation'];  
     
    SHC_7T_name = [Subjects{ii} '_SHC_7T_b2000.mat'];
    SHC_7T = importdata(SHC_7T_name);  
    SCH_7T_b2000_harmonized = zeros(size(SHC_7T));
    b2000_7T_harmonized_inversed = zeros(200,200,132,length(bvec_7T_b2000(1,:)));
    
    for k = 1:132
        for j = 1:200
            for i = 1:200
                if brain_mask(i,j,k) == 1
                    SCH_7T_b2000_harmonized(i,j,k,1) = SHC_7T(i,j,k,1).*T_b2000_R0(i,j,k);
                    SCH_7T_b2000_harmonized(i,j,k,2:6) = SHC_7T(i,j,k,2:6).*T_b2000_R2(i,j,k);    
                    SCH_7T_b2000_harmonized(i,j,k,7:15) = SHC_7T(i,j,k,7:15).*T_b2000_R4(i,j,k);
                    SCH_7T_b2000_harmonized(i,j,k,16:28) = SHC_7T(i,j,k,16:28).*T_b2000_R6(i,j,k);
                    SCH_7T_b2000_harmonized(i,j,k,29:45) = SHC_7T(i,j,k,29:45).*T_b2000_R8(i,j,k);
                    
                    F_N(1,1) = squeeze(SCH_7T_b2000_harmonized(i,j,k,1));
                    F_N(2:4,1) = [0;0;0];
                    F_N(5:9,1) = squeeze(SCH_7T_b2000_harmonized(i,j,k,2:6));
                    F_N(10:16,1) = [0;0;0;0;0;0;0];
                    F_N(17:25,1) = squeeze(SCH_7T_b2000_harmonized(i,j,k,7:15));
                    F_N(26:36,1) = [0;0;0;0;0;0;0;0;0;0;0];
                    F_N(37:49,1) = squeeze(SCH_7T_b2000_harmonized(i,j,k,16:28));
                    F_N(50:64,1) = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
                    F_N(65:81,1) = squeeze(SCH_7T_b2000_harmonized(i,j,k,29:45));
                       
                    b2000_7T_harmonized_inversed(i,j,k,:) = inverseSHT(F_N, bvec_7T_b2000_angle, 'real');
                end
            end
        end
    end
    
    save_name = [Subjects{ii} '_SCH_7T_b2000_harmonized.mat'];
    save(save_name, 'SCH_7T_b2000_harmonized');
    save_name2 = [Subjects{ii} '_b2000_7T_harmonized_inversed.mat'];
    save(save_name2, 'b2000_7T_harmonized_inversed');
    
    
    % b1000
    index_7T_b1000 = find(bval_7T<1020 & bval_7T>900);  % bval_7T<2100 & bval_7T>1900  % bval_7T<1020 & bval_7T>900 % bval_7T<80
    bvec_7T_b1000 = bvec_7T(:,index_7T_b1000);
    [azimuth,elevation,~] = cart2sph( bvec_7T_b1000(1,:),bvec_7T_b1000(2,:),bvec_7T_b1000(3,:) );
    bvec_7T_b1000_angle = [azimuth',pi/2 - elevation'];

    SHC_7T_name = [Subjects{ii} '_SHC_7T_b1000.mat'];
    SHC_7T = importdata(SHC_7T_name);  
    SCH_7T_b1000_harmonized = zeros(size(SHC_7T));
    b1000_7T_harmonized_inversed = zeros(200,200,132,length(bvec_7T_b1000(1,:)));
    
    for k = 1:132
        for j = 1:200
            for i = 1:200
                if brain_mask(i,j,k) == 1
                    SCH_7T_b1000_harmonized(i,j,k,1) = SHC_7T(i,j,k,1).*T_b1000_R0(i,j,k);
                    SCH_7T_b1000_harmonized(i,j,k,2:6) = SHC_7T(i,j,k,2:6).*T_b1000_R2(i,j,k);
                    SCH_7T_b1000_harmonized(i,j,k,7:15) = SHC_7T(i,j,k,7:15).*T_b1000_R4(i,j,k);
                    SCH_7T_b1000_harmonized(i,j,k,16:28) = SHC_7T(i,j,k,16:28).*T_b1000_R6(i,j,k);
                    SCH_7T_b1000_harmonized(i,j,k,29:45) = SHC_7T(i,j,k,29:45).*T_b1000_R8(i,j,k);
                    
                    F_N(1,1) = squeeze(SCH_7T_b1000_harmonized(i,j,k,1));
                    F_N(2:4,1) = [0;0;0];
                    F_N(5:9,1) = squeeze(SCH_7T_b1000_harmonized(i,j,k,2:6));
                    F_N(10:16,1) = [0;0;0;0;0;0;0];
                    F_N(17:25,1) = squeeze(SCH_7T_b1000_harmonized(i,j,k,7:15));
                    F_N(26:36,1) = [0;0;0;0;0;0;0;0;0;0;0];
                    F_N(37:49,1) = squeeze(SCH_7T_b1000_harmonized(i,j,k,16:28));
                    F_N(50:64,1) = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
                    F_N(65:81,1) = squeeze(SCH_7T_b1000_harmonized(i,j,k,29:45));
                    
                    b1000_7T_harmonized_inversed(i,j,k,:) = inverseSHT(F_N, bvec_7T_b1000_angle, 'real');
                end
            end
        end
    end
    
    save_name = [Subjects{ii} '_SCH_7T_b1000_harmonized.mat'];
    save(save_name, 'SCH_7T_b1000_harmonized');
    save_name2 = [Subjects{ii} '_b1000_7T_harmonized_inversed.mat'];
    save(save_name2, 'b1000_7T_harmonized_inversed');    
end



