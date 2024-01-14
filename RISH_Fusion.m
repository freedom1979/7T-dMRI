Subjects = {'126426','130114','130518','134627','135124',...
            '146735','165436','167440','177140','180533',...
            '193845','239136','360030','385046','401422',...
            '463040','550439','644246','654552','757764',...
            '765864','878877','905147','943862','971160',...
            '995174'};
        
Num_of_Subs = length(Subjects);
for iii = 2:26
    mask = load_untouch_nii([Subjects{iii} '_mask.nii']);
    brain_mask = mask.img;
    
    bval_3T = load([Subjects{iii} '_3T_DWI_dir95_LR.bval']);
    bvec_3T = load([Subjects{iii} '_3T_DWI_dir95_LR.bvec']);
    index_3T_b2000 = find(bval_3T<2100 & bval_3T>1group900);
    bvec_3T_b2000 = bvec_3T(:,index_3T_b2000);
    index_3T_b1000 = find(bval_3T<1080 & bval_3T>900);
    bvec_3T_b1000 = bvec_3T(:,index_3T_b1000);
    index_3T_b0 = find(bval_3T<100);
    bval_3T_b0 = bval_3T(index_3T_b0);
    bvec_3T_b0 = bvec_3T(:,index_3T_b0);
    bval_3T_b0_b1000_b2000 = bval_3T([index_3T_b0 index_3T_b1000 index_3T_b2000]);
    bvec_3T_b0_b1000_b2000 = bvec_3T(:,[index_3T_b0 index_3T_b1000 index_3T_b2000]);
    save([Subjects{iii} '_bval_3T_b0_b1000_b2000_rish.bval'],'bval_3T_b0_b1000_b2000','-ascii');
    save([Subjects{iii} '_bvec_3T_b0_b1000_b2000_rish.bvec'],'bvec_3T_b0_b1000_b2000','-ascii');
    
    [~,N_gradient] = size(bval_3T_b0_b1000_b2000);
    nii_3T = load_untouch_nii([Subjects{iii} '_3T_new.nii']); 
    nii_3T_b0 = nii_3T.img(:,:,:,index_3T_b0);
    

    load([Subjects{iii} '_RISH_3T_b1000.mat'])
    load([Subjects{iii} '_RISH_3T_b2000.mat'])
    load([Subjects{iii} '_RISH_7T_b1000.mat'])
    load([Subjects{iii} '_RISH_7T_b2000.mat'])
    
    load([Subjects{iii} '_SHC_7T_b1000.mat'])
    load([Subjects{iii} '_SHC_7T_b2000.mat'])
    load([Subjects{iii} '_SHC_3T_b1000.mat'])
    load([Subjects{iii} '_SHC_3T_b2000.mat'])
    
    
    t3 = RISH_3T_b1000(:,:,:,1); t7 = RISH_7T_b1000(:,:,:,1);
    R3 = SHC_3T_b1000(:,:,:,1); R7 = SHC_7T_b1000(:,:,:,1);
    R3(t3 < t7) = 0; R7(t3 > t7) = 0;
    Fusion_3T_7T_b1000(:,:,:,1) = R3 + R7;
    
    t3 = RISH_3T_b1000(:,:,:,2); t7 = RISH_7T_b1000(:,:,:,2);
    R3 = SHC_3T_b1000(:,:,:,2:6); R7 = SHC_7T_b1000(:,:,:,2:6);
    t3 = repmat(t3,[1 1 1 5]); t7 = repmat(t7,[1 1 1 5]);
    R3(t3 < t7) = 0; R7(t3 > t7) = 0;
    Fusion_3T_7T_b1000(:,:,:,2:6) = R3 + R7;
    
    t3 = RISH_3T_b1000(:,:,:,3); t7 = RISH_7T_b1000(:,:,:,3);
    R3 = SHC_3T_b1000(:,:,:,7:15); R7 = SHC_7T_b1000(:,:,:,7:15);
    t3 = repmat(t3,[1 1 1 9]); t7 = repmat(t7,[1 1 1 9]);
    R3(t3 < t7) = 0; R7(t3 > t7) = 0;
    Fusion_3T_7T_b1000(:,:,:,7:15) = R3 + R7;
    
    t3 = RISH_3T_b1000(:,:,:,4); t7 = RISH_7T_b1000(:,:,:,4);
    R3 = SHC_3T_b1000(:,:,:,16:28); R7 = SHC_7T_b1000(:,:,:,16:28);
    t3 = repmat(t3,[1 1 1 13]); t7 = repmat(t7,[1 1 1 13]);
    R3(t3 < t7) = 0; R7(t3 > t7) = 0;
    Fusion_3T_7T_b1000(:,:,:,16:28) = R3 + R7;
    
    t3 = RISH_3T_b1000(:,:,:,5); t7 = RISH_7T_b1000(:,:,:,5);
    R3 = SHC_3T_b1000(:,:,:,29:45); R7 = SHC_7T_b1000(:,:,:,29:45);
    t3 = repmat(t3,[1 1 1 17]); t7 = repmat(t7,[1 1 1 17]);
    R3(t3 < t7) = 0; R7(t3 > t7) = 0;
    Fusion_3T_7T_b1000(:,:,:,29:45) = R3 + R7;
    
    t3 = RISH_3T_b2000(:,:,:,1); t7 = RISH_7T_b2000(:,:,:,1);
    R3 = SHC_3T_b2000(:,:,:,1); R7 = SHC_7T_b2000(:,:,:,1);
    R3(t3 < t7) = 0; R7(t3 > t7) = 0;
    Fusion_3T_7T_b2000(:,:,:,1) = R3 + R7;
    
    t3 = RISH_3T_b2000(:,:,:,2); t7 = RISH_7T_b2000(:,:,:,2);
    R3 = SHC_3T_b2000(:,:,:,2:6); R7 = SHC_7T_b2000(:,:,:,2:6);
    t3 = repmat(t3,[1 1 1 5]); t7 = repmat(t7,[1 1 1 5]);
    R3(t3 < t7) = 0; R7(t3 > t7) = 0;
    Fusion_3T_7T_b2000(:,:,:,2:6) = R3 + R7;
    
    t3 = RISH_3T_b2000(:,:,:,3); t7 = RISH_7T_b2000(:,:,:,3);
    R3 = SHC_3T_b2000(:,:,:,7:15); R7 = SHC_7T_b2000(:,:,:,7:15);
    t3 = repmat(t3,[1 1 1 9]); t7 = repmat(t7,[1 1 1 9]);
    R3(t3 < t7) = 0; R7(t3 > t7) = 0;
    Fusion_3T_7T_b2000(:,:,:,7:15) = R3 + R7;
    
    t3 = RISH_3T_b2000(:,:,:,4); t7 = RISH_7T_b2000(:,:,:,4);
    R3 = SHC_3T_b2000(:,:,:,16:28); R7 = SHC_7T_b2000(:,:,:,16:28);
    t3 = repmat(t3,[1 1 1 13]); t7 = repmat(t7,[1 1 1 13]);
    R3(t3 < t7) = 0; R7(t3 > t7) = 0;
    Fusion_3T_7T_b2000(:,:,:,16:28) = R3 + R7;
    
    t3 = RISH_3T_b2000(:,:,:,5); t7 = RISH_7T_b2000(:,:,:,5);
    R3 = SHC_3T_b2000(:,:,:,29:45); R7 = SHC_7T_b2000(:,:,:,29:45);
    t3 = repmat(t3,[1 1 1 17]); t7 = repmat(t7,[1 1 1 17]);
    R3(t3 < t7) = 0; R7(t3 > t7) = 0;
    Fusion_3T_7T_b2000(:,:,:,29:45) = R3 + R7;
    
    
    
    [azimuth,elevation,~] = cart2sph( bvec_3T_b1000(1,:),bvec_3T_b1000(2,:),bvec_3T_b1000(3,:) );
    bvec_3T_b1000_angle = [azimuth', pi/2-elevation']; %
    [azimuth,elevation,~] = cart2sph( bvec_3T_b2000(1,:),bvec_3T_b2000(2,:),bvec_3T_b2000(3,:) );
    bvec_3T_b2000_angle = [azimuth', pi/2-elevation']; %
    
    Fusion_3T_7T_b1000_inversed = zeros(200,200,132,length(bvec_3T_b1000(1,:)));
    Fusion_3T_7T_b2000_inversed = zeros(200,200,132,length(bvec_3T_b2000(1,:)));
    for k = 1:132
        for j = 1:200
            for i = 1:200
                if brain_mask(i,j,k) == 1
                    F_N1(1,1) = squeeze(Fusion_3T_7T_b1000(i,j,k,1));
                    F_N1(2:4,1) = [0;0;0];
                    F_N1(5:9,1) = squeeze(Fusion_3T_7T_b1000(i,j,k,2:6));
                    F_N1(10:16,1) = [0;0;0;0;0;0;0];
                    F_N1(17:25,1) = squeeze(Fusion_3T_7T_b1000(i,j,k,7:15));
                    F_N1(26:36,1) = [0;0;0;0;0;0;0;0;0;0;0];
                    F_N1(37:49,1) = squeeze(Fusion_3T_7T_b1000(i,j,k,16:28));
                    F_N1(50:64,1) = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
                    F_N1(65:81,1) = squeeze(Fusion_3T_7T_b1000(i,j,k,29:45));
                    
                    Fusion_3T_7T_b1000_inversed(i,j,k,:) = inverseSHT(F_N1, bvec_3T_b1000_angle, 'real');
                    
                    
                    F_N2(1,1) = squeeze(Fusion_3T_7T_b2000(i,j,k,1));
                    F_N2(2:4,1) = [0;0;0];
                    F_N2(5:9,1) = squeeze(Fusion_3T_7T_b2000(i,j,k,2:6));
                    F_N2(10:16,1) = [0;0;0;0;0;0;0];
                    F_N2(17:25,1) = squeeze(Fusion_3T_7T_b2000(i,j,k,7:15));
                    F_N2(26:36,1) = [0;0;0;0;0;0;0;0;0;0;0];
                    F_N2(37:49,1) = squeeze(Fusion_3T_7T_b2000(i,j,k,16:28));
                    F_N2(50:64,1) = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
                    F_N2(65:81,1) = squeeze(Fusion_3T_7T_b2000(i,j,k,29:45));
                    
                    Fusion_3T_7T_b2000_inversed(i,j,k,:) = inverseSHT(F_N2, bvec_3T_b2000_angle, 'real');
                end
            end
        end
    end
    
    Fusion_3T_7T_inversed = cat(4, nii_3T_b0, Fusion_3T_7T_b1000_inversed, Fusion_3T_7T_b2000_inversed);
    nii_3T.hdr.dime.dim = [4,200,200,132,N_gradient,1,1,1];
    nii_3T.img = Fusion_3T_7T_inversed;
    save_untouch_nii(nii_3T, [Subjects{iii} '_Fusion_3T_7T_inversed_RISH.nii']);   
end


%%
%{
mask = load_untouch_nii('126426_mask.nii');     %%%%%%%%%%%%%%%%%%
brain_mask = mask.img;

bval_3T = load('126426_3T_DWI_dir95_LR.bval');  %%%%%%%%%%%%%%%%%%
bvec_3T = load('126426_3T_DWI_dir95_LR.bvec');  %%%%%%%%%%%%%%%%%%
index_3T_b2000 = find(bval_3T<2100 & bval_3T>1900);
bvec_3T_b2000 = bvec_3T(:,index_3T_b2000);
index_3T_b1000 = find(bval_3T<1080 & bval_3T>900);
bvec_3T_b1000 = bvec_3T(:,index_3T_b1000);
index_3T_b0 = find(bval_3T<100);
bval_3T_b0 = bval_3T(index_3T_b0);
bvec_3T_b0 = bvec_3T(:,index_3T_b0);
bval_3T_b0_b1000_b2000 = bval_3T([index_3T_b0 index_3T_b1000 index_3T_b2000]);
bvec_3T_b0_b1000_b2000 = bvec_3T(:,[index_3T_b0 index_3T_b1000 index_3T_b2000]);
save bval_3T_b0_b1000_b2000.bval -ascii bval_3T_b0_b1000_b2000
save bvec_3T_b0_b1000_b2000.bvec -ascii bvec_3T_b0_b1000_b2000


[~,N_gradient] = size(bval_3T_b0_b1000_b2000);
nii_3T = load_untouch_nii('126426_3T_new.nii');   %%%%%%%%%%%%%%%%%%
nii_3T_b0 = nii_3T.img(:,:,:,index_3T_b0);

load 126426_RISH_3T_b1000.mat;
load 126426_RISH_3T_b2000.mat;
load 126426_RISH_7T_b1000.mat;
load 126426_RISH_7T_b2000.mat;

load 126426_SHC_7T_b1000.mat;                    %%%%%%%%%%%%%%%%%%
load 126426_SHC_7T_b2000.mat;                    %%%%%%%%%%%%%%%%%%
load 126426_SHC_3T_b1000.mat;                     %%%%%%%%%%%%%%%%%%
load 126426_SHC_3T_b2000.mat                      %%%%%%%%%%%%%%%%%%

% fusion rules
% Fusion_3T_7T_b1000 = max(SHC_3T_b1000,SHC_7T_b1000);  
% Fusion_3T_7T_b2000 = max(SHC_3T_b2000,SHC_7T_b2000);  
% Fusion_3T_7T_b1000 = (SHC_3T_b1000+SHC_7T_b1000)/2;
% Fusion_3T_7T_b2000 = (SHC_3T_b2000+SHC_7T_b2000)/2;

t3 = RISH_3T_b1000(:,:,:,1); t7 = RISH_7T_b1000(:,:,:,1);
R3 = SHC_3T_b1000(:,:,:,1); R7 = SHC_7T_b1000(:,:,:,1);
R3(t3 < t7) = 0; R7(t3 > t7) = 0;
Fusion_3T_7T_b1000(:,:,:,1) = R3 + R7;

t3 = RISH_3T_b1000(:,:,:,2); t7 = RISH_7T_b1000(:,:,:,2);
R3 = SHC_3T_b1000(:,:,:,2:6); R7 = SHC_7T_b1000(:,:,:,2:6);
t3 = repmat(t3,[1 1 1 5]); t7 = repmat(t7,[1 1 1 5]); 
R3(t3 < t7) = 0; R7(t3 > t7) = 0;
Fusion_3T_7T_b1000(:,:,:,2:6) = R3 + R7;

t3 = RISH_3T_b1000(:,:,:,3); t7 = RISH_7T_b1000(:,:,:,3);
R3 = SHC_3T_b1000(:,:,:,7:15); R7 = SHC_7T_b1000(:,:,:,7:15);
t3 = repmat(t3,[1 1 1 9]); t7 = repmat(t7,[1 1 1 9]); 
R3(t3 < t7) = 0; R7(t3 > t7) = 0;
Fusion_3T_7T_b1000(:,:,:,7:15) = R3 + R7;

t3 = RISH_3T_b1000(:,:,:,4); t7 = RISH_7T_b1000(:,:,:,4);
R3 = SHC_3T_b1000(:,:,:,16:28); R7 = SHC_7T_b1000(:,:,:,16:28);
t3 = repmat(t3,[1 1 1 13]); t7 = repmat(t7,[1 1 1 13]); 
R3(t3 < t7) = 0; R7(t3 > t7) = 0;
Fusion_3T_7T_b1000(:,:,:,16:28) = R3 + R7;

t3 = RISH_3T_b1000(:,:,:,5); t7 = RISH_7T_b1000(:,:,:,5);
R3 = SHC_3T_b1000(:,:,:,29:45); R7 = SHC_7T_b1000(:,:,:,29:45);
t3 = repmat(t3,[1 1 1 17]); t7 = repmat(t7,[1 1 1 17]); 
R3(t3 < t7) = 0; R7(t3 > t7) = 0;
Fusion_3T_7T_b1000(:,:,:,29:45) = R3 + R7;

t3 = RISH_3T_b2000(:,:,:,1); t7 = RISH_7T_b2000(:,:,:,1);
R3 = SHC_3T_b2000(:,:,:,1); R7 = SHC_7T_b2000(:,:,:,1);
R3(t3 < t7) = 0; R7(t3 > t7) = 0;
Fusion_3T_7T_b2000(:,:,:,1) = R3 + R7;

t3 = RISH_3T_b2000(:,:,:,2); t7 = RISH_7T_b2000(:,:,:,2);
R3 = SHC_3T_b2000(:,:,:,2:6); R7 = SHC_7T_b2000(:,:,:,2:6);
t3 = repmat(t3,[1 1 1 5]); t7 = repmat(t7,[1 1 1 5]); 
R3(t3 < t7) = 0; R7(t3 > t7) = 0;
Fusion_3T_7T_b2000(:,:,:,2:6) = R3 + R7;

t3 = RISH_3T_b2000(:,:,:,3); t7 = RISH_7T_b2000(:,:,:,3);
R3 = SHC_3T_b2000(:,:,:,7:15); R7 = SHC_7T_b2000(:,:,:,7:15);
t3 = repmat(t3,[1 1 1 9]); t7 = repmat(t7,[1 1 1 9]); 
R3(t3 < t7) = 0; R7(t3 > t7) = 0;
Fusion_3T_7T_b2000(:,:,:,7:15) = R3 + R7;

t3 = RISH_3T_b2000(:,:,:,4); t7 = RISH_7T_b2000(:,:,:,4);
R3 = SHC_3T_b2000(:,:,:,16:28); R7 = SHC_7T_b2000(:,:,:,16:28);
t3 = repmat(t3,[1 1 1 13]); t7 = repmat(t7,[1 1 1 13]); 
R3(t3 < t7) = 0; R7(t3 > t7) = 0;
Fusion_3T_7T_b2000(:,:,:,16:28) = R3 + R7;

t3 = RISH_3T_b2000(:,:,:,5); t7 = RISH_7T_b2000(:,:,:,5);
R3 = SHC_3T_b2000(:,:,:,29:45); R7 = SHC_7T_b2000(:,:,:,29:45);
t3 = repmat(t3,[1 1 1 17]); t7 = repmat(t7,[1 1 1 17]); 
R3(t3 < t7) = 0; R7(t3 > t7) = 0;
Fusion_3T_7T_b2000(:,:,:,29:45) = R3 + R7;



[azimuth,elevation,~] = cart2sph( bvec_3T_b1000(1,:),bvec_3T_b1000(2,:),bvec_3T_b1000(3,:) );
bvec_3T_b1000_angle = [azimuth', pi/2-elevation']; %
[azimuth,elevation,~] = cart2sph( bvec_3T_b2000(1,:),bvec_3T_b2000(2,:),bvec_3T_b2000(3,:) );
bvec_3T_b2000_angle = [azimuth', pi/2-elevation']; %

Fusion_3T_7T_b1000_inversed = zeros(200,200,132,length(bvec_3T_b1000(1,:)));
Fusion_3T_7T_b2000_inversed = zeros(200,200,132,length(bvec_3T_b2000(1,:)));
for k = 1:132
    for j = 1:200
        for i = 1:200
            if brain_mask(i,j,k) == 1               
                F_N1(1,1) = squeeze(Fusion_3T_7T_b1000(i,j,k,1));
                F_N1(2:4,1) = [0;0;0];
                F_N1(5:9,1) = squeeze(Fusion_3T_7T_b1000(i,j,k,2:6));
                F_N1(10:16,1) = [0;0;0;0;0;0;0];
                F_N1(17:25,1) = squeeze(Fusion_3T_7T_b1000(i,j,k,7:15));
                F_N1(26:36,1) = [0;0;0;0;0;0;0;0;0;0;0];
                F_N1(37:49,1) = squeeze(Fusion_3T_7T_b1000(i,j,k,16:28));
                F_N1(50:64,1) = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
                F_N1(65:81,1) = squeeze(Fusion_3T_7T_b1000(i,j,k,29:45));
                
                Fusion_3T_7T_b1000_inversed(i,j,k,:) = inverseSHT(F_N1, bvec_3T_b1000_angle, 'real');
                
                
                F_N2(1,1) = squeeze(Fusion_3T_7T_b2000(i,j,k,1));
                F_N2(2:4,1) = [0;0;0];
                F_N2(5:9,1) = squeeze(Fusion_3T_7T_b2000(i,j,k,2:6));
                F_N2(10:16,1) = [0;0;0;0;0;0;0];
                F_N2(17:25,1) = squeeze(Fusion_3T_7T_b2000(i,j,k,7:15));
                F_N2(26:36,1) = [0;0;0;0;0;0;0;0;0;0;0];
                F_N2(37:49,1) = squeeze(Fusion_3T_7T_b2000(i,j,k,16:28));
                F_N2(50:64,1) = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
                F_N2(65:81,1) = squeeze(Fusion_3T_7T_b2000(i,j,k,29:45));

                Fusion_3T_7T_b2000_inversed(i,j,k,:) = inverseSHT(F_N2, bvec_3T_b2000_angle, 'real');
            end
        end
    end
end

Fusion_3T_7T_inversed = cat(4, nii_3T_b0, Fusion_3T_7T_b1000_inversed, Fusion_3T_7T_b2000_inversed);
nii_3T.hdr.dime.dim = [4,200,200,132,N_gradient,1,1,1];
nii_3T.img = Fusion_3T_7T_inversed;
save_untouch_nii(nii_3T, 'Fusion_3T_7T_inversed_RISH.nii');   
%}
