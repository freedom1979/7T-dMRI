%{ 
load S765864_b2000_7T_harmonized_bymyself.mat;
load S765864_b1000_7T_harmonized_bymyself.mat;


bval_3T = load('765864_3T_DWI_dir95_LR.bval');
bvec_3T = load('765864_3T_DWI_dir95_LR.bvec');
nii_3T_new = load_untouch_nii('765864_3T_new.nii');
CCC = nii_3T_new.img;

bval_7T = load('765864_7T_DWI_dir72_AP.bval');
bvec_7T = load('765864_7T_DWI_dir72_AP.bvec');
index_7T_b2000 = find(bval_7T<2100 & bval_7T>1900);
bvec_7T_b2000 = bvec_7T(:,index_7T_b2000);
bval_7T_b2000 = bval_7T(index_7T_b2000);
index_7T_b1000 = find(bval_7T<1020 & bval_7T>900);
bvec_7T_b1000 = bvec_7T(:,index_7T_b1000);
bval_7T_b1000 = bval_7T(index_7T_b1000);

bval_3T_b0_b1000_b2000 = [bval_3T bval_7T_b1000 bval_7T_b2000];
bvec_3T_b0_b1000_b2000 = [bvec_3T bvec_7T_b1000 bvec_7T_b2000];
save bval_3T_b0_b1000_b2000.bval -ascii bval_3T_b0_b1000_b2000
save bvec_3T_b0_b1000_b2000.bvec -ascii bvec_3T_b0_b1000_b2000
[~,N_gradient] = size(bval_3T_b0_b1000_b2000);


Fusion_3T_7T_inversed = cat(4, CCC, S765864_b1000_7T_harmonized_bymyself, S765864_b2000_7T_harmonized_bymyself);
nii_3T_new.hdr.dime.dim = [4,200,200,132,N_gradient,1,1,1];
nii_3T_new.img = Fusion_3T_7T_inversed;
save_untouch_nii(nii_3T_new, 'Fusion_3T_7T_DIRECT.nii');   
%}


%% 
Subjects = {'126426','130114','130518','134627','135124',...
            '146735','165436','167440','177140','180533',...
            '193845','239136','360030','385046','401422',...
            '463040','550439','644246','654552','757764',...
            '765864','878877','905147','943862','971160',...
            '995174'};
        
Num_of_Subs = length(Subjects);
for iii = 2:2
    b1000_7T_harmonized_bymyself = importdata([Subjects{iii} 'b1000_7T_harmonized_bymyself.mat']);
    b2000_7T_harmonized_bymyself = importdata([Subjects{iii} 'b2000_7T_harmonized_bymyself.mat']);
    
    bval_3T = load([Subjects{iii} '_3T_DWI_dir95_LR.bval']);
    bvec_3T = load([Subjects{iii} '_3T_DWI_dir95_LR.bvec']);
    index_3T_b2000 = find(bval_3T<2100 & bval_3T>1900);
    bvec_3T_b2000 = bvec_3T(:,index_3T_b2000);
    bval_3T_b2000 = bval_3T(index_3T_b2000);
    index_3T_b1000 = find(bval_3T<1080 & bval_3T>900);
    bvec_3T_b1000 = bvec_3T(:,index_3T_b1000);
    bval_3T_b1000 = bval_3T(index_3T_b1000);
    
    nii_3T_new = load_untouch_nii([Subjects{iii} '_3T_new.nii']);
    index_3T_b0 = find(bval_3T<100);
    bvec_3T_b0 = bvec_3T(:,index_3T_b0);
    bval_3T_b0 = bval_3T(index_3T_b0);
    nii_3T_b0 = nii_3T_new.img(:,:,:,index_3T_b0);
    nii_3T_b1000 = nii_3T_new.img(:,:,:,index_3T_b1000);
    nii_3T_b2000 = nii_3T_new.img(:,:,:,index_3T_b2000);
    
    bval_3T_b0_b1000_b2000 = [bval_3T_b0 bval_3T_b1000 bval_3T_b2000];
    bvec_3T_b0_b1000_b2000 = [bvec_3T_b0 bvec_3T_b1000 bvec_3T_b2000];
    save([Subjects{iii} '_bval_3T_b0_b1000_b2000_direct.bval'],'bval_3T_b0_b1000_b2000','-ascii');
    save([Subjects{iii} '_bvec_3T_b0_b1000_b2000_direct.bvec'],'bvec_3T_b0_b1000_b2000','-ascii');
    [~,N_gradient] = size(bval_3T_b0_b1000_b2000);
    dwi_3T_b0_b1000_b2000 = cat(4, nii_3T_b0, nii_3T_b1000, nii_3T_b2000);
    nii_3T_new.hdr.dime.dim = [4,200,200,132,N_gradient,1,1,1];
    nii_3T_new.img = dwi_3T_b0_b1000_b2000;
    save_untouch_nii(nii_3T_new, [Subjects{iii} '_dwi_3T_b0_b1000_b2000.nii']);
    
    
    bval_7T = load([Subjects{iii} '_7T_DWI_dir72_AP.bval']);
    bvec_7T = load([Subjects{iii} '_7T_DWI_dir72_AP.bvec']);
    index_7T_b2000 = find(bval_7T<2100 & bval_7T>1900);
    bvec_7T_b2000 = bvec_7T(:,index_7T_b2000);
    bval_7T_b2000 = bval_7T(index_7T_b2000);
    index_7T_b1000 = find(bval_7T<1020 & bval_7T>900);
    bvec_7T_b1000 = bvec_7T(:,index_7T_b1000);
    bval_7T_b1000 = bval_7T(index_7T_b1000);
    
    
    bval_3T_7T_b0_b1000_b2000 = [bval_3T_b0_b1000_b2000 bval_7T_b1000 bval_7T_b2000];
    bvec_3T_7T_b0_b1000_b2000 = [bvec_3T_b0_b1000_b2000 bvec_7T_b1000 bvec_7T_b2000];
    save([Subjects{iii} '_bval_3T_7T_b0_b1000_b2000_direct.bval'],'bval_3T_7T_b0_b1000_b2000','-ascii');
    save([Subjects{iii} '_bvec_3T_7T_b0_b1000_b2000_direct.bvec'],'bvec_3T_7T_b0_b1000_b2000','-ascii');
    [~,N_gradient] = size(bval_3T_7T_b0_b1000_b2000);
    
    
    dwi_3T_7T_b0_b1000_b2000 = cat(4, dwi_3T_b0_b1000_b2000, b1000_7T_harmonized_bymyself, b2000_7T_harmonized_bymyself);
    nii_3T_new.hdr.dime.dim = [4,200,200,132,N_gradient,1,1,1];
    nii_3T_new.img = dwi_3T_7T_b0_b1000_b2000;
    save_untouch_nii(nii_3T_new, [Subjects{iii} '_dwi_3T_7T_b0_b1000_b2000_direct.nii']);
end


%%
%{
load 126426b1000_7T_harmonized_bymyself.mat;
load 126426b2000_7T_harmonized_bymyself.mat;


bval_3T = load('126426_3T_DWI_dir95_LR.bval');
bvec_3T = load('126426_3T_DWI_dir95_LR.bvec');
index_3T_b2000 = find(bval_3T<2100 & bval_3T>1900);
bvec_3T_b2000 = bvec_3T(:,index_3T_b2000);
bval_3T_b2000 = bval_3T(index_3T_b2000);
index_3T_b1000 = find(bval_3T<1080 & bval_3T>900);
bvec_3T_b1000 = bvec_3T(:,index_3T_b1000);
bval_3T_b1000 = bval_3T(index_3T_b1000);

nii_3T_new = load_untouch_nii('126426_3T_new.nii');
index_3T_b0 = find(bval_3T<100);
bvec_3T_b0 = bvec_3T(:,index_3T_b0);
bval_3T_b0 = bval_3T(index_3T_b0);
nii_3T_b0 = nii_3T_new.img(:,:,:,index_3T_b0);
nii_3T_b1000 = nii_3T_new.img(:,:,:,index_3T_b1000);
nii_3T_b2000 = nii_3T_new.img(:,:,:,index_3T_b2000);

bval_3T_b0_b1000_b2000 = [bval_3T_b0 bval_3T_b1000 bval_3T_b2000];
bvec_3T_b0_b1000_b2000 = [bvec_3T_b0 bvec_3T_b1000 bvec_3T_b2000];
save bval_3T_b0_b1000_b2000.bval -ascii bval_3T_b0_b1000_b2000
save bvec_3T_b0_b1000_b2000.bvec -ascii bvec_3T_b0_b1000_b2000
[~,N_gradient] = size(bval_3T_b0_b1000_b2000);
dwi_3T_b0_b1000_b2000 = cat(4, nii_3T_b0, nii_3T_b1000, nii_3T_b2000);
nii_3T_new.hdr.dime.dim = [4,200,200,132,N_gradient,1,1,1];
nii_3T_new.img = dwi_3T_b0_b1000_b2000;
save_untouch_nii(nii_3T_new, 'dwi_3T_b0_b1000_b2000.nii');  


bval_7T = load('126426_7T_DWI_dir72_AP.bval');
bvec_7T = load('126426_7T_DWI_dir72_AP.bvec');
index_7T_b2000 = find(bval_7T<2100 & bval_7T>1900);
bvec_7T_b2000 = bvec_7T(:,index_7T_b2000);
bval_7T_b2000 = bval_7T(index_7T_b2000);
index_7T_b1000 = find(bval_7T<1020 & bval_7T>900);
bvec_7T_b1000 = bvec_7T(:,index_7T_b1000);
bval_7T_b1000 = bval_7T(index_7T_b1000);


bval_3T_7T_b0_b1000_b2000 = [bval_3T_b0_b1000_b2000 bval_7T_b1000 bval_7T_b2000];
bvec_3T_7T_b0_b1000_b2000 = [bvec_3T_b0_b1000_b2000 bvec_7T_b1000 bvec_7T_b2000];
save bval_3T_7T_b0_b1000_b2000.bval -ascii bval_3T_7T_b0_b1000_b2000
save bvec_3T_7T_b0_b1000_b2000.bvec -ascii bvec_3T_7T_b0_b1000_b2000
[~,N_gradient] = size(bval_3T_7T_b0_b1000_b2000);


dwi_3T_7T_b0_b1000_b2000 = cat(4, dwi_3T_b0_b1000_b2000, b1000_7T_harmonized_bymyself, b2000_7T_harmonized_bymyself);
nii_3T_new.hdr.dime.dim = [4,200,200,132,N_gradient,1,1,1];
nii_3T_new.img = dwi_3T_7T_b0_b1000_b2000;
save_untouch_nii(nii_3T_new, 'dwi_3T_7T_b0_b1000_b2000_direct.nii');   
%}
