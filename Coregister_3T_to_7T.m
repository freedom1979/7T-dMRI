%% coregister 3T to 7T space

Subjects = {'126426','130114','130518','134627','135124',...
            '146735','165436','167440','177140','180533',...
            '193845','239136','360030','385046','401422',...
            '463040','550439','644246','654552','757764',...
            '765864','878877','905147','943862','971160',...
            '995174'};

        
Num_of_Subs = length(Subjects);
for i = 1:1
    nii_3T = load_untouch_nii([Subjects{i} '_3T_DWI_dir95_final.nii']);
    bval_3T = load([Subjects{i} '_3T_DWI_dir95_LR.bval']);
    bval_3T_0 = find(bval_3T<70);
    B0_3T = nii_3T.img(:,:,:,bval_3T_0);
    B0_3T_Mean = mean(B0_3T, 4);
    [m_3T,n_3T,s_3T] = size(B0_3T_Mean);
    
    nii_7T = load_untouch_nii([Subjects{i} '_7T_DWI_dir72_final.nii']);
    bval_7T = load([Subjects{i} '_7T_DWI_dir72_AP.bval']);
    bval_7T_0 = find(bval_7T<80);
    B0_7T = nii_7T.img(:,:,:,bval_7T_0); 
    B0_7T_Mean = mean(B0_7T, 4);
    [m_7T,n_7T,s_7T] = size(B0_7T_Mean);
    
    nii_3T.hdr.dime.dim = [3,m_3T,n_3T,s_3T,1,1,1,1];
    nii_3T.original.hdr.dime.dim = [3,m_3T,n_3T,s_3T,1,1,1,1];
    nii_3T.img = B0_3T_Mean;
    save_untouch_nii(nii_3T, 'B0_3T.nii');
    nii_7T.hdr.dime.dim = [3,m_7T,n_7T,s_7T,1,1,1,1];
    nii_7T.original.hdr.dime.dim = [3,m_7T,n_7T,s_7T,1,1,1,1];
    nii_7T.img = B0_7T_Mean;
    save_untouch_nii(nii_7T, 'B0_7T.nii');
    
    !d:/tools/diffusionkit/reg_aladin -ref B0_7T.nii -flo B0_3T.nii -res B0_3T_r.nii -aff B0_3T_r.txt

    img_data = zeros(200,200,132,95);
 
    nii_3T_s = nii_3T;  % nii_3T different from the origin
    nii_3T = load_untouch_nii([Subjects{i} '_3T_DWI_dir95_final.nii']);
    
    for j = 1:95
        cmd1 = ['BD' num2str(j) '_3T = nii_3T.img(:,:,:,' num2str(j) ');']; eval(cmd1);
        cmd2 = ['[m_3T,n_3T,s_3T] = size(BD' num2str(j) '_3T);']; eval(cmd2);
        nii_3T_s.hdr.dime.dim = [3,m_3T,n_3T,s_3T,1,1,1,1];
        nii_3T_s.original.hdr.dime.dim = [3,m_3T,n_3T,s_3T,1,1,1,1];
        cmd3 = ['nii_3T_s.img = BD' num2str(j) '_3T;']; eval(cmd3);
        cmd4 = ['save_untouch_nii(nii_3T_s,' '''nii_3T_s.nii'');']; eval(cmd4);
        
        cmd5 = ['!d:/tools/diffusionkit/reg_aladin -ref B0_3T_r.nii -flo ' 'nii_3T_s.nii -res BD_3T_r.nii -aff BD_3T_r.txt'];
        eval(cmd5);
        
        BD_3T_r = load_untouch_nii('BD_3T_r.nii');
        img_data(:,:,:,j) = BD_3T_r.img;
    end
    
    nii_3T_new = load_untouch_nii([Subjects{i} '_3T_DWI_dir95_final.nii']);
    nii_3T_new.hdr.dime.dim = [4,200,200,132,95,1,1,1,1];
    nii_3T_new.original.hdr.dime.dim = [4,200,200,132,95,1,1,1,1];
    nii_3T_new.img = img_data;
    save_untouch_nii(nii_3T_new, [Subjects{i} '_3T_new.nii']);
end



%%
%%
%%

%{
%% 3T_B0 co-registered to 7T_B0
nii_3T = load_nii('126426_3T_DWI_dir95_final.nii');
bval_3T = load('126426_3T_DWI_dir95_LR.bval');
bval_3T_0 = find(bval_3T<70);
B0_3T = nii_3T.img(:,:,:,bval_3T_0);
B0_3T_Mean = mean(B0_3T, 4);
[m_3T,n_3T,s_3T] = size(B0_3T_Mean);

nii_7T = load_nii('126426_7T_DWI_dir72_final.nii');
bval_7T = load('126426_7T_DWI_dir72_AP.bval');
bval_7T_0 = find(bval_3T<80);
B0_7T = nii_7T.img(:,:,:,bval_7T_0); 
B0_7T_Mean = mean(B0_7T, 4);
[m_7T,n_7T,s_7T] = size(B0_7T_Mean);

nii_3T.hdr.dime.dim = [3,m_3T,n_3T,s_3T,1,1,1,1];
nii_3T.original.hdr.dime.dim = [3,m_3T,n_3T,s_3T,1,1,1,1];
nii_3T.img = B0_3T_Mean;
save_nii(nii_3T, 'B0_3T.nii'); 
nii_7T.hdr.dime.dim = [3,m_7T,n_7T,s_7T,1,1,1,1];
nii_7T.original.hdr.dime.dim = [3,m_7T,n_7T,s_7T,1,1,1,1];
nii_7T.img = B0_7T_Mean;
save_nii(nii_7T, 'B0_7T.nii'); 

!d:/tools/diffusionkit/reg_aladin -ref B0_7T.nii -flo B0_3T.nii -res B0_3T_r.nii -aff B0_3T_r.txt


%% whole 3T volume registered to 7T voulume
% split, co-register
img_data = zeros(200,200,132,95);
B0_3T_r = load_nii('B0_3T_r.nii');
img_data(:,:,:,1) = B0_3T_r.img;

nii_3T = load_nii('126426_3T_DWI_dir95_final.nii');
nii_3T_s = nii_3T;
for i = 2:95
    cmd1 = ['BD' num2str(i) '_3T = nii_3T.img(:,:,:,' num2str(i) ');']; eval(cmd1);
    cmd2 = ['[m_3T,n_3T,s_3T] = size(BD' num2str(i) '_3T);']; eval(cmd2);
    nii_3T_s.hdr.dime.dim = [3,m_3T,n_3T,s_3T,1,1,1,1];
    nii_3T_s.original.hdr.dime.dim = [3,m_3T,n_3T,s_3T,1,1,1,1];
    cmd3 = ['nii_3T_s.img = BD' num2str(i) '_3T;']; eval(cmd3);
    cmd4 = ['save_nii(nii_3T_s,''BD' num2str(i) '_3T.nii'');']; eval(cmd4);
    
    cmd5 = ['!d:/tools/diffusionkit/reg_aladin -ref B0_3T_r.nii -flo BD' num2str(i) '_3T.nii -res BD_3T_r.nii -aff BD_3T_r.txt'];
    eval(cmd5); 
    
    BD_3T_r = load_nii('BD_3T_r.nii');
    img_data(:,:,:,i) = BD_3T_r.img;
end

% save 3T_new first
nii_3T_new = load_nii('126426_3T_DWI_dir95_final.nii');
nii_3T_new.hdr.dime.dim = [4,200,200,132,95,1,1,1,1];
nii_3T_new.original.hdr.dime.dim = [4,200,200,132,95,1,1,1,1];
nii_3T_new.img = img_data;
save_nii(nii_3T_new, '126426_3T_new.nii'); 
%}


