Subjects = {'126426','130114','130518','134627','135124',...
            '146735','165436','167440','177140','180533',...
            '193845','239136','360030','385046','401422',...
            '463040','550439','644246','654552','757764',...
            '765864','878877','905147','943862','971160',...
            '995174'};
        
Num_of_Subs = length(Subjects);
for iii = 19:26
    mask_name = [Subjects{iii} '_mask.nii'];
    mask = load_untouch_nii(mask_name);
    brain_mask = mask.img;
    
    bval_7T = load([Subjects{iii} '_7T_DWI_dir72_AP.bval']);
    bvec_7T = load([Subjects{iii} '_7T_DWI_dir72_AP.bvec']);
    nii_7T = load_untouch_nii([Subjects{iii} '_7T_DWI_dir72_final.nii']);
    
    bval_3T = load([Subjects{iii} '_3T_DWI_dir95_LR.bval']);
    bvec_3T = load([Subjects{iii} '_3T_DWI_dir95_LR.bvec']);
    nii_3T_new = load_untouch_nii([Subjects{iii} '_3T_new.nii']);

    index_3T_b2000 = find(bval_3T<2100 & bval_3T>1900);
    bvec_3T_b2000 = bvec_3T(:,index_3T_b2000);
    nii_3T_b2000 = nii_3T_new.img(:,:,:,index_3T_b2000);
    bval_3T_b2000 = bval_3T(index_3T_b2000);
    
    index_3T_b1000 = find(bval_3T<1080 & bval_3T>900);
    bvec_3T_b1000 = bvec_3T(:,index_3T_b1000);
    nii_3T_b1000 = nii_3T_new.img(:,:,:,index_3T_b1000);
    bval_3T_b1000 = bval_3T(index_3T_b1000);
    
    index_3T_b0 = find(bval_3T<100);
    nii_3T_b0 = nii_3T_new.img(:,:,:,index_3T_b0);
    [~,~,~,N_3T_b0] = size(nii_3T_b0);
    nii_3T_b0_Mean = mean(nii_3T_b0,4);
    
    index_7T_b2000 = find(bval_7T<2100 & bval_7T>1900);
    bvec_7T_b2000 = bvec_7T(:,index_7T_b2000);
    nii_7T_b2000 = nii_7T.img(:,:,:,index_7T_b2000);
    bval_7T_b2000 = bval_7T(index_7T_b2000);
    
    index_7T_b1000 = find(bval_7T<1020 & bval_7T>900);
    bvec_7T_b1000 = bvec_7T(:,index_7T_b1000);
    nii_7T_b1000 = nii_7T.img(:,:,:,index_7T_b1000);
    bval_7T_b1000 = bval_7T(index_7T_b1000);
    
    index_7T_b0 = find(bval_7T<100);
    nii_7T_b0 = nii_7T.img(:,:,:,index_7T_b0);
    [~,~,~,N_7T_b0] = size(nii_7T_b0);
    nii_7T_b0_Mean = mean(nii_7T_b0,4);
    
    [~,~,~,N_7T_b2000] = size(nii_7T_b2000);
    b2000_7T_harmonized_bymyself = zeros(200,200,132,N_7T_b2000);
    
    nii_7T_b2000 = double(nii_7T_b2000);
    nii_7T_b1000 = double(nii_7T_b1000);
    nii_7T_b0_Mean = double(nii_7T_b0_Mean);
    nii_3T_b0_Mean = double(nii_3T_b0_Mean);
    
    for vol = 1:N_7T_b2000
        D_b2000 = (-1)/bval_7T_b2000(vol).*log(nii_7T_b2000(:,:,:,vol)./nii_7T_b0_Mean);
        D_b2000(brain_mask ~= 1) = 0;
        Q = nii_3T_b0_Mean.*exp(0-bval_7T_b2000(vol).*D_b2000); 
        Q(brain_mask ~= 1) = 0;
        b2000_7T_harmonized_bymyself(:,:,:,vol) = Q;
    end
    Sum_Inter = 0;
    N_Inter = 0;
    for vol = 1:N_7T_b2000
        for k = 2:132-1
            for j = 2:200-1
                for i = 2:200-1
                    if brain_mask(i,j,k) == 1
                        if b2000_7T_harmonized_bymyself(i,j,k,vol) == Inf
                            if b2000_7T_harmonized_bymyself(i-1,j,k,vol) ~= Inf
                                Sum_Inter = Sum_Inter + b2000_7T_harmonized_bymyself(i-1,j,k,vol);
                                N_Inter = N_Inter + 1;
                            end
                            if b2000_7T_harmonized_bymyself(i+1,j,k,vol) ~= Inf
                                Sum_Inter = Sum_Inter + b2000_7T_harmonized_bymyself(i+1,j,k,vol);
                                N_Inter = N_Inter + 1;
                            end
                            if b2000_7T_harmonized_bymyself(i,j-1,k,vol) ~= Inf
                                Sum_Inter = Sum_Inter + b2000_7T_harmonized_bymyself(i,j-1,k,vol);
                                N_Inter = N_Inter + 1;
                            end
                            if b2000_7T_harmonized_bymyself(i,j+1,k,vol) ~= Inf
                                Sum_Inter = Sum_Inter + b2000_7T_harmonized_bymyself(i,j+1,k,vol);
                                N_Inter = N_Inter + 1;
                            end
                            if b2000_7T_harmonized_bymyself(i,j,k-1,vol) ~= Inf
                                Sum_Inter = Sum_Inter + b2000_7T_harmonized_bymyself(i,j,k-1,vol);
                                N_Inter = N_Inter + 1;
                            end
                            if b2000_7T_harmonized_bymyself(i,j,k+1,vol) ~= Inf
                                Sum_Inter = Sum_Inter + b2000_7T_harmonized_bymyself(i,j,k+1,vol);
                                N_Inter = N_Inter + 1;
                            end
                            
                            if N_Inter == 0
                                disp('Error: Cannot Interplate 7T_b2000 Inf values!');
                                pause;
                            end
                            
                            b2000_7T_harmonized_bymyself(i,j,k,vol) = Sum_Inter/N_Inter;
                        end
                    end
                end
            end
        end
    end
    
    save_name = [Subjects{iii} 'b2000_7T_harmonized_bymyself.mat'];
    save(save_name, 'b2000_7T_harmonized_bymyself');
    
    [~,~,~,N_7T_b1000] = size(nii_7T_b1000);
    b1000_7T_harmonized_bymyself = zeros(200,200,132,N_7T_b1000);
    for vol = 1:N_7T_b1000
        D_b1000 = (-1)/bval_7T_b1000(vol).*log(nii_7T_b1000(:,:,:,vol)./nii_7T_b0_Mean);
        D_b1000(brain_mask ~= 1) = 0;
        Q = nii_3T_b0_Mean.*exp(0-bval_7T_b1000(vol).*D_b1000); %
        Q(brain_mask ~= 1) = 0;
        b1000_7T_harmonized_bymyself(:,:,:,vol) = Q;
    end
    Sum_Inter = 0;
    N_Inter = 0;
    for vol = 1:N_7T_b1000
        for k = 2:132-1
            for j = 2:200-1
                for i = 2:200-1
                    if brain_mask(i,j,k) == 1
                        if b1000_7T_harmonized_bymyself(i,j,k,vol) == Inf
                            if b1000_7T_harmonized_bymyself(i-1,j,k,vol) ~= Inf
                                Sum_Inter = Sum_Inter + b1000_7T_harmonized_bymyself(i-1,j,k,vol);
                                N_Inter = N_Inter + 1;
                            end
                            if b1000_7T_harmonized_bymyself(i+1,j,k,vol) ~= Inf
                                Sum_Inter = Sum_Inter + b1000_7T_harmonized_bymyself(i+1,j,k,vol);
                                N_Inter = N_Inter + 1;
                            end
                            if b1000_7T_harmonized_bymyself(i,j-1,k,vol) ~= Inf
                                Sum_Inter = Sum_Inter + b1000_7T_harmonized_bymyself(i,j-1,k,vol);
                                N_Inter = N_Inter + 1;
                            end
                            if b1000_7T_harmonized_bymyself(i,j+1,k,vol) ~= Inf
                                Sum_Inter = Sum_Inter + b1000_7T_harmonized_bymyself(i,j+1,k,vol);
                                N_Inter = N_Inter + 1;
                            end
                            if b1000_7T_harmonized_bymyself(i,j,k-1,vol) ~= Inf
                                Sum_Inter = Sum_Inter + b1000_7T_harmonized_bymyself(i,j,k-1,vol);
                                N_Inter = N_Inter + 1;
                            end
                            if b1000_7T_harmonized_bymyself(i,j,k+1,vol) ~= Inf
                                Sum_Inter = Sum_Inter + b1000_7T_harmonized_bymyself(i,j,k+1,vol);
                                N_Inter = N_Inter + 1;
                            end
                            
                            if N_Inter == 0
                                disp('Error: Cannot Interplate 7T_b2000 Inf values!');
                                pause;
                            end
                            
                            b1000_7T_harmonized_bymyself(i,j,k,vol) = Sum_Inter/N_Inter;
                        end
                    end
                end
            end
        end
    end
    
    save_name = [Subjects{iii} 'b1000_7T_harmonized_bymyself.mat'];
    save(save_name, 'b1000_7T_harmonized_bymyself');
end        
