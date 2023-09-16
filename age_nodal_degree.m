%% 节点度的年龄相关性变化
load('E:\01niujinpeng\graduate\ISFC\article_fourth\struc_predict\AAL_out0.1-0.7\DegreeCentrality\DegreeCentrality.mat');
 
covariate(:,1) = xlsread('E:\01niujinpeng\graduate\behavior.xlsx','struc_predict','G2:G423');
covariate(:,2) = xlsread('E:\01niujinpeng\graduate\behavior.xlsx','struc_predict','E2:E423');
covariate(:,3:5) = xlsread('E:\01niujinpeng\graduate\behavior.xlsx','struc_predict','H2:J423');
DC = zeros(422,90);
for i = 1:90
    DC(:,i) = zscore(aDc(:,i));
end
[beta,ResMS,T,P] = dong_multi_regress(DC,covariate);
beta_age = beta(2,:);
p_age = P(2,:);
t_age = T(2,:);
clear beta P T ResMS
a_age = find(p_age < 0.05/90);
b_age = find(p_age >= 0.05/90);
t_age(b_age(:)) = 0;
% make images
aal_cortical = spm_vol('H:\Cam_CAN\struc_predict\AAL_out0.1-0.7\AAL_cortical.nii');
data = spm_read_vols(aal_cortical);
for i = 1:78
    site = find(data==i);
    data(site()) = t_age(i,1);
    clear site
end
aal_cortical.fname = 'H:\Cam_CAN\struc_predict\AAL_out0.1-0.7\DegreeCentrality\image_out\age_DC_cortical.nii';
spm_write_vol(aal_cortical,data);
clear i data aal_cortical

aal_all = spm_vol('H:\Cam_CAN\struc_predict\AAL_out0.1-0.7\AAL_all.nii');
data = spm_read_vols(aal_all);
for i = 79:90
    site = find(data==i);
    data(site()) = t_age(i,1);
    clear site
end
for i = 1:78
    site = find(data==i);
    data(site()) = 0;
    clear site
end
aal_all.fname = 'H:\Cam_CAN\struc_predict\AAL_out0.1-0.7\DegreeCentrality\image_out\age_DC_subcortical.nii';
spm_write_vol(aal_all,data);
