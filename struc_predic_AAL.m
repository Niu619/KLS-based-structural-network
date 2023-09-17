%% 
% PLS regression
response_var_file = 'E:\01niujinpeng\graduate\ISFC\article_fourth\struc_predict\AAL_out0.1-0.7\DegreeCentrality\PLSinput\fluid_intelligence.csv';
predictor_var_file = 'E:\01niujinpeng\graduate\ISFC\article_fourth\struc_predict\AAL_out0.1-0.7\DegreeCentrality\PLSinput\aDC.csv';
output_dir = 'E:\01niujinpeng\graduate\ISFC\article_fourth\struc_predict\AAL_out0.1-0.7\DegreeCentrality\PLSoutput';
% response variable names 
importdata(response_var_file);
subject_id=ans.data(:,1);
ResponseVarNames=ans.textdata(1,:);
ResponseVarNames(1)=[];
fluid_intell=ans.data(:,2);
clear ans

%import predictor variables
indata=importdata(predictor_var_file);
DC_data=indata.data;
DC_data(1,:)=[];
ROI_name=indata.textdata;
ROI_name=ROI_name(2:length(ROI_name));
ROIindex=1:length(ROI_name);
clear indata

%DO PLS in 2 dimensions (with 2 components) 
X = DC_data';
X = zscore(X);
Y = zscore(fluid_intell);
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

figure,
plot(1:dim,cumsum(100*PCTVAR(2,1:dim)),'-o','LineWidth',1.5,'Color',[140/255,0,0]);
set(gca,'Fontsize',14)
xlabel('Number of PLS components','FontSize',14);
ylabel('Percent Variance Explained in Y','FontSize',14);

[R1,p1]=corr([XS(:,1),XS(:,2)],fluid_intell);

%align PLS components with desired direction%
if R1(1,1)<0
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end

[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=ROI_name(x1);
ROIindex1=ROIindex(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids=ROI_name(x2);
ROIindex2=ROIindex(x2);

%print out results
csvwrite(fullfile(output_dir,'PLS1_scores.csv'),XS(:,1));
csvwrite(fullfile(output_dir,'PLS2_scores.csv'),XS(:,2));
% permutation
rep=5000;
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
temp=cumsum(100*PCTVAR(2,1:dim));
Rsquared = temp(dim);
    for j=1:rep
        %j
        order=randperm(size(Y,1));
        Yp=Y(order,:);

        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp,dim);

        temp=cumsum(100*PCTVAR(2,1:dim));
        Rsq(j) = temp(dim);
    end
R(dim)=Rsquared
p(dim)=length(find(Rsq>=Rsquared))/rep

% PLS1,PLS2 and fluid permutation
rep=5000;
dim=2;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
null_PLS1 = zeros(length(rep),1);
null_PLS2 = zeros(length(rep),1);
for m=1:rep
    null_order = randperm(size(XS,1));    
    pls = XS(null_order,:);
    [null_PLS1(m,1),~] = corr(pls(:,1),fluid_intell);
    [null_PLS2(m,1),~] = corr(pls(:,2),fluid_intell);
end
p_pls1 = length(find(null_PLS1>=R1(1,1)))/rep
p_pls2 = length(find(null_PLS2>=R1(2,1)))/rep
%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];

%start bootstrap number of bootstrap iterations
bootnum=5000;
disp('Bootstrapping - could take a while')
for i=1:bootnum
    myresample = randsample(size(X,1),size(X,1));
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled regions
    Yr=Y(myresample,:); % define Y for resampled regions
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
    
    temp=stats.W(:,1);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
       newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run
    
    temp=stats.W(:,2);%extract PLS2 weights
    newW=temp(x2); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
       newW=-1*newW;
    end
    PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run    
end

%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');

%get bootstrap weights
temp1=PLS1w./PLS1sw';
temp2=PLS2w./PLS2sw';

%order bootstrap weights (Z) and names of regions
[Z1 ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
ROIindex1=ROIindex1(ind1);
[Z2 ind2]=sort(temp2,'descend');
PLS2=PLS2ids(ind2);
ROIindex2=ROIindex2(ind2);

% correction
a1 = 1 - normcdf(abs(Z1));
h1 = fdr(a1,0.01);
index1 = find(h1==1);

a2 = 1 - normcdf(abs(Z2));
h2 = fdr(a2,0.01);
index2 = find(h2==1);

%print out results
fid1 = fopen(fullfile(output_dir,'PLS1_ROIWeights.csv'),'w');
for i=1:length(ROI_name)
  fprintf(fid1,'%s, %d, %f\n', PLS1{i}, ROIindex1(i), Z1(i));
end
fclose(fid1);

fid2 = fopen(fullfile(output_dir,'PLS2_ROIWeights.csv'),'w');
for i=1:length(ROI_name)
  fprintf(fid2,'%s, %d, %f\n', PLS2{i},ROIindex2(i), Z2(i));
end
fclose(fid2);

% make images
% PLS1
PLS1 = importdata('H:\Cam_CAN\struc_predict\AAL_out0.1-0.7\DegreeCentrality\PLSoutput\PLS1_ROIWeights.csv');
ROI_index = PLS1.textdata;
PLS1 = PLS1.data(:,4);
z_PLS1 = mapminmax(PLS1',0,1)';

aal_cortical = spm_vol('H:\Cam_CAN\struc_predict\AAL_out0.1-0.7\AAL_cortical.nii');
data = spm_read_vols(aal_cortical);
for i = 1:78
    site = find(data==i);
    data(site()) = z_PLS1(i,1);
    clear site
end
aal_cortical.fname = 'H:\Cam_CAN\struc_predict\AAL_out0.1-0.7\DegreeCentrality\image_out\PLS1cortical.nii';
spm_write_vol(aal_cortical,data);
clear i data aal_cortical

aal_all = spm_vol('H:\Cam_CAN\struc_predict\AAL_out0.1-0.7\AAL_all.nii');
data = spm_read_vols(aal_all);
for i = 79:90
    site = find(data==i);
    data(site()) = z_PLS1(i,1);
    clear site
end
for i = 1:78
    site = find(data==i);
    data(site()) = 0;
    clear site
end
aal_all.fname = 'H:\Cam_CAN\struc_predict\AAL_out0.1-0.7\DegreeCentrality\image_out\PLS1subcortical.nii';
spm_write_vol(aal_all,data);

% PLS2
PLS2 = importdata('H:\Cam_CAN\struc_predict\AAL_out0.1-0.7\DegreeCentrality\PLSoutput\PLS2_ROIWeights.csv');
ROI_index = PLS2.textdata;
PLS2 = PLS2.data(:,4);
z_PLS2 = mapminmax(PLS2',0,1)';

aal_cortical = spm_vol('H:\Cam_CAN\struc_predict\AAL_out0.1-0.7\AAL_cortical.nii');
data = spm_read_vols(aal_cortical);
for i = 1:78
    site = find(data==i);
    data(site()) = z_PLS2(i,1);
    clear site
end
aal_cortical.fname = 'H:\Cam_CAN\struc_predict\AAL_out0.1-0.7\DegreeCentrality\image_out\PLS2cortical.nii';
spm_write_vol(aal_cortical,data);
clear i data aal_cortical

aal_all = spm_vol('H:\Cam_CAN\struc_predict\AAL_out0.1-0.7\AAL_all.nii');
data = spm_read_vols(aal_all);
for i = 79:90
    site = find(data==i);
    data(site()) = z_PLS2(i,1);
    clear site
end
for i = 1:78
    site = find(data==i);
    data(site()) = 0;
    clear site
end
aal_all.fname = 'H:\Cam_CAN\struc_predict\AAL_out0.1-0.7\DegreeCentrality\image_out\PLS2subcortical.nii';
spm_write_vol(aal_all,data);




