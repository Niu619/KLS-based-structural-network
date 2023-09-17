function gretna_preprocessing_extractOtherSurf_CAT12(Data_path, Para)

%==========================================================================
% This function is used to perform tissue segmentation of images (typically
% structural MRI images) with the CAT12 toolbox.
%
%
% Syntax: function gretna_preprocessing_Surf_others_CAT12(Data_path, Para)
%
% Inputs:
%       Data_path:
%                   The directory & filename of a .txt file that contains
%                   the directory of those files to be processed (can be
%                   obtained by gretna_gen_data_path.m).
%   Para (optional):
%                   Para.Nthreads:
%                     The number of threads for parallel calculation.
%                   Para.OtherSurfMetrics_GI:
%                     'yes'-- Extract additional GI (gyrification index).
%                     'no' -- Not extract additional GI (gyrification index).
%                   Para.OtherSurfMetrics_FD:
%                     'yes'-- Extract additional FD (fractal dimension).
%                     'no' -- Not extract additional FD (fractal dimension).
%                   Para.OtherSurfMetrics_SD:
%                     'yes'-- Extract additional SD (sulcus depth).
%                     'no' -- Not extract additional SD (sulcus depth).
%
% Jinhui WANG, HZNU, Hangzhou, 2017/01/16, jinhui.wang.1982@gmail.com
%==========================================================================

if nargin == 1
    Para.Nthreads = spm_input('Number of Threads',1,'e',[],1);
    Para.OtherSurfMetrics_GI = spm_input('Whether gyrification index','-1','y/n',[1,0],2);
    Para.OtherSurfMetrics_FD = spm_input('Whether fractal dimension','-1','y/n',[1,0],2);
    Para.OtherSurfMetrics_SD = spm_input('Whether sulcus depth','-1','y/n',[1,0],2);
end
close

load gretna_extractOtherSurf_CAT12.mat
batch_othersurfmeasures = matlabbatch;
batch_othersurfmeasures{1}.spm.tools.cat.stools.surfextract.nproc = Para.Nthreads;
batch_othersurfmeasures{1}.spm.tools.cat.stools.surfextract.GI    = Para.OtherSurfMetrics_GI;
batch_othersurfmeasures{1}.spm.tools.cat.stools.surfextract.FD    = Para.OtherSurfMetrics_FD;
batch_othersurfmeasures{1}.spm.tools.cat.stools.surfextract.SD    = Para.OtherSurfMetrics_SD;

fid = fopen(Data_path);
Dir_data = textscan(fid, '%s');
fclose(fid);

Num_subs = size(Dir_data{1},1);
Sour = cell(Num_subs,1);

for isub = 1:Num_subs
    cd(fullfile([Dir_data{1}{isub}]));
    
    Gii_central = ls('lh.central.*.gii');
    Sour(isub,1) = {fullfile([Dir_data{1}{isub}], Gii_central)};
end

batch_othersurfmeasures{1}.spm.tools.cat.stools.surfextract.data_surf = Sour;
spm_jobman('run',batch_othersurfmeasures);

return