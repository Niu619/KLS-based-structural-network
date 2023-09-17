function gretna_preprocessing_Surf_resample_CAT12(Data_path, File_filter, Para)

%==========================================================================
% This function is used to perform tissue segmentation of images (typically
% structural MRI images) with the CAT12 toolbox.
%
%
% Syntax: function gretna_preprocessing_Surf_resample_CAT12(Data_path, File_filter, Para)
%
% Inputs:
%       Data_path:
%                   The directory & filename of a .txt file that contains
%                   the directory of those files to be processed (can be
%                   obtained by gretna_gen_data_path.m).
%       File_filter:
%                   The prefix of those surface files for resample.
%   Para (optional):
%                   Para.Nthreads:
%                     The number of threads for parallel calculation.
%                   Para.SmoothSize:
%                     Smoohting size (default 15 mm). Note, 15 mm is
%                     suggested for cortical thickness and 25 mm for
%                     others.
%
% Jinhui WANG, HZNU, Hangzhou, 2017/01/16, jinhui.wang.1982@gmail.com
%==========================================================================

if nargin == 2
    Para.Nthreads = spm_input('Number of Threads',1,'e',[],1);
    Para.SmoothSize = spm_input('Smoothing Size','+1','e',15); % default 15 mm
end
close

load gretna_resample_CAT12.mat
batch_resamplesmooth = matlabbatch;
batch_resamplesmooth{1}.spm.tools.cat.stools.surfresamp.nproc = Para.Nthreads;
batch_resamplesmooth{1}.spm.tools.cat.stools.surfresamp.fwhm_surf = Para.SmoothSize;

fid = fopen(Data_path);
Dir_data = textscan(fid, '%s');
fclose(fid);

Num_subs = size(Dir_data{1},1);

data_surf = [];
for isub = 1:Num_subs
    cd (fullfile([Dir_data{1}{isub}]));
    
    surf = ls(['lh.' File_filter '*']);
    data_surf = char(data_surf, [repmat([pwd filesep],size(surf,1),1) surf]);
end

batch_resamplesmooth{1}.spm.tools.cat.stools.surfresamp.data_surf = cellstr(data_surf(2:end,:));

spm_jobman('run',batch_resamplesmooth)

return