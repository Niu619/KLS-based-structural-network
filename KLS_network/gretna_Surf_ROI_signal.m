function gretna_Surf_ROI_signal(Data_path, File_filter_LH, File_filter_RH, Path_template_LH, Path_template_RH, Ind_ROI_LH, Ind_ROI_RH, Output_path, Level, Num_thres)

%==========================================================================
% This function is used to extract signals within each ROI defined in a
% brain surface atlas. NB. brain surface atlases and surface results
% (e.g., cortical thickness) are typically sorted in two files corresponding
% to the two brain hemishperes.
%
% Syntax: function gretna_Surf_ROI_signal(Data_path, File_filter_LH, File_filter_RH, Path_template_LH, Path_template_RH, Ind_ROI_LH, Ind_ROI_RH, Output_path, Level, Num_thres)
%
% Inputs:
%         Data_path:
%                   The directory & filename of a .txt file that contains
%                   the directory of those files to be copied (can be
%                   obtained by gretna_gen_data_path.m).
%    File_filter_LH:
%                   The prefix of those files to be processed
%                   (LEFT hemisphere).
%    File_filter_RH:
%                   The prefix of those files to be processed
%                   (RIGHT hemisphere).
%  Path_template_LH:
%                   The directory & filename of a surface brain atlas
%                   (LEFT hemisphere).
%  Path_template_RH:
%                   The directory & filename of a surface brain atlas
%                   (RIGHT hemisphere).
%        Ind_ROI_LH:
%                   The index of ROIs in the LEFT hemisphere.
%        Ind_ROI_RH:
%                   The index of ROIs in the RIGHT hemisphere.
%                   NOTE that the order of extracted data is the same as
%                   the order of [RoiIndex_LH RoiIndex_RH].
%       Output_path:
%                   The directory where the resultant files are sorted.
%             Level:
%                   'all':  extract all signals in each ROI;
%                   'mean': extract mean signal in each ROI;
%                   'both': extract both all and mean signals in each ROI.
%         Num_thres:
%                   Positive integer threshold for regional number of
%                   vertices. Regions with vertices less than the threshold
%                   are excluded (Default = 128).
%
% Outputs:
% Signal_LR_xxx.mat:
%                   1 * M (# of ROIs) cell array for each subject with each
%                   cell containing a data array. NB. the outputted files
%                   are sorted in the same order as that in Data_path.
% mean_Signal_LR_xxx.mat:
%                   1 * M (# of ROIs) data array for each subject.NB. the
%                   outputted files are sorted in the same order as that in
%                   Data_path.
%
% Jinhui WANG, CCBD, HZNU, Hangzhou, 2016/12/14, Jinhui.Wang.1982@gmail.com
% Yinzhi Lee,  CCBD, HZNU, Hangzhou, 2016/12/14, Yuyedaren@gmail.com
% Ningkai Wang,IBRR, SCNU, Guangzhou,2020/6/25,  Ningkai.Wang.1993@gmail.com
%==========================================================================

[~, Name_template_LH, Ext_template_LH] = fileparts(Path_template_LH);
[~, ~,                Ext_template_RH] = fileparts(Path_template_RH);

if ~strcmpi(Ext_template_LH, Ext_template_RH)
    error('File type is not the same for atlases between the two hemispheres!');
end

if ~isdir(Output_path)
    mkdir(Output_path)
end
if  nargin == 9
    Num_thres = 128;
end

%% Determine template type
switch Ext_template_LH
    case '.annot'
        [~, Label_LH, Colortable_LH] = read_annotation(Path_template_LH, 0);
        Label_region_LH = Colortable_LH.table(2:end, 5);
        
        [~, Label_RH, colortable_RH] = read_annotation(Path_template_RH, 0);
        Label_region_RH = colortable_RH.table(2:end, 5);
        
    case '.gii'
        Colortable_LH   = gifti(Path_template_LH);
        Label_LH        = Colortable_LH.cdata;
        Label_region_LH = unique(Label_LH);
        Label_region_LH = Label_region_LH(Label_region_LH ~= 0);
        
        colortable_RH   = gifti(Path_template_RH);
        Label_RH        = colortable_RH.cdata;
        Label_region_RH = unique(Label_RH);
        Label_region_RH = Label_region_RH(Label_region_RH ~= 0);
        
    otherwise
        error('Not supported file type %s\n', Ext_template_LH);
end

Num_label_LH = length(Label_region_LH);

fid      = fopen(Data_path);
Dir_data = textscan(fid, '%s');
fclose(fid);

Num_subs    = size(Dir_data{1},1);
Num_regs_LH = length(Ind_ROI_LH);
Num_regs_RH = length(Ind_ROI_RH);

%% Extracting regional signal for each subject
Num_ver_reg       = zeros(Num_subs, Num_regs_LH + Num_regs_RH);
Removed_subs_reg  = zeros(Num_subs, Num_regs_LH + Num_regs_RH);
Name_sub          = cell(Num_subs,1);
Cell_removed_regs = cell(Num_subs,1);
getname           = @(x) inputname(1);

for i_sub = 1:Num_subs
    fprintf('Extracting regional signal for %s\n', Dir_data{1}{i_sub});
    
    File_name_LH = dir([Dir_data{1}{i_sub} filesep File_filter_LH '*.gii']);
    File_name_RH = dir([Dir_data{1}{i_sub} filesep File_filter_RH '*.gii']);
    
    if isempty(File_name_LH) || isempty(File_name_RH)
        error('There are no such specified .gii files in %s ', Dir_data{1}{i_sub});
    end
    
    Num_imgs_LH = size(File_name_LH, 1);
    Num_imgs_RH = size(File_name_RH, 1);
    
    if Num_imgs_LH ~= Num_imgs_RH
        error('The number of images is not equal between the two hemispheres for %s \n', [Dir_data{1}{i_sub}]);
    end
    
    Sig_roi           = cell(Num_imgs_LH, Num_regs_LH + Num_regs_RH);
    Mean_sig_roi      = zeros(Num_imgs_LH, Num_regs_LH + Num_regs_RH);
    Num_vers_nonan_LH = zeros(1, Num_regs_LH);
    Num_vers_nonan_RH = zeros(1, Num_regs_RH);

    for i_img = 1:Num_imgs_LH
        
        Sub          = gifti([Dir_data{1}{i_sub} filesep File_name_LH(i_img).name]);
        Data_LH      = double(Sub.cdata);
        nonanMask_LH = ~isnan(Data_LH);
        
        for i_roi    = 1:Num_regs_LH
            Ind_roi  = Label_region_LH(Ind_ROI_LH(i_roi));
            Index    = (Label_LH == Ind_roi);
            Num_vers_nonan_LH(1, i_roi) = sum(Index & nonanMask_LH);
            Sig_roi{i_img, i_roi}       = Data_LH(Index & nonanMask_LH);
            Mean_sig_roi(i_img, i_roi)  = nanmean(Data_LH(Index));
        end
        
        Sub          = gifti([Dir_data{1}{i_sub} filesep File_name_RH(i_img).name]);
        Data_RH      = double(Sub.cdata);
        nonanMask_RH = ~isnan(Data_RH);
        
        for i_roi   = 1:Num_regs_RH
            Ind_roi = Label_region_RH(Ind_ROI_RH(i_roi));
            Index   = (Label_RH == Ind_roi);
            Num_vers_nonan_RH(1, i_roi)            = sum(Index & nonanMask_RH);
            Sig_roi{i_img, i_roi+Num_regs_LH}      = Data_RH(Index & nonanMask_RH);
            Mean_sig_roi(i_img, i_roi+Num_regs_LH) = nanmean(Data_RH(Index));
        end

    end

    Num_ver_reg(i_sub,:)         = cellfun('size', Sig_roi(1,:), 1);
    Removed_reg                  = Num_ver_reg(i_sub,:) < Num_thres;
    Removed_subs_reg(i_sub, :)   = Removed_reg;
    Sig_roi(:,Removed_reg)       = [];
    Mean_sig_roi(:,Removed_reg)  = [];
    
    Name_sub{i_sub, 1} = ['sub_' num2str(i_sub,'%04d') '.mat'];
    
    switch lower(Level)
        case 'all'
            save([Output_path filesep 'Signal_LR_'      Name_sub{i_sub, 1}], getname(Sig_roi));
        case 'mean'
            save([Output_path filesep 'mean_Signal_LR_' Name_sub{i_sub, 1}], getname(Mean_sig_roi));
        case 'both'
            save([Output_path filesep 'Signal_LR_'      Name_sub{i_sub, 1}], getname(Sig_roi));
            save([Output_path filesep 'mean_Signal_LR_' Name_sub{i_sub, 1}], getname(Mean_sig_roi));
    end
    
    fprintf('Extracting regional signal for %s ...... is done\n', [Dir_data{1}{i_sub}]);
    
end

%% Record number of vertices in each region
Ind_regs             = 1:Num_regs_LH + Num_regs_RH;
Ind_regs             = num2cell(Ind_regs');
Num_survived_ver     = num2cell([Num_vers_nonan_LH(:); Num_vers_nonan_RH(:)]);

%% Generate an excel file containing important atlas informations
Altas_name = 'Unknown';

Supported_atlas = gretna_label('surf');
Supported_atlas = Supported_atlas.abbr;

Num_ver_ind_LH   = zeros(Num_regs_LH, 1);
for i_roi        = 1:Num_regs_LH
    Ind_roi      = Label_region_LH(Ind_ROI_LH(i_roi));
    Num_ver_ind_LH(i_roi, 1) = sum(Label_LH == Ind_roi);
end

Num_ver_ind_RH   = zeros(Num_regs_RH, 1);
for i_roi        = 1:Num_regs_RH
    Ind_roi      = Label_region_RH(Ind_ROI_RH(i_roi));
    Num_ver_ind_RH(i_roi,1) = sum(Label_RH == Ind_roi);
end

for i_surf = 1:length(Supported_atlas)
    if ~isempty(strfind(Name_template_LH, Supported_atlas{i_surf}))
        Altas_name     = Supported_atlas{i_surf};
        Region_list    = gretna_label(Altas_name);
        Region_list_LH = Region_list.abbr(Ind_ROI_LH, 1);
        Region_list_RH = Region_list.abbr(Num_label_LH + Ind_ROI_RH,1);
        Region_list    = [Region_list_LH; Region_list_RH];
    end
end

if ~exist('Region_list','var')
    Region_list = cell(Num_regs_LH + Num_regs_RH, 1);
end

% Write .csv file
Output_switch  = mean(Removed_subs_reg, 1);
Cell_all_regs  = cat(2, Ind_regs, Region_list, Num_survived_ver);

if mean(bsxfun(@eq, Output_switch, abs(fix(Output_switch)))) == 1
    Num_rem_regs          = sum(~Removed_reg);
    Cell_ind_rem          = num2cell(1:Num_rem_regs);
    Cell_ori_ind_rem      = Cell_all_regs(~Removed_reg,1);
    Cell_reserved         = Cell_all_regs(~Removed_reg,2:3);
    Cell_reserved         = cat(2, Cell_ind_rem', Cell_ori_ind_rem, Cell_reserved);
    Ind_removed           = 1:sum(Removed_reg);
    Cell_ori_ind_removed  = num2cell(Ind_removed');
    Cell_removed          = Cell_all_regs(Removed_reg, :);
    
    Output_label_csv                              = cell(Num_regs_LH + Num_regs_RH+3, 13);
    Output_label_csv(4:end, 1:2)                  = Cell_all_regs(:,[1, 2]);
    Output_label_csv(4:end, 3)                    = num2cell([Num_ver_ind_LH; Num_ver_ind_RH]);
    Output_label_csv(4:3+Num_rem_regs, 5:8)       = Cell_reserved;
    Output_label_csv(4:3+sum(Removed_reg), 10)    = Cell_ori_ind_removed;
    Output_label_csv(4:3+sum(Removed_reg), 11:13) = Cell_removed;
    
    Output_label_csv(3, [1, 5, 10]) = {'Index'};
    Output_label_csv(3, [2, 7, 12]) = {'Label'};
    Output_label_csv(3, [3, 8, 13]) = {'Vertices'};
    Output_label_csv(3, [6, 11])    = {'Original_Index'};
    Output_label_csv(2, [1, 5, 10]) = {'Original', 'Reserved', 'Removed'};
else
    Output_label_csv       = cell(Num_regs_LH+Num_regs_RH+3, 3);
    Output_label_csv(3, 1) = {'Index'};
    Output_label_csv(3, 2) = {'Label'};
    Output_label_csv(3, 3) = {'Vertices'};
    Output_label_csv(2, 1) = {'Original'};
    
end

for i_sub = 1:Num_subs
    Cell_removed_regs{i_sub, 1} = find(Removed_subs_reg(i_sub, :)==1);
end

Output_label_csv(1,1) = {['The threshold you set is: ' num2str(Num_thres)]};
Output_subs_csv       = table(Dir_data{1}, Name_sub, sum(~Removed_subs_reg, 2), sum(Removed_subs_reg, 2), Cell_removed_regs,...
    'VariableNames', {'Location_file', 'ID_sub', 'Num_reserved_regs', 'Num_removed_regs', 'Ind_removed_regs'});

writetable(Output_subs_csv, [Output_path filesep sprintf('Atlas_%s_SubjectModification.csv', Altas_name)], 'WriteVariableNames', 1);
writetable(cell2table(Output_label_csv), [Output_path filesep sprintf('Atlas_%s_InfoSummary.csv', Altas_name)], 'WriteVariableNames', 0);

fprintf('\nPlease check the two .csv files in %s for important info!\n', Output_path);

return