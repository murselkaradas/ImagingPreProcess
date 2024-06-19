clc; 
clear all;
%%
set(0, 'DefaultLineLineWidth', 3);
set(0,'defaultTextFontSize', 10);
set(0,'defaultAxesFontSize',10);
set(0,'DefaultFigureColormap',jet)
set(0,'defaultfigurecolor',[1 1 1])
set(0,'DefaultAxesTitleFontWeight', 'normal')
set(0,'DefaultTextInterpreter','none')

addpath(genpath('/gpfs/scratch/karadm01/ImagingPreProcess'))
%%

folderpath = '/gpfs/scratch/karadm01/MultiScaleMapping/240613';
roifile ='20240613_RoiSet';
WSfieldname = 'OMP02_230207_HighConcWS';

fieldname = 'JG11571_240613_field1_2xDMD2P_DMDstim_00001';
VoyeurH5_file = '11571_1_01_D2024_6_13T17_33_31_odor';

get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false,'calculate_diff_image',true,...
   'isWS',false, 'inh_realign',false,'pre',1000,'post',2000,'usealigned_tiff',true)


fieldname = 'JG11571_240613_field1_2xDMD2P_DMDstim_00002';
VoyeurH5_file = '11571_1_02_D2024_6_13T18_22_23_odor';
get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false,'calculate_diff_image',true,...
   'isWS',false, 'inh_realign',false,'pre',1000,'post',2000,'usealigned_tiff',true)


fieldname = 'JG11571_240613_field1_2xSLM2P_cellseries_00001';
VoyeurH5_file = '11571_1_01_D2024_6_13T20_46_10_odor';
get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false,'calculate_diff_image',true,...
   'isWS',false, 'inh_realign',false,'pre',1000,'post',2000,'usealigned_tiff',true)


fieldname = 'JG11571_240613_field1_2xSLM2P_cellseries_00002';
VoyeurH5_file = '11571_1_02_D2024_6_13T21_23_27_odor';
get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false,'calculate_diff_image',true,...
   'isWS',false, 'inh_realign',false,'pre',1000,'post',2000,'usealigned_tiff',true)


fieldname = 'JG11571_240613_field1_2xSLM2P_cellseries_00003';
VoyeurH5_file = '11571_1_03_D2024_6_13T21_47_7_odor';
get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false,'calculate_diff_image',true,...
   'isWS',false, 'inh_realign',false,'pre',1000,'post',2000,'usealigned_tiff',true)


fieldname = 'JG11571_240613_field1_2xSLM2P_cellseries_00004';
VoyeurH5_file = '11571_1_04_D2024_6_13T22_39_1_odor';
get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false,'calculate_diff_image',true,...
   'isWS',false, 'inh_realign',false,'pre',1000,'post',2000,'usealigned_tiff',true)
