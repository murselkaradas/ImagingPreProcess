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

addpath(genpath('/gpfs/data/rinberglab/Mursel/ImagingPreProcess/'))
%%

folderpath = '/gpfs/data/rinberglab/Mursel/240520';
roifile ='RoiSet';


fieldname = 'JG41851_240520_field1_00001';
VoyeurH5_file = 'JG41851_240520_field1_1_01_D2024_5_20T17_57_7_odor';
WSfieldname = 'OMP02_230207_HighConcWS';

get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false,'calculate_diff_image',true,...
   'isWS',false, 'inh_realign',false,'pre',1000,'post',2000,'usealigned_tiff',false)

%%
%fieldname = 'JG44524_230614_glomfield4_seq_00001';
%VoyeurH5_file = 'JG44524_230614_glomfield4_seq_1_01';
%get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false)

%fieldname = 'JG44524_230614_glomfield4_seq_00002';
%VoyeurH5_file = 'JG44524_230614_glomfield4_seq_1_02';
%get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false)
