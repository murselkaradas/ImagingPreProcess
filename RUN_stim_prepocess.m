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

%addpath(genpath('/gpfs/scratch/karadm01/ImagingPreProcess/'))
%%

folderpath = 'E:\SingleGlom\Mouse00691\Imaging_2900V_231019';
roifile ='SC69B_MCImaging_ROI09_2900V_231019_RoiSet';


fieldname = 'SC69B_MCImaging_ROI09_2900V_231019';
VoyeurH5_file = 'mouse0691_sess01_2900V_D231019';
WSfieldname = 'OMP02_230207_HighConcWS';

get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false,'calculate_diff_image',true,'inh_realign',true)

%%
%fieldname = 'JG44524_230614_glomfield4_seq_00001';
%VoyeurH5_file = 'JG44524_230614_glomfield4_seq_1_01';
%get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false)

%fieldname = 'JG44524_230614_glomfield4_seq_00002';
%VoyeurH5_file = 'JG44524_230614_glomfield4_seq_1_02';
%get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false)
