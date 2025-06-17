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

addpath(genpath('/gpfs/scratch/karadm01/ImagingPreProcess/'))
%%

folderpath = '/gpfs/scratch/karadm01/2Pdata/SC92/250520/';
roifile ='RoiSet';

fieldname = 'SC92_250520_Field12SEQ';
VoyeurH5_file = 'SC92_250520_Field12SEQ_1_01_D2025_5_20T17_35_45_odor';
WSfieldname = 'OMP02_230207_HighConcWS';

get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false,'calculate_diff_image',true,...
   'isWS',false, 'inh_realign',false,'pre',1000,'post',2000,'usealigned_tiff',false,'stim_cell',[1,2,3,4])






