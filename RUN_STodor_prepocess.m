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

addpath(genpath('/gpfs/scratch/karadm01/ImagingPreProcess//'))
%%

folderpath = '/gpfs/scratch/karadm01/2Pdata/STolfa/SK6362/250703';
roifile ='RoiSet_6362_TestifWorkst';


fieldname = 'SK6362_250516_STOlfaTestField1';
VoyeurH5_file = 'SK6362_250516_STOlfaTestField1_1_01_D2025_7_3T10_36_5_odor';
WSfieldname = 'SK6362_250516_STOlfaTestField1_0001-0094';

get_prepocessed_STodordata('inh_realign',true,'tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file, 'WSfieldname',WSfieldname,'isOdor',true,'usealigned_tiff', false, 'calculate_diff_image', true,'dfflim', [-0.2,1.0])

