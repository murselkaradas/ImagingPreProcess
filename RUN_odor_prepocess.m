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

folderpath = '/gpfs/scratch/karadm01/2Pdata/M72/9280/250605';
roifile ='RoiSet';


fieldname = '9280_250605_OdorMapping';
VoyeurH5_file = '9280_250605_OdorMapping_1_01_D2025_6_5T14_49_15_odor';
WSfieldname = '9280_250605_OdorMappingWS_0001-0211';

get_prepocessed_odordata('inh_realign',true,'tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file, 'WSfieldname',WSfieldname,'isOdor',true,'usealigned_tiff', false, 'calculate_diff_image', true,'dfflim', [-0.1,0.2])

