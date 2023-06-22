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

folderpath = '/gpfs/scratch/karadm01/2Pdata/OMP02/230207';
roifile ='OMP02_230207';


fieldname = 'OMP02_230207_HighConc';
VoyeurH5_file = 'OMP02_230207_HighConc_1_01';
WSfieldname = 'OMP02_230207_HighConcWS';
get_prepocessed_odordata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file, 'WSfieldname',WSfieldname,'isOdor',true,'usealigned_tiff', false)


fieldname = 'OMP02_230207_MidConc';
VoyeurH5_file = 'OMP02_230207_MidConc_1_01';
WSfieldname = 'OMP02_230207_MidConcWS';
get_prepocessed_odordata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file, 'WSfieldname',WSfieldname,'isOdor',true,'usealigned_tiff', false)


fieldname = 'OMP02_230207_LowConc';
VoyeurH5_file = 'OMP02_230207_LowConc_1_01';
WSfieldname = 'OMP02_230207_LowConcWS';
get_prepocessed_odordata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file, 'WSfieldname',WSfieldname,'isOdor',true,'usealigned_tiff', false)
