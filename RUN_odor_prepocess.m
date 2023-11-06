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

folderpath = 'G:/.shortcut-targets-by-id/1Z9x_DTUNRkFooVBWp3g69Z9A7D8lPZJB/sampleData/130723';
roifile ='mouse0049_230713_triangle';


fieldname = 'mouse0049_230713_triangle';
VoyeurH5_file = 'mouse0049_1_01_D2023';
WSfieldname = 'mouse0049_230713_triangleWS';
get_prepocessed_odordata('inh_realign',false,'tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file, 'WSfieldname',WSfieldname,'isOdor',true,'usealigned_tiff', true)

