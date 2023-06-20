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


%%
addpath '/gpfs/scratch/karadm01/2Panalysis'
addpath '/gpfs/scratch/karadm01/2Panalysis/2P_splitSVD'
addpath '/gpfs/scratch/karadm01/WaveSurfer-1.0.2'
addpath(genpath('/gpfs/scratch/karadm01/breathmetrics-master'))
addpath '/gpfs/scratch/karadm01/SessionAnalysis'

addpath 'P:/2P/ImagingPreProcess/UtilFuncs'
addpath 'P:/2P/ImagingPreProcess/WaveSurfer-1.0.2'
addpath(genpath('P:/2P/ImagingPreProcess/breathmetrics-master'))
addpath 'P:/2P/ImagingPreProcess/SessionAnalysisScripts'


%%

folderpath = 'P:/2P/2Pdata\JG44524\230614';
roifile ='JG44524_230614_glomfield4_stim';


fieldname = 'JG44524_230614_glomfield4_stim_00001';
VoyeurH5_file = 'JG44524_230614_glomfield4_stim_1_01';
get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false,'calculate_diff_image',true)

%%
fieldname = 'JG44524_230614_glomfield4_seq_00001';
VoyeurH5_file = 'JG44524_230614_glomfield4_seq_1_01';
get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false)

fieldname = 'JG44524_230614_glomfield4_seq_00002';
VoyeurH5_file = 'JG44524_230614_glomfield4_seq_1_02';
get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false)
