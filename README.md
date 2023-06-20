# ImagingPreProcess
 2P imaging preprocess functions I used in RinbergLab

It relies on the following packages:

NORMCorre: This package is used for motion correction. https://github.com/flatironinstitute/NoRMCorre

Wavesurfer: This package is required for blanked odor recording. Each trial has its own TIF file. For continuous recording, this package is not necessary. https://github.com/JaneliaSciComp/Wavesurfer

Breathmetrics: This package is optional  and  used for sniff processing. It is useful if you want to realign inhalation onset. Voyeur's inhalation onset can be slightly inaccurate. https://github.com/zelanolab/breathmetrics


## Odor preprocessing 
### Odor Imaging :  Continuous acquisition

**Must have**: ScanImage Tiff Files  and Voyeur h5 File
*Optional*: Breathmetric

### Odor Imaging :  Blanked Recording
Acquisition is discrete

**Must have**:  ScanImage Tiff Files  and Voyeur h5 File and Wavesurfer H5 

*Optional* : Breathmetric


### Only Stim

```matlab
addpath(genpath('/gpfs/scratch/karadm01/ImagingPreProcess/'))
%%
tiffpath = '/gpfs/data/rinberglab/Jon/JG44524/230614';
roifile ='JG44524_230614_glomfield4_stim';

fieldname = 'JG44524_230614_glomfield4_stim_00001';
VoyeurH5_file = 'JG44524_230614_glomfield4_stim_1_01';
get_prepocessed_stimdata('tiffpath',folderpath, 'fieldname', fieldname, 'roiname', roifile,'VoyeurH5_name', VoyeurH5_file,'isOdor',false,'calculate_diff_image',true)
```


All above code will generate .MAT files
```matlab
%% Data format, I used it to share data with Jon and Saeed. My Python scripts are also written based on this format
Session.OdorResponse = {};
Session.F= Fluo_cell_Kalman';
Session.blockTrials = {};
Session.fieldname = fieldname;

for i = 1: size(OdorInfo.odors,1)
    Session.OdorResponse{i} = permute(dffKalman(:,:,OdorInfo.odorTrials{i}),[2,1,3]);
    Session.blockTrials{i} = ones(length(OdorInfo.odorTrials{i}),1);
end
% df images [-0.5s 1s] wrt inhalation onset
if calculate_diff_image
   Session.diff_image = img_df_percond;
end
% Odor/stim conditions
Session.UniqueConds = cellstr(OdorInfo.odors);
% Odor/stim trial numbers
Session.OdorTrials = OdorInfo.odorTrials;
%sniff per trial
Session.Sniffs = Sniff_trial';
Session.CellMask = cellMask_vec;

Session.Infos.fps = fps;
Session.Infos.ImgFormat = img_format;
Session.Infos.imgwithROIs = img;
Session.Infos.pre_inh = pre_inh;
Session.Infos.post_inh = post_inh;
Session.Infos.TrialsRead = trials_read;
%  Voyeur data
Session.VoyeurData = Data;
save(strcat(fieldname,'_S_v73.mat'), 'Session','-v7.3')
```

MAT file has following parameters
| Parameter name | Description |
|----------------|-------------|
| ```OdorResponse``` | dF/F responses, 1X N unique odor/stim condition cell array, every cell is 3D array, Ntime x Nroi x Ntrial   |
| ```F``` | Fluorescence data, NtotalFrame x Nroi|
| ```overlap_pre```| size of overlapping region in each direction before upsampling  |
| ```fieldname```    | fieldname is retrieved from tiffs |
| ```UniqueConds ``` | Unique trial condition presented |
| ```OdorTrials``` | Odor trials for each unique conditions | 
| ```Sniffs``` | Ntime x Ntrials sniff data, [-pre  post] wrt inhalation onset |
| ```CellMask``` | Cell masks for each ROI, NxNy  x Nroi |
| ```VoyeurData``` | Voyeur data with prepocessed infos |
| ```Infos``` | it has extra information regarding pre, post inhalation duration, fps, image format and trial read infos 
