function   get_prepocessed_stimdata(varargin)

% This functions read beahivor/imaging data and do preprocessing. It requires at least 
% tiff stacks, Voyeur generated H5 file and ROIs in either .zip or roi format
% It generates two different .mat file. The MAT file end by  _S_v73.mat is smaller and 
% compact. It has necessary data to be used later. Depends on user input it also generates heatmap
% and dFF traces for all ROIs
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'fieldname'       tiff stack name, use proper identifier
%                       we will use all tiffs with fieldname_
%     'VoyeurH5_name'   Voeyur H5 file   it will read  as dir([strcat(VoyeurH5_name,'*.h5')]);           
%     'roiname'         ROI file either in zip or roi format. It will read as dir([strcat(roiname,'*.zip')]);
%     'minDuration'     min ripple duration. Keeping this input nomenclature for backwards
%                       compatibility
%     'pre'             pre inhalation duration in ms (default = 2000)
%     'post'            post inhalation duration in ms (default = 4000)
%     'img_format'      image format for acquisition pixel x pixel (default = [512, 512])
%     'stim_cell'       stim cell  ROI IDs for plotting (default = [1,2,3])
%     'dfflim'          ylim for dff plots (default = [-0.3,0.6])
%     'inh_realign'     inhalation realign using breatmetric (default = false)'
%     'usealigned_tiff' do you have aligned tiff in aligned folder (default = yes)
%     'isWS'            use Wavesurfer (WS) available WS recording to determine frame numbers, (default = false)
%                       it is crucial for blanked recording. Behavior box drops frame occasinally. WS more reliable.
%     'isplot   '       do you want to plot dFF heatmap and traces? 
%     'kalman_gain   '  the strength of the filter [0 to 1]. Larger gain values means more
%                       aggressive filtering in time so a smoother function with a lower 
%                       peak. Gain values above 0.5 will weight the predicted value of the 
%                       pixel higher than the observed value
%     'calculate_diff_image   '  Generaate difference image. it is half second before 
%                       stim/inh. onset and 1 second after. 
%     'stimcorrection'  Stim Correction type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SET DEFAULT FREE PARAMETERS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   See also: get_prepocessed_odordata
p = inputParser;
addParameter(p, 'tiffpath', pwd, @isstr);
addParameter(p, 'fieldname', pwd, @isstr);
addParameter(p, 'roiname', pwd, @isstr);
addParameter(p, 'VoyeurH5_name', pwd, @isstr);

addParameter(p,'pre',2000, @isnumeric);
addParameter(p,'post',4000, @isnumeric);
addParameter(p,'Odor_duration',1, @isnumeric);

addParameter(p, 'img_format', [512,512], @isnumeric);

addParameter(p, 'stim_cell', [1,2,3], @isnumeric);
addParameter(p, 'dfflim', [-0.3,0.6], @isnumeric);

addParameter(p, 'inh_realign', false, @islogical);
addParameter(p, 'usealigned_tiff', true, @islogical);
addParameter(p, 'isOdor', true, @islogical);
addParameter(p, 'isplot', true, @islogical);

addParameter(p, 'isWS', false, @islogical);
addParameter(p, 'WSfieldname', pwd, @isstr);

addParameter(p, 'kalman_gain', 0.5,@isnumeric);
addParameter(p, 'calculate_diff_image', false, @islogical);

addParameter(p,'correct_stim', true, @islogical);
addParameter(p, 'stimcorrection', "filloutliers",@isstr);
addParameter(p, 'external_trigger_clicked', true,@islogical);
parse(p,varargin{:});

tiffpath = p.Results.tiffpath;
fieldname = p.Results.fieldname;
roiname = p.Results.roiname;
VoyeurH5_name = p.Results.VoyeurH5_name;
pre = p.Results.pre;
post = p.Results.post;
Odor_duration = p.Results.Odor_duration;
img_format = p.Results.img_format;
inh_realign = p.Results.inh_realign;
usealigned_tiff = p.Results.usealigned_tiff;
stim_cell = p.Results.stim_cell;

isOdor = p.Results.isOdor;
isplot = p.Results.isplot;

dfflim = p.Results.dfflim;
WSfieldname = p.Results.WSfieldname;
kalman_gain = p.Results.kalman_gain;
isWS = p.Results.isWS;
calculate_diff_image = p.Results.calculate_diff_image;
stimcorrection = p.Results.stimcorrection;
correct_stim = p.Results.correct_stim;
external_trigger_clicked = p.Results.external_trigger_clicked;
p.Results
 if img_format(1)==256
    fps =58.20;
else
    fps  = 30;
end
%%
cd(tiffpath);
H5Name = dir([strcat(VoyeurH5_name,'*.h5')]);
path_h5=H5Name.folder;
h5_name= H5Name.name;
%%
RoiName = dir([strcat(roiname,'*.zip')]);
pathroi = fullfile(RoiName.folder,RoiName.name);
[cellMask1,cellMask_sep]=create_ROImask_manual2(img_format,pathroi);%
cellMask_vec=reshape(cellMask_sep,[],size(cellMask_sep,3));
cellMask_vec=cellMask_vec./sum(cellMask_vec);
num_cell=size(cellMask_vec,2);
%%
if inh_realign && isWS
    WsH5name = dir([strcat(WSfieldname ,'*.h5')]);
    path_wsh5=WsH5name.folder;
    wsh5_name= WsH5name.name;

    wsdata = loadDataFile(wsh5_name);
    wsNtrials = size(fieldnames(wsdata),1)- 1; % First field is header
    AIchannelnames = deblank(string(cell2mat(wsdata.header.AIChannelNames)));
    FrameTriggerChannel = find(strcmp(AIchannelnames, 'FrameTrigger'));
    PMTgateChannel = find(strcmp(AIchannelnames, 'FrameSync_DMD'));
    FVChannel = find(strcmp(AIchannelnames, 'FinalValve'));
    PIDChannel = find(strcmp(AIchannelnames, 'PID'));
    ShutterChannel = find(strcmp(AIchannelnames, 'LaserShutter'));

    wsrate = wsdata.header.AcquisitionSampleRate/1e3;

    wsFrameNumbers = zeros(wsNtrials,1);
    ws_stimframe = zeros(wsNtrials,1);
    n=20;
    for i = 1:wsNtrials
        an = eval(sprintf('wsdata.sweep_%04d.analogScans',i));
        framediff= find(diff(an(:,1))>2);  %% Frame trigger is channel 1
        wsFrameNumbers(i) = length(find(diff(framediff/20)>(0.8*1e3/fps)))+1;
    end
    [sniff,sniff_smooth,frame_trigger_trial,frametrigger,Data,~]=Read_Trial_Info(h5_name,path_h5,pre,post,true,fps,wsFrameNumbers);
elseif inh_realign
    
    [sniff,sniff_smooth,frame_trigger_trial,frametrigger,Data,~]=Read_Trial_Info(h5_name,path_h5,'pre', 2000,'post', 4000,'fps',30, 'inh_detect',inh_realign);
else
    [sniff,frametrigger,Data,~]=read_sniff_frametrigger_trialinfo(h5_name,path_h5,pre,post,false,fps);
    sniff_smooth = sniff;

end
odor_duration = mean(Data.fvdur)/1e3;
figure(1);imagesc(sniff_smooth,[min(sniff_smooth(:))/1.5,max(sniff_smooth(:))/1.5]);
xlim([1500,3000])
colorbar
savefig(strcat(fieldname, '_sniff', '.fig'))
saveas(figure(1),strcat(fieldname, '_sniff', '.png'))
clear n an


%% In Case ScanImageTiffReader does not work use following 
%% Read TIFF files and generate fluorescence signal per ROI
if usealigned_tiff
    tiff_names = dir([strcat('aligned/',fieldname,'_*.tif')]);
else
    tiff_names = dir([strcat(fieldname,'_*.tif')]);
end

tiff_filenames = {tiff_names.name};
tiff_foldername = {tiff_names.folder};
filenum = size(tiff_names,1);
datamean = zeros(img_format);
data_time = [];
Fluo_cell = [];
Nframestart= 0 ;
datafull = [];
if filenum == 0
     error('Error. \nNo tiff file found.')
else
    disp(strcat(num2str(filenum),' tiff files will be loaded'));
    for i = 1: filenum
        i
        tiff_loaded = loadtiff(fullfile(tiff_foldername{i},tiff_filenames{i}));
        if calculate_diff_image
            datafull = cat(3, datafull, tiff_loaded);
        end
        data= double(tiff_loaded);
        Nframe = size(data,3);
        datamean = datamean + mean(data,3);
        %data_time = [data_time, mean(mean(data,1),1)];
        Fluo_cell =[Fluo_cell, double(cellMask_vec')*double(reshape(data,[img_format(1)*img_format(2),Nframe]))];
        
    end
end
%% STIM artefact correction

if correct_stim
    disp('Filtering traces without Stim Artefacts frames');
    stimframe_perroi= [];
    Voyeur_stim_frame = Data.stim_frame;
    Voyeur_stim_frame(Voyeur_stim_frame ==0) = [];
    Nstack = size(Fluo_cell,2);
    %in case, tiff file is shorter than Voyeur stim
    stim_frame_tifflimit = Voyeur_stim_frame((Voyeur_stim_frame<Nstack)  );
    stim_frame_tifflimit = stim_frame_tifflimit (stim_frame_tifflimit>1);
    for k = 1:num_cell
        stimframe_perroi(k,:) = findstimframe(squeeze(Fluo_cell(k,:,:)),stim_frame_tifflimit,1);
        Fluo_substracted{k} = Fluo_cell(k,:);
    end
    stimframes_allroi = findmostrepeatedstims(stimframe_perroi);
    stim_frame_mod = mode(stimframes_allroi,1);
    if stimcorrection == "filloutliers"
        Fluo_raw = Fluo_cell;
        Nframes = size(Fluo_cell,2);
        for cellid =1:num_cell
            for stimid = 1: length(stim_frame_mod)
                stimbound = stim_frame_mod(stimid)-round(fps/2):stim_frame_mod(stimid)+round(3*fps/2);
                if stimbound(1)>1 && stimbound(end)<Nframes
                    dff_stimbound = Fluo_raw(cellid,stimbound);
                    outlier_removed = filloutliers(dff_stimbound,"pchip","percentiles",[3 97]);      
                    Fluo_raw(cellid,stimbound(2:end-1)) = outlier_removed(2:end-1); %edges can be problematic if outliers are there
                end
            end
            fluotemp = reshape(Kalman_Stack_Filter(Fluo_raw(cellid,:) , kalman_gain, 0.5),1,[]);
            Fluo_cell_Kalman(cellid,:)  = fluotemp;
            Fluo_cell_stimcorrected(cellid,:) = Fluo_raw(cellid,:);
        end
        disp('FillOutliers is used to correct stim frames.');
    else
        stims = cat(1, stim_frame_mod-1, stim_frame_mod, stim_frame_mod+1);
        stim_frames = stims(:);
        Fluo_cell_Kalman = Fluo_cell;
        Fluo_cell_stimcorrected = Fluo_cell;
        for j =1:num_cell
            Fluo_substracted{j}(stim_frames) = [];
            Fluo_raw = reshape(Fluo_substracted{j},1,1,[]);
            fluotemp = reshape(Kalman_Stack_Filter(Fluo_raw , kalman_gain, 0.5),1,[]);
            fluotemp1 = insertstimframeback(fluotemp,Fluo_cell(j,:), stim_frames,[]);
            Fluo_cell_Kalman(j,:)  = fluotemp1;
            Fluo_cell_stimcorrected(j,:) = insertstimframeback(squeeze(Fluo_raw)', Fluo_cell(j,:), stim_frames,[]);
        end
        disp('Kalman Filtered without stim artefacts.Stim Frames are inserted back');

    end
else

    Fluo_cell_Kalman = Fluo_cell;
    Fluo_cell_stimcorrected = Fluo_cell;
    disp('No stim artefact correction is done');
end
%%
clear data1 data  Fluo_substracted Fluo_raw fluotemp option

%% SAVE AVG tiff
img = repmat(imadjust(mat2gray(datamean)),1,1,3)*0.8;
opt = 1;
img(:,:,2) = img(:,:,2)+double(logical(cellMask1))*0.1;
figure(22);imagesc(img);CenterFromRoiMasks(cellMask1,1:num_cell,opt);axis square
savefig(strcat(fieldname, 'ROI', '.fig'))
saveas(figure(22),strcat(fieldname, 'ROI', '.png'))
options.overwrite = true;
saveastiff(int16(datamean),strcat(fieldname,'_AVG.tif'), options)
disp('Avg tiff file has been saved');

if ~external_trigger_clicked
   Nremove = sum(Data.raw_frame_triggers(1:end-100) ==0);
   if Nremove > 0
       fprintf('%d Tiffs were not included processing\n',Nremove );
       Fluo_cell_Kalman = Fluo_cell_Kalman(:,Nremove+1:end);
       datafull = datafull(:,:,Nremove+1:end);
   end
end
%% Find Included Trials
Nframe = size(frametrigger,1);
pre_inh=round(pre*fps/1e3);
post_inh=round(post*fps/1e3);
if isOdor
    inh_onset = Data.inh_onset;
    disp('Using Inh. onset time');
else
    % for no odor blocks, there is a delay to ensure mechanical shutter is
    % completely opened, usually it is around 50 ms.
    if isfield(Data, 'laserontime_1st')
        Data.laserontime= Data.laserontime_1st;
    end
    inh_onset = Data.inh_onset + mode(Data.laserontime- Data.inh_onset);
    disp('Using Laser onset time');
end
% inh_onset = Data.inh_onset;
up_sampling_fac=(1000/fps);

%% Kalman  Filters result
FKalman = [];
dffKalman = [];
meanAllKalman = [];
ind =1;
img_df = [];
for tr=1:length(inh_onset)
    [~,inh_frame(tr)]=min(abs(frametrigger-double(inh_onset(tr))));%frame in which inh_onset is included
    time_diff = frametrigger(inh_frame(tr))-double(inh_onset(tr));
    if ~isempty(inh_frame(tr))
        sv_frame_range=inh_frame(tr)-pre_inh:inh_frame(tr)+post_inh-1;
        inh_diff = abs(frametrigger(inh_frame(tr))-double(inh_onset(tr)));
        if max(sv_frame_range) < size(Fluo_cell,2) && (inh_frame(tr)>pre_inh) && inh_diff<1e3
            fcellKalman=Fluo_cell_Kalman(:,sv_frame_range);
            baseline_frame=pre_inh-32:pre_inh-2; % be conservative
            dffKalman(:,:,ind)=(fcellKalman-mean(fcellKalman(:,baseline_frame),2))./mean(fcellKalman(:,baseline_frame),2);
            FKalman(:,:,ind) = fcellKalman;
            trials_read(tr)=true;
            meanAllKalman(:,ind) = (mean(fcellKalman(:,baseline_frame),2));
            
            if calculate_diff_image
                img_trial = double(datafull(:,:,sv_frame_range));
                img_baseline = (mean(img_trial(:,:,baseline_frame),3));
                img_df(:,:,:,ind) = (img_trial(:,:,pre_inh-round(fps/2):pre_inh+round(fps)-1) -img_baseline);
            end
            ind = ind+1;
        else
            trials_read(tr)=false;
            fprintf('trial %d inh_frame(tr) was not included in tiff stack\n',tr)

        end
    else
        trials_read(tr)=false;
        fprintf('trial %d inh_frame(tr) was not included in tiff stack\n',tr)
    end
    
end
disp('dFF of each trial is calculated!')

%%
if isOdor
    OdorInfo = HDF5_getOdors(path_h5,h5_name,trials_read,10);
else
    OdorInfo = HDF5_getStimID(path_h5,h5_name,'trials_read', trials_read);
end

if calculate_diff_image
    fig55 = figure2('diff_image');
    p=numSubplots(size(OdorInfo.odors,1));
    index = (reshape(1:p(2)*p(1),p(2),p(1)).');
    img_df_percond = [];
    std_df = double(std(single(img_df(:))));
    for i = 1: size(OdorInfo.odors,1)
        img_df_percond(:,:,:,i) =  int16(mean(img_df(:,:,:,OdorInfo.odorTrials{i}),4));
        ii = index(i);
        subplot(p(1),p(2),ii)
        imagesc(imgaussfilt(mean(mean(img_df(:,:,:,OdorInfo.odorTrials{i}),4),3),1));axis square
        caxis([-1*std_df/2,std_df/2]);
        colormap(NegativeEnhancingColormap(256, [-1*std_df/2,std_df/2], ...
            [0 0 1], [1 0 0],[1 1 1], 1.25));
        title(OdorInfo.odors{i})
    end
    saveas(fig55,strcat(fieldname, 'diff_images', '.png'))
    %clear datafull img_trial img_baseline img_df
end
%%
if isplot
    Sniff_trial = -2*sniff(trials_read,:)./peak2peak(sniff(trials_read,:),2);
    framelim = -pre_inh:post_inh-1;
    cellid = 1:num_cell;
    sniff_ds = downsample(Sniff_trial', round(1000/fps));
    p=numSubplots(size(OdorInfo.odors,1));
    % p = [4,1];
    index = (reshape(1:p(2)*p(1),p(2),p(1)).');
    fig4 = figure2('map');
    fig5 = figure2('dff');
    fig6 = figure2('dffstim_cell');
    for i = 1: size(OdorInfo.odors,1)
        figure(fig4.Number)
        ii = index(i);
        subplot(p(1),p(2),ii)
        dff1 = mean(dffKalman(:,:,OdorInfo.odorTrials{i}),3);
        imagesc(framelim, cellid, dff1)
        title(OdorInfo.odors{i})
        caxis(dfflim)
        colormap(NegativeEnhancingColormap(256, dfflim, ...
            [0 0 1], [1 0 0],[1 1 1], 1.25));
        colorbar
        ylabel('Cell ID')
        xlabel('Frame number')
        hold on
        plot([1 1]*0, ylim, '--b','LineWidth',1.5)
        plot([1 1]*odor_duration*fps, ylim, '--b','LineWidth',1.5)
        xlim([-fps 2*fps])
        hold off
        
        figure(fig5.Number)
        subplot(p(1),p(2),ii)
        plot(framelim, dff1(:,:)', 'LineWidth',2)
        hold on
        ylim(dfflim)
        plot([1 1]*0, ylim, '--b', 'LineWidth',2)
        plot([1 1]*odor_duration*fps, ylim, '--b','LineWidth',2)
        plot(downsample((-pre:post-1),round(1000/fps))./round(1000/fps), mean(sniff_ds(:,OdorInfo.odorTrials{i}),2)*0.25-0.3, 'k', 'LineWidth',2)
        hold off
        title(OdorInfo.odors{i})
        xlabel('#')
        ylabel('\DeltaF/F_0')
        xlim([-fps 2*fps])
        
        figure(fig6.Number)
        subplot(p(1),p(2),ii)
        plot(framelim, dff1(stim_cell,:)', 'LineWidth',2)
        hold on
        ylim(dfflim)
        plot([1 1]*0, ylim, '--b', 'LineWidth',2)
        plot([1 1]*odor_duration*fps, ylim, '--b','LineWidth',2)
        plot(downsample((-pre:post-1),round(1000/fps))./round(1000/fps), mean(sniff_ds(:,OdorInfo.odorTrials{i}),2)*0.3-0.25, 'k', 'LineWidth',2)
        hold off
        title(OdorInfo.odors{i})
        xlabel('#')
        ylabel('\DeltaF/F_0')
        xlim([-fps 2*fps])

    end
    %savefig(fig6, strcat(fieldname, '_DFFTrace_stim_cellKalman', '.fig'))
    saveas(fig6, strcat(fieldname, '_DFFTrace_stim_cellKalman', '.png'))
    %savefig(fig5,strcat(fieldname, '_DFFTrace_StimKalman', '.fig'))
    saveas(fig5,strcat(fieldname, '_DFFTrace_Kalman', '.png'))
    %savefig(fig4,strcat(fieldname, '_DFFMAPKalman', '.fig'))
    saveas(fig4,strcat(fieldname, '_DFFMAPKalman', '.png'))
    %  SAVE all workspace
    clear options opt option j k i fig1 fig2 fig3 fig4 fig5 fig6
    close all

    %

end
%% SAVE everything in workspace, it is usefull in some cases
save(strcat(fieldname, '.mat'));

%% Data format, I used it to share data with Jon and Saeed. My Python scripts are written based on this format
Session.OdorResponse = {};
Session.F= Fluo_cell_Kalman';
Session.blockTrials = {};
Session.fieldname = fieldname;
for i = 1: size(OdorInfo.odors,1)
    Session.OdorResponse{i} = permute(dffKalman(:,:,OdorInfo.odorTrials{i}),[2,1,3]);
    Session.blockTrials{i} = ones(length(OdorInfo.odorTrials{i}),1);
end
if calculate_diff_image
   Session.diff_image = img_df_percond;
end
Session.UniqueConds = cellstr(OdorInfo.odors);
Session.OdorTrials = OdorInfo.odorTrials;
Session.Sniffs = Sniff_trial';
Session.CellMask = cellMask_vec;

Session.Infos.fps = fps;
Session.Infos.ImgFormat = img_format;
Session.Infos.imgwithROIs = img;
Session.Infos.pre_inh = pre_inh;
Session.Infos.post_inh = post_inh;
Session.Infos.TrialsRead = trials_read;
Session.datamean = datamean; 
Session.data_time = data_time;
Session.inhframe = inh_frame;
Session.VoyeurData = Data;
save(strcat(fieldname,'_S_v73.mat'), 'Session','-v7.3')
disp("Saving DONE!!")
end