function   get_prepocessed_odordata(varargin)

% This functions read behavior/imaging data and do necessary preprocessing. It requires at least 
% tiff stacks, Voyeur generated H5 file and ROIs in either .zip or roi format
% It generates two different .mat file. The MAT file end by  _S_v73.mat is smaller and more 
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
%     'inh_realign'     inhalation realign using breatmetric (default = true)'
%     'usealigned_tiff' do you have aligned tiff in aligned folder (default = yes)
%     'isWS'            use Wavesurfer (WS) available WS recording to determine frame numbers, (default = false)
%                       it is crucial for blanked recording. Behavior box drops frame occasinally. WS more reliable.
%     'isplot   '       do you want to plot dFF heatmap and traces? 
%     'kalman_gain   '  the strength of the filter [0 to 1]. Larger gain values means more
%                       aggressive filtering in time so a smoother function with a lower 
%                       peak. Gain values above 0.5 will weight the predicted value of the 
%                       pixel higher than the observed value
%     'calculate_diff_image   '  Generate difference image. it is half second before 
%                       stim/inh. onset and 1 second after. 
%     'stimcorrection'  Stim Correction type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SET DEFAULT FREE PARAMETERS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   See also: get_prepocessed_stimdata
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

addParameter(p, 'inh_realign', true, @islogical);
addParameter(p, 'usealigned_tiff', true, @islogical);
addParameter(p, 'isOdor', true, @islogical);
addParameter(p, 'isplot', true, @islogical);

addParameter(p, 'isWS', true, @islogical);
addParameter(p, 'WSfieldname', pwd, @isstr);

addParameter(p, 'kalman_gain', 0.5,@isnumeric);
addParameter(p, 'calculate_diff_image', false, @islogical);


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
isWS = p.Results.isWS;
dfflim = p.Results.dfflim;
WSfieldname = p.Results.WSfieldname;
kalman_gain = p.Results.kalman_gain;
calculate_diff_image = p.Results.calculate_diff_image;

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
    
    [sniff,sniff_smooth,frame_trigger_trial,frametrigger,Data,~]=Read_Trial_Info(h5_name,path_h5,pre,post,true,fps,[]);
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
Fluo_cell = [];
Nframestart= 0 ;
datafull = [];
Fluo_cell_stimcorrected = [];
Fluo_cell_Kalman =[];
if filenum == 0
     error('Error. \n No tiff file found.')
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
        Fluo_trial = double(cellMask_vec')*double(reshape(data,[img_format(1)*img_format(2),Nframe]));
        datamean = datamean + mean(data,3)/filenum;
        Fluo_cell =[Fluo_cell, Fluo_trial];
        Fluo_raw = zeros(1,num_cell,size(Fluo_trial,2));
        Fluo_raw(1,:,:) = Fluo_trial;
        Fluo_cell_Kalman = [Fluo_cell_Kalman,squeeze(Kalman_Stack_Filter(double(Fluo_raw),kalman_gain,0.5))];
        Fluo_cell_stimcorrected=[Fluo_cell_stimcorrected Fluo_cell];
    end
end
clear  data  Fluo_substracted Fluo_raw fluotemp option tiff_loaded

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
    inh_onset = Data.inh_onset + mode(Data.laserontime- Data.inh_onset);
    disp('Using Laser onset time');
end
up_sampling_fac=(1000/fps);

ind = 1;
FKalman = [];
dffKalman = [];
meanAllKalman = [];
Inh_frame_all = [];
for tr=1:length(inh_onset)
    Frame_offset = sum(Data.Voyeur_frame_numbers(1:tr-1));
    inh = inh_onset(tr);
	%% Consider to change here use either first or second line. I mostly ysed first.
    inh_frame=find(frame_trigger_trial{tr}<inh, 1, 'last' ) + Frame_offset;%frame in which inh_onset is included
	% [~,inh_frame]=min(abs(frametrigger-double(inh_onset(tr))))+ Frame_offset;
	
    if ~isempty(inh_frame)
        Inh_frame_all(tr)= inh_frame;
        sv_frame_range=inh_frame-pre_inh:inh_frame+post_inh-1;
        inh_diff = frametrigger(inh_frame)-double(inh);
        if max(sv_frame_range) < size(Fluo_cell,2) && (inh_frame>pre_inh) && abs(inh_diff)<1e3
            fcellKalman=Fluo_cell_Kalman(:,sv_frame_range); %upsample imaging data
            baseline_frame=pre_inh-32:pre_inh-2; %
            dffKalman(:,:,ind)=(fcellKalman-mean(fcellKalman(:,baseline_frame),2))./mean(fcellKalman(:,baseline_frame),2);
            FKalman(:,:,ind) = fcellKalman;
            trials_read(tr) = true;
            meanAllKalman(:,ind) = (mean(fcellKalman(:,baseline_frame),2));
            if calculate_diff_image
                img_trial = double(datafull(:,:,sv_frame_range));
                img_baseline = (mean(img_trial(:,:,baseline_frame),3));
                img_df(:,:,:,ind) = (img_trial(:,:,pre_inh-round(fps/2):pre_inh+round(fps)-1) -img_baseline);
            end
            ind = ind+1;
        else
            trials_read(tr)=false;
            fprintf('trial %d inh_frame was not included in tiff stack',tr)

        end
    else
        trials_read(tr)=false;
        fprintf('trial %d inh_frame was not included in tiff stack',tr)
    end
    
end
disp('dFF of each trial is calculated!')

%%
if isOdor
    OdorInfo = HDF5_getOdors(path_h5,h5_name,trials_read,10);
else
    OdorInfo = HDF5_getStimID(path_h5,h5_name,trials_read,10);
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
    clear datafull img_trial img_baseline img_df
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
    for j = 1:length(stim_cell)
        cellid = stim_cell(j);
        options.x_axis = framelim;
        options.color_area = [128 193 219]./255;
        options.color_line = 'r';
        options.alpha = 0.5;
        options.line_width =2;
        options.error = 'sem';
        fig7 = figure2('dffstim_cell');
        for i = 1: size(OdorInfo.odors,1)
            ii = index(i);
            dff_singlecell = squeeze(dffKalman(cellid,:,OdorInfo.odorTrials{i}));

            figure(fig7.Number)
            subplot(p(1),p(2),ii)
            plotareaerrorbar(dff_singlecell(:,:)', gcf, options)
        %     plot(framelim, dff_singlecell(:,1:6)', 'LineWidth',2)
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
        %savefig(fig7, strcat(fieldname, '_DFFTrace_stim_cellIDsem_',num2str(cellid),'.fig'))
        saveas(fig7, strcat(fieldname, '_DFFTrace_stim_cellIDsem_',num2str(cellid),'.png'))

    end
    close all

end

%%% Read PID data
Ntrials = sum(trials_read)
PIDmeas = zeros(Ntrials,7*wsrate*1e3);
kk = 1
trials_read_2 = trials_read(2:end)
for i = 1:wsNtrials
    if trials_read(i)
        an = eval(sprintf('wsdata.sweep_%04d.analogScans',i));
        framediff= find(diff(an(:,1))>2);  %% Frame trigger is channel 1
        FVdiff= find(diff(an(:,FVChannel))>2);  %% Frame trigger is channel 1
        if isempty(FVdiff)
            FVdiff1 = 65e3;
        elseif ~isempty(FVdiff)
           [M,I]= (min(abs(FVdiff)));
           FVdiff1 = FVdiff(I);
        else
           FVdiff1 = FVdiff(1);
        end
        startPID = FVdiff1 - (wsdata.header.AcquisitionSampleRate)*3;
        endPID = FVdiff1 +(wsdata.header.AcquisitionSampleRate)*4;
        if startPID<0
            PIDmeas(kk,:) = vertcat(an(1,PIDChannel)*ones(abs(startPID)+1,1),an(1:endPID-1,PIDChannel));
        else
            PIDmeas(kk,:) = an(startPID:endPID-1,PIDChannel);
        end
        kk = kk + 1
    end
end

fig3 = figure2('PID');
for i = 1: size(OdorInfo.odors,1)
    ii = index(i);
    figure(fig3.Number)
    subplot(p(1),p(2),ii)
    pid_plot = downsample(PIDmeas(OdorInfo.odorTrials{i},:)',20);
    plot(pid_plot, 'LineWidth',2)
    title(OdorInfo.odors{i})
    xlabel('#')
    ylabel('PID_voltage')
    legend()
end
saveas(fig3,strcat(fieldname, '_PID', '.png'))

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
    if exist('PIDmeas','var')
        Session.PID{i} = downsample(PIDmeas(OdorInfo.odorTrials{i},:)',20);
    end
end
if calculate_diff_image
   Session.diff_image = img_df_percond;
end
Session.UniqueConds = cellstr(OdorInfo.odors);
Session.OdorTrials = OdorInfo.odorTrials;
Session.Sniffs = Sniff_trial';
Session.CellMask = cellMask_vec;
Session.InhFrames = Inh_frame_all;
Session.Infos.OdorDuration = OdorDuration;

Session.Infos.fps = fps;
Session.Infos.ImgFormat = img_format;
Session.Infos.imgwithROIs = img;
Session.Infos.pre_inh = pre_inh;
Session.Infos.post_inh = post_inh;
Session.Infos.TrialsRead = trials_read;
Session.VoyeurData = Data;
save(strcat(fieldname,'_S_v73.mat'), 'Session','-v7.3')
disp("Saving DONE!!")
end