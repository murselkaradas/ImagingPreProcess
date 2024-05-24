function [Sniff,Sniff_smooth,frame_trigger_trial,frame_trigger,Voyeurdata,Sniff_time]=Read_Trial_Info(h5,path_used,varargin)

% This functions read behavior/imaging data and do preprocessing for frames and sniffs trials. It requires at least 
% Voyeur generated H5 file 
%   [Sniff,Sniff_smooth,frame_trigger_trial,frame_trigger,data,Sniff_time]=Read_Trial_Info(h5,path_used,varargin)
%   Read_Trial_Info reads h5 file and extract sniff traces, inhalation
%   timing, frametrigger, and other trial information
%   Input:
%       h5: h5 file name
%       path_used: path where h5 file is located
%       varargin:
%           pre: time before inhalation onset (default: 1000ms)
%           post: time after inhalation onset (default: 1000ms)
%           inh_detect: true/false, if true, it detects inhalation onset
%           using sniff trace (default: true)
%           fps: frame per second (default: 30)
%           pmtblank_dur: pmt blank duration (default: 16ms)
%           wsFrameNumbers: frame numbers from wavesurfer (default: [])
%   Output:
%       Sniff: sniff traces for each trial
%       Sniff_smooth: sniff traces after baseline correction
%       frame_trigger_trial: frametrigger for each trial
%       frame_trigger: frametrigger for all trials
%       Voyeurdata: trial information
%       Sniff_time: time for each sniff trace

% Mursel Karadas 2022

%   See also: get_prepocessed_odordata
p = inputParser;
addParameter(p, 'pre', 1000, @isnumeric);
addParameter(p, 'post', 1000, @isnumeric);
addParameter(p, 'inh_detect', true,@islogical);
addParameter(p, 'fps',30, @isnumeric);
addParameter(p, 'pmtblank_dur', 16, @isnumeric);
addParameter(p,'wsFrameNumbers',[],@isnumeric);

parse(p,varargin{:});
pre = p.Results.pre;
post = p.Results.post;
inh_detect = p.Results.inh_detect;
fps = p.Results.fps;
pmtblank_dur = p.Results.pmtblank_dur;
wsFrameNumbers = p.Results.wsFrameNumbers;

p.Results

path_orig=pwd;
cd(path_used)

try
    Voyeurdata=h5read(h5,'/Trials');
catch
    error('%s not found in %s',h5, path_used)
end
%DMTS task has different field names
if isfield(Voyeurdata, 'inh_onset_1st')
    Voyeurdata.inh_onset = Voyeurdata.inh_onset_1st;
    Voyeurdata.fvOnTime = Voyeurdata.fvOnTime_1st;
end

% check if there is any wsFrameNumbers coming from wavesurfer
if isempty(wsFrameNumbers)
    Blanked_recording = false;
else
    Blanked_recording = true;
end
%Need to check if this manipulation is ok for 2p Voyeurdata
inh_onset=double(Voyeurdata.inh_onset);
num_trial=length(inh_onset);

h_info=h5info(h5);
h_info=h_info.Groups;
% Keys is Number_of_trials x 1 
Keys={h_info.Name};%h5 file from 2p rig doesn't record trial0 analog signal
fvOnTime=Voyeurdata.fvOnTime;
%%
%obtain frametrigger
record_onset=zeros(size(inh_onset));
frame_trigger= [];
frame_tol = 1e3/fps*1.2;  
delta_t = 1e3/fps;
frame_fv = zeros(num_trial,1);
for i=1:num_trial
    %i
    data_Events=h5read(h5,strcat(Keys{i},'/Events'));

    packet_sent_time=data_Events.packet_sent_time;
    sniff_samples=data_Events.sniff_samples;    %array for the duration of each packet
    record_onset(i)=packet_sent_time(1)-sniff_samples(1); %onset time for Sniff{trial(i)}
    trialframes = cell2mat(h5read(h5,strcat(Keys{i},'/frame_triggers')));
    if ~isempty(trialframes)
        frame_trigger = [frame_trigger; trialframes];
        [~,frame_fv(i)] = min(abs( frame_trigger -double(fvOnTime(i))));
    end
end
Blanked_Duration = 1e3;  % To identify each trial since voyeur start trial by FV opening.

%% Frame Trigger Correction
% This section of the code is dedicated to detecting and correcting missing frame triggers.

% Frame triggers can occasionally be missed due to various reasons. This code identifies such instances.

% In some cases, frame triggers are not missed, but are instead replaced by a later time (typically delayed by 300 - 400ms). 
% This results in the same frame trigger occurring twice.

% The code also handles these duplicates by removing them, ensuring each frame trigger is unique.

% Check for duplicate frame triggers
if length(frame_trigger)~=length(unique(frame_trigger))
    fprintf('duplicated frametrigger \n')
end

% Remove duplicates and sort the frame triggers
frametrigger2 = sort(unique(frame_trigger));
if any(frametrigger2 ==0)
    frametrigger2(frametrigger2==0) = [];
    fprintf('There is invalid frame trigger times \n')
end
if Blanked_recording
    Frame_endindices = [find(diff(frametrigger2)>Blanked_Duration); length(frametrigger2)];
    if length(wsFrameNumbers) >1
        VoyeurFrameNumbers = [Frame_endindices(1); diff(Frame_endindices)];
        fprintf(['VoyeurFrameTotal: ',num2str(sum(VoyeurFrameNumbers))]);
        fprintf(['WSFrameTotal: ',num2str(sum(wsFrameNumbers))]);
        [ProblematicTrials ,~,~]= find((VoyeurFrameNumbers(:)-wsFrameNumbers(:)) ~=0);
        if isempty(ProblematicTrials(:))
            fprintf('No frame loss in Voyeur \n');
        else
            fprintf(['Trials numbers :  ',num2str(ProblematicTrials'), '   have missing frames.','\n']);
            Frametrigger_temp = frametrigger2;
         for i = 1:length(ProblematicTrials)
          if ProblematicTrials(i)==1
                    ProblematicFrameBounds = [1 Frame_endindices(ProblematicTrials(i))];
          else
                ProblematicFrameBounds = [Frame_endindices(ProblematicTrials(i)-1)+1 Frame_endindices(ProblematicTrials(i))];
                frameidx=find(diff(Frametrigger_temp(ProblematicFrameBounds(1): ProblematicFrameBounds(2)))>frame_tol,1);
                 if ~isempty(frameidx)
                    frameidx =  frameidx + ProblematicFrameBounds(1) -1;
                 end
                fill_trial = 0;
                while ~isempty(frameidx)
                    fill=round((Frametrigger_temp(frameidx+1)-Frametrigger_temp(frameidx))/delta_t); %missing frames
                    errorvals = Frametrigger_temp(frameidx:frameidx+1);
                    Frametrigger_temp=[Frametrigger_temp(1:frameidx);...
                        round(Frametrigger_temp(frameidx)+delta_t:delta_t:Frametrigger_temp(frameidx)+delta_t*(fill-1))';...
                        Frametrigger_temp((frameidx+1):end)]; % make sure to add the right number of frames
                    fprintf(['We replaced ',num2str(errorvals'),' with ',...
                        num2str(Frametrigger_temp(frameidx:frameidx+fill)'),'\n']);
                    ProblematicFrameBounds  = [ProblematicFrameBounds(1) ProblematicFrameBounds(2)+fill-1];
                    frameidx=find(diff(Frametrigger_temp(ProblematicFrameBounds(1): ProblematicFrameBounds(2)))>frame_tol,1);
                    fill_trial = fill_trial + fill-1;
                    if ~isempty(frameidx)
                        frameidx =  frameidx + ProblematicFrameBounds(1) -1;
                    end
                    if length(errorvals') == length(Frametrigger_temp(frameidx:frameidx+fill)')
                        frameidx = [];
                    end
                end
                Frame_endindices = [find(diff(Frametrigger_temp)>1e3); length(Frametrigger_temp)];
                VoyeurFrameNumbers = [Frame_endindices(1); diff(Frame_endindices)];
                if VoyeurFrameNumbers(ProblematicTrials(i)) == wsFrameNumbers(ProblematicTrials(i))
                   fprintf(['Trials numbers :  ',num2str(ProblematicTrials(i)'), '  matched to wavesurfer.','\n']);
                end
            end
            frametrigger2 = Frametrigger_temp;
    end
        end
        %% check again especially for first trial since Voyeur can drop many frames at the beginning
        if(VoyeurFrameNumbers(1)-wsFrameNumbers(1)) ~=0
            Nmiss = wsFrameNumbers(1) -VoyeurFrameNumbers(1);
            frametrigger2 = vertcat(linspace(frametrigger2(1)-delta_t-Nmiss*delta_t, frametrigger2(1)-delta_t, Nmiss)', frametrigger2);
            Frame_endindices = [find(diff(frametrigger2)>Blanked_Duration); length(frametrigger2)];
            VoyeurFrameNumbers = [Frame_endindices(1); diff(Frame_endindices)];

        end

    else
        if ~sum( ~diff(frametrigger2)<1.2*1e3/fps &  ~diff(frametrigger2)>0.8*1e3/fps)
            fprintf('No problem in frames \n');
        end
    end
    
    %% Assign frame times to each trials by checking starttrial infos
    frame_trigger_trial = {};
    frames = frametrigger2;
    N_tiff = length(VoyeurFrameNumbers);
    Voyeurdata.Voyeur_frame_numbers = [];
    for i = 1: N_tiff
        Nf = VoyeurFrameNumbers(end-i+1);
        frames_i = frames(end-Nf+1:end);
        for j = 1:num_trial
            [val(j),frame_diff(j)] = min(abs(frames_i -double(Voyeurdata.starttrial(j))));
        end
        [~,trial_id] = min(val); 
        frame_trigger_trial{trial_id} = frames_i;
        frames(end-Nf+1:end) = [];
        Voyeurdata.Voyeur_frame_numbers(trial_id) = size(frames_i,1);
    end

else
    frame_trigger_trial = {};
    if nnz(diff(frametrigger2)>frame_tol)
        frameidx=find(diff(frametrigger2)>frame_tol,1);
        while ~isempty(frameidx)
            fill=round((frametrigger2(frameidx+1)-frametrigger2(frameidx))/delta_t); %missing frames
            errorvals = frametrigger2(frameidx:frameidx+1);
            frametrigger2=[frametrigger2(1:frameidx);...
                round(frametrigger2(frameidx)+delta_t:delta_t:frametrigger2(frameidx)+delta_t*(fill-1))';...
                frametrigger2((frameidx+1):end)]; % make sure to add the right number of frames
            fprintf(['We replaced ',num2str(errorvals'),' with ',...
                num2str(frametrigger2(frameidx:frameidx+fill)'),'\n']);
            %         frameidx=find(diff(frametrigger2)>35,1); Caused infinite loop in
            %         14393_190914_field4
            frameidx=find(diff(frametrigger2)>frame_tol,1);
        end
    end
end

frame_trigger = frametrigger2;

%%
%Inhalation timing within sniff traces of each trial
tot_pre_post = pre + post;
inh_onset_local=inh_onset-record_onset;
Sniff=zeros(num_trial,pre+post);
Sniff_time=zeros(num_trial,pre+post);

%check missing sniff packets
packet_sent_time=[];
sniff_samples=[];
trial_index = [];
trial_subind = [];
for i=1:num_trial
    %check missing sniff packet
    events = h5read(h5,strcat(Keys{i},'/Events'));
    packet_sent_time = [packet_sent_time;events.packet_sent_time];
    sniff_samples = [sniff_samples;events.sniff_samples];
    trial_index = [trial_index;i*ones(length(events.sniff_samples),1)];
    trial_subind = [trial_subind;[1:length(events.sniff_samples)]'];
end

lost_packet = find(diff(packet_sent_time)~=sniff_samples(2:end))+1;
sniff_all=[];
for i=1:num_trial
    st=zeros(1,pre+post);
    try
        sniffcell=h5read(h5,strcat(Keys{i},'/sniff'));
        if ~isempty(lost_packet)
            pos_lost=lost_packet(trial_index(lost_packet)==i);%position of lost packet in time from start
            
            %             subind=trial_subind(lost_packet);
            %             subind(trial_index(lost_packet)~=i)=[];
            if ~isempty(pos_lost)
                for k=1:length(pos_lost)
                    j = pos_lost(k);
                    fprintf('trial %d, %d-th packet lost',i,trial_subind(j))
                    sniffcell{trial_subind(j)}=[int16(zeros((packet_sent_time(j)-packet_sent_time(j-1))-length(sniffcell{trial_subind(j)}),1));sniffcell{trial_subind(j)}];
                end
                
            end
        end
        sniff0=cell2mat(sniffcell);
        sniff_all=[sniff_all;sniff0];
        %pad with 10000 zeros before and after sniff
        sniff=[zeros(2*tot_pre_post,1);sniff0;zeros(2*tot_pre_post,1)];
        sniff_pos=[(-2*tot_pre_post+1):0,1:length(sniff),length(sniff)+1:length(sniff)+(2*tot_pre_post)];
        
        sniff_range=inh_onset_local(i)-pre:inh_onset_local(i)+post-1;
        
        st=sniff(ismember(sniff_pos,sniff_range));
        
        Sniff(i,:)=st;
        Sniff_time(i,:)=inh_onset(i)-pre:inh_onset(i)+post-1;
    catch
        sprintf('Error in trial %d',i)
    end
    
end

% what is -29999  and 30000
sniffall_time=[(-5*tot_pre_post+1):0,1:length(sniff_all),length(sniff_all)+1:length(sniff_all)+(5*tot_pre_post)];
sniff_all=[int16(zeros((5*tot_pre_post),1));sniff_all;int16(zeros((5*tot_pre_post),1))];

%%
%find real inhalation onset based on the threshold crossing of sniff signal
num_trial=size(Sniff,1);
inh_onset=Voyeurdata.inh_onset;
fvOnTime=Voyeurdata.fvOnTime;

%Check sniff alignment to fvOnTime or inh_onset
Sniff(Sniff_time==0)=NaN;
Sniff_time(Sniff_time==0) = NaN; %% Sniff_time =0 is problematic acqusition and cause wrong estiamtion for inh and fv bin, MK 21/07/10
row_NaN = isnan(mean(Sniff_time,2)); % We will add this to inh and fv bin;
locs_NaN = find(row_NaN ==1);
[inh_bin,~]=find(Sniff_time'==double(inh_onset)',length(inh_onset),'first');
inh_bin = insert_element(inh_bin,locs_NaN, zeros(length(locs_NaN),1));
% inh_bin=findfirst(Sniff_time==double(inh_onset),2);
[bin_fv,trial_id]=find(Sniff_time'==double(fvOnTime)',length(fvOnTime),'first');
% fv_bin=findfirst(Sniff_time==double(fvOnTime),2);
fv_bin = 1950*ones(length(fvOnTime),1);
fv_bin(trial_id) = bin_fv;
fv_bin = insert_element(fv_bin,locs_NaN, zeros(length(locs_NaN),1));

%% Find Stim trigger
stim_time = inh_onset + Voyeurdata.pulseOnsetDelay_1;
for i = 1:num_trial-1
    %stim_frame(i)=find(frame_trigger<stim_time(i+1), 1, 'last');%frame in which inh_onset is included
    %if Voyeurdata.amplitude_1(i) >0
        [Voyeurdata.stim_diff(i), ind] = min(abs(double(frame_trigger)+fps/2.0 -double(stim_time(i+1)- pmtblank_dur/2.0)));
        Voyeurdata.stim_frame(i) = ind;
    %end
end
%% Need to change this it causes problem in inh detection
inh_inSniff=inh_onset-record_onset(1);
Snifftemp=zeros(size(Sniff));
for i=1:num_trial
    if nnz(ismember(sniffall_time,inh_inSniff(i)-pre:inh_inSniff(i)+post-1))>0
        
        Snifftemp(i,:)=sniff_all(ismember(sniffall_time,inh_inSniff(i)-pre:inh_inSniff(i)+post-1));
        Sniff_time(i,:)=inh_onset(i)-pre:inh_onset(i)+post-1;
    end
end
inh_bin2= inh_bin;
%% Realign inhalation traces using breathmetrics
if inh_detect
    for kk = 2:num_trial
        kk
        respiratoryTrace = Snifftemp(kk,:);
        bmObj = breathmetrics(-1*respiratoryTrace', 1e3, 'rodentAirflow');
        bmObj.estimateAllFeatures(1,'simple', 0, 1);
        [~, index] = min(abs(bmObj.inhaleOnsets-2001));
        inh_bin2(kk) = bmObj.inhaleOnsets(index);
        %fig = bmObj.plotFeatures({'onset'});
    end

    Voyeurdata.inh_onset_voyeur=Voyeurdata.inh_onset;
    inh_onset=Voyeurdata.inh_onset+int32(inh_bin2(:)-inh_bin(:));
    Voyeurdata.inh_onset=inh_onset;
    
else
    Voyeurdata.inh_onset_voyeur=Voyeurdata.inh_onset;
end
%%
inh_inSniff=inh_onset-record_onset(1);
Sniff2=zeros(size(Sniff));
for i=2:num_trial
    if nnz(ismember(sniffall_time,inh_inSniff(i)-pre:inh_inSniff(i)+post-1))>0
        
        Sniff2(i,:)=sniff_all(ismember(sniffall_time,inh_inSniff(i)-pre:inh_inSniff(i)+post-1));
        Sniff_time(i,:)=inh_onset(i)-pre:inh_onset(i)+post-1;
    end
end
Sniff=single(Sniff2);
Voyeurdata.frametrigger = frame_trigger;
Voyeurdata.fv_bin = fv_bin;
Sniff_smooth = Sniff;
for i=2:num_trial
    bmObj = breathmetrics(-1*Sniff(i,:)', 1e3, 'rodentAirflow');
    bmObj.estimateAllFeatures(1,'simple', 0, 1);
    Sniff_smooth(i,:) = -1*bmObj.baselineCorrectedRespiration;
    [~, index] = min(abs(bmObj.inhaleOnsets-2001));
    Voyeurdata.pre_inhs{i} = bmObj.inhaleOnsets(1:index-1);
    Voyeurdata.post_inhs{i} = bmObj.inhaleOnsets(index:end);
end
cd(path_orig);


function [v] = insert_element(vin,locs, locsval)
   v = vin;
   for i =1: length(locs)
       v = vertcat(v(1:locs(i)-1), locsval(i), v(locs(i):end));
   end    
