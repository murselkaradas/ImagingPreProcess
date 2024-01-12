function [Sniff,Sniff_smooth,data,Sniff_time]=ReadVoyeurInfo(h5,path_used,pre,post,inh_detect)
%This function is to extract sniff trace of individual trials
%Modified for use in 2p imaging
%Sniff(trial x ms): sniff traces of each trial
% Sniff_smooth: smoothed version of sniff
%Sniff_time(trial x ms): time in data file for corresponding time bins in
%Sniff
%
% Mursel Karadas 2022

%Duration before inahlation to be analyzed
if ~exist('pre','var')
    pre=1000;
end

%Duration after inhalation to be analyzed
if ~exist('post','var')
    post=1000;
end

if ~exist('inh_detect','var')
    %Whether or not detect inhalation by threshold crossing
    inh_detect=true;
end


path_orig=pwd;
cd(path_used)

try
    data=h5read(h5,'/Trials');
catch
    error('%s not found in %s',h5, path_used)
end
%Need to check if this manipulation is ok for 2p data
inh_onset=double(data.inh_onset);
num_trial=length(inh_onset);

h_info=h5info(h5);
h_info=h_info.Groups;
% Keys is Number_of_trials x 1 
Keys={h_info.Name};%h5 file from 2p rig doesn't record trial0 analog signal
fvOnTime=data.fvOnTime;
%%
%obtain frametrigger
record_onset=zeros(size(inh_onset));
for i=1:num_trial
    %i
    data_Events=h5read(h5,strcat(Keys{i},'/Events'));
    packet_sent_time=data_Events.packet_sent_time;
    sniff_samples=data_Events.sniff_samples;    %array for the duration of each packet
    record_onset(i)=packet_sent_time(1)-sniff_samples(1); %onset time for Sniff{trial(i)}
end

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
inh_onset=data.inh_onset;
fvOnTime=data.fvOnTime;

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
%%
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

    data.inh_onset_voyeur=data.inh_onset;
    inh_onset=data.inh_onset+int32(inh_bin2(:)-inh_bin(:));
    data.inh_onset=inh_onset;
    
else
    data.inh_onset_voyeur=data.inh_onset;
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
data.fv_bin = fv_bin;
Sniff_smooth = Sniff;
for i=2:num_trial
    bmObj = breathmetrics(-1*Sniff(i,:)', 1e3, 'rodentAirflow');
    bmObj.estimateAllFeatures(1,'simple', 0, 1);
    Sniff_smooth(i,:) = -1*bmObj.baselineCorrectedRespiration;
    [~, index] = min(abs(bmObj.inhaleOnsets-2001));
    data.pre_inhs{i} = bmObj.inhaleOnsets(1:index-1);
    data.post_inhs{i} = bmObj.inhaleOnsets(index:end);
end
cd(path_orig);


function [v] = insert_element(vin,locs, locsval)
   v = vin;
   for i =1: length(locs)
       v = vertcat(v(1:locs(i)-1), locsval(i), v(locs(i):end));
   end    
