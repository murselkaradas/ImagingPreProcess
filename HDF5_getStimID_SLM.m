function [odorInfo] = HDF5_getStimID_SLM(fpathH5,fnameH5,varargin)

p = inputParser;
addParameter(p, 'latfactor', 10, @isnumeric);
addParameter(p, 'sortflow', false, @islogical);
addParameter(p, 'trials_read', [],@islogical);

parse(p,varargin{:});
latfactor = p.Results.latfactor;
sortflow = p.Results.sortflow;
trials_read = p.Results.trials_read;
p.Results

H5=h5read(fullfile(fpathH5,fnameH5),'/Trials');

Ntrial = length(H5.trialNumber);
%Check trials read exist or not
if isempty(trials_read)
    trials_read = logical(ones(1,Ntrial));
end
PulseDur = H5.pulseOnDur_1/1e3;
try
    AOMpower = H5.AOMpower/1000;
catch
    warning('No AOM field in H5 file. Using Laser2 Amplitude_2');
    AOMpower = H5.amplitude_2/1000;
end
PulseLat = H5.pulseOnsetDelay_1;
try 
    PatternID = deblank(string((H5.SLM_pattern)'));
catch
    warning('No patternid field in H5 file. Using stimtype');
    PatternID =  deblank(string((H5.stimtype)'));
end

OffDur = H5.pulseOffDur_1;

% Combined condition to check for empty PatternID  
condition_check = (PatternID == "") ;

PulseDur(condition_check) = 0;
AOMpower(condition_check) = 0;
PulseLat(condition_check) = 0;
OffDur(condition_check) = 0;
PatternID(condition_check) = 'Blank';

cond1 = size(unique(PulseDur(PatternID ~= 'Blank')),1)>1;
cond2 = size(unique(OffDur(PatternID ~= 'Blank')),1)>1;
cond3 = size(unique(PulseLat(PatternID ~= 'Blank')),1)>1;
try
    cond4 = size(unique(AOMpower(PatternID ~= 'Blank')),1)>1;
catch
    cond4 = false;
end
cond4 = size(unique(AOMpower(PatternID ~= 'Blank')),1)>1;
ID  = strings(Ntrial,1);
for i = 1:Ntrial
    if PatternID(i) =='Blank'
        ID(i) = strcat(ID(i) ,'id:',PatternID(i));
    else
        ID(i) = strcat(ID(i),'id:',PatternID(i));
        if cond1
            ID(i) = strcat(ID(i),'/PD:',num2str(PulseDur(i)));
        end
        if cond2
            ID(i) = strcat(ID(i),'/Off:',num2str(OffDur(i)));
        end
        if cond3
            ID(i) = strcat(ID(i),'/L:',num2str(PulseLat(i)));
        end
        if cond4
            ID(i) = strcat(ID(i),'/P:',num2str(AOMpower(i),'%.1f'));
        end
    end
end
%ID = strcat('ID:',PatternID,'/ONDur:',num2str(PulseDur),'/OFFDur:',num2str(OffDur),'/Lat:',num2str(PulseLat),'/Power:',num2str(AOMpower,'%.1f'));
IDs_ = ID(trials_read);
ID_unique = unique(IDs_);
for idx = 1:length(ID_unique)
    IDsTrials{idx} = find(IDs_==ID_unique(idx));
end

odorInfo.odors = ID_unique;
odorInfo.odorTrials = IDsTrials;
end