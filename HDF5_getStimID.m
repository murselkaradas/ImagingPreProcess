function [odorInfo] = HDF5_getStimID(fpathH5,fnameH5,varargin)
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

% Try to get PatternID
try 
    PatternID = deblank(string((H5.patternid)'));
catch
    warning('No patternid field in H5 file. Using stimtype');
    PatternID =  deblank(string((H5.stimtype)'));
end

% Try to get SLMID
try
    SLMID = deblank(string((H5.SLM_pattern)'));
    has_SLMID = true;
catch
    has_SLMID = false;
end

PatternID(PatternID =="None") = "";
[common1, pat_unique] = findUniquePatternParts(PatternID);
DMD_pat = strcat("P:",pat_unique);
DMD_pat(AOMpower==0) = "";
SLMID(SLMID =="None") = "";
[SLM_common, slm_unique] = findUniquePatternParts(SLMID);
OffDur = H5.pulseOffDur_1;
PulseDur(DMD_pat == "") = 0;
AOMpower(DMD_pat == "") = 0;
PulseLat(DMD_pat == "") = 0;
OffDur(DMD_pat == "") = 0;
DMD_pat(DMD_pat == "") = 'Blank';

cond1 = size(unique(PulseDur(DMD_pat ~= 'Blank')),1)>1;
cond2 = size(unique(OffDur(DMD_pat ~= 'Blank')),1)>1;
cond3 = size(unique(PulseLat(DMD_pat ~= 'Blank')),1)>1;
try
    cond4 = size(unique(AOMpower(DMD_pat ~= 'Blank')),1)>1;
catch
    cond4 = false;
end

ID = strings(Ntrial,1);
for i = 1:Ntrial
    if DMD_pat(i) == 'Blank'
        ID(i) = strcat(ID(i), 'DMD:', DMD_pat(i));
    else
        ID(i) = strcat(ID(i), 'DMD:', DMD_pat(i));
        if cond1
            ID(i) = strcat(ID(i), '/PD:', num2str(PulseDur(i)));
        end
        if cond2
            ID(i) = strcat(ID(i), '/Off:', num2str(OffDur(i)));
        end
        if cond3
            ID(i) = strcat(ID(i), '/L:', num2str(PulseLat(i)));
        end
        if cond4
            ID(i) = strcat(ID(i), '/P:', num2str(AOMpower(i), '%.1f'));
        end
    end
    if has_SLMID && ~isempty(SLM_common)
            ID(i) = strcat(ID(i), '/SLM:', slm_unique(i));
    end
end

IDs_ = ID(trials_read);
ID_unique = unique(IDs_);
for idx = 1:length(ID_unique)
    IDsTrials{idx} = find(IDs_ == ID_unique(idx));
end

odorInfo.odors = ID_unique;
odorInfo.odorTrials = IDsTrials;
end

function [commonPart, uniqueParts] = findUniquePatternParts(patterns)
    % Remove file extensions
    patterns = regexprep(patterns, '\..*$', '');
    uniqueParts = patterns;
    % Filter out empty or NaN patterns
    validPatterns = ~cellfun(@isempty, patterns) & ~cellfun(@(x) all(isnan(x)), patterns);
    patterns2 = patterns(validPatterns);
    
    % If no valid patterns are left, return empty results
    if isempty(patterns2)
        commonPart = '';
        uniqueParts = {};
        return;
    end
    
    % Find the longest common substring
    commonPart = '';
    minLength = min(cellfun(@length, patterns2));
    for i = 1:minLength
        allChars = cellfun(@(x) x(i), patterns2, 'UniformOutput', false);
        if all(strcmp(allChars{1}, allChars))
            commonPart = [commonPart allChars{1}];
        else
            break;
        end
    end
    
    % Extract unique parts
    uniqueParts(validPatterns) = cellfun(@(x) strrep(x, commonPart, ''), patterns2, 'UniformOutput', false);

end