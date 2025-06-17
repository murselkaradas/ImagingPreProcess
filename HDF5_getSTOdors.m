function [odorInfo] = HDF5_getSTOdors(fpathH5,fnameH5,trials_read)

H5=h5read(fullfile(fpathH5,fnameH5),'/Trials');

Ntrial = length(H5.trialNumber);
%Check trials read exist or not
if ~exist('trials_read','var')
    trials_read = logical(ones(1,Ntrial));
end



% Define letters for each olfactory channel
letters = 'A':'H';  % A to H for 8 channels

% Number of trials (assuming H5 is structured with trials)
num_trials = length(H5.olfa_1_odor);  

% Initialize cell array to store the olfa_odor string for each trial
olfa_odor = strings(num_trials, 1);
% Loop through each trial
for t = 1:num_trials
    trial_str = "";  % Initialize empty string for this trial
    
    % Loop through all 8 olfactory channels
    for i = 1:8
        % Construct field names dynamically
        odor_field = sprintf('olfa_%d_odor', i);
        flow_field = sprintf('olfa_%d_flow', i);
        
        % Assign letter
        letter = letters(i);
        odor_tmp = (deblank(string(H5.(odor_field)')));
        % Read values for this trial
        odor = odor_tmp(t);
        flow = H5.(flow_field)(t);
        
        % Append to string if odor is not 'none'
        if odor ~= "None"
            trial_str = trial_str + letter + num2str(flow);
        end
    end
    
    % Store the resulting string for this trial
    olfa_odor(t) = trial_str;
end

olfa_odor(olfa_odor =='') = 'Blank';

olfa_odors = olfa_odor(trials_read);
odors = unique(olfa_odors);

% Separate "Blank" from other odors
is_blank = odors == "Blank";
blank_odor = odors(is_blank);
odors(is_blank) = [];

% Determine the max number of odors in any combination
max_odors = max(cellfun(@(x) sum(isletter(x)), odors));

% Initialize sorted array
sorted_odors = [];

% Loop through 1 to max_odors (single odors, double odors, ..., max combo)
for num_odor = 1:max_odors
    % Select odors with exactly 'num_odor' unique letters
    selected_odors = odors(cellfun(@(x) sum(isletter(x)) == num_odor, odors));
    
    % Sort alphabetically
    selected_odors = sort(selected_odors);
    
    % Append to sorted list
    sorted_odors = [sorted_odors; selected_odors];
end

% Append "Blank" at the end
sorted_odors = [sorted_odors; blank_odor];
for idx = 1:length(sorted_odors)
    odorTrials{idx} = find(olfa_odors==sorted_odors(idx));
end

odorInfo.odors = sorted_odors;
odorInfo.odorTrials = odorTrials;
end
