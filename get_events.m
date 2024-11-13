% 1. get all event frames per subject
clear all; close all; clc;

% load data
load test_data.mat;

% make directory
mkdir output;

% Note: these are randomly assigned ages for subjects in this test data;
subjage = randi([6 85], 1, 100);

nsubj = length(subjts); % number of subjects

folder = cd; % path to save
outputfolder = strcat(folder, '\output'); % directory of output

% get events for each subject
for i = 1:nsubj
    disp(i);
    
    try
        age = subjage(i);
        % Note: here we are negating the motion since the frames to retain 
        % were set as 1. Ideally this would be the framewise displacement 
        % values at each timepoint
        smt = subjmot{i};
        
        sts = subjts{i};

        output_tmp = fcn_get_events(sts, 100, smt, 0.15, 3, 0.05, 'circshift');

        mthrframes = sum(smt < 0.15);

        output.events{i}.sage = age;
        output.events{i}.edge_output = output_tmp;
        output.events{i}.mthrframes = mthrframes;
    catch
        smt = subjmot{i};
        mthrframes = sum(smt < 0.15);

        output.events{i}.sage = age;
        output.events{i}.edge_output = [];
        output.events{i}.mthrframes = mthrframes;
    end
end

outputname = fullfile(outputfolder, sprintf('events.mat'));
save(outputname, 'output', 'subjage');