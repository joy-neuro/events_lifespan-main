% 2. sample subjects randomly, getting equal numbers of events 
clear all; close all; clc;

% add data path and load data
load test_data.mat;
load output/events.mat; % load event data

% iteration of kmeans
niter = 100;

% maximum iteration of kmeans for convergence
maxiter = 1000;

% variables
nag = 5; % number of age groups; nag = 7;
nrand = 2; % number of subjects randomly selected in each age group with replacement; nrand = 20;
nevnt = 10; % number of events for each subject
nsubj = length(subjts); % number of subjects

nnn = size(subjts{1}, 2); % number of nodes
utmask = triu(ones(nnn), 1) > 0; % upper triangle mask for matrix
mmm = nnn*(nnn-1)/2; % number of edges

folder = cd; % path to save
outputfolder = strcat(folder, '\output'); % directory of output
[u, v] = find(triu(ones(nnn), 1)); % edge indices

% set randomization seed
s1 = rng(1);

for k = 2
%for k = 2:20
    % keeping randomization indices
    randsid = zeros(niter, nag, nrand);
    randevnt = zeros(niter, nag, nrand, nevnt);

    for iter = 1:niter
        X = ['sampling iter...', num2str(iter)];
        disp(X);
        efc_p = nan(mmm, nrand*nevnt*nag);

        count = 1;
        for i = 1:nag
           % age maximum/minimum for subjects within each age bin 
           amin = 8*i - 3; % amin = 10*i - 5;
           amax = 8*i + 7; % amax = 10*i + 6;

           aid = subjage < amax & subjage > amin;
           sid = find(aid);

           rand_sid = randsample(length(sid), nrand);
           randsid(iter, i, :) = rand_sid;

           for j = 1:nrand
               ep = output.events{sid(rand_sid(j))}.edge_output.edge.tspeaks;

               % randomly select nevnt frames 
               sep = size(ep, 1);

               rand_sep = randsample(sep, nevnt);

               randevnt(iter, i, j, :) = rand_sep;
               %%

               for l = 1:nevnt    
                   ep_tmp = ep(rand_sep(l), :);
                   p = ep_tmp(:, u).*ep_tmp(:, v);
                   efc_p(:, count) = p;
                   ep_iter(:, count, iter) = ep_tmp;
                   count = count + 1;
               end
           end
        end
       [cip(:, iter), pC(:, :, iter)] = kmeans(efc_p', k, 'maxiter', maxiter, 'distance', 'correlation');
    end

   outputname = fullfile(outputfolder, sprintf('events_kmeans_%i.mat', k));
   save(outputname);
end