%% kmeans centroid evnts
clear all; close all; clc;

load test_data.mat;
niter = 100; %100
kmax = 20;
maxiter = 1000;

nag = 5; % number of age groups; nag = 7;
nrand = 2; % number of randomly selected subjects per age group; nrand = 20;
nevnt = 10; % number of events per subject
nsubj = length(subjts); % number of subjects

nnn = size(subjts{1}, 2); % number of nodes
nutmask = triu(ones(nnn), 1) > 0; % nodal upper triangle mask
mmm = nnn*(nnn-1)/2;

folder = cd; % path to save
outputfolder = strcat(folder, '\output'); % directory of output
[u, v] = find(triu(ones(nnn), 1));

%% run kmeans on each subsample 
%for k = 2:kmax
for k = 2
    tic
    filename = ['\output\events_kmeans_' num2str(k) '.mat'];
    inputfile = strcat(folder, filename);
    load(inputfile);
    
    % calculate matchpair similarity of iterations
    psim = zeros(niter, niter);
    psim_gc_align = zeros(niter, 1);

    % calc matchpair similarity of iterations (min diag)
    psim_min = zeros(niter, niter);

    % matchpairs idx
    psim_mm = zeros(k, 2, niter, niter);

    for m = 1:niter
       X = ['similarity...', num2str(m)];
       disp(X);
       for n = m+1:niter
           cp1 = zeros(mmm, k); cp2 = cp1;
           
           pm = ep_iter(:, :, m);
           pn = ep_iter(:, :, n);
           
           % efc_p for calculating similarity across iterations
           pefc_m = pm(u, :).*pm(v, :);
           pefc_n = pn(u, :).*pn(v, :);
           
           for p = 1:k
               % similiarity across iterations for kmeans, etspeaks
               idx = cip(:, m) == p;
               cp1(:, p) = nanmean(pefc_m(:, idx), 2);
               idx = cip(:, n) == p;
               cp2(:, p) = nanmean(pefc_n(:, idx), 2);
           end
           
           %kmeans p
           rho = corr(cp1, cp2);
           mm = matchpairs(rho, 0, 'max');

           Rmincost = rho(mm(:, 1), mm(:, 2));
           cmask = logical(eye(k));
           psim_mm(:, :, m, n) = mm;

           psim(m, n) = min(Rmincost(cmask));
       end
    end
    
    % create symmetric matrix
    psim = psim + psim';

    % find best partition
    [max_psim, bc_psim] = max(mean(psim)); % best centroid!
    cicon_p = cip(:, bc_psim);

    % align all 100 iterations to that centroid and get grand centroid
    efc_psim = zeros(mmm, k, niter);
    a = psim_mm(:, :, :, bc_psim);
    for b = 1:bc_psim
        ep = ep_iter(:, :, b);
        efc_p = ep(u, :).*ep(v, :);   
        if isequal(a(:, 1, b), a(:, 2, b))
            for p = 1:k
                idx = cip(:, b) == p;
                efc_psim(:, p, b) = nanmean(efc_p(:, idx), 2);
            end
            psim_gc_align(b) = 1;
        else
            for p = 1:k
               idx = cip(:, b) == a(p, 1, b);
               efc_psim(:, p, b) = nanmean(efc_p(:, idx), 2);
            end
            psim_gc_align(b) = 0;
        end
    end
    
    a = psim_mm(:, :, bc_psim, :);
    for b = bc_psim+1:niter
        ep = ep_iter(:, :, b);
        efc_p = ep(u, :).*ep(v, :);
        if isequal(a(:, 1, 1, b), a(:, 2, 1, b))
            for p = 1:k
                idx = cip(:, b) == p;
                efc_psim(:, p, b) = nanmean(efc_p(:, idx), 2);
            end
            psim_gc_align(b) = 1;
        else
            for p = 1:k
               idx = cip(:, b) == a(p, 1, 1, b);
               efc_psim(:, p, b) = nanmean(efc_p(:, idx), 2);
            end
            psim_gc_align(b) = 0;
        end
    end
    
    % get grand centroid - psim
    psim_gc = zeros(size(efc_p, 1), k);
    for p = 1:k
       psim_gc(:, p) = nanmean(squeeze(efc_psim(:, p, :)), 2);
    end
    
    toc
    clear output psim_mm u v cp1 cp2; 

    % save run
    outputname = fullfile(outputfolder, sprintf(['events_kmeans_best_centroid_k_%i.mat'], k));
    
    save(outputname); 
end

% plot best aligned centroid
load hcp100.mat;
[gx,gy,jdx] = fcn_plot_blocks(lab);

for i = 1:k
    mat = zeros(nnn);
    mat(triu(ones(nnn), 1) > 0) = psim_gc(:, i);
    mat = mat + mat';
    subplot(1,k,i); imagesc(mat(jdx,jdx));%,[-1.5,1.5]);
    axis square; colormap(jet); 
    title(sprintf('centroid %i',i)); hold on; plot(gx,gy,'k');
end
