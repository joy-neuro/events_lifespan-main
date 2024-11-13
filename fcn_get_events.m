function output = fcn_get_events(ts,numrand,motion,motion_thr,dilateframes,q,nullmodel,eventtype)
% extracts events from time series
%
%   output = fcn_get_events(ts,numrand,motion,motion_thr,dilateframes,q);
%
%   inputs:
%
%           ts  -  [time x node] regional time series
%      numrand  -  number of randomizations to perform using permutation
%                  null model (default).
%       motion  -  [time x 1] vector of motion characteristics.
%    motionthr  -  threshold for considering a frame "high motion"
% dilateframes  -  number of frames to remove before/after high motion
%                  frames.
%            q  -  false discover rate used as part of statistical testing.
%    nullmodel  -  string for determining which null model to to use. if
%                  empty defaults to 'permutation'. other option is
%                  'circshift'.
%    eventtype  -  default is to calculate events based on both edge and
%                  node  RSS amplitudes. if you only want nodes, you can
%                  specify 'node'.

% Rick Betzel, IU 2021
%   updated 3/11/21 - added optional input for detecting events based on
%                     node activity rather than edge co-activity.


% zscore regional time series
ts = zscore(ts);

if ~exist('nullmodel','var') || isempty(nullmodel)
    nullmodel = 'permutation';
end

if ~exist('eventtype','var') || isempty(nullmodel)
    eventtype = 'edge';
end

% number of time points, nodes
[t,n] = size(ts);

% calculate rss

if ~strcmp(eventtype,'node')
    
    % calculate edge time series
    ets = fcn_edgets(ts);
    
    % get rss
    rss.edge = sum(ets.^2,2).^0.5;
    
end

% get rss
rss.node = sum(ts.^2,2).^0.5;

% find high-motion frames
idx = find(motion >= motion_thr);

% relmot is a logical array 1 = good frame, 0 = bad. this loop dilates
% high-motion frames by zeroing out frames before and after that might also
% be contaminated by motion.
relmot = ones(size(ts,1),1);
for i = 1:length(idx)
    jdx = (idx(i) - dilateframes):(idx(i) + dilateframes);
    jdx(jdx < 1) = [];
    jdx(jdx > t) = [];
    relmot(jdx) = 0;
end

% randomly permute order of time series (can replace with circular shift if
% you want).

for irand = 1:numrand
    
    % switch null model
    switch nullmodel
        case 'circshift'
            tsr = fcn_circshift(ts);
        case 'permutation'
            tsr = fcn_tsperm(ts);
    end
    
    if ~strcmp(eventtype,'node')
        
        % calculate edge time series
        etsr = fcn_edgets(tsr);
        
        % get rss
        rredge = sum(etsr.^2,2).^0.5;
        
        % store rss
        rssrand.edge(:,irand) = rredge;
        
    end
    
    % get rss
    rrnode = sum(tsr.^2,2).^0.5;
    
    % store rss
    rssrand.node(:,irand) = rrnode;
    
end
%%

% repeat statistical analysis for edge and node amplitudes
names = fieldnames(rss);
for i = 1:length(names)
    
    % calculate non-parametric p-values at each frame using pooled, global
    % distribution of rss values
    r = rss.(names{i});
    rrand = rssrand.(names{i});
    
    % calculate non-parametric p-values at each frame using pooled, global
    % distribution of rss values
    p = zeros(length(r),2);
    rrand_vec = rrand(:);
    for t = 1:length(r)
        p(t,1) = mean(rrand_vec >= r(t));
        p(t,2) = mean(rrand_vec <= r(t));
    end
    
    % calculate adjusted p-value
    padj = fcn_linear_step_up(p(:),q);
    
    % find frames that pass
    sig = p < padj;
    
    % mask those frames
    mask = bsxfun(@times,sig,[1,-1]);
    
    % loop over significant frames with rss > and < than expected and calculate
    % some statistics
    peak_mask = zeros(size(mask));
    labels = zeros(size(mask));
    count = 0;
    for j = 1:2
        vals = [];
        
        % group significant frames into contiguous blocks
        jdx = find(sig(:,j));
        dff = (jdx') - (1:length(jdx));
        unq = unique(dff);
        
        % for each block, find the frame with extreme rss (max if >, min if <)
        for k = 1:length(unq)
            kdx = jdx(dff == unq(k));
            count = count + 1;
            if j == 1
                rk = r(kdx);
                [~,pkindx] = max(rk);
            elseif j == 2
                rk = -r(kdx);
                [~,pkindx] = max(rk);
            end
            peak_mask(kdx(pkindx),j) = count;
            labels(kdx,j) = count;
        end
    end
    
    % get z-scored activity from peaks and troughs
    tspeaks = ts(peak_mask(:,1) > 0,:);
    tstroughs = ts(peak_mask(:,2) > 0,:);
    
    I = dummyvar(labels(:,1) + 1);
    I = bsxfun(@rdivide,I,sum(I));
    I = I(:,2:end);
    tsmeanpeaks = I'*ts;
    
    I = dummyvar(labels(:,2) + 1);
    I = bsxfun(@rdivide,I,sum(I));
    I = I(:,2:end);
    tsmeantroughs = I'*ts;
    
    % get motion scores for those frames
    motionpeaks = relmot(peak_mask(:,1) > 0);
    motiontroughs = relmot(peak_mask(:,2) > 0);
    
    %%
    
    % arrange everything into a structure
    output.(names{i}).r = r;
    output.(names{i}).rrand = rrand;
    output.(names{i}).tspeaks = tspeaks;
    output.(names{i}).tstroughs = tstroughs;
    output.(names{i}).motion = motion;
    output.(names{i}).sig = sig;
    output.(names{i}).peak_mask = peak_mask;
    output.(names{i}).padj = padj;
    output.(names{i}).p = p;
    output.(names{i}).motionpeaks = motionpeaks;
    output.(names{i}).motiontroughs = motiontroughs;
    output.(names{i}).labels = labels;
    output.(names{i}).tsmeanpeaks = tsmeanpeaks;
    output.(names{i}).tsmeantroughs = tsmeantroughs;
    
end
output.ts = ts;
output.relmot = relmot;
output.q = q;
output.dilateframes = dilateframes;
output.nullmodel = nullmodel;
output.motion_thr = motion_thr;

%% extra functions

function [a,u,v] = fcn_edgets(ts,zflag)
if nargin == 1
    zflag = true;
end
[~,n] = size(ts);               % number samples/nodes
if zflag
    z = zscore(ts);                 % z-score
else
    z = ts;
end
[u,v] = find(triu(ones(n),1));  % get edges
a = z(:,u).*z(:,v);             % edge ts products

return

function [padj] = fcn_linear_step_up(p,q)
%FCN_LINEAR_STEP_UP         computes adjusted p-values
%
%   PADJ = FCN_LINEAR_STEP_SUP(P,Q) takes a vector of values, P, and an
%          accepted false discovery rate, Q, and computes an adjusted
%          p-value, PADJ.
%
%   INPUTS:     P,      vector of p values
%               Q,      false discovery rate
%
%   OUTPUTS: PADJ,      adjusted p value
%
%   Richard Betzel, Indiana University, 2012
%

%modification history
%02.18.2012 - original

p(isnan(p)) = [];
[r,c] = size(p);
if r > c
    p = p';
end

n = length(p);
t = (1:n)*(q/n);
p = sort(p,'ascend');

ind = find((t - p) > 0,1,'last');
padj = p(ind);

if isempty(padj)
    padj = 0;
end

return

function tsr = fcn_circshift(ts,offsets)
if nargin == 1
    offsets = -size(ts,1):size(ts,1);
end
tsr = zeros(size(ts));
for i = 1:size(ts,2)
    tsr(:,i) = circshift(ts(:,i),offsets(randi(length(offsets))));
end
return

function tsr = fcn_tsperm(ts)
[~,r] = sort(rand(size(ts)));
r = bsxfun(@plus,r,(0:(size(ts,2) - 1))*size(ts,1));
tsr = ts(r);
return