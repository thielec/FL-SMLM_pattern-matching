%% Example code for pattern matching
%
% The code loads two files and uses the single molecule TCSPC curves as
% reference data. Based on the reference data, the single molecules in multiple
% files are classified.
% 
% The loaded files (*_TNT.mat) are generated with TrackNTrace using the lifetime
% fitting plugin with the option to export the TCSPC data.
%
% In this example, confocal dSTORM measurments of cells labeled with either
% Alexa 647 or Atto 655 are used. The raw data can be found at:
% https://projects.gwdg.de/projects/cfl-smlm/repository
%
% Christoph Thiele, 2021

%% Load the analysed data from the Tubulin-Alexa647 and Clathrin-Atto655 sample
% Folder with the mat files
resfolder = '.\tntres\';

% Files for reference patterns
tntres_ref  = [...
    load([resfolder 'COS7_beta-tubulin-Alexa647_D2O+MEA_power6+OD0p6_1_2021-m09-d03-14h49_TNT.mat']),...
    load([resfolder 'COS7_clathrin-Atto655_D2O+MEA_power6+OD0p6_3_2021-m09-d03-14h49_TNT.mat'])...
    ];

% Files to classify
testfiles = {'COS7_beta-tubulin-Alexa647_D2O+MEA_power6+OD0p6_2_2021-m09-d03-14h49_TNT.mat','COS7_beta-tubulin-Alexa647_D2O+MEA_power6+OD0p6_3_2021-m09-d03-14h49_TNT.mat','COS7_beta-tubulin-Alexa647_D2O+MEA_power6+OD0p6_4_2021-m09-d03-14h49_TNT.mat';...
             'COS7_clathrin-Atto655_D2O+MEA_power6+OD0p6_1_2021-m09-d03-14h49_TNT.mat','COS7_clathrin-Atto655_D2O+MEA_power6+OD0p6_2_2021-m09-d03-14h49_TNT.mat','COS7_clathrin-Atto655_D2O+MEA_power6+OD0p6_4_2021-m09-d03-14h49_TNT.mat'};

% Ground truth (1 == Alexa647, 2 == Atto 655)
testind = 1*contains(testfiles,'647') + 2*contains(testfiles,'655');

%% Define filter

% Define filter to remove bad localisations (background, multiple emitters),
% same syntax as in TrackNTrace
filter_str = 'Frame>500 & abs(sigma-1.4)<.4 & nphoton>50';

% Generates filter function
getFilterFun = @(paramsNames)str2func(['@(posData)true(size(posData,1),1)&' regexprep(vectorize(['( ' filter_str ' )']),...
    strcat('(?<!\w)',matlab.lang.makeValidName(paramsNames),'(?!\w)'),... % Replace spaces with underscore and makes sure the match is not within a name.
    cellfun(@(n)sprintf('posData(:,%i)',n),num2cell(1:numel(paramsNames)),'UniformOutput',false),...
    'ignorecase')]);


%% Generate reference decays
rnum = numel(tntres_ref); % number of references

% Variabel for the reference patterns
tcspc_ref = nan(max(arrayfun(@(tntres)size(tntres.postprocOptions.TCSPC,2),tntres_ref)),rnum); % this only works well if all TCSPC curves have the same length (+/-1) and and peak position

for ridx = 1:rnum  
    % For each reference: find all molecules that pass the filter and sum their
    % TCSPC histograms
    assert(isfield(tntres_ref(ridx).postprocOptions,'TCSPC')); % Make sure TCSPC data is exported
    
    % find index of Track-ID (should usully be 1)
    data_track_ind = find(contains(tntres_ref(ridx).postprocOptions.outParamDescription,'Track-ID','IgnoreCase',true),1);
    
    % get Track-IDs for molecules passing the filter
    filter_ind = feval(getFilterFun(tntres_ref(ridx).postprocOptions.outParamDescription),tntres_ref(ridx).postprocData);
    filter_TrackIDs = unique(tntres_ref(ridx).postprocData(filter_ind,data_track_ind));
    % save summed TCSPC histograms of selected molecules
    tcspc_ref(1:size(tntres_ref(ridx).postprocOptions.TCSPC,2),ridx) = sum(tntres_ref(ridx).postprocOptions.TCSPC(filter_TrackIDs,:),1);
end
    
%% Prepare patterns liklihood measure Q for each TCSPC
% Avoid zeros in (too) sparse reference data
tcspc_ref(tcspc_ref==0|isnan(tcspc_ref)) = 0.01;
% calculate log of normalised reference pattern
% (Enderlein&Sauer 2001, eq 3)
ln_p_ref = log(tcspc_ref./sum(tcspc_ref,1));


%% Classify test data
tnum = numel(testfiles);
tdata = {};
for tidx = 1:tnum
    tntres = load([resfolder testfiles{tidx}]);
    
    data_track_ind = find(contains(tntres.postprocOptions.outParamDescription,'Track-ID','IgnoreCase',true),1);
    assert(~isempty(data_track_ind));
    data_tau_ind = find(contains(tntres.postprocOptions.outParamDescription,'lt-tau','IgnoreCase',true),1);
    assert(~isempty(data_tau_ind));
    
    % apply filter
    filter_ind = find(feval(getFilterFun(tntres.postprocOptions.outParamDescription),tntres.postprocData));
    
    % save index and lifetime of selected molecules
    [~,idx] = unique(tntres.postprocData(filter_ind,data_track_ind));
    data = struct();
    data.Track_ID = tntres.postprocData(filter_ind(idx),data_track_ind);
    data.lt_tau = tntres.postprocData(filter_ind(idx),data_tau_ind);
    
    assert(isfield(tntres.postprocOptions,'TCSPC')); % Make sure TCSPC data is exported
    % Calculate the likelihood matrix lambda
    Lambda = - tntres.postprocOptions.TCSPC(data.Track_ID,1:size(ln_p_ref,1)) * ln_p_ref;
    % Find the position of the minima
    [~,PM_ind] = min(Lambda,[],2);
    
    % Calculate the posterior proberbility for the nth species
    % f_n = exp(-Lambda)./(sum(exp(-Lambda),2)); % P(S1)/(P(S1)+P(S2))
    f_n = exp(-Lambda-mean(-Lambda,2))./(sum(exp(-Lambda-mean(-Lambda,2)),2)); % Avoids nans due to overflowing floats
    
    data.PM_ind = PM_ind;
    for qidx = 1:size(f_n,2)
        data.(sprintf('f%i',qidx)) = f_n(:,qidx);
    end
    
    % save the ground truth
    data.true_ind = ones(size(PM_ind)) * testind(tidx);
    %%
    tdata{tidx} = struct2table(data);
end

tdata = vertcat(tdata{:});


%% Plot the results
% make figure if it does not exist
if ~exist('fig','var') || isempty(fig) || ~isgraphics(fig)
    %%
    fig = figure('Color',[1 1 1],'Units','centimeters');
    fig.Position(3:4) = [18,9];
    movegui(fig,'onscreen');
else
    clf(fig);
end
% make subplots
ax_patterns = subplot(1,3,1,'Parent',fig);
ax_roc_PM = subplot(1,3,2,'Parent',fig);
ax_roc_LT = subplot(1,3,3,'Parent',fig);

% plot patterns
semilogy(ax_patterns,tcspc_ref./sum(tcspc_ref,1));
ax_patterns.XLabel.String = 'time bin';
ax_patterns.YLabel.String = 'probability';
lgd = legend(ax_patterns,{'pattern Alexa 647','pattern Atto 655'});
lgd.ItemTokenSize = [9 18];
%% Plot ROC curve for PM

[X1,Y1,T1,AUC] = perfcurve(tdata.true_ind,tdata.f1,1);
plot(ax_roc_PM,X1,Y1,'DisplayName',sprintf('Alexa 647: AUC %.2f',AUC),'LineWidth',1);
ax_roc_PM.XLabel.String = 'false positive rate';
ax_roc_PM.YLabel.String = 'true positive rate';

ax_roc_PM.NextPlot = 'add';
[X2,Y2,T2,AUC] = perfcurve(tdata.true_ind,tdata.f2,2);
plot(ax_roc_PM,X2,Y2,'DisplayName',sprintf('Atto 655: AUC %.2f',AUC),'LineWidth',1);
ax_roc_PM.NextPlot = 'replace';
ax_roc_PM.XTick = 0:0.2:1;
ax_roc_PM.NextPlot = 'replace';

lgd = legend(ax_roc_PM,'Location','southeast');
lgd.Title.String = 'pattern matching';
lgd.ItemTokenSize = [9 18];
%% ROC lifetime

[X,Y,T,AUC] = perfcurve(tdata.true_ind(~isnan(tdata.lt_tau)),-tdata.lt_tau(~isnan(tdata.lt_tau)),1);
plot(ax_roc_LT,X,Y,'DisplayName',sprintf('Alexa 647: AUC %.2f',AUC),'LineWidth',1);
ax_roc_LT.XLabel.String = 'false positive rate';
ax_roc_LT.YLabel.String = 'true positive rate';

ax_roc_LT.NextPlot = 'add';
[X,Y,T,AUC] = perfcurve(tdata.true_ind(~isnan(tdata.lt_tau)),tdata.lt_tau(~isnan(tdata.lt_tau)),2);
plot(ax_roc_LT,X,Y,'DisplayName',sprintf('Atto 655: AUC %.2f',AUC),'LineWidth',1);
ax_roc_LT.NextPlot = 'replace';
ax_roc_LT.XTick = 0:0.2:1;

lgd = legend(ax_roc_LT,'Location','southeast');
lgd.Title.String = 'lifetime classification';
lgd.ItemTokenSize = [9 18];