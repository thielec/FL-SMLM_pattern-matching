%% Example code for an MLE grid search based on pattern matching
%
% The code loads the single molecule TCSPC histograms from a TrackNTrace results
% file and fits the lifetime with an MLE grid search.
% Alternatively, simulated decays can be fitted. 
% 
% The loaded file (*_TNT.mat) is generated with TrackNTrace using the lifetime
% fitting plugin with the option to export the TCSPC data.
%
% In this example, a confocal dSTORM measurment of a cell labeled with Alexa 647
% is used. The raw data can be found at:
% https://projects.gwdg.de/projects/cfl-smlm/repository
%
% Christoph Thiele, 2021


%% Get experimental single molecule TCSPC curves
tntres = load('tntres\COS7_beta-tubulin-Alexa647_D2O+MEA_power6+OD0p6_1_2021-m09-d03-14h49_TNT.mat');

% Define filter to remove bad localisations (background, multiple emitters),
% same syntax as in TrackNTrace
filter_str = 'Frame>500 & abs(sigma-1.4)<.4 & nphoton>50';

% Generates filter function
getFilterFun = @(paramsNames)str2func(['@(posData)true(size(posData,1),1)&' regexprep(vectorize(['( ' filter_str ' )']),...
    strcat('(?<!\w)',matlab.lang.makeValidName(paramsNames),'(?!\w)'),... % Replace spaces with underscore and makes sure the match is not within a name.
    cellfun(@(n)sprintf('posData(:,%i)',n),num2cell(1:numel(paramsNames)),'UniformOutput',false),...
    'ignorecase')]);

% find index of Track-ID (should usully be 1)
data_track_ind = find(contains(tntres.postprocOptions.outParamDescription,'Track-ID','IgnoreCase',true),1);

% get Track-IDs for molecules passing the filter
filter_ind = feval(getFilterFun(tntres.postprocOptions.outParamDescription),tntres.postprocData);
[filter_TrackIDs,uind] = unique(tntres.postprocData(filter_ind,data_track_ind));

% save fitted single molecule lifetimes
data_tau_ind = find(contains(tntres.postprocOptions.outParamDescription,'lt-tau','IgnoreCase',true),1);
filter_ind = find(filter_ind);
fullMLE_tau = tntres.postprocData(filter_ind(uind),data_tau_ind);

% save single molecule TCSPCs
tcspc_sm = tntres.postprocOptions.TCSPC(filter_TrackIDs,:);
tcspc_dt = mean(diff(tntres.postprocOptions.TCSPC_t)); % TCSPC time resolution in ns

% find tail and get time axis    
tailfit_cutoff = 0.2;
% find peak position    
[~,maxpos] = max(sum(tcspc_sm,1));
% index of first time bin in the tail
tailstart = maxpos + ceil(tailfit_cutoff/tcspc_dt);
% make time axis
tcspc_t = tcspc_dt*(0:(size(tcspc_sm,2)-tailstart));
% Duration of the decay curve
tcspc_T = diff(tcspc_t([1 end]));
% Normalised monoexponetial decay with background
pfun_monoexp = @(tau) exp(-tcspc_t(:)./tau); 
pfun_monoexpBG = @(tau,b)b./numel(tcspc_t(:))+(1-b).*pfun_monoexp(tau)./sum(pfun_monoexp(tau),1);

%% ALTERNATIVE Simulate TCSPC decays
% n_decays = 1e4;  % number of molecules
% n_photons = 1e3; % average number of photons
% fullMLE_tau = 2+0.5*randn(n_decays,1); % normal distributed lifetimes 2+/-0.2 ns
% b_decay = max(0,0.1+0.02*randn(n_decays,1)); % background of the decay
% 
% tcspc_dt = 0.016;% TCSPC time resolution
% tcspc_T = 25; % Duration of the decay curve
% tcspc_t = 0:tcspc_dt:tcspc_T;
% 
% % Normalised monoexponetial decay with background
% pfun_monoexp = @(tau) exp(-tcspc_t(:)./tau); 
% pfun_monoexpBG = @(tau,b)b./numel(tcspc_t(:))+(1-b).*pfun_monoexp(tau)./sum(pfun_monoexp(tau),1);
%
% % Generate TCSPC curves with Poissonian noise
% tcspc_sm = poissrnd(n_photons*pfun_monoexpBG(fullMLE_tau(:)',b_decay(:)')');
% tailstart = 1;
%% Calculate patterns

% define grid
grid_lts = linspace(0.01,5,500); % lifetime values
grid_bs = linspace(0.0,0.6,60);  % background values

% get all combinations
[grid_ltmat,grid_bmat] = meshgrid(grid_lts,grid_bs);
grid_ltmat = grid_ltmat(:)';
grid_bmat = grid_bmat(:)';

grid_logdecays = log(pfun_monoexpBG(grid_ltmat,grid_bmat));

%% Find best matching pattern
% Calculate Lambda and directly find the index of the minimum value
[~,grid_ind] = min(-tcspc_sm(:,tailstart:end) * grid_logdecays,[],2);

% get the correponding lifetime and background values
gridMLE_tau = grid_ltmat(grid_ind);
gridMLE_b = grid_bmat(grid_ind);
gridMLE_N = sum(tcspc_sm(:,tailstart:end),2);

%% Compare results
ax = cla;
histogram2(ax,fullMLE_tau(:),gridMLE_tau(:),grid_lts,grid_lts,'DisplayStyle','tile');
ax.DataAspectRatio = [1 1 1];
ax.XLabel.String = 'lifetime full MLE fit (ns)';
ax.YLabel.String = 'lifetime MLE grid search (ns)';
line(ax,ax.XLim,ax.XLim,'Color',[0 0 0]);