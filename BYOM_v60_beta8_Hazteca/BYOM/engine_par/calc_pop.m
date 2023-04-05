function calc_pop(par_out,X0mat,t,c,Th,opt_conf,opt_pop)

% Usage: calc_pop(par_out,X0mat,t,c,Th,opt_conf,opt_pop)
%
% This function allows for a calculation of the intrinsic rate of increase
% using the Euler-Lotka equation. I removed discrete reproduction here as I
% think this has the potential to lead to nonsense (e.g., calibrating a
% model assuming continuous reproduction, and then calculating rgr with
% discrete repro events). Further, note that this function does not
% calculate negative rgr's; it stops at rgr=0, which is bad enough.
%
% This function needs some work, and its call may be changed in the future.
% It now makes use of calc_epx_helper to perform the rgr calculation.
%
% par_out : fitted parameter structure
% X0mat   : matrix with scenarios and initial values
% t       : time vector for population calculations
% c       : concentration range for population calculations
% Th      : hatching time
% opt_conf: structure with options for confidence intervals
% opt_pop : structure with options for population calculations
%
% Author     : Tjalling Jager
% Date       : August 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM.

global glo glo2

% Wrap the glo and glo2 structures into a single wrapper (this is an
% adaptation for the use of the parallel toolbox)
WRAP.glo  = glo;
WRAP.glo2 = glo2;

savenm1 = glo.basenm;

fscen       = opt_pop.fscen; % scenarios (food levels)
plt_fly     = opt_pop.plt_fly; % switch for plotting on the fly
use_par_out = opt_conf.use_par_out; % set to 1 to use par_out, as entered in this function, rather than from saved set

pop_txt{1} = glo.leglab1; % use the legend text for the x-axis
pop_txt{2} = glo.leglab2;

t       = t(:); % make sure time is a column vector
t2      = mean([t(1:end-1) t(2:end)],2); % new averaged time vector
plotmax = -1; % initialise variable to capture max growth rate

% create wrapper to allow use of calc_epx_helper
WRAP2.t  = t;
WRAP2.t2 = t2;
WRAP2.Th = Th; % hatching time

if c(1) == -1 % if user wants to use only the scenarios in X0mat ...
    cpop = X0mat(1,:);
else % a scenario range (e.g., a large concentration vector) is provided
    cpop = c;
end

rgr       = nan(length(cpop),length(fscen)); % initialise matrix to catch the pop. growth rate
[figh,ft] = make_fig(1,2); % make figure of correct size

subplot(1,2,1) % first subplot is for absolute rgr, plotted on the fly
h = gca; % handle to current axes
set(h,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
hold on
xlabel([pop_txt{1} pop_txt{2}],ft.name,ft.label)
ylabel('intrinsic rate of increase (1/d)',ft.name,ft.label)

plotcol = 'krbgcm'; % 6 colors should be enough ...

legtxt = cell(length(fscen),1); % pre-define cell array to collect legend entries
for i = 1:length(fscen) % run through scenarios (food levels)
    
    if i == 1
        rgr_init = 1; % initial value for RGR for fzero (should not matter too much)
    elseif isfinite(rgr(1,i-1))
        rgr_init = rgr(1,i-1); % for 2nd round and further, use first one from previous scenario
    end
    
    if length(fscen)>1
        % Here, different food levels are used, which means that parameter f is overwritten.
        par_out.f(1) = fscen(i); % take the next scenario, replace f in par_out
    end
    
    for j = 1:length(cpop) % run through all treatments
        if c(1) == -1
            X0mat_tmp    = X0mat(:,j); % take exact initial states for THIS treatment!
        else
            X0mat_tmp    = X0mat(:,1); % take first set of initial states for ALL treatments!
            X0mat_tmp(1) = c(j); % only replace the concentration
        end
        
        % continuous repro
        WRAP2.rgr_init = rgr_init; % put it in the wrapper as well
        [~,rgr_tmp] = calc_epx_helper(1,1,t,par_out,X0mat_tmp,glo,1,[],[],WRAP2);
        rgr(j,i)    = rgr_tmp;
        rgr_init    = rgr_tmp; % update initial guess (new one will not be far from old one)
        
        if plt_fly == 1 % plot results on the fly
            if c(1) == -1 % if user wants to use only the scenarios in X0mat ...
                plot(h,cpop(j),rgr(j,i),'ko')
            else
                plot(h,cpop(j),rgr(j,i),'k.')
            end
            drawnow
        end
    end
    
    if c(1) == -1 % if user wants to use only the scenarios in X0mat ...
        plot(h,cpop,rgr(:,i),[plotcol(i),'--'],'LineWidth',1);
    else
        hf(i) = plot(h,cpop,rgr(:,i),[plotcol(i),'-'],'LineWidth',1.5);
        legtxt{i} = ['f = ',num2str(fscen(i))];
    end
    plotmax = max(plotmax,max(rgr(:,i)));
    drawnow
    
end

ylim(h,[0 plotmax*1.05]) % restrict y axis
if c(1) ~= -1 && length(fscen)>1
    h_leg = legend(hf,legtxt); % create a legend if helpful
    set(h_leg,ft.name,ft.legend); % use standard font formats
end

% Use second panel for the relative population growth rate (relative to
% first treatment, which should be the 'control'). Should not be needed to
% remember the axes handle, and use it for plotting, since this is done
% once (so little risk of the user switching to another figure window).
figure(figh) % make correct figure current
subplot(1,2,2)
set(gca,'LineWidth',1,ft.name,ft.ticks) % adapt axis formatting
hold on
xlabel([pop_txt{1} pop_txt{2}],ft.name,ft.label)
ylabel('relative rate of increase (-)',ft.name,ft.label)
for i = 1:length(fscen)
    if c(1) == -1 % if user wants to use only the scenarios in X0mat ...
        plot(cpop,rgr(:,i)/rgr(1,i),[plotcol(i),'--'],'LineWidth',1);
    else
        hf(i) = plot(cpop,rgr(:,i)/rgr(1,i),[plotcol(i),'-'],'LineWidth',1.5);
    end
end
ylim([0 1.05]) % restrict y axis
drawnow
if c(1) ~= -1 && length(fscen)>1
    h_leg = legend(hf,legtxt); % create a legend if helpful
    set(h_leg,ft.name,ft.legend); % use standard font formats
end

%% Go for confidence intervals on RGR

type_conf = opt_conf.type; 

if type_conf > 0
    [rnd,par] = load_rnd(opt_conf); % loading and selecting the sample is handled in a separate function
    if numel(rnd) == 1 % that means that no sample was found
        type_conf = -1; % no need to produce an error, just do LC50s without CI
    end
    
    drawnow % empty plot buffer if there's something in it
    
    % Start/check parallel pool
    if glo2.n_cores > 0
        poolobj = gcp('nocreate'); % get info on current pool, but don't create one just yet
        if isempty(poolobj) % if there is no parallel pool ...
            parpool('local',glo2.n_cores) % create a local one with specified number of cores
        end
    end
    
    if use_par_out == 1
        par = par_out; % then we'll use the input par, rather than the saved one
        % Note: par_plot must already be structured in the main script, such
        % that the fitted parameters match the ones in the saved set, etc.
        % Note: when par_plot is NOT entered, it will have been made equal
        % to par from the saved set.
    end
    
    n_sets   = size(rnd,1); % number of MCMC samples
    pmat     = packunpack(1,par,0,WRAP); % transform structure *from saved set* into a regular matrix
    % it is better to use the saved par, as there may be differences in the
    % log-setting of parameters between the saved set and the optimised
    % par_out matrix (especially when using the alllog option in
    % calc_slice).
    
    par_comp(par,par_out,'f') % compare par from input with the one from the MAT file (ignore f!)
    ind_fit    = (pmat(:,2)==1); % indices to fitted parameters
    ind_logfit = (pmat(:,5)==0 & pmat(:,2)==1); % indices to pars on log scale that are also fitted!
    
    disp('Parallel toolbox is used; no progress is shown.')
    
    RGRlo_all = nan(length(cpop),length(fscen));
    RGRhi_all = nan(length(cpop),length(fscen));
    
    for i = 1:length(fscen) % run through scenarios (food levels)
        
        rgr_init = rgr(1,i); % initial value for RGR for fzero, taken from best parameter run
        RGR_coll = zeros(n_sets,length(cpop)); % initialise vector to collect RGR values
        
        % some trickery to get parfor running ... First create a pmat_coll with
        % all sets from the sample!
        pmat_coll = cell(n_sets,1); % cell array to construct pmat for each sample
        for k = 1:n_sets % run through all sets in the sample
            pmat(ind_fit,1) = rnd(k,:); % replace values in pmat with the k-th random sample from the MCMC
            % put parameters that need to be fitted on log scale back on normal
            % scale (as call_deri requires normal scale, in contrast to transfer.m)
            if sum(ind_logfit)>0
                pmat(ind_logfit,1) = 10.^(pmat(ind_logfit,1));
            end
            % Note: pmat is still on normal scale here, but the sample in rnd
            % contains the value on a log scale, if a parameter is fitted on log
            % scale. Call_deri needs normal scale structure par_k.
            
            % this has to be done after the tranformation to normal scale!
            pmat_coll{k} = pmat; % collect it in a huge cell array
        end
        
        Lcpop   = length(cpop); % to get parfor running ...
        glo_tmp = glo; % copy the global to get parfor running ...

        parfor k = 1:n_sets % run through all sets in the sample
            
            par_k = packunpack(2,0,pmat_coll{k},WRAP); % transform parameter matrix into a structure
            
            if length(fscen)>1
                % different food levels are used, which means that parameter f is overwritten.
                par_k.f(1) = fscen(i); % take the next scenario, replace f in par_out
            end
            
            WRAP3          = WRAP2; % this is needed to please parfor
            WRAP3.rgr_init = rgr_init; % put it in the wrapper as well
            
            for j = 1:Lcpop % run through all treatments
                if c(1) == -1
                    X0mat_tmp    = X0mat(:,j); % take exact initial states for THIS treatment!
                else
                    X0mat_tmp    = X0mat(:,1); % take first set of initial states for ALL treatments!
                    X0mat_tmp(1) = c(j); % only replace the concentration
                end
                
                % continuous repro
                [~,rgr_tmp]    = calc_epx_helper(1,1,t,par_k,X0mat_tmp,glo_tmp,1,[],[],WRAP3);
                RGR_coll(k,j)  = rgr_tmp;
                WRAP3.rgr_init = rgr_tmp; % update initial guess (new one will not be far from old one)
                
            end
            
        end
        
        % No need to exclude NaNs for prctile and min/max as they are
        % ignored (by default).
        switch type_conf % depending on the type of sample, take percentiles or min-max
            case 1
                RGRlo = prctile(RGR_coll,2.5,1);  % 2.5 percentile of LCx
                RGRhi = prctile(RGR_coll,97.5,1); % 97.5 percentile of LCx
            case 2
                RGRlo = min(RGR_coll,[],1); % take minimum
                RGRhi = max(RGR_coll,[],1); % take maximum
            case 3
                RGRlo = min(RGR_coll,[],1); % take minimum
                RGRhi = max(RGR_coll,[],1); % take maximum
        end
        RGRlo_all(:,i) = RGRlo(:); % make sure it is a column and collect it
        RGRhi_all(:,i) = RGRhi(:); % make sure it is a column and collect it
        plotmax = max(plotmax,max(RGRhi_all(:,i))); % update plotmax
        
        figure(figh) % make sure the correct figure is current to add the CIs in
        subplot(1,2,1)
        if c(1) == -1 % if user wants to use only the scenarios in X0mat ...
            plot(cpop,RGRlo_all(:,i),[plotcol(i),':'],'LineWidth',1);
            plot(cpop,RGRhi_all(:,i),[plotcol(i),':'],'LineWidth',1);
        else
            plot(cpop,RGRlo_all(:,i),[plotcol(i),':'],'LineWidth',1.5);
            plot(cpop,RGRhi_all(:,i),[plotcol(i),':'],'LineWidth',1.5);
        end
        ylim([0 plotmax*1.05]) % restrict y axis
        subplot(1,2,2)
        if c(1) == -1 % if user wants to use only the scenarios in X0mat ...
            plot(cpop,RGRlo_all(:,i)/RGRlo_all(1,i),[plotcol(i),':'],'LineWidth',1);
            plot(cpop,RGRhi_all(:,i)/RGRhi_all(1,i),[plotcol(i),':'],'LineWidth',1);
        else
            plot(cpop,RGRlo_all(:,i)/RGRlo_all(1,i),[plotcol(i),':'],'LineWidth',1.5);
            plot(cpop,RGRhi_all(:,i)/RGRhi_all(1,i),[plotcol(i),':'],'LineWidth',1.5);
        end
        drawnow
    end
end

% Add a title to the plot window
axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],'Box','off','Visible','off','Units','normalized', 'clipping','off');
switch type_conf
    case 1
        h_txt = text(0.5, 1,'Intervals: 95% of model predictions from posterior','HorizontalAlignment','center','VerticalAlignment', 'top');
    case 2
        h_txt = text(0.5, 1,'Intervals: predictions from likelihood region (likreg)','HorizontalAlignment','center','VerticalAlignment', 'top');
    case 3
        h_txt = text(0.5, 1,'Intervals: predictions from likelihood region (parspace)','HorizontalAlignment','center','VerticalAlignment', 'top');
    otherwise 
        h_txt = text(0.5, 1,['Plotted from: ',savenm1, ' (',date,')'],'Interpreter','none','HorizontalAlignment','center','VerticalAlignment', 'top');
end
set(h_txt,ft.name,ft.text); % use standard formatting for this header

if glo.saveplt > 0 % if we want to save the plot
    savenm = ['population_',savenm1];%
    save_plot(figh,savenm,h_txt);
end
