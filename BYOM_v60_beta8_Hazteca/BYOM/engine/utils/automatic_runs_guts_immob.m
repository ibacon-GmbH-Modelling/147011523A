function [par_out,best_sel] = automatic_runs_guts_immob(fit_tox,par,skip_sg,SEL,CONFIG,opt_optim,opt_plot,varargin)

% This is a helpful function to run a GUTS-immobility analysis in parts,
% based on the setting in fit_tox. Furthermore, it allows to run several
% death mechanisms sequentially (SD and IT, and even mixed), providing a
% summary comparison of the AICs. In general, this function will be used in
% conjunction with the parameter-space explorer (opt_optim.type = 4), but
% most of it will also work with other optimisation routines.
%
% Note: this should also work in case there are multiple data sets per
% state. However, startgrid_guts will only use the time vector of the first
% one for its derivation of search ranges for parspace.
% 
% Select what to fit with fit_tox (this is a 2-element vector). First 
% element of fit_tox is which part of the data set to use:
%   fit_tox(1) = -1  control survival (c=0) only
%   fit_tox(1) = 0   not used for GUTS
%   fit_tox(1) = 1   fit tox parameter, but when fitting, keep hb fixed; run through all
%               elements in SEL sequentially and provide a table at the end
%   fit_tox(1) = 2   fit tox parameter, but when fitting, also fit hb; run through all
%               elements in SEL sequentially and provide a table at the end
% 
% Second element of fit_tox is whether to fit or only to plot:
%   fit_tox(2) = 0   don't fit; for standard optimisations, plot results for
%               parameter values in [par], for parspace optimisations, use saved mat file.
%   fit_tox(2) = 1   fit parameters
% 
% Author     : Tjalling Jager 
% Date       : November 2021
% Web support: http://www.debtox.info/byom.html

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo glo2 DATA X0mat

% first, a quick check whether there could be problems with the input
if size(DATA,1) > 1
    if isfield(glo,'names_sep') && ~isempty(glo.names_sep) && fit_tox(1)==-1 && any(strcmp(glo.names_sep,'hb'))
        error('Function automatic_run_guts_immob cannot (yet) work properly with hb differing per data set!')
    end
end
if isfield(glo,'int_scen') && ~isempty(glo.int_scen) && ismember(0.1,glo.int_scen)
    if fit_tox(1) == -1
        warning('Note that automatic_run_guts_immob cannot (yet) identify solvent controls with ID=0.1! The hb is fitted on the control only.')
    end
end

basenm_rem = glo.basenm; % remember basename as we will modify it!
best_sel   = [NaN NaN]; % this will be used to return the best death mechanism

% Remember settings so we can turn back everything at the end
glo_rem       = glo;
glo2_rem      = glo2;
X0mat_rem     = X0mat;
DATA_rem      = DATA;

opt_optim.fit = 1; % fit the parameters (1), or don't (0)
% Note: this option is placed here, rather than in the main script. Reason
% is that, perhaps counter intuitively, opt_optim.fit needs to be 1 for the
% parameter-space explorer, even when NOT fitting, otherwise the saved set
% is not loaded. The second element of fit_tox takes over as switch for
% fitting or not fitting.
opt_optim.ps_saved = 0; % use saved set for parameter-space explorer (1) or not (0);
% Note: this is set to fit parameters. Users have to set fit_tox(2) to 0 to
% use the saved set! This setting then automatically sets this
% opt_optim.ps_saved to 1 when needed.

if isempty(skip_sg) % if this variable is not input ...
    skip_sg = 0; % use startgrid for search ranges
end

id_ctrl = 0; % control is treatment that has identifier 0

switch fit_tox(1)
                
    case -1 % for fitting the survival control response only ...
        % This piece of code selects only the treatment with identifier 0
        % (zero), removes all data apart from survival, and turns fitting
        % off for all parameters but background hazard.
        
        % Some GUTS calculations can use feedbacks, with the matrix CONFIG.
        % The exact setting feedback is not relevant for the
        % control, but it must be defined to prevent errors.
        glo.damconfig = CONFIG(1); % use FIRST entry in CONFIG matrix
        glo.sel       = SEL(1); % use FIRST entry in SEL matrix
        
        % Specifically for the immobility package: remove all data other
        % than those for the dead animals. However, startgrid_guts will
        % only use the time vector of the first one for its derivation of
        % search ranges.
        for i_sets = 1:size(DATA,1)
            DATA(i_sets,[glo.locC glo.locD glo.loc_h glo.loc_i glo.loc_id]) = {0}; % remove all but the survival data from the complete data set
            % Note: this also works in case there are multiple data sets per state. 
        end
        X0mat = X0mat(:,X0mat(1,:)==id_ctrl); % only keep controls
        
        data_tst = DATA{glo.loc_d}(:,2:end);
        data_tst(:,data_tst(1,:)~=id_ctrl) = [];
        data_tst = data_tst(3:end,:);
        
        pmat      = packunpack(1,par,[]); % turn parameter structure into a matrix
        pmat(:,2) = 0; % don't fit ANY parameters now
        par       = packunpack(2,[],pmat); % turn parameter matrix back into a structure
        
        if fit_tox(2) == 0 % don't fit, just use values in par or saved set
            if opt_optim.type ~= 4 % for simplex fitting
                opt_optim.fit = 0; % fit the parameters (1), or don't (0)
            else
                opt_optim.ps_saved = 1; % use saved set for parameter-space explorer (1) or not (0);
                opt_optim.fit      = 1; % fit the parameters (1), or don't (0)
                % Fit needs to be 1 otherwise the saved set is not loaded.
            end
        else
            par.hb(2) = 1;    % but DO fit the background hazard rate
            if opt_optim.type == 4 && sum(data_tst(:)) > 0 % for the parspace explorer ...
                par = startgrid_guts_immob(par); % replace par by estimates based on data set (as used by openGUTS)
            else
                par.hb(1) = 0.01; % and start from a reasonable default (if parspace explorer is NOT used)
            end
        end
        
        glo2.ctot = 0; % make sure that the total concentration vector to analyse only has a zero
        
        opt_plot.statsup = glo.locD; % vector with states to suppress in plotting fits (no need for damage)
        % No need for annotations or legends as we fit hb only
        opt_plot.annot   = 0; % annotations in multiplot for fits: 1) box with parameter estimates 2) single legend
        opt_plot.legsup  = 1; % set to 1 to suppress legends on fits
        glo.basenm       = [basenm_rem,'_hb']; % modify basename to make a separate PS.mat file
        
        if sum(data_tst(:)) == 0
            disp('There are no dead animals in the control, so can stop here')
            par.hb(1) = par.hb(3); % set to minimum allowed
            par_out   = par;
        else
            % optimise and plot (fitted parameters in par_out)
            par_out = calc_optim(par,opt_optim); % start the optimisation
            calc_and_plot(par_out,opt_plot); % calculate model lines and plot them
        end
        
        print_par(par_out,-2) % this prints out the optimised parameter values in a
        % formatted way so they can be directly copied into the main script; the
        % -2 means that only fitted parameters are displayed
        
        if opt_optim.type ~= 4 && ~isempty(varargin)
            opt_prof = varargin{1};
            calc_proflik(par_out,'all',opt_prof,opt_optim);  % calculate a profile
        end
        
    case 0 % this one is in the DEBtox version, but not in GUTS
        
        error('This option (fit_tox=0) is not available in GUTS')
         
    otherwise % for fitting tox parameters ...
        
        if fit_tox == 1 % then we fix hb to its value in par
            par.hb(2) = 0; % turn fitting for the background hazard off
        end
        par_rem = par; % remember the old par as we'll modify it in the loop below
        
        if fit_tox(2) == 0 % don't fit, just use values in par or saved set
            if opt_optim.type ~= 4 % for simplex fitting
                opt_optim.fit = 0; % fit the parameters (1), or don't (0)
            else
                opt_optim.ps_saved = 1; % use saved set for parameter-space explorer (1) or not (0);
                opt_optim.fit      = 1; % fit the parameters (1), or don't (0)
                % Fit needs to be 1 otherwise the saved set is not loaded.
            end
        else
            pmat = packunpack(1,par,[]); % turn output parameter structure into a matrix
        end
        
        % GUTS-immobility calculations use a damage configuration, with the
        % matrix FEEDB. 
        if isempty(CONFIG) % if this variable is not in the input ...
            CONFIG = 1; % define as configuration 1
            config_chk = 0;
        else
            config_chk = 1;
        end
        
        % Run through settings for the death mechanism as given in SEL (SD,
        % IT and/or mixed), and configurations (if used) as given in
        % CONFIG. The basename, as used for naming the mat file and the
        % output plots, is adapted to identify the settings.
        MLL_coll = zeros(length(SEL),length(CONFIG)); % initialise matrix to catch MLLs
        par_coll = cell(length(SEL),length(CONFIG));  % initialise matrix to catch parameter sets
        for i = 1:length(SEL) % run through all MoAs
            for j = 1:length(CONFIG) % run through all feedback configuration
                
                par           = par_rem; % reset par before modifying it
                glo.sel       = SEL(i); % change global for death mechanism
                glo.damconfig = CONFIG(j); % change global for feedback configuration
                if config_chk == 0 % no feedbacks, so no need to add to basename
                    glo.basenm = [basenm_rem,'_sel',sprintf('%d',glo.sel)]; % create a new basenm for each death mechanism
                else
                    glo.basenm = [basenm_rem,'_sel',sprintf('%d',glo.sel),'_config',sprintf('%d',glo.damconfig)];
                end
                
                switch glo.sel % make sure that right parameters are fitted
                    % This is also included in <startgrid_guts>, but would
                    % still be needed for opt_optim.type ~= 4.
                    case 1 % for SD ...
                        if isfield(par,'Fs') % then we have a single spread factor
                            par.Fs(2) = 0; % never fit the threshold spread!
                        else
                            par.Fsi(2) = 0; % never fit the threshold spread!
                            par.Fsd(2) = 0; % never fit the threshold spread!
                        end
                    case 2 % for IT ...
                        par.bii(2) = 0; % never fit the killing rate!
                        par.bid(2) = 0; % never fit the killing rate!
                        par.bir(2) = 0; % never fit the killing rate!
                end
                % automatic creation of search ranges for GUTS tox parameters
                if fit_tox(2) == 1 && opt_optim.type == 4 && skip_sg == 0 % only when using parspace explorer for fitting
                    par = startgrid_guts_immob(par); % calculate the search ranges and settings from data set
                end
                
                % optimise and plot (fitted parameters in par_out)
                [par_tmp,MLL] = calc_optim(par,opt_optim); % start the optimisation
                calc_and_plot(par_tmp,opt_plot); % calculate model lines and plot them
                drawnow
                MLL_coll(i,j) = MLL; % collect the MLL so we can compare them later
                par_coll{i,j} = par_tmp; % collect the par so we can output the best one
            end
        end
        
        if fit_tox(2) == 0 % don't fit, just use values in par or saved set
            AIC_coll = nan(length(SEL),length(CONFIG)); % use NaNs matrix to catch MLLs
            % If you did not fit, AIC is meaningless
            if length(SEL) == 1 && length(CONFIG) == 1
                best_sel = [1 1]; % indices for best configuration
                par_out  = par_coll{1,1}; % output the best par
            end
            % If a script is run to display a saved set, we can safely say
            % it's the best (that may prevent an unnecessary error).
        else
            AIC_coll  = 2*sum(pmat(:,2))+2*MLL_coll; % calculate AIC from MLL
        end
        AIC_coll  = AIC_coll - min(AIC_coll(:)); % calculate delta AIC relative to best one
        loss_coll = exp(-AIC_coll/2); % relative probability to minimise information loss
        
        % display on screen and in diary
        diary ('results.out') % collect output in the diary "results.out"
        disp(' ')
        disp('Death mech. damage config. MLL    delta-AIC     prob.    best')
        disp('====================================================================')
        for i = 1:length(SEL) % run through all MoAs
            for j = 1:length(CONFIG) % run through all feedback configuration
                str1 = sprintf('%d',SEL(i)); % string for the MoA configuration
                str2 = sprintf('%d',CONFIG(j)); % string for the feedbacks configuration
                fprintf('%s            %s     %#10.2f %#10.2f %#10.2f ',str1,str2,MLL_coll(i,j),AIC_coll(i,j),loss_coll(i,j))
                if AIC_coll(i,j) == 0 % if this is the best one ...
                    fprintf('%s','     *') % give it a star
                    best_sel = [i j]; % indices for best configuration
                    par_out  = par_coll{i,j}; % output the best par
                end
                fprintf('\n')
            end
            if config_chk == 1
                disp('====================================================================')
            end
        end
        if config_chk == 0
            disp('====================================================================')
        end
        diary off  % turn off the diary function

end

% Return settings to their original value
if fit_tox(1) ~= 1 % when we fit tox data, it is good to keep the modified glo
    glo   = glo_rem;
end
glo2      = glo2_rem;
X0mat     = X0mat_rem;
DATA      = DATA_rem;
