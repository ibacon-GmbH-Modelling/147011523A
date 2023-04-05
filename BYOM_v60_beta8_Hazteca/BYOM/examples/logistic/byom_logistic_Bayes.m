%% BYOM, byom_logistic_Bayes.m, a quick fitting example 
%
% * Author: Tjalling Jager 
% * Date: November 2021
% * Web support: <http://www.debtox.info/byom.html>
% * Back to index <walkthrough_byom.html>
%
% BYOM is a General framework for simulating model systems in terms of
% ordinary differential equations (ODEs). The model itself needs to be
% specified in <derivatives.html derivatives.m>, and <call_deri.html
% call_deri.m> may need to be modified to the particular problem as well.
% The files in the engine directory are needed for fitting and plotting.
% Results are shown on screen but also saved to a log file (results.out).
%
% *The model:* The standard logistic growth model, with parameters for the
% population growth rate (_r_) and the carrying capacity (_K_). Optionally,
% a constant predation rate (_p_) can be included. Model system:
%
% $$ \frac{dN}{dt} = r N \left(1-\frac{N}{K} \right) -p $$
%
% *This script:* Data set for algal growth from Schanz & Zahler (1981).
% Fitting with simplex optimisation, followed by slice sampling to make
% this a Bayesian analysis.
% 
%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

%% Initial things
% Make sure that this script is in a directory somewhere *below* the BYOM
% folder.

clear, clear global % clear the workspace and globals
global DATA W X0mat % make the data set and initial states global variables
global glo          % allow for global parameters in structure glo
diary off           % turn off the diary function (if it is accidentaly on)
% set(0,'DefaultFigureWindowStyle','docked'); % collect all figure into one window with tab controls
set(0,'DefaultFigureWindowStyle','normal'); % separate figure windows

pathdefine(0) % set path to the BYOM/engine directory (option 1 uses parallel toolbox)
glo.basenm  = mfilename; % remember the filename for THIS file for the plots
glo.saveplt = 0; % save all plots as (1) Matlab figures, (2) JPEG file or (3) PDF (see all_options.txt)

%% The data set
% Data are entered in matrix form, time in rows, scenarios (exposure
% concentrations) in columns. First column are the exposure times, first
% row are the concentrations or scenario numbers. The number in the top
% left of the matrix indicates how to calculate the likelihood:
%
% * -1 for multinomial likelihood (for survival data)
% * 0  for log-transform the data, then normal likelihood
% * 0.5 for square-root transform the data, then normal likelihood
% * 1  for no transformation of the data, then normal likelihood

% observed population size (in ml KMnO4)
DATA{1} = [1 1
    5	0.53
    7	1.29
    8	1.65
    9	2.25
    10	2.31
    11	2.91
    12	3.62
    13	4.11
    14	4.39
    16	4.67
    19	4.93];  

% Removing 4 highest points produces identifiability issues

% if weight factors are not specified, ones are assumed in start_calc.m

%% Initial values for the state variables
% Initial states, scenarios in columns, states in rows. First row are the
% 'names' of all scenarios.

X0mat = [1      % the scenario(s) to run
         1];  % initial values state 1: population size at t=0 (overwritten by parameter N0)

%% Initial values for the model parameters
% Model parameters are part of a 'structure' for easy reference. 

% syntax: par.name = [startvalue fit(0/1) minval maxval];
par.r    = [0.5  1 0 5];  % population growth rate, d-1
par.K    = [5    1 0 20]; % carrying capacity, ml
par.N0   = [0.1  1 0 1];  % initial population size, ml
par.p    = [0    0 0 1];  % predation rate, ml d-1

%% Priors for Bayesian analyses
% Optionally, prior distributions can be specified for parameters, see the
% file calc_prior.m for the definition of the distributions. You must use
% the exact same names for the prior parameters as used in the _par_
% structure. If you do not specify _pri_ in your scripts, uniform priors
% are assumed, defined by the min-max range in _par_ above.
 
% % UNCOMMENT LINE(S) TO TRY PRIORS
% % First element in pri is the choice of distribution.
% pri.K   = [2 3 7 5];    % triangular with min, max and center
% pri.N0  = [3 0.1 0.05]; % normal with mean and sd
% % Note that the prior is always defined on normal scale, also when the
% % parameter will be fitted on log scale.

%% Time vector and labels for plots
% Specify what to plot. If time vector glo.t is not specified, a default is
% constructed, based on the data set.

% specify the y-axis labels for each state variable
glo.ylab{1} = 'population size (ml KMnO4)';
% specify the x-axis label (same for all states)
glo.xlab    = 'time (days)';
glo.leglab1 = 'Scen. '; % legend label before the 'scenario' number
glo.leglab2 = '';   % legend label after the 'scenario' number

prelim_checks % script to perform some preliminary checks and set things up
% Note: prelim_checks also fills all the options (opt_...) with defauls, so
% modify options after this call, if needed.

%% Calculations and plotting
% Here, the functions are called that will do the calculation and the
% plotting. 
% 
% Options for the plotting can be set using opt_plot (see prelim_checks.m).
% Options for the optimsation routine can be set using opt_optim. Options
% for the ODE solver are part of the global glo. For the demo, the
% iterations were turned off (opt_optim.it = 0).

opt_optim.fit = 1; % fit the parameters (1), or don't (0)
opt_optim.it  = 1; % show iterations of the optimisation (1, default) or not (0)
glo.useode    = 0; % calculate model using ODE solver (1) or analytical solution (0)

% optimise and plot (fitted parameters in par_out)
par_out = calc_optim(par,opt_optim); % start the optimisation
calc_and_plot(par_out,opt_plot);     % calculate model lines and plot them

%% Slice sampler
% The slice sampler can be used for a Bayesian analysis as it provides a
% sample from the posterior distribution. A .mat file is saved which
% contains the sample, to use later for e.g., intervals on model
% predictions. The output includes the Markov chain and marginal
% distributions for each fitted parameter. The function calc_conf can be
% used to put confidence intervals on the model lines. 
%
% Notes: MCMC is also slow process for complex models. If a prior is
% specified, it is plotted in the marginal posterior plot as a red
% distribution. Two types of credible intervals are generated:
%
% # Highest probability regions. This captures the 95% of the parameter
% values with the highest posterior probability. it is calculated by
% lowering the horizotal dotted line in the figures for each parameters,
% until the area between the cut-off points captures 95% of the density.
% This can lead to unequal probabilities in the tails.
% # Quantiles. The region between the 2.5 and 97.5% quantiles is taken.
% This leads to 2.5% probability in each tail. This is shown in the plots
% by the broken vertical lines.
% 
% Options for the slice sampling can be set using opt_slice (see
% prelim_checks.m). Options for confidence bounds on model curves can be
% set in opt_conf (see prelim_checks).
% 
% Taking on parameters on log scale helps, as the slice sampler seems quite
% bad at obtaining a representative sample when the parameters have very
% different ranges. Thinning will be needed to reduce the autocorrelation.
% Note the plot that is produced with details of the sample!

% UNCOMMENT LINE(S) TO CALCULATE
opt_slice.thin     = 50; % thinning of the sample (keep one in every 'thin' samples)
opt_slice.burn     = 100; % number of burn-in samples (0 is no burn in)
opt_slice.alllog   = 1; % set to 1 to put all parameters on log-scale before taking the sample
calc_slice(par_out,1000,opt_slice); % second argument number of samples (-1 to re-use saved sample from previous runs)

%% Plot results with confidence intervals
% The following code can be used to make a standard plot (the same as for
% the fits), but with confidence intervals. Options for confidence bounds
% on model curves can be set using opt_conf (see prelim_checks).
% 
% Use opt_conf.type to tell calc_conf which sample to use: 
% -1) Skips CIs (zero does the same, and an empty opt_conf as well).
% 1) Bayesian MCMC sample (default); CI on predictions are 95% ranges on 
% the model curves from the sample 
% 2) parameter sets from a joint likelihood region using the shooting 
% method (limited sets can be used), which will yield (asymptotically) 95% 
% CIs on predictions
% 3) as option 2, but using the parameter-space explorer

opt_conf.type    = 1; % make intervals from 1) slice sampler, 2) likelihood region shooting, 3) parspace explorer
opt_conf.lim_set = 0; % use limited set of n_lim points (1) or outer hull (2, likelihood methods only) to create CIs
opt_conf.sens    = 0; % type of analysis 0) no sensitivities 1) corr. with state, 2) corr. with state/control, 3) corr. with relative change of state over time

out_conf = calc_conf(par_out,opt_conf);   % calculate confidence intervals on model curves
calc_and_plot(par_out,opt_plot,out_conf); % call the plotting routine again to plot fits with CIs
