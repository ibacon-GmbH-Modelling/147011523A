function [DATA,W,LabelTable] = data_AZT_Hazteca_calibration_nosurv(TR,opt,study)

% Handy function to load the data sets from standard text files, and
% prepare them for BYOM analyses.
%
% This file: calibration data set for AZT in Daphnia magna. Data received
% as text files from Marie Trijau (Ibacon) on 25/6/2021 by email.
% Specifications:
% - Repro data are live offspring only. Since there is a direct effect on
%   offspring this is the best way to include mortality in the analysis
%   without additional survival module for the neonates in the brood pouch.
% - -1 in repro data specifies first egg appearance, and moults without
%   neonate release.
% - A NaN for repro data is entered on the time point AFTER mother dies, so
%   any repro is still accounted for. At this moment, make_repro_ind does 
%   NOT compensate this situation with weight factor or averaged-expected
%   repro.
% 
% Modifications made by Tjalling Jager:
% - The file for length data included some quotes before or after the NaNs.
%   These were removed by loading in Excel and filling the cells with new NaN
%   text.
%
% * Author: Tjalling Jager
% * Date: July 2021
% * Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo

% survivors
S = load(fullfile('data_AZT_hyalella','Data_surv_Hazteca_AZT_calibration.txt'));
S(1,1) = -1; % make sure to use multinomial likelihood for survival

% body length (mm) on each observation time (d)
L = load(fullfile('data_AZT_hyalella','Data_length_Hazteca_AZT_calibration.txt'));
L(1,1) = TR; % make sure to use transformation as specified by input
L(2,2:end) = 2.665;
LW = ones(size(L)-1); % these are individual replicates, so weights are ones
% NOTE: the initial lengths are the same 10 numbers for each treatment.

% Reproduction, as neonates released, as observed on each observation time (d)
R = load(fullfile('data_AZT_hyalella','Data_repro_Hazteca_AZT_calibration.txt'));
R(1,1) = TR; % make sure to use transformation as specified by input

% number of females in the replicates
F = load(fullfile('data_AZT_hyalella','Data_females_Hazteca_AZT_calibration.txt'));

[data,w] = makerepro_sex(R,S,F);
R = data;

if opt == 0 % if we only want to check the repro data, we can stop here
    DATA       = [];
    W          = [];
    LabelTable = [];
    return
end

% and define a scenario for the exposure treatments (ug/L??)
Cw = load(fullfile('data_AZT_hyalella','Data_exposure_Hazteca_AZT_calibration.txt'));
Cw_type = 4; % linear interpolation pulses   

% next, remove the controls from Cw; it is probably faster to define them
% separately since there is no need to run through them in steps
[~,loc_ctrl] = ismember([0 0.1],Cw(1,2:end));
Cw(:,loc_ctrl+1) = [];

% Define exposure scenarios for controls (for id=0.1 this is no exposure,
% as it is the solvent control)
Cw0 = [0         0 0.1 
       0         0 0 
       Cw(end,1) 0 0 ];

Cw0_type = 2; % block pulses   

% Create a table with nicer custom labels for the legends
Scenario = [0 0.1 1:4]'; % scenario identifiers that get a label; % scenario identifiers that get a label
Label = {'control';'solvent control';'20';'40';'80';'160'};

% Modify the scenario identifiers using the study number provided
S(1,2:end)   = S(1,2:end)   + (study-1)*100;
L(1,2:end)   = L(1,2:end)   + (study-1)*100;
R(1,2:end)   = R(1,2:end)   + (study-1)*100;
Cw0(1,2:end) = Cw0(1,2:end) + (study-1)*100;
Cw(1,2:end)  = Cw(1,2:end)  + (study-1)*100;
Scenario     = Scenario     + (study-1)*100;

% Create a table with nicer custom labels for the legends
LabelTable = table(Scenario,Label); % create a Matlab table for the labels
glo.LabelTable = [LabelTable]; % temporarily defined for make_scen to provide labels

glo.scen_plot = 0; % don't make a plot for the solvent scenario
make_scen(Cw0_type,Cw0); % type 2 creates block pulses (fine for controls and constant exposure)

glo.scen_plot = 1; % but do make a plot for the exposure scenarios 
make_scen(Cw_type,Cw); % type 4 creates linear interpolation, type 2 creates block pulses

S=[0]; % Exclude now survival data!!

DATA{1,1} = [0]; % there are never data for scaled damage (state 1)
DATA{1,2} = L;   % length data
W{1,2}    = LW;  % length weights data
DATA{1,3} = R;   % reproduction data
W{1,3}    = w;  % reproduction weights data
DATA{1,4} = S;   % survival data

