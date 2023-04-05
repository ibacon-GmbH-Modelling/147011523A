function [DATA,W,LabelTable] = data_AZT_Hazteca_validation1_nosurv(TR,opt,study)

% Handy function to load the data sets from standard text files, and
% prepare them for BYOM analyses.
%
% This file: validation data set for AZT in Daphnia magna. Data received
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
% - Repro is quite erratic with individuals releasing on neonates on
% consecutive days. Also in the controls.
% 
% Modifications made by Tjalling Jager:
% - The file for length data was unicode text rather than tab-delimited.
% - Added a -1 at several positions in T2-T6, at places where I suspect a 
%   moult without neonate release.
%
% * Author: Tjalling Jager
% * Date: July 2021
% * Web support: <http://www.debtox.info/byom.html>

%  Copyright (c) 2012-2021, Tjalling Jager, all rights reserved.
%  This source code is licensed under the MIT-style license found in the
%  LICENSE.txt file in the root directory of BYOM. 

global glo

% survivors
S = load(fullfile('data_AZT_hyalella','PS1-B_Data_surv3_Hazteca_AZT_validation.txt'));
S(1,1) = -1; % make sure to use multinomial likelihood for survival

% body length (mm) on each observation time (d)
L = load(fullfile('data_AZT_hyalella','PS1-B_Data_length_Hazteca_AZT_validation.txt'));
L(1,1) = TR; % make sure to use transformation as specified by input
L(2,2:end) = 2.373;
LW = ones(size(L)-1); % these are individual replicates, so weights are ones
% NOTE: the initial lengths are the same 10 numbers for each treatment.

% Reproduction, as neonates released, as observed on each observation time (d)
R = load(fullfile('data_AZT_hyalella','PS1-B_Data_repro_Hazteca_AZT_validation.txt'));
R(1,1) = TR; % make sure to use transformation as specified by input

% number of females in the replicates
F = load(fullfile('data_AZT_hyalella','PS1-B_Data_females_Hazteca_AZT_validation.txt'));

% create intermediate arrays to slice data only for the reproduction
Si=S;
Fi=F;
Ri=R;
%Fi(:,5)=[];Fi(:,26)=[];Si(:,5)=[];Si(:,26)=[];Ri(:,5)=[];Ri(:,26)=[];
%Ri(isnan(Ri))=0;
% take only the first 3 replicate
S1 = [Si(:,1:4) Si(:,12:14) Si(:,22:24) Si(:,32:34) Si(:,42:44)];
F1 = [Fi(:,1:4) Fi(:,12:14) Fi(:,22:24) Fi(:,32:34) Fi(:,42:44)];
R1 = [Ri(:,1:4) Ri(:,12:14) Ri(:,22:24) Ri(:,32:34) Ri(:,42:44)];

%take only the last replicates and exclude those that do not have females
%measured
S2 = [Si(:,1) Si(:,6:11) Si(:,15:21) Si(:,25:26) Si(:,28:31) Si(:,35:41) Si(:,45:51)];
F2 = [Fi(:,1) Fi(:,6:11) Fi(:,15:21) Fi(:,25:26) Fi(:,28:31) Fi(:,35:41) Fi(:,45:51)];
R2 = [Ri(:,1) Ri(:,6:11) Ri(:,15:21) Ri(:,25:26) Ri(:,28:31) Ri(:,35:41) Ri(:,45:51)];
% remove time 14d as they are all NaN for this time measurement
S2(6,:) = []; F2(2,:) = []; R2(6,:) = [];

% compute reproduction separately for the 2 datasets
[data1,w1] = makerepro_sex(R1,S1,F1);
R1 = data1;

[data2,w2] = makerepro_sex(R2,S2,F2);
R2 = data2;

% join again together the reproduction data
R2 = [R2(1:5,:);(NaN(1, size(R2,2))); R2(6:end,:)];
Rt = [R1(:,1:4) R2(:,2:7) R1(:,5:7) R2(:,8:14) R1(:,8:10) R2(:,15:20) R1(:,11:13) R2(:,21:27) R1(:,14:16) R2(:,28:34)];
w2 = [w2(1:4,:);(NaN(1, size(w2,2))); w2(5:end,:)];
wt = [w1(:,1:3) w2(:,1:6) w1(:,4:6) w2(:,7:13) w1(:,7:9) w2(:,14:19) w1(:,10:12) w2(:,20:26) w1(:,13:15) w2(:,27:33)];

if opt == 0 % if we only want to check the repro data, we can stop here
    DATA       = [];
    W          = [];
    LabelTable = [];
    return
end

% and define a scenario for the exposure treatments (ug/L??)
Cw = load(fullfile('data_AZT_hyalella','PS1-B_Data_exposure_Hazteca_AZT_validation.txt'));
Cw_type = 4; % block pulses   

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
Scenario = [Cw0(1,2:end) Cw(1,2:end)]'; % scenario identifiers that get a label
Label = {'control';'solvent control';'PS1-B 1';'PS1-B 2';'PS1-B 3'};%;'wide 1';'wide 2';'wide 3'};

% Modify the scenario identifiers using the study number provided
%S(1,2:end)   = S(1,2:end)   + (study-1)*100;
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

DATA{1,1} = [0]; % there are never data for scaled damage (state 1)
%DATA{2,1} = [0];
DATA{1,2} = L;   % length data
%DATA{2,2} = [0];
W{1,2}    = LW;  % length weights data
%W{2,2}    = [0];
DATA{1,3} = Rt;   % reproduction data
%DATA{2,3} = R2;
W{1,3}    = wt;  % reproduction weights data
%W{2,3}    = w2;
%DATA{1,4} = S;   % survival data
%DATA{2,4} = [0];


