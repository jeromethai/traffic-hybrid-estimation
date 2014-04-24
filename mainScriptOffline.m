%% load data

clear all
close all
clc

%load data/data0301_7_9.mat
% load data/data0301_7_19.mat
%load data/data0301_17_19.mat

%load data/data0305_7_9.mat
load data/data0305_7_19.mat
%load data/data0305_17_19.mat


PLOT_FIG = true;
% ALGORITHM = 'EnKFmode';
% ALGORITHM = 'EKFmode';
% ALGORITHM = 'RIMMSimple';
% ALGORITHM = 'RIMM';
ALGORITHM = 'RIMMtest';
% ALGORITHM = 'EKFmodeLL';
 
%% run algorithm
THRESHOLD = 1;

[rho, vel, output] = pCTM(route, VFF, RHO_MAX, RHO_C, DT, NUM_ENSEMBLES, ...
    PLOT_FIG, ALGORITHM, THRESHOLD); 