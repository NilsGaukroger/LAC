% 46320 LAC Course
% Assignment 3: Post-processing
% Authors: Nikolaos Stamatopoulos, Nils Joseph Gaukroger, Stefanos Masalas,
% Eduardo Fredrich de Miranda
% 28th October 2021

close all; clear variables; clc

%% Add functions folder to path
addpath('functions\')

%% Set default plot settings
disp('Setting default plot parameters...');
set(    0,          'defaulttextInterpreter', 'tex');
set(groot, 'defaultAxesTickLabelInterpreter', 'tex');
set(groot,        'defaultLegendInterpreter', 'tex');
set(    0,             'defaultAxesFontSize',    15);
set(    0,            'DefaultLineLineWidth',     2);
disp('Default plot parameters set.');

%% Figure saving settings
save_var = false; % true: saves figures, false: doesn't
local    = './plots/report_3/';
locs     = {local};
for i = 1:length(locs) % if any directory doesn't exist don't attempt to save there
    if (not(isfolder(locs{i})))
        save_var = false;
    end
end

%% Load data from postProcessing_stab
load('postProcessing_stab.mat');

%% Create table for gains (cases C1-C6)
caseName        = ["C1";"C2";"C3";"C4";"C5";"C6"];
generatorTorque = ["Constant Power";"Constant Power";"Constant Power";...
    "Constant Torque";"Constant Torque";"Constant Torque"];
genTorque_switch = [1;1;1;0;0;0];
partialLoad_omega = 0.05 * ones(6,1);
partialLoad_zeta  = 0.7  * ones(6,1);
fullLoad_omega  = [0.05; 0.01; 0.10; 0.05; 0.01; 0.10]; % natural frequency [Hz]
fullLoad_zeta   = [ 0.7;  0.7;  0.7;  0.7;  0.7;  0.7]; % damping ratio [-]

gains = table(caseName,generatorTorque,partialLoad_omega,partialLoad_zeta,fullLoad_omega,fullLoad_zeta);

%% Add details for C7
gains = [gains;{"C7","Constant Power", 0.05, 0.7, 0.03, 0.7}];
genTorque_switch = [genTorque_switch; 1];

%% Import gains from *_ctrl_tuning.txt
% File import settings
order = 2; % Order of aerodynamic gain scheduling fit (1 = linear, 2 = quadratic)

% Preallocation
[gains.("Region 1: K"),gains.("Region 2: I"),gains.("Region 2: Kp"),gains.("Region 2: Ki"),...
    gains.("Region 3: Kp"),gains.("Region 3: Ki"),...
    gains.("Region 3: K1"),gains.("Region 3: K2")] = deal(NaN(size(gains,1),1));

for i = 1:size(gains,1)
    filename = strcat("./your_model/results_redesign/cont/", gains.caseName(i), "/redesign_cont_HS2_ctrl_tuning.txt");
    if isfile(filename) % check if file exists
        [gains.("Region 1: K")(i),gains.("Region 2: I")(i),gains.("Region 2: Kp")(i),gains.("Region 2: Ki")(i),...
            gains.("Region 3: Kp")(i),gains.("Region 3: Ki")(i),...
            gains.("Region 3: K1")(i),gains.("Region 3: K2")(i)] = import_gains(filename,order);
    else
        fprintf('Controller tuning .txt for case C%d does not exist\n',i)
    end
end

%% Write gains.dat files
for i = 1:size(gains,1)
    if any(isnan(gains{i,7:end}))
        fprintf('Missing gains for C%d\n',i)
    else
        write_gains(gains,genTorque_switch,redesign,i)
    end
end

%% Import HAWC2 results
data = cell(size(gains,1),1);
for i = 1:length(data)
    filename = strcat('your_model\results_redesign\cont\',gains.caseName(i),'\redesign_cont.dat');
    if isfile(filename) % check if file exists
        data{i} = readtable(filename);
    else
        fprintf('HAWC2 results .dat file for case C%d does not exist\n',i)
    end
end

%% Plot HAWC2 results
idx  = [15,4,10,100]; scale = [1,1,1,1e-3];
ylabels = ["Wind speed [m/s]", "Pitch angle [deg]", "Rotational speed [rad/s]", "Electrical power [kW]"];

figure
for j = 1:4
    subplot(2,2,j)
    for i = 3:-1:1
        plot(data{i}{:,1},data{i}{:,idx(j)}.*scale(j))
        hold on
    end
    hold off
    if j > 2
        xlabel('Time [s]');
    end
    ylabel(ylabels(j))
    legend('C3','C2','C1','Location','NW')
    box on
    grid on
end
if save_var
    saveFig(locs,'C1-C3',"png");
end

figure
for j = 1:4
    subplot(2,2,j)
    for i = 6:-1:4
        plot(data{i}{:,1},data{i}{:,idx(j)}.*scale(j))
        hold on
    end
    hold off
    if j > 2
        xlabel('Time [s]');
    end
    ylabel(ylabels(j))
    legend('C6','C5','C4','Location','NW')
    box on
    grid on
end
if save_var
    saveFig(locs,'C4-C6',"png");
end

% C1 vs. C4
x = [1,4];
idx = [99,100];
ylabels = ["LSS torque [Nm]", "Electrical power [kW]"];

figure
for j = 1:2
    subplot(2,2,j)
    for i = 1:length(x)
        plot(data{x(i)}{:,1},data{x(i)}{:,idx(j)}.*scale(j))
        hold on
    end
    hold off
    if j > 2
        xlabel('Time [s]');
    end
    ylabel(ylabels(j))
    legend('C1','C4','Location','NW')
    box on
    grid on
    xlim([387 inf])
end
if save_var
    saveFig(locs,'C1vsC4',"png");
end

% Low-pass filtered tower fore-aft acceleration
x = [1,3,6];
idx = [96,3,94,95];
ylabels = ["Low-pass filtered tower fore-aft acceleration [m/s^2]", "Pitch angle [deg]", "Flag for mechanical brake [-]", "Flag for emergency pitch stop [-]"];

figure
for j = 1:4
    subplot(2,2,j)
    for i = 1:length(x)
        plot(data{x(i)}{:,1},data{x(i)}{:,idx(j)}.*scale(j))
        hold on
    end
    hold off
    if j > 2
        xlabel('Time [s]');
    end
    ylabel(ylabels(j))
    legend('C1','C3','C6','Location','NW')
    box on
    grid on
end
if save_var
    saveFig(locs,'C1vsC3',"png");
end

%% C7 tuning
