% 46320 LAC Course
% Assignment 3: Post-processing
% Authors: Nikolaos Stamatopoulos, Nils Joseph Gaukroger, Stefanos Masalas,
% Eduardo Fredrich de Miranda
% 28th October 2021

close all; clear variables; clc

%% Add functions folder to path
addpath('functions\')

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
gains = [gains;{"C7","Constant Power", 0.05, 0.7, 0.065, 0.7}];
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
        fprintf('File %s does not exist.\n',filename)
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
data = cell(6,1);
for i = 1:length(data)
    filename = strcat('your_model\results_redesign\cont\',caseName(i),'\redesign_cont.dat');
    if isfile(filename) % check if file exists
        data{i} = readtable(filename);
    else
        fprintf('File %s does not exist.\n',filename)
    end
end

%% Plot HAWC2 results
idx  = [15,3,10,100]; scale = [1,1,1,1e-3];
ylabels = ["Wind speed [m/s]", "Pitch angle [deg]", "Rotational speed [rad/s]", "Electrical power [kW]"];

figure
for j = 1:4
    subplot(2,2,j)
    for i = 1:3
        plot(data{i}{:,1},data{i}{:,idx(j)}.*scale(j))
        hold on
    end
    hold off
    xlabel('Time [s]'); ylabel(ylabels(j))
    legend('C1','C2','C3')
    box on
    grid minor
end

figure
for j = 1:4
    subplot(2,2,j)
    for i = 4:6
        plot(data{i}{:,1},data{i}{:,idx(j)}.*scale(j))
        hold on
    end
    hold off
    xlabel('Time [s]'); ylabel(ylabels(j))
    legend('C4','C5','C6')
    box on
    grid minor
end

x = [2,5];

figure
for j = 1:4
    subplot(2,2,j)
    for i = 1:length(x)
        plot(data{x(i)}{:,1},data{x(i)}{:,idx(j)}.*scale(j))
        hold on
    end
    hold off
    xlabel('Time [s]'); ylabel(ylabels(j))
    legend('C2','C5')
    box on
    grid minor
end