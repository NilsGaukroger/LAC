% 46320 LAC Course
% Assignment 3: Post-processing
% Authors: Nikolaos Stamatopoulos, Nils Joseph Gaukroger, Stefanos Masalas,
% Eduardo Fredrich de Miranda
% 28th October 2021

close all; clear variables; clc

%% Add functions folder to path
addpath('..\..\functions\')

%% Set default plot settings
disp('Setting default plot parameters...');
set(    0,          'defaulttextInterpreter', 'tex');
set(groot, 'defaultAxesTickLabelInterpreter', 'tex');
set(groot,        'defaultLegendInterpreter', 'tex');
set(    0,             'defaultAxesFontSize',    15);
set(    0,            'DefaultLineLineWidth',     2);
disp('Default plot parameters set.');

%% Figure saving settings
save_var = true; % true: saves figures, false: doesn't
local    = '../figs/';
locs     = {local};
for i = 1:length(locs) % if any directory doesn't exist don't attempt to save there
    if (not(isfolder(locs{i})))
        save_var = false;
    end
end

%% Load data from postProcessing_stab
load('../../mat/postProcessing_stab.mat');

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
    filename = strcat("../res/redesign/", gains.caseName(i), "/redesign_cont_HS2_ctrl_tuning.txt");
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
    filename = strcat('..\res\redesign\',gains.caseName(i),'\redesign_cont.dat');
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

%% Low-pass filtered tower fore-aft acceleration
x = [1,3,6];
idx = [96,3,94,95];
ylabels = {["Low-pass filtered tower"; "fore-aft acceleration [m/s^2]"], "Pitch angle [deg]", "Flag for mechanical brake [-]", "Flag for emergency pitch stop [-]"};

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
    if j == 1
        yline(1.5,'--','LineWidth',1.5)
    end
    ylabel(ylabels{j})
    legend('C1','C3','C6','Location','NW')
    box on
    grid on
end
if save_var
    saveFig(locs,'C1vsC3',"png");
end

%% Import C7 HAWC2 results
filepath = '..\res\redesign\C7\';
folders = dir('..\res\redesign\C7\');
folders = folders(3:end);
folders = folders([folders.isdir] == 1);
C7      = cell(length(folders),2);
for i = 1:length(folders)
    filename = strcat(filepath, folders(i).name, '\redesign_cont.dat');
    C7{i,1} = [str2double(folders(i).name(6:10)), str2double(folders(i).name(16:end))];
    C7{i,2} = readtable(filename);
end

%% Plot C7 HAWC2 results
figure
plot(data{1}{:,1},data{1}{:,10}); hold on
comb = [0.055, 0.7 ;
        0.05 , 0.65;
        0.055, 0.65];
letter = ['a','b','c'];
leg = cell(size(comb,1)+1,1);
leg{1} = "C1:   \omega = 0.050, \zeta = 0.70";
for j = 1:size(comb,1)
    for i = 1:length(folders)
        if C7{i,1} == comb(j,:)
            plot(C7{i,2}{:,1},C7{i,2}{:,10})
        end
    end
    leg{j+1} = sprintf('C7%s: \\omega = %.3f, \\zeta = %.2f',letter(j),comb(j,1),comb(j,2));
end
yline(0.774,'--','LineWidth',1.5) % target rotational speed
hold off
xlabel('Time [s]'); ylabel('Rotational speed [rad/s]')
legend(leg,'FontSize',18)
box on
grid on
xlim([427 468])
if save_var
    saveFig(locs,'C1vsC7_rotationalSpeed_1Step',"png");
end

%%
comb = [0.055, 0.65];
leg = cell(size(comb,1)+1,1);
leg{1} = "C1:   \omega = 0.050, \zeta = 0.70";
bounds = [427, 468;
          674, 714;
          961, 1000];

figure
for k = 1:3
    subplot(1,3,k)
plot(data{1}{:,1},data{1}{:,10}); hold on
for j = 1:size(comb,1)
    for i = 1:length(folders)
        if C7{i,1} == comb(j,:)
            plot(C7{i,2}{:,1},C7{i,2}{:,10})
        end
    end
    leg{j+1} = sprintf('C7%s: \\omega = %.3f, \\zeta = %.2f',letter(3),comb(j,1),comb(j,2));
end
yline(0.774,'--','LineWidth',1.5) % target rotational speed
hold off
xlabel('Time [s]');
if k == 1
    ylabel('Rotational speed [rad/s]')
end
legend(leg)
box on
grid on
xlim([bounds(k,1) bounds(k,2)])
end
if save_var
    saveFig(locs,'C1vsC7_rotationalSpeed_3Steps',"png");
end

%%
figure
plot(data{1}{:,1},data{1}{:,4}); hold on
comb = [0.055, 0.7 ;
        0.05 , 0.65;
        0.055, 0.65];
leg = cell(size(comb,1)+1,1);
leg{1} = "C1:   \omega = 0.050, \zeta = 0.70";
for j = 1:size(comb,1)
    for i = 1:length(folders)
        if C7{i,1} == comb(j,:)
            plot(C7{i,2}{:,1},C7{i,2}{:,4})
        end
    end
    leg{j+1} = sprintf('C7%s: \\omega = %.3f, \\zeta = %.2f',letter(j),comb(j,1),comb(j,2));
end
yline(6.864350,'--','LineWidth',1.5) % target pitch angle
hold off
xlabel('Time [s]'); ylabel('Pitch angle [deg]')
legend(leg,'FontSize',18)
box on
grid on
xlim([427 468])
if save_var
    saveFig(locs,'C1vsC7_pitchAngle_1Step',"png");
end

%%
comb = [0.055, 0.65];
leg = cell(size(comb,1)+1,1);
leg{1} = "C1:   \omega = 0.050, \zeta = 0.70";
bounds = [427, 468;
          674, 714;
          961, 1000];

figure
for k = 1:3
    subplot(1,3,k)
plot(data{1}{:,1},data{1}{:,4}); hold on
for j = 1:size(comb,1)
    for i = 1:length(folders)
        if C7{i,1} == comb(j,:)
            plot(C7{i,2}{:,1},C7{i,2}{:,4})
        end
    end
    leg{j+1} = sprintf('C7%s: \\omega = %.3f, \\zeta = %.2f',letter(3),comb(j,1),comb(j,2));
end
hold off
xlabel('Time [s]');
if k == 1
    ylabel('Pitch angle [deg]')
end
if k == 3
    legend(leg,'Location','SE')
else
    legend(leg)
end
box on
grid on
xlim([bounds(k,1) bounds(k,2)])
end
if save_var
    saveFig(locs,'C1vsC7_pitchAngle_3Steps',"png");
end