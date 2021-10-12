% 46320 LAC Course
% Assignment 2: Post-processing
% Authors: Nikolaos Stamatopoulos, Nils Joseph Gaukroger, Stefanos Masalas,
% Eduardo Fredrich de Miranda
% 12th October 2021

close all; clear variables; clc

%% Load data from postProcessing_struct
load('postProcessing_struct.mat');

%% Modal identification
n_modes = 12; % number of modes

% DTU 10MW
DTU.modes.struc = {'1st Twr FA', '1st Twr SS', '1st BW flap', '1st SYM flap',...
              '1st FW flap', '1st BW edge', '1st FW edge', '2nd BW flap',...
              '2nd FW flap', '2nd SYM flap', '1st COL edge','3rd BW flap'};
DTU.modes.ael = {'1st Twr SS', '1st Twr FA', '1st BW flap', '1st SYM flap',...
              '1st FW flap', '1st BW edge', '1st FW edge', '2nd BW flap',...
              '2nd FW flap', '2nd SYM flap', '1st COL edge', '3rd BW flap'};
% Redesign
redesign.modes.struc = {'1st Twr FA', '1st Twr SS', '1st BW flap', '1st SYM flap',...
               '1st FW flap', '1st BW edge', '1st FW edge', '2nd BW flap',...
               '2nd FW flap', '2nd SYM flap', '1st COL edge','3rd BW flap'};
redesign.modes.ael = {'1st Twr FA', '1st Twr SS', '1st BW flap', '1st SYM flap',...
              '1st FW flap', '1st BW edge', '1st FW edge', '2nd BW flap',...
              '2nd FW flap', '2nd SYM flap', '1st COL edge', '3rd BW flap'};

%% Set default plot settings
disp('Setting default plot parameters...');
set(    0,          'defaulttextInterpreter', 'tex');
set(groot, 'defaultAxesTickLabelInterpreter', 'tex');
set(groot,        'defaultLegendInterpreter', 'tex');
set(    0,             'defaultAxesFontSize',    15);
set(    0,            'DefaultLineLineWidth',   1.5);
disp('Default plot parameters set.');

%% Figure saving settings
save_var = false; % true: saves figures, false: doesn't
overleaf = 'C:\Users\nilsg\Dropbox\Apps\Overleaf\LAC Assignment 2\figures\stab\'; % overleaf not working currently
local    = './plots/report_2/stab/';
locs     = {overleaf,local};
for i = 1:length(locs) % if any directory doesn't exist don't attempt to save there
    if (not(isfolder(locs{i})))
        save_var = false;
    end
end

%% Add functions folder to path
addpath('functions\')

%% Import .cmb files
turbine   = {'results_dtu10mw/stab/','results_redesign/stab/'};
names     = {'DTU 10MW','Redesign'};
type      = {'st','ael'}; % structural 'st' or aeroelastic 'ael' Campbell diagram
file      = {'DTU_10MW_','redesign_'};
cmb.camp     = cell(2,3); % array for storing campbell diagram inputs files, rows = turbine, cols = elasticity (+ turbine name)
cmb.damp     = cell(2,3); % same for damping
cmb.realPart = cell(2,3); % dame for real parts of aeroelastic eigenvalues
modeNames = {DTU.modes.struc, DTU.modes.ael; redesign.modes.struc, redesign.modes.ael};
idx = cell(2);

for i = 1:size(cmb.camp,1)
    for j = 1:size(cmb.camp,2)
        if j == 1
            [cmb.camp{i,j},cmb.damp{i,j}] = deal(names{i});
        elseif j == 2
            [cmb.camp{i,j},cmb.damp{i,j},idx{i,j}] = import_cmb(n_modes,strcat("your_model/",turbine{i},file{i},type{j-1},".cmb"),...
                modeNames{i,j-1},true);
        elseif j == 3
            [cmb.camp{i,j},cmb.damp{i,j},idx{i,j},cmb.realPart{i,j}] = import_cmb(n_modes,strcat("your_model/",turbine{i},file{i},type{j-1},".cmb"),...
                modeNames{i,j-1},true);
        end
    end
end

fn = fieldnames(cmb);
for i = 1:length(fn)
    cmb.(fn{i}) = cell2table(cmb.(fn{i}));
    cmb.(fn{i}).Properties.VariableNames = {'Turbine','Structural','Aeroelastic'};
end

%% Plot Campbell diagrams
analyses = ["Structural", "Aeroelastic"];
markers  = ['o','+','*','x','|','s','d','^','v','>','<','p'];
leg      = cell(1,n_modes);

figure
counter = 1;
for i = 1:size(cmb.camp,1) % turbine loop
    for j = 2:size(cmb.camp,2) % structural or aeroelastic
        subplot(2,2,counter)
        for k = 2:n_modes+1
            plot(cmb.camp.(j){i}.("Wind Speed"),cmb.camp.(j){i}{:,k},'Marker',markers(k-1)); hold on
        end
        hold off
        title(strcat(names{i}," ",analyses{j-1}))
        if i == 2
            xlabel('Wind speed [m/s]')
        end
        if j == 2
            ylabel('Natural frequency [Hz]')
        end
        xlim([4 25])
        grid minor
        for l = 1:n_modes
            leg{l} = strcat(num2str(l),". ",cmb.camp.(j){i}.Properties.VariableNames{l+1});
        end
        legend(leg)
        counter = counter + 1;
    end
end
sgtitle('Campbell diagrams')

%% Plot damping diagrams
figure
counter = 1;
for i = 1:size(cmb.camp,1) % turbine loop
    for j = 2:size(cmb.camp,2) % structural or aeroelastic
        subplot(2,2,counter)
        for k = 2:n_modes+1
            plot(cmb.damp.(j){i}.("Wind Speed"),cmb.damp.(j){i}{:,k},'Marker',markers(k-1)); hold on
        end
        hold off
        title(strcat(names{i}," ",analyses{j-1}))
        if i == 2
            xlabel('Wind speed [m/s]')
        end
        if j == 2
            ylabel('Damping [% critical]')
        end
        xlim([4 25])
        grid minor
        for l = 1:n_modes
            leg{l} = strcat(num2str(l),". ",cmb.damp.(j){i}.Properties.VariableNames{l+1});
        end
        legend(leg)
        counter = counter + 1;
    end
end
sgtitle('Damping diagrams')