% 46320 LAC Course
% Assignment 2: Post-processing
% Authors: Nikolaos Stamatopoulos, Nils Joseph Gaukroger, Stefanos Masalas,
% Eduardo Fredrich de Miranda
% 1st October 2021

close all; clear variables; clc

%% Set default plot settings
disp('Setting default plot parameters...');
set(    0,          'defaulttextInterpreter', 'tex');
set(groot, 'defaultAxesTickLabelInterpreter', 'tex');
set(groot,        'defaultLegendInterpreter', 'tex');
set(    0,             'defaultAxesFontSize',      15);
set(    0,            'DefaultLineLineWidth',       1.5);
disp('Default plot parameters set.');

%% Add functions folder to path
addpath('functions\')

%% Load data from Structural_scaling
load('struct.mat');

%% Create new variables
elas = {'st_flex','st_rigid'}; % filenames of different elasticities
% Array containing new variable and variables from which it is derived
% NB: order as such 'EI_xx' = f('E','I_x') etc.
vars = {'EI_xx','E','I_x';
    'EI_yy','E','I_y';
    'GJ','G','I_p';
    'EA','E','A'};

for i = 1:2
    for j = 1:length(elas)
        for k = 1:size(vars,1)
            if i == 1
                DTU.(elas{j}).(vars{k,1}) =...
                    DTU.(elas{j}).(vars{k,2}) .* DTU(i).(elas{j}).(vars{k,3});
            elseif i == 2
                redesign.(elas{j}).(vars{k,1}) =...
                    redesign.(elas{j}).(vars{k,2}) .* redesign.(elas{j}).(vars{k,3});
            end
        end
        % create non-dimensional span column
        if i == 1
            DTU.(elas{j}).s_S = DTU.(elas{j}).r / DTU.bladeLength;
        elseif i == 2
            redesign.(elas{j}).s_S = redesign.(elas{j}).r / redesign.bladeLength;
        end
    end
end

%% Stiffness distributions against blade span
ylabels = {'Flexural rigidity, EI_{xx} [Nm^2]','Flexural rigidity, EI_{yy} [Nm^2]',...
    'Torsional rigidity, GJ [Nm^2]','Axial rigidity, EA [N]'};
figure
for i = 1:4
    subplot(2,2,i)
    semilogy(DTU.st_flex.s_S,DTU.st_flex.(vars{i,1})); hold on
    plot(redesign.st_flex.s_S,redesign.st_flex.(vars{i,1})); hold off
    xlabel('Non-dimensionalised blade span [-]'); ylabel(ylabels{i})
    legend('DTU 10 MW','Redesign','Location','SW')
    grid minor
end

%% HAWC post-processing
turbine = {'results_dtu10mw/struct/','results_redesign/struct/'};
elas    = {'flex','rigid'};
file    = {'/DTU_10MW_struct_','/redesign_struct_'};
pwr     = cell(2,2); % array for storing .pwr files, rows = turbine, cols = elasticity

for i = 1:size(pwr,1)
    for j = 1:size(pwr,2)
        pwr{i,j} = import_pwr(strcat("your_model/",turbine{i},elas{j},file{i},elas{j},".pwr"));
    end
end

%% Calculation of deflection
for i = 1:size(pwr,1)
    pwr{i,1}.x_defl = pwr{i,1}.('Tip x') - pwr{i,2}.('Tip x');
    pwr{i,1}.y_defl = pwr{i,1}.('Tip y') - pwr{i,2}.('Tip y');
    pwr{i,1}.z_defl = pwr{i,1}.('Tip z') - pwr{i,2}.('Tip z');
end

%% Plotting
vars = {'P','T','Speed','Pitch','x_defl','y_defl','z_defl'};
ylabels = {'Power [kW]','Thrust [kN]','Rotor speed [rpm]','Pitch angle [deg]',...
    'Flapwise tip deflection [m]','Edgewise tip deflection [m]','Spanwise tip deflection [m]'};

% Operational parameters
figure
for i = 1:4
    subplot(2,2,i)
    plot(pwr{1,1}.V,pwr{1,1}.(vars{i})); hold on
    plot(pwr{2,1}.V,pwr{2,1}.(vars{i})); hold off
    ylabel(ylabels(i));
    legend('DTU 10MW','Redesign','Location','best')
    xlim([4 25])
    if i > 2
        xlabel('Wind speed [m/s]')
    end
    grid minor
end

% Tip deflections
figure
for i = 1:3
    subplot(3,1,i)
    plot(pwr{1,1}.V,pwr{1,1}.(vars{i+4})); hold on
    plot(pwr{2,1}.V,pwr{2,1}.(vars{i+4})); hold off
    ylabel(ylabels(i+4));
    legend('DTU 10MW','Redesign','Location','best')
    xlim([4 25])
    if i > 2
        xlabel('Wind speed [m/s]')
    end
    grid minor
end