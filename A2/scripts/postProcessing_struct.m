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
set(    0,             'defaultAxesFontSize',    15);
set(    0,            'DefaultLineLineWidth',   1.5);
disp('Default plot parameters set.');

%% Add functions folder to path
addpath('..\..\functions\')

%% Figure saving settings
save_var = true; % true: saves figures, false: doesn't
overleaf = 'C:\Users\nilsg\Dropbox\Apps\Overleaf\LAC Assignment 2\figures\struct\'; % overleaf not working currently
local    = '../figs/v1/struct/';
locs     = {local};
for i = 1:length(locs) % if any directory doesn't exist don't attempt to save there
    if (not(isfolder(locs{i})))
        save_var = false;
    end
end

%% Load data from Structural_scaling
load('..\..\mat\Structural_scaling.mat');

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
    semilogy(redesign.st_flex.s_S,redesign.st_flex.(vars{i,1})); hold off
    if i > 2
        xlabel('Non-dimensionalised blade span [-]')
    end
    ylabel(ylabels{i})
    legend('DTU 10 MW','Redesign','Location','SW')
    grid minor
end
if save_var
    saveFigasPDF(locs,'stiffness');
end

figure
for i = 1:3
    semilogy(redesign.st_flex.s_S,redesign.st_flex.(vars{i,1})); hold on
%     semilogy(DTU.st_flex.s_S,DTU.st_flex.(vars{i,1}),'--')
end
hold off
xlabel('Non-dimensionalised blade span [-]')
ylabel('Rigidity [Nm^2]')
legend('Flapwise','Edgewise','Torsional')
grid minor
if save_var
    saveFigasPDF(locs,'stiffness_comp');
end

%% HAWC post-processing
turbine = {'dtu10mw/struct/','redesign_v1/struct/'};
names   = {'DTU 10MW','Redesign'};
elas    = {'flex','rigid'};
file    = {'/DTU_10MW_struct_','/redesign_struct_'};
pwr     = cell(2,3); % array for storing .pwr files, rows = turbine, cols = elasticity (+ turbine name)

for i = 1:size(pwr,1)
    for j = 1:size(pwr,2)
        if j == 1
            pwr{i,j} = names{i};
        else
            pwr{i,j} = import_pwr(strcat("../res/",turbine{i},elas{j-1},file{i},elas{j-1},".pwr"));
        end
    end
end

pwr = cell2table(pwr); % convert to table
pwr.Properties.VariableNames = {'Turbine',elas{1},elas{2}};

%% Calculation of deflection
for i = 1:size(pwr,1)
    pwr.flex{i}.x_defl = pwr.flex{i}.('Tip x') - pwr.rigid{i}.('Tip x');
    pwr.flex{i}.y_defl = pwr.flex{i}.('Tip y') - pwr.rigid{i}.('Tip y');
    pwr.flex{i}.z_defl = pwr.flex{i}.('Tip z') - pwr.rigid{i}.('Tip z');
end

%% Plotting
vars = {'P','T','Speed','Pitch','Tip y','x_defl','z_defl'};
ylabels = {'Power [kW]','Thrust [kN]','Rotor speed [rpm]','Pitch angle [deg]',...
    {'Flapwise'; 'tip position [m]'},{'Edgewise'; 'tip deflection [m]'},...
    {'Spanwise'; 'tip deflection [m]'}};

% Operational parameters
figure
for i = 1:4
    subplot(2,2,i)
    plot(pwr.flex{1}.V,pwr.flex{1}.(vars{i}),'marker','o','MarkerSize',4); hold on
    plot(pwr.flex{2}.V,pwr.flex{2}.(vars{i}),'marker','o','MarkerSize',4); hold off
    ylabel(ylabels(i));
    legend('DTU 10MW','Redesign','Location','best')
    xlim([4 25])
    if i == 1 || 3
        legend('DTU 10MW','Redesign','Location','SE')
    end
    if i > 2
        xlabel('Wind speed [m/s]')
    end
    grid minor
end
if save_var
    saveFigasPDF(locs,'operationalParameters');
end

% Tip deflections
figure
for i = 1:3
    subplot(3,1,i)
    plot(pwr.flex{1}.V,pwr.flex{1}.(vars{i+4}),'marker','o','MarkerSize',4); hold on
    plot(pwr.flex{2}.V,pwr.flex{2}.(vars{i+4}),'marker','o','MarkerSize',4); hold off
    ylabel(ylabels{i+4});
    legend('DTU 10MW','Redesign','Location','NE');
    xlim([4 25])
    if i == 1
        h = yline(0,'--','Planar Rotor','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','LineWidth',1.5,'color','k','FontSize',12);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    if i > 2
        xlabel('Wind speed [m/s]')
    end
    grid minor
end
if save_var
    saveFigasPDF(locs,'deflections');
end

%% Save outputs
save('..\..\mat\postProcessing_struct.mat','DTU','redesign','pwr');