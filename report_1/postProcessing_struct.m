% 46320 LAC Course
% Assignment 2: Post-processing
% Authors: Nikolaos Stamatopoulos, Nils Joseph Gaukroger, Stefanos Masalas,
% Eduardo Fredrich de Miranda
% 1st October 2021

close all; clear variables; clc

%% Load data from Structural_scaling
load('struct.mat');

%% Create new variables
elas = {'st_flex','st_rigid'};
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
DTU_folder_flex = 'your_model/results_dtu10mw/struct/flex/';
DTU_pwr_fileName_flex = 'DTU_10MW_struct_flex.pwr';
redesign_folder_flex = 'your_model/results_redesign/struct/flex/';
redesign_pwr_fileName_flex = 'redesign_struct_flex.pwr';

DTU_pwr_rigid = readtable([DTU_folder_flex DTU_pwr_fileName_flex],'FileType','text');
DTU_pwr_rigid.Properties.VariableNames = {'V','P','T','Cp','Ct',...
    'Pitch Q','Flap M','Edge M','Pitch','Speed','Tip x','Tip y','Tip z',...
    'J_rot','J_DT'};
DTU_pwr_flex = readtable([DTU_folder_rigid DTU_pwr_fileName_rigid],'FileType','text');
DTU_pwr_flex.Properties.VariableNames = DTU_pwr_rigid.Properties.VariableNames;
redesign_pwr_rigid = readtable([redesign_folder_rigid redesign_pwr_fileName_rigid],'FileType','text');
redesign_pwr_flex = readtable([redesign_folder_flex redesign_pwr_fileName_flex],'FileType','text');
redesign_pwr_flex.Properties.VariableNames = DTU_pwr_flex.Properties.VariableNames;

%% 
vars = {'P','T','Speed','Pitch','Tip y','Tip x','Tip z'};
ylabels = {'Power [kW]','Thrust [kN]','Rotor speed [rpm]','Pitch angle [deg]',...
    'Flapwise tip deflection [m]','Edgewise tip deflection [m]','Tip-z [m]'};
figure
for i = 1:4
    subplot(2,2,i)
    plot(DTU_pwr_flex.V,DTU_pwr_flex.(vars{i})); hold on
    plot(redesign_pwr_flex.V,redesign_pwr_flex.(vars{i})); hold off
    ylabel(ylabels(i));
    legend('DTU 10MW','Redesign','Location','best')
    if i > 2
        xlabel('Wind speed [m/s]')
    end
    grid minor
end
figure
for i = 1:3
    subplot(3,1,i)
    plot(DTU_pwr_flex.V,DTU_pwr_flex.(vars{i+4})); hold on
    plot(redesign_pwr_flex.V,redesign_pwr_flex.(vars{i+4})); hold off
    ylabel(ylabels(i+4));
    legend('DTU 10MW','Redesign','Location','best')
    grid minor
end