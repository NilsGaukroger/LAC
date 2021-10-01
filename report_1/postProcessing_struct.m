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

% %% HAWC post-processing
% DTU_folder = 'your_model/results_dtu10mw/flexible/';
% DTU_pwr_fileName = 'DTU_10MW_flexible_hawc2s.pwr';
% redesign_folder = 'your_model/results_redesign/struct/flex/';
% redesign_pwr_fileName = 'redesign_struct_flex.pwr';
% 
% DTU_pwr = readtable([DTU_folder DTU_pwr_fileName],'FileType','text');
% redesign_pwr = readtable([redesign_folder redesign_pwr_fileName],'FileType','text');