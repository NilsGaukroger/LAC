% 46320 LAC Course
% Assignment 4: Post-processing
% Authors: Nikolaos Stamatopoulos, Nils Joseph Gaukroger, Stefanos Masalas,
% Eduardo Fredrich de Miranda
% 16th November 2021

close all; clear variables; clc

%% Add functions folder to path
addpath('..\..\functions\')

%% Stats directories
stat_dir(1) = "../res/dtu10mw/turb/"; % stats file, turbine 1
stat_dir(2) = "../res/redesign/turb/";  % stats file, turbine 2
labels = ["DTU 10 MW, TC A", "Redesign, TC B"]; % legend labels

%% Values and indices
wind_idxs = [15,15];
%                   Label              idx_t1 idx_t2
plot_vals = {"Pitch angle [deg]",        4,    4;
             "Rotor speed [rad/s]",      10,   10;
             "Thrust [kN]",              13,   13;
             "AoA @ 2/3R m",             NaN,  NaN;
             "Cl @ 2/3R m",              NaN,  NaN;
             "Tower-base FA [kNm]",      19,   19;
             "Tower-base SS [kNm]",      20,   20;
             "Yaw-bearing pitch [kNm]",  22,   22;        
             "Yaw-bearing roll [kNm]",   23,   23;
             "Shaft torsion [kNm]",      27,   27;
             "OoP BRM [kNm]",            28,   28;
             "IP BRM [kNm]",             29,   29;
             "Generator torque [Nm]",    72,   72;
             "Electrical power [W]",     102,  102; 
             "Tower clearance [m]",      110,  110;
             };
plot_vals = cell2table(plot_vals,'VariableNames',["Label","idx_t1","idx_t2"]);

%% Import stats files
stat = ["mean","max","min"];
for j = 1:2
    for i = 1:length(stat)
        stat_file = strcat(stat_dir(j),'stats_', stat(i), '.txt');
        if j == 1
            if i == 1
                [~, DTU.channels, DTU.mean] = load_stats(stat_file,j);
            elseif i == 2
                [~, ~, DTU.max] = load_stats(stat_file,j);
            elseif i == 3
                [~, ~, DTU.min] = load_stats(stat_file,j);
            end
        elseif j == 2
            if i == 1
                [~, redesign.channels, redesign.mean] = load_stats(stat_file,j);
            elseif i == 2
                [~, ~, redesign.max] = load_stats(stat_file,j);
            elseif i == 3
                [~, ~, redesign.min] = load_stats(stat_file,j);
            end
        end
    end
end

%% Plot
for i = 1:length(DTU.channels)
    if DTU.channels(i) == wind_idxs(2)
        continue
    else
        figure
        plot(reshape(DTU.mean(:,:,4),[],1),reshape(DTU.mean(:,:,i),[],1),'ok','MarkerSize',3)
        hold on
        plot(reshape(redesign.mean(:,:,4),[],1),reshape(redesign.mean(:,:,i),[],1),'or','MarkerSize',3)
%         plot(reshape(DTU.mean(:,:,4),[],1),reshape(DTU.max(:,:,i),[],1),'xk','MarkerSize',5)
%         plot(reshape(DTU.mean(:,:,4),[],1),reshape(DTU.min(:,:,i),[],1),'xk','MarkerSize',5)
        hold off
        xlabel('Wind speed [m/s]')
        ylabel(plot_vals.Label(find(plot_vals.idx_t1 == DTU.channels(i))))
        xlim([4 25])
        legend('DTU','Redesign','Location','best')
        grid on
        box on
    end
end