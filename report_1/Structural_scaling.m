% 46320 LAC Course
% Assignment 2
% Authors: Nikolaos Stamatopoulos, Nils Joseph Gaukroger, Stefanos Masalas,
% Eduardo Fredrich de Miranda
% 30th September 2021

close all; clear variables; clc

%% Add functions folder to path
addpath('functions');

%% Name input and output files
% st.dat file
original_fileName  = 'your_model/dtu10mw/DTU_10MW_RWT_Blade_st.dat';
new_fileName       = 'your_model/data/redesign_Blade_st.dat';

% c2def.txt file
original_c2def_fileName = 'your_model/dtu10mw/DTU_10MW_c2def.txt';
new_c2def_fileName = 'your_model/redesign_c2def.txt';

%% Run structuralScaling
[DTU,redesign] = structuralScaling(original_fileName,new_fileName,...
    original_c2def_fileName,new_c2def_fileName);

%% Save data for post-processing
save('struct.mat','DTU','redesign');