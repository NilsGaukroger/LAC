% 46320 LAC Course
% Assignment 2
% Stefanos Masalas & Nils Joseph Gaukroger
% 30th September 2021

close all; clear variables; clc

%% Name input and output files
% st.dat file
original_fileName  = 'your_model/dtu10mw/DTU_10MW_RWT_Blade_st.dat';
new_fileName       = 'your_model/data/redesign_Blade_st.dat';

% c2def.txt file
original_c2def_fileName = 'your_model/dtu10mw/DTU_10MW_c2def.txt';
new_c2def_fileName = 'your_model/redesign_c2def.txt';

%% Load aero_design outputs
load('aero_design.mat');

%% Import original DTU 10MW structural file
% Define import options
opts = detectImportOptions(original_fileName);
opts.VariableNames = {'r','m','x_cg','y_cg','ri_x','ri_y','x_sh','y_sh',...
    'E','G','I_x','I_y','I_p','k_x','k_y','A','pitch','x_e','y_e'}; % set variable names
opts.DataLines     = [6 56]; % only read flexible data (ignoring headers)
% in some places the delimiter is spaces not tabs, so readtable should recognise both
opts.Delimiter     = {'\t',' '};

% Import flexible data as table
f_original_st_flex  = readtable(original_fileName,opts);

% Import rigid data as table
opts.DataLines      = [59 109]; % only read rigid data (ignoring headers)
f_original_st_rigid = readtable(original_fileName,opts);

%% Set scaling factor
s = rotor.bladeLength / DTU.bladeLength;

%% Preallocate new (scaled) structural tables
f_new_st_flex  = f_original_st_flex;
f_new_st_rigid = f_original_st_rigid;

%% Scale variables by power of s and save in tables
vars{1} = {'r','x_cg','y_cg','ri_x','ri_y','x_sh','y_sh','x_e','y_e'}; % s^1
vars{2} = {'m','A'};                                                   % s^2
% (NB 'm' = mass per unit length)
vars{4} = {'I_x','I_y','I_p'};                                         % s^4

for i = 1:length(vars) % for each s^i (each cell in vars)
    for j = 1:length(vars{i}) % for each variable with vars{i}
        % Variable in scaled table = old variable * scaling factor ^ i
        f_new_st_flex.(vars{i}{j}) = f_original_st_flex.(vars{i}{j}) * s^i;
        f_new_st_rigid.(vars{i}{j}) = f_original_st_rigid.(vars{i}{j}) * s^i;
    end
end

%% Copy flexible headings from original structural file and paste in new
data   = importdata(original_fileName);  % import original file
fileID = fopen(new_fileName,'w');        % open new file with write permissions
fprintf(fileID,'%s\n',data.textdata{:}); % write original headings to file
fclose(fileID);                          % close file

%% Append new flexible data to file
writetable(f_new_st_flex,new_fileName,'WriteMode','append','Delimiter','tab');

%% Copy rigid headings from original structural file and paste in new
% Extract original headings
linenum        = 57;                           % line number where headings start
fileID         = fopen(original_fileName,'r'); % open original file with read permissions
rigid_headings = textscan(fileID,'%s',2,'delimiter','\n','headerlines',linenum-1);
% read 2 lines starting from linenum and save to rigid_headings
fclose(fileID);                                % close file

% Append headings to new file
fileID         = fopen(new_fileName,'a');    % open new file with append permissions
fprintf(fileID,'%s\n',rigid_headings{1}{:}); % write original headings to file
fclose(fileID);                              % close file

%% Append new flexible data to file
writetable(f_new_st_rigid,new_fileName,'WriteMode','append','Delimiter','tab');

%% Scale DTU 10 MW c2def coordinates
% Define import options
opts_c2def    = detectImportOptions(original_c2def_fileName);
opts_c2def.VariableNames = {'sec','idx','x','y','z','twist'}; % set variable names

f_original_c2 = readtable(original_c2def_fileName,opts_c2def); % read original data as table
f_new_c2      = f_original_c2; % preallocate new table

c2_vars = {'x','y','z'}; % specify variables to be scaled
for i = 1:length(c2_vars)
    f_new_c2.(c2_vars{i}) = f_original_c2.(c2_vars{i}) * s; % scale and save in new table
end

% Extract twist from aero_design output
for i = 1:length(HAWC_in.beta)
    f_new_c2.twist{i} = append(string(-rad2deg(HAWC_in.beta(i))),';');
    % NB: convert from radians, change sign, and append semicolon
end

% Write new data to file
writetable(f_new_c2,new_c2def_fileName,'Delimiter','tab','WriteVariableNames',0);

%% Save data for post-processing
save('struct.mat');