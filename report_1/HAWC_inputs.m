% 46320 LAC Course
% Assignment 1
% Nils Joseph Gaukroger
% 30th September 2021

close all; clear variables; clc

%% Options
elasticity = 'rigid'; %'flexible' or 'rigid' operation in htc file 4
mode = '1pt'; % '1pt', 'tsr', or 'ws'
ctsr = 1; % chosen tsr index for blade geometry
auto = 0;

%% Load inputs
load('aero_design.mat');
DTU.c2def = readtable('your_model/dtu10mw/DTU_10MW_c2def.txt');
DTU.c2def.Properties.VariableNames = {'sec','idx','x','y','z','beta'};

%% Scale coordinates
scale = rotor.R / DTU.R;

HAWC_in.c2def   = DTU.c2def;
HAWC_in.c2def.x = DTU.c2def.x * scale;
HAWC_in.c2def.y = DTU.c2def.y * scale;
HAWC_in.c2def.z = HAWC_in.r;
HAWC_in.c2def.beta = -rad2deg(HAWC_in.beta);

%% Check for max rpm
TS_max    = 90; % maximum tip speed [m/s]
omega_max = (TS_max / rotor.R);
N_max     = omega_max * (60/(2*pi));
U0_max    = (omega_max * rotor.R) / HAWC_in.tsr_opt;

if U0_max < rotor.V_rated
    disp('Warning: maximum tip speed exceeded within partially loaded region')
end

%% Single point operation
U0 = 8; % wind speed [m/s]
omega = (HAWC_in.tsr_opt * U0) / rotor.R;
N     = omega * (60/(2*pi));
pitch = 0;

% Write data to file
fileID = fopen('your_model/data/operation_rigid_1pt.dat','w');
fprintf(fileID,'%3s %17s %21s %21s\n',num2str(length(U0)),'wind speed [m/s]','pitch [deg]','rot. speed [rpm]');
fprintf(fileID,' %.6f  %.6f  %.6f\n',U0,pitch,N);
fclose(fileID);

%% Multiple TSR operation
n_tsr = 10;
tsr = linspace(5,10,n_tsr);

U0 = 8:0.001:(8+(n_tsr-1)*0.001); % wind speed [m/s]
omega = (tsr .* U0) / rotor.R;
N     = omega * (60/(2*pi));
pitch = zeros(1,n_tsr);

% Write data to file
fileID = fopen('your_model/data/operation_rigid_tsr.dat','w');
fprintf(fileID,'%3s %17s %21s %21s\n',num2str(length(U0)),'wind speed [m/s]','pitch [deg]','rot. speed [rpm]');
for i = 1:n_tsr
    fprintf(fileID,' %.6f  %.6f  %.6f\n',U0(i),pitch(i),N(i));
end
fclose(fileID);

%% Multiple wind speed operation
n_U0 = 10;
U0 = linspace(4,rotor.V_rated,n_U0);

omega = (HAWC_in.tsr_opt .* U0) / rotor.R;
N     = omega * (60/(2*pi));
pitch = zeros(1,n_tsr);

% Write data to file
fileID = fopen('your_model/data/operation_rigid_ws.dat','w');
fprintf(fileID,'%3s %17s %21s %21s\n',num2str(length(U0)),'wind speed [m/s]','pitch [deg]','rot. speed [rpm]');
for i = 1:n_tsr
    fprintf(fileID,' %.6f  %.6f  %.6f\n',U0(i),pitch(i),N(i));
end
fclose(fileID);

%% Write new .ae file
fileID = fopen('your_model/data/redesign_ae.dat','w');
fprintf(fileID,'1\n');
fprintf(fileID,'1 %d\n',length(HAWC_in.r));
for i = 1:length(HAWC_in.r)
    fprintf(fileID,' %.14f  %.14f  %.14f 1 ;\n',HAWC_in.r(i),HAWC_in.c(i),HAWC_in.that(i)*100);
end
fclose(fileID);

%% Print for .htc file
if auto == 1
    nsec = 27;
    n = linspace(1,nsec,nsec);
    z = (n-1)/(nsec-1).*(rotor.R-rotor.r(1)); %z position along blade

    text = fileread(['your_model/dtu10mw/DTU_10MW_' elasticity '_hawc2s.htc']); %load .htc to text

    %find indices for start and end of c2_def section to replace
    nsecline=sprintf('nsec	%d;',nsec);
    B=strfind(text, nsecline)+length(nsecline);
    text(B:B+length(nsecline)-1);
    countline = 0;
    for i=B+2:length(text)
        if text(i) == newline
            countline = countline+1;
        end
        if countline == nsec
            Bf = 1 + i;
            break
        end
    end

    A=sscanf(text(B+2:Bf),'        sec	%i %f %f %f %f;\n', [5 nsec])'; %read values to matrix
    prebent=A(:,1:3); %gets numbering and x and y coordinates
    prebent(:,2:3) = prebent(:,2:3)*(HAWC_in.r(end)/(DTU.R-2.8)); %reescaling prebend with blade length
    printhtcmat = [prebent HAWC_in.r HAWC_in.c2def.beta]; %matrix to be printed in the new htc
    text(B+2:Bf)=[]; %delete old values
    text=[text(1:B+1) sprintf('        sec	%d %f %f %f %f;\n',printhtcmat') text(B+2:end)]; %insert new values
    text=strrep(text,'DTU_10MW_RWT_ae.dat','redesign_ae.dat'); %changes ae file to the redesigned
    % text=strrep(text,'DTU_10MW_RWT_Blade_st.dat','redesign_Blade_st.dat'); %changes st file to the redesigned
    text=strrep(text,'operation_rigid.dat',['operation_' elasticity '_' mode '.dat']); %changes operational data file to the redesigned

    newhtcfileID = fopen(['your_model/redesign_aero_' elasticity '.htc'],'w'); %creates new htc file
    fprintf(newhtcfileID,'%c',text); %writes to new htc file

    fclose('all'); %saves and closes all files
else
    % Alternative printing
    fprintf('nsec %d;\n',length(HAWC_in.r))
    for i = 1:length(HAWC_in.r)
        fprintf('sec %d %.8f %.8f %.8f %.8f;\n',i,HAWC_in.c2def.x(i),HAWC_in.c2def.y(i),HAWC_in.c2def.z(i),HAWC_in.c2def.beta(i));
    end
end

%% Max gen speed for structural
gearratio = 50;
genspeed = N(end) * gearratio;
fprintf('genspeed 0 %.2f; RPM (0 %.1f)\n', genspeed, N(end));