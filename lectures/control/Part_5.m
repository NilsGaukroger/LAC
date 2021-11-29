% =========================================================================
%           Mahmood Mirzaei                                               %
%           mmir@dtu.dk                                                   %
%           Department of Wind Energy                                     %
%           Technical University of Denmark                               %
%           25-08-2016
% =========================================================================
clc; fclose all; clear all; warning off; profile on;beep off;
global SimParams Controller;

%% -------------------------------------------------------------------------
%                          Initialization
% -------------------------------------------------------------------------

SimParams.Simulation_TEND   = 400;          % Simulation lenght in seconds

WT.Model = 'WT0';                           % Model complexity,  WT0:Rotor,   WT1:Rotor+DT,  WT2:Rotor+DT+Tower fore-aft

Controller.Type             = 'PI';%'gs-PI';
SimParams.Open_Loop_Pitch   = 1;            	% pitch angle in degrees
SimParams.Open_Gen_Torque   = 1;               % generator reaction torque in N.m


% Controller.PiParams.KII     = 13013100;  % NOTE: student in the lecture need to play around the value
Controller.PiParams.KII     = 13013100*2;


%%

% -------------------------------------------------------------------------
%                   Choose the wind speed profile here:
% -------------------------------------------------------------------------


% Here you can choose the type of wind speed you'd like to use!
wind_profile_options = 8;   
% 1: for step wind speed, use WSP_Profile_Generator to produce wind steps!
% 2: for stochastic wind speed, mean wind speed: 8  m/s
% 3: for stochastic wind speed, mean wind speed: 12 m/s
% 4: for stochastic wind speed, mean wind speed: 15 m/s
% 5: for stochastic wind speed, mean wind speed: 18 m/s
% 6: for stochastic wind speed, mean wind speed: 15 m/s ,no wind shear.
% 7: EOG
% 8: for step wind speed below rated, use WSP_Profile_Generator to produce wind steps!
%%

main_script;
