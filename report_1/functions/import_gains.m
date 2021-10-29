function [K,I,Kp2,Ki2,Kp3,Ki3,gs_K1,gs_K2] = import_gains(filename,order)
%% Import data from .txt
data = importdata(filename);

%% Assign gains by region
% Region 1: Optimal CP tracking (T = K\omega^2)
K  = str2double(data.textdata{2,3});  % Optimal CP tracking K factor [Nm/(rad/s)^2]
% Region 2: PI generator torque control
I  = str2double(data.textdata{4,3});  % Combined rotor-generator inertia? [kg*m^2]
Kp2 = str2double(data.textdata{5,3});  % Proportional gain of torque controller [Nm/(rad/s)]
Ki2 = str2double(data.textdata{6,3});  % Integral gain of torque controller [Nm/rad]
% Region 3: PI blade pitch control (w/ gain scheduling)
Kp3 = str2double(data.textdata{8,3});  % Proportional gain of pitch controller [rad/(rad/s)]
Ki3 = str2double(data.textdata{9,3});  % Integral gain of pitch controller [rad/rad]
gs_K1 = str2double(data.textdata{10,3}); % Coefficient of linear term in aerodynamic gain scheduling [deg]
if order == 2 % if using a quadratic fit for the aerodynamic gain scheduling
    gs_K2 = str2double(data.textdata{10,7}); % Coefficient of quadratic term in aerodynamic gain scheduling [deg^2]
end

end