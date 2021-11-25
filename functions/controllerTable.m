function [gains,genTorque_switch] = controllerTable(full_load,operation)
%% Create table for gains (cases C1-C2)
caseName          = ["C1";"C2"];
generatorTorque   = ["Constant Power";"Constant Power"];
genTorque_switch  = [1;1];
partialLoad_omega = 0.05 * ones(2,1);
partialLoad_zeta  = 0.7  * ones(2,1);
fullLoad_omega    = [0.05; 0.055]; % natural frequency [Hz]
fullLoad_zeta     = [ 0.7;  0.65]; % damping ratio [-]

gains = table(caseName,generatorTorque,partialLoad_omega,partialLoad_zeta,fullLoad_omega,fullLoad_zeta);

%% Add details for C7
operationName = ["Constant Torque","Constant Power"];
gains = [gains;{"C3",operationName(operation+1), 0.05, 0.7, full_load(1), full_load(2)}];
genTorque_switch = [genTorque_switch; 1];

end