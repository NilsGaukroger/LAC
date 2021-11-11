% 46320 LAC Course
% Assignments
% Authors: Nikolaos Stamatopoulos, Nils Joseph Gaukroger, Stefanos Masalas,
% Eduardo Fredrich de Miranda
% 10th November 2021

close all; clear variables; clc

%% Run all scripts
% Assignment 1: Aerodynamic rotor design for IEC class IIIb
run('A1/scripts/aero_design')
% Assignment 1: Writing inputs for HAWC
run('A1/scripts/HAWC_inputs')
% Assignment 1: Aerodynamic post-processing
run('A1/scripts/postProcessing_aero')
% Assignment 2: Scaling of structural parameters
run('A2/scripts/Structural_scaling')
% Assignment 2: Structural post-processing
run('A2/scripts/postProcessing_struct')
% Assignment 2: Stability post-processing
run('A2/scripts/postProcessing_stab')
% Assignment 3: Controller design post-processing
run('A2/scripts/postProcessing_cont')