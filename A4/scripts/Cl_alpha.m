% 46320 LAC Course
% Assignment 1
% Authors: Nikolaos Stamatopoulos, Nils Joseph Gaukroger, Stefanos Masalas,
% Eduardo Fredrich de Miranda
% 18th September 2021

close all; clear variables; clc

%% Add functions folder to path
addpath('../../functions');

load('../../mat/aero_design.mat')

span = 79.7;
r = span + redesign.r_hub;

t_c = interp1(redesign.r,(redesign.t./result.c)*100,r);

[p1,p2,cl_max,alpha_max,x_desi] = desPolys(aerofoil,0.3,1);

figure
plot(redesign.r/redesign.R,(redesign.t./result.c)*100,'Linewidth',2)
line([0 r/redesign.R r/redesign.R],[t_c t_c 0],'LineStyle','--','Color','k','LineWidth',1)
xlabel('Non-dimensional radius [-]'); ylabel('t/c [%]')
% xline(r/redesign.R,'--')
grid minor

%% Aerofoil polar
n = 1; % ignore 60.0% and cylinder
limits = [0, 90; 0, 2]; % [xlim, ylim]
plot_polars(aerofoil,limits,n,1,x_desi);
line([0 alpha_max alpha_max],[cl_max cl_max 0],'LineStyle','--','Color','k','LineWidth',1)
% yline(cl_max,'--','label','C_{l,max}','LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom')
ax = gca();