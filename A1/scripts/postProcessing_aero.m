% 46320 LAC Course
% Assignment 1: Post-processing
% Authors: Nikolaos Stamatopoulos, Nils Joseph Gaukroger, Stefanos Masalas,
% Eduardo Fredrich de Miranda
% 26th September 2021

close all; clear variables; clc
addpath('..\..\functions\');

%% Set default plot settings
disp('Setting default plot parameters...');
set(    0,          'defaulttextInterpreter', 'tex');
set(groot, 'defaultAxesTickLabelInterpreter', 'tex');
set(groot,        'defaultLegendInterpreter', 'tex');
set(    0,             'defaultAxesFontSize',      15);
set(    0,            'DefaultLineLineWidth',       1.5);
disp('Default plot parameters set.');
% set(    0,           'DefaultFigureColormap',feval('turbo'));

%% Load MATLAB results
load('aero_design.mat');

%% Load HAWC results
turbine = {'dtu10mw/','redesign/'};
simul   = {'tsr_opt/','tsr/','ws/'};
file    = {'DTU_10MW_aero_rigid','redesign_aero_rigid'};

u   = cell(length(turbine),length(simul));
pwr = cell(length(turbine),length(simul));
idx = [1,4,4];

for i = 1:length(turbine)
    for j = 1:length(simul)
        filepath = strcat('../res/', turbine{i}, simul{j});
        u{i,j} = import_u(filepath,file{i},idx(j));
        if j > 1
            pwr{i,j} = import_pwr(strcat(filepath,file{i},'.pwr'));
        end
    end
end

%% Calculated values
% HAWC_out.tsr_pwr.tsr = ((HAWC_out.tsr_pwr.Speed * ((2*pi)/60)) * redesign.R)...
%     ./ HAWC_out.tsr_pwr.V;
% chord = interp1(HAWC_in.r,HAWC_in.c,HAWC_out.tsr_u.s);
% chord(end) = 1e-2;
% HAWC_out.tsr_u.that = thickness(HAWC_out.tsr_u.s,p,t_max,redesign.R) ./ chord;
% HAWC_out.tsr_u.that(HAWC_out.tsr_u.that > 1) = 1;

% r_R = 0.98;
% a = HAWC_out.tsr_u.that(find((HAWC_out.tsr_u.s/rotor.bladeLength) > r_R,1));
% HAWC_out.tsr_u.that((HAWC_out.tsr_u.s/rotor.bladeLength) > r_R) = a;
% 
% %% Part 2: Required plots
% % Side-by-side plots of the power and thrust coefficients at design pitch
% % vs. TSR
% 
% figure
% subplot(2,1,1)
% plot(HAWC_out.tsr_pwr.tsr,HAWC_out.tsr_pwr.Cp);
% xlabel('TSR [-]'); ylabel('C_P [-]')
% grid minor
% subplot(2,1,2)
% plot(HAWC_out.tsr_pwr.tsr,HAWC_out.tsr_pwr.Ct);
% xlabel('TSR [-]'); ylabel('C_T [-]')
% grid minor
% 
% % Side-by-side plots of the actual lift coefficient and the design lift
% % coefficient versus relative thickness (left plot) and versus radius
% % (right plot) for design pitch and TSR.
% 
% figure
% subplot(2,1,1)
% plot(HAWC_out.tsr_u.that(1:end-9)*100,x_des(HAWC_out.tsr_u.that(1:end-9)*100,p1(1,:),p2(1,:))); hold on
% plot(HAWC_out.tsr_u.that(1:end-9)*100,HAWC_out.tsr_u.CL0(1:end-9),'color',[0.8500, 0.3250, 0.0980]); hold off
% xlabel('Relative thickness, t/c [%]'); ylabel('c_l [-]');
% xlim([24.1 100]);
% legend('Design','Actual','Location','best')
% grid minor
% subplot(2,1,2)
% plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),x_des(HAWC_out.tsr_u.that*100,p1(1,:),p2(1,:))); hold on
% plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),HAWC_out.tsr_u.CL0,'color',[0.8500, 0.3250, 0.0980]); hold off
% xlabel('Non-dimensional radius [-]'); ylabel('c_l [-]');
% legend('Design','Actual','Location','best')
% grid minor
% 
% figure
% subplot(2,1,1)
% plot(HAWC_out.tsr_u.that(1:end-9)*100,x_des(HAWC_out.tsr_u.that(1:end-9)*100,p1(2,:),p2(2,:))); hold on
% plot(HAWC_out.tsr_u.that(1:end-9)*100,rad2deg(HAWC_out.tsr_u.ALPHA0(1:end-9)),'color',[0.8500, 0.3250, 0.0980]); hold off
% xlabel('Relative thickness, t/c [%]'); ylabel('\alpha [deg]');
% xlim([24.1 100]);
% legend('Design','Actual','Location','best')
% grid minor
% subplot(2,1,2)
% % plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),x_des(HAWC_out.tsr_u.that*100,p1(2,:),p2(2,:))); hold on
% plot(redesign.r/redesign.R,result2.alpha); hold on
% plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),rad2deg(HAWC_out.tsr_u.ALPHA0),'color',[0.8500, 0.3250, 0.0980]); hold off
% xlabel('Non-dimensional radius [-]'); ylabel('\alpha [deg]');
% legend('Design','Actual','Location','best')
% grid minor
% 
% figure
% subplot(2,1,1)
% plot(HAWC_out.tsr_u.that(1:end-9)*100,x_des(HAWC_out.tsr_u.that(1:end-9)*100,p1(3,:),p2(3,:))); hold on
% plot(HAWC_out.tsr_u.that(1:end-9)*100,HAWC_out.tsr_u.CL0(1:end-9)./HAWC_out.tsr_u.CD0(1:end-9),'color',[0.8500, 0.3250, 0.0980]); hold off
% xlabel('Relative thickness, t/c [%]'); ylabel('c_l/c_d [-]');
% xlim([24.1 100]);
% legend('Design','Actual','Location','best')
% grid minor
% subplot(2,1,2)
% plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),x_des(HAWC_out.tsr_u.that*100,p1(3,:),p2(3,:))); hold on
% plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),HAWC_out.tsr_u.CL0./HAWC_out.tsr_u.CD0,'color',[0.8500, 0.3250, 0.0980]); hold off
% xlabel('Non-dimensional radius [-]'); ylabel('c_l/c_d [-]');
% legend('Design','Actual','Location','best')
% grid minor
% 
% % Plots of the actual lift coefficient, lift-drag ratio, AoA, axial
% % induction (a), local CT, and local CP vs. radius at design TSR.
% 
% figure
% subplot(3,2,1)
% plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),HAWC_out.tsr_u.CL0); hold on
% plot(DTU.tsr_u.s/DTU.tsr_u.s(end),DTU.tsr_u.CL0); hold off
% legend('Redesign','DTU 10MW','Location','best');
% % xlabel('Non-dimensional blade length [-]');
% ylabel('c_l [-]')
% grid minor
% subplot(3,2,2)
% plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),HAWC_out.tsr_u.CL0./HAWC_out.tsr_u.CD0); hold on
% plot(DTU.tsr_u.s/DTU.tsr_u.s(end),DTU.tsr_u.CL0./DTU.tsr_u.CD0); hold off
% % legend('Redesign','DTU 10MW','Location','best');
% % xlabel('Non-dimensional blade length [-]');
% ylabel('c_l/c_d [-]')
% grid minor
% subplot(3,2,3)
% plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),rad2deg(HAWC_out.tsr_u.ALPHA0)); hold on
% plot(DTU.tsr_u.s/DTU.tsr_u.s(end),rad2deg(DTU.tsr_u.ALPHA0)); hold off
% % legend('Redesign','DTU 10MW','Location','best');
% % xlabel('Non-dimensional blade length [-]');
% ylabel('\alpha [deg]')
% grid minor
% subplot(3,2,4)
% plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),HAWC_out.tsr_u.A); hold on
% plot(DTU.tsr_u.s/DTU.tsr_u.s(end),DTU.tsr_u.A)
% yline(1/3,'--'); hold off
% % legend('Redesign','DTU 10MW','Location','best');
% % xlabel('Non-dimensional blade length [-]');
% ylabel('Axial induction, a [-]')
% grid minor
% subplot(3,2,5)
% plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),HAWC_out.tsr_u.CP); hold on
% plot(DTU.tsr_u.s/DTU.tsr_u.s(end),DTU.tsr_u.CP)
% yline(16/27,'--'); hold off
% xlabel('Non-dimensional blade length [-]'); ylabel('C_{P,loc} [-]')
% % legend('Redesign','DTU 10MW','Location','best');
% grid minor
% subplot(3,2,6)
% plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),HAWC_out.tsr_u.CT); hold on
% plot(DTU.tsr_u.s/DTU.tsr_u.s(end),DTU.tsr_u.CT); hold off
% xlabel('Non-dimensional blade length [-]'); ylabel('C_{T,loc} [-]')
% % legend('Redesign','DTU 10MW','Location','best');
% grid minor
% 
% % Side-by-side plots of the operational rotor speed and pitch angles vs.
% % wind speed
% 
% figure
% subplot(2,1,1)
% plot(HAWC_out.ws_pwr.V,HAWC_out.ws_pwr.Speed);
% % xlabel('Wind speed [m/s]');
% xlim([4 10.7196])
% ylabel('Rotor speed [rpm]')
% grid minor
% subplot(2,1,2)
% plot(HAWC_out.ws_pwr.V,HAWC_out.ws_pwr.Pitch);
% xlim([4 10.7196])
% xlabel('Wind speed [m/s]'); ylabel('Pitch angle [deg]')
% grid minor
% 
% % Side-by-side plots of the aerodynamic power and its coefficient and the
% % thrust and its coefficient vs. wind speed
% 
% figure
% subplot(2,1,1)
% plot(HAWC_out.ws_pwr.V,HAWC_out.ws_pwr.P/1e3); hold on
% ylabel('Aerodynamic power [MW]')
% yyaxis right
% yline(0.4621,'--','MATLAB','LabelHorizontalAlignment','left')
% plot(HAWC_out.ws_pwr.V,HAWC_out.ws_pwr.Cp); hold off
% % xlabel('Wind speed [m/s]');
% ylabel('C_P [-]')
% ylim([0 (16/27)])
% xlim([4 10.7196])
% grid minor
% subplot(2,1,2)
% plot(HAWC_out.ws_pwr.V,HAWC_out.ws_pwr.T); hold on
% ylabel('Thrust [kN]')
% yyaxis right
% plot(HAWC_out.ws_pwr.V,HAWC_out.ws_pwr.Ct); hold off
% xlabel('Wind speed [m/s]'); ylabel('C_T [-]')
% ylim([0 1])
% xlim([4 10.7196])
% grid minor
% 
% %% Compare CP and CT with DTU 10MW
% fprintf('Redesign CP = %.3f\n',result2.CP);
% fprintf('DTU 10 MW CP = 0.478\n');
% fprintf('Redesign max. thrust = %0.1f kN\n',HAWC_out.ws_pwr.T(end));
% fprintf('DTU 10 MW max. thrust = 1507.4 kN\n');
% 
% %% Max power for .htc operational_data block
% fprintf('maxpow %.1f\n', HAWC_out.ws_pwr.P(end));