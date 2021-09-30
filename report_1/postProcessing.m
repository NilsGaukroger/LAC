% 46320 LAC Course
% Assignment 1: Post-processing
% Nils Joseph Gaukroger
% 26th September 2021

close all; clear variables; clc

%% Load MATLAB results
load('aero_design.mat');

%% Load HAWC results
% Retrieve file names
DTU_folder = 'results_dtu10mw\rigid\';
DTU_name   = 'DTU_10MW_rigid_hawc2s_';
tsr_folder = 'results_redesign\tsr\';
ws_folder  = 'results_redesign\ws\';
rd_name    = 'redesign_rigid_hawc2s_';
DTU_files_tsr.u = dir([DTU_folder, DTU_name, 'u*']);
files_tsr.defl = dir([tsr_folder, rd_name, 'd*']);
files_tsr.fext = dir([tsr_folder, rd_name, 'f*']);
files_tsr.u    = dir([tsr_folder, rd_name, 'u*']);
files_tsr.pwr  = dir([tsr_folder, '\*.pwr']);
files_ws.defl = dir([ws_folder, rd_name, 'd*']);
files_ws.fext = dir([ws_folder, rd_name, 'f*']);
files_ws.u    = dir([ws_folder, rd_name, 'u*']);
files_ws.pwr  = dir([ws_folder, '\*.pwr']);

DTU_N_tsr_u = {DTU_files_tsr.u.name};
N_tsr_u = {files_tsr.u.name};
N_tsr_pwr = {files_tsr.pwr.name};
N_ws_u = {files_ws.u.name};
N_ws_pwr = {files_ws.pwr.name};

% Parse results as tables
u_vars = ["s","A","AP","PHI0","ALPHA0","U0","FX0","FY0","M0","UX0","UY0","UZ0",...
    "Twist","X_AC0","Y_AC0","Z_AC0","CL0","CD0","CM0","CLp0","CDp0","CMp0",...
    "F0","F'","CL_FS0","CLFS'","V_a","V_t","Tors.","vx","vy","chord",...
    "CT","CP","angle","v_1","v_2","v_3"];
HAWC_out.tsr_u = readtable([tsr_folder, N_tsr_u{4}],'FileType','text');
HAWC_out.tsr_u.Properties.VariableNames = u_vars;

pwr_vars = ["V","P","T","Cp","Ct","Pitch Q","Flap M","Edge M","Pitch","Speed",...
    "Tip x","Tip y","Tip z","J_rot","J_DT"];
HAWC_out.tsr_pwr = readtable([tsr_folder, N_tsr_pwr{1}],'FileType','text');
HAWC_out.tsr_pwr.Properties.VariableNames = pwr_vars;

HAWC_out.ws_u = readtable([ws_folder, N_ws_u{4}],'FileType','text');
HAWC_out.ws_u.Properties.VariableNames = u_vars;

HAWC_out.ws_pwr = readtable([ws_folder, N_ws_pwr{1}],'FileType','text');
HAWC_out.ws_pwr.Properties.VariableNames = pwr_vars;

DTU.tsr_u = readtable([DTU_folder, DTU_N_tsr_u{4}],'FileType','text');
DTU.tsr_u.Properties.VariableNames = u_vars;

%% Calculated values
HAWC_out.tsr_pwr.tsr = ((HAWC_out.tsr_pwr.Speed * ((2*pi)/60)) * rotor.R)...
    ./ HAWC_out.tsr_pwr.V;
chord = interp1(HAWC_in.r,HAWC_in.c,HAWC_out.tsr_u.s);
chord(end) = 1e-2;
HAWC_out.tsr_u.that = thickness(HAWC_out.tsr_u.s,p,t_max,rotor.R) ./ chord;
% HAWC_out.tsr_u.that(HAWC_out.tsr_u.that > 1) = 1;

% r_R = 0.98;
% a = HAWC_out.tsr_u.that(find((HAWC_out.tsr_u.s/rotor.bladeLength) > r_R,1));
% HAWC_out.tsr_u.that((HAWC_out.tsr_u.s/rotor.bladeLength) > r_R) = a;

%% Part 2: Required plots
% Side-by-side plots of the power and thrust coefficients at design pitch
% vs. TSR

figure
subplot(2,1,1)
plot(HAWC_out.tsr_pwr.tsr,HAWC_out.tsr_pwr.Cp);
xlabel('TSR [-]'); ylabel('C_P [-]')
grid minor
subplot(2,1,2)
plot(HAWC_out.tsr_pwr.tsr,HAWC_out.tsr_pwr.Ct);
xlabel('TSR [-]'); ylabel('C_T [-]')
grid minor

% Side-by-side plots of the actual lift coefficient and the design lift
% coefficient versus relative thickness (left plot) and versus radius
% (right plot) for design pitch and TSR.

figure
subplot(2,1,1)
plot(HAWC_out.tsr_u.that(1:end-9)*100,x_des(HAWC_out.tsr_u.that(1:end-9)*100,p1(1,:),p2(1,:))); hold on
plot(HAWC_out.tsr_u.that(1:end-9)*100,HAWC_out.tsr_u.CL0(1:end-9),'color',[0.8500, 0.3250, 0.0980]); hold off
xlabel('Relative thickness, t/c [%]'); ylabel('c_l [-]');
xlim([24.1 100]);
legend('Design','Actual','Location','best')
grid minor
subplot(2,1,2)
plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),x_des(HAWC_out.tsr_u.that*100,p1(1,:),p2(1,:))); hold on
plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),HAWC_out.tsr_u.CL0,'color',[0.8500, 0.3250, 0.0980]); hold off
xlabel('Non-dimensional radius [-]'); ylabel('c_l [-]');
legend('Design','Actual','Location','best')
grid minor

figure
subplot(2,1,1)
plot(HAWC_out.tsr_u.that(1:end-9)*100,x_des(HAWC_out.tsr_u.that(1:end-9)*100,p1(2,:),p2(2,:))); hold on
plot(HAWC_out.tsr_u.that(1:end-9)*100,rad2deg(HAWC_out.tsr_u.ALPHA0(1:end-9)),'color',[0.8500, 0.3250, 0.0980]); hold off
xlabel('Relative thickness, t/c [%]'); ylabel('\alpha [deg]');
xlim([24.1 100]);
legend('Design','Actual','Location','best')
grid minor
subplot(2,1,2)
plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),x_des(HAWC_out.tsr_u.that*100,p1(2,:),p2(2,:))); hold on
plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),rad2deg(HAWC_out.tsr_u.ALPHA0),'color',[0.8500, 0.3250, 0.0980]); hold off
xlabel('Non-dimensional radius [-]'); ylabel('\alpha [deg]');
legend('Design','Actual','Location','best')
grid minor

figure
subplot(2,1,1)
plot(HAWC_out.tsr_u.that(1:end-9)*100,x_des(HAWC_out.tsr_u.that(1:end-9)*100,p1(3,:),p2(3,:))); hold on
plot(HAWC_out.tsr_u.that(1:end-9)*100,HAWC_out.tsr_u.CL0(1:end-9)./HAWC_out.tsr_u.CD0(1:end-9),'color',[0.8500, 0.3250, 0.0980]); hold off
xlabel('Relative thickness, t/c [%]'); ylabel('c_l/c_d [-]');
xlim([24.1 100]);
legend('Design','Actual','Location','best')
grid minor
subplot(2,1,2)
plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),x_des(HAWC_out.tsr_u.that*100,p1(3,:),p2(3,:))); hold on
plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),HAWC_out.tsr_u.CL0./HAWC_out.tsr_u.CD0,'color',[0.8500, 0.3250, 0.0980]); hold off
xlabel('Non-dimensional radius [-]'); ylabel('c_l/c_d [-]');
legend('Design','Actual','Location','best')
grid minor

% Plots of the actual lift coefficient, lift-drag ratio, AoA, axial
% induction (a), local CT, and local CP vs. radius at design TSR.

figure
subplot(3,2,1)
plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),HAWC_out.tsr_u.CL0); hold on
plot(DTU.tsr_u.s/DTU.tsr_u.s(end),DTU.tsr_u.CL0); hold off
legend('Redesign','DTU 10MW','Location','best');
% xlabel('Non-dimensional blade length [-]');
ylabel('c_l [-]')
grid minor
subplot(3,2,2)
plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),HAWC_out.tsr_u.CL0./HAWC_out.tsr_u.CD0); hold on
plot(DTU.tsr_u.s/DTU.tsr_u.s(end),DTU.tsr_u.CL0./DTU.tsr_u.CD0); hold off
% legend('Redesign','DTU 10MW','Location','best');
% xlabel('Non-dimensional blade length [-]');
ylabel('c_l/c_d [-]')
grid minor
subplot(3,2,3)
plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),rad2deg(HAWC_out.tsr_u.ALPHA0)); hold on
plot(DTU.tsr_u.s/DTU.tsr_u.s(end),rad2deg(DTU.tsr_u.ALPHA0)); hold off
% legend('Redesign','DTU 10MW','Location','best');
% xlabel('Non-dimensional blade length [-]');
ylabel('\alpha [deg]')
grid minor
subplot(3,2,4)
plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),HAWC_out.tsr_u.A); hold on
plot(DTU.tsr_u.s/DTU.tsr_u.s(end),DTU.tsr_u.A)
yline(1/3,'--'); hold off
% legend('Redesign','DTU 10MW','Location','best');
% xlabel('Non-dimensional blade length [-]');
ylabel('Axial induction, a [-]')
grid minor
subplot(3,2,5)
plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),HAWC_out.tsr_u.CP); hold on
plot(DTU.tsr_u.s/DTU.tsr_u.s(end),DTU.tsr_u.CP)
yline(16/27,'--'); hold off
xlabel('Non-dimensional blade length [-]'); ylabel('C_{P,loc} [-]')
% legend('Redesign','DTU 10MW','Location','best');
grid minor
subplot(3,2,6)
plot(HAWC_out.tsr_u.s/HAWC_out.tsr_u.s(end),HAWC_out.tsr_u.CT); hold on
plot(DTU.tsr_u.s/DTU.tsr_u.s(end),DTU.tsr_u.CT); hold off
xlabel('Non-dimensional blade length [-]'); ylabel('C_{T,loc} [-]')
% legend('Redesign','DTU 10MW','Location','best');
grid minor

% Side-by-side plots of the operational rotor speed and pitch angles vs.
% wind speed

figure
subplot(2,1,1)
plot(HAWC_out.ws_pwr.V,HAWC_out.ws_pwr.Speed);
% xlabel('Wind speed [m/s]');
xlim([4 10.7196])
ylabel('Rotor speed [rpm]')
grid minor
subplot(2,1,2)
plot(HAWC_out.ws_pwr.V,HAWC_out.ws_pwr.Pitch);
xlim([4 10.7196])
xlabel('Wind speed [m/s]'); ylabel('Pitch angle [deg]')
grid minor

% Side-by-side plots of the aerodynamic power and its coefficient and the
% thrust and its coefficient vs. wind speed

figure
subplot(2,1,1)
plot(HAWC_out.ws_pwr.V,HAWC_out.ws_pwr.P/1e3); hold on
ylabel('Aerodynamic power [MW]')
yyaxis right
yline(0.4621,'--','MATLAB','LabelHorizontalAlignment','left')
plot(HAWC_out.ws_pwr.V,HAWC_out.ws_pwr.Cp); hold off
% xlabel('Wind speed [m/s]');
ylabel('C_P [-]')
ylim([0 (16/27)])
xlim([4 10.7196])
grid minor
subplot(2,1,2)
plot(HAWC_out.ws_pwr.V,HAWC_out.ws_pwr.T); hold on
ylabel('Thrust [kN]')
yyaxis right
plot(HAWC_out.ws_pwr.V,HAWC_out.ws_pwr.Ct); hold off
xlabel('Wind speed [m/s]'); ylabel('C_T [-]')
ylim([0 1])
xlim([4 10.7196])
grid minor