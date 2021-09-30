% 46320 LAC Course
% Assignment 1
% Nils Joseph Gaukroger
% 18th September 2021

close all; clear variables; clc

%% Design polynomials from DTU 10MW reports
dis = 101; % number of points in DTU.r discretisation
[DTU.c,DTU.that,DTU.beta,DTU.t,DTU.r] = DTU10MW_des(1,dis);

%% New rotor radius
DTU.R = 89.1660;
[rotor.R,rotor.V_rated,~] = rotorScaling(11.4,0.16,DTU.R,0.14);

%% Aerofoil data
aerofoil = polars();

%% New design polynomials
x = [0.3, 0.35, 0.47, 0.5]; % cl,max - cl,des
[p1,p2,cl_max,x_desi] = desPolys(aerofoil,x,4);

%% cl vs alpha
n = 4; % ignore 60.0% and cylinder
limits = [0, 90; 0, 2]; % [xlim, ylim]
plot_polars(aerofoil,limits,n,1,x_desi);

%% cl vs cd
limits = [0, 90]; % limits for alpha
plot_polars(aerofoil,limits,n,2,x_desi);

%% New absolute thickness
p = newThickness(DTU.t,DTU.r,DTU.R,rotor.R,1); % output as coeffs of fitted polynomial
 
%% Create geometries by minimising residuals
rotor.B = 3;   % number of blades [-]
rotor.a = 1/3; % axial induction [-]
% tsr     = 5:1:10; % tsr(s) [-] (NB: minimum 5)
tsr     = 6.61; % optimal tsr [-]

% Constraints
t_max    = 5.38; % maximum absolute thickness [m]
c_max    = (rotor.R/DTU.R) * 6.17447367; % scaled max chord of DTU 10 MW [m]
beta_max = (rotor.R/DTU.R) * 14.5;   % maximum twist [deg]
% in future, maybe limit twist up to same non-dimensional radius as DTU
% 10MW

% Radial discretisation
spacing = 0.2; % increment for radial discretisation [m]
rotor.r_hub = 2.8; % hub radius [m]
rotor.r = linspace(rotor.r_hub,rotor.R,475); % blade span [m]
rotor.bladeLength = rotor.R - rotor.r_hub;

% Preallocation
[rotor.t, result.c, result.phi, result.alpha, result.beta,...
    result.cl, result.cd, result.ap, result.cp, result.ct]...
    = deal(NaN(length(tsr),length(rotor.r))); % spanwise values
[result.CP, result.CT] = deal(NaN(length(tsr),1)); % global values

% Least-squares parameters
lb = [0, 0];     % lower bounds
ub = [inf, 1];   % upper bounds
opts = optimset('display','off'); % suppress lsqnonlin messages

if length(tsr) == 1
    fprintf('Creating design for TSR = %.2f\n',tsr)
else
    fprintf('Creating designs for %d TSRs between %.2f and %.2f\n',length(tsr),tsr(1),tsr(end))
end
for j = 1:length(tsr)
    if length(tsr) ~= 1
        fprintf('TSR = %.1f\n',tsr(j))
    end
    x0 = [3, 0.001]; % initial guess
    for i = 1:length(rotor.r)
        x0 = lsqnonlin(@(x)residuals(x,rotor,tsr(j),i,p,p1,p2,t_max),x0,lb,ub,opts);
        [~,values] = residuals(x0,rotor,tsr(j),i,p,p1,p2,t_max);
        rotor.t(j,i)     = values(1);
        result.c(j,i)     = values(2);
        result.phi(j,i)   = values(3);
        result.alpha(j,i) = values(4);
        result.beta(j,i)  = values(5);
        result.cl(j,i)    = values(6);
        result.cd(j,i)    = values(7);
        result.ap(j,i)    = values(8);
        result.cp(j,i)    = values(9);
        result.ct(j,i)    = values(10);
    end
    result.CP(j,1) = (2/rotor.R^2) * trapz(rotor.r(2:end),rotor.r(2:end).*result.cp(j,2:end));
%     result_tsr.CT(j,1) % can't remember the equation for this right now
end

%% Apply general constraints
% cap chord at c_max
result.c(result.c > c_max) = c_max;
result.c(1:27) = result.c(1); %----->>> SOS added that to make it a bit flat at the beginning

% cap twist at beta_max
result.beta(result.beta > deg2rad(beta_max)) = deg2rad(beta_max);

% set lower limit on relative thickness of 24.1%
result.c(result.c > (rotor.t / 0.241)) =...
    rotor.t(result.c > (rotor.t / 0.241)) / 0.241;

% remove flick after r/R > 0.98 for twist
r_R = 0.98;
result.beta = flattenTip(result.beta,result.beta,rotor,r_R);

% remove flick after r/R > 0.98 for relative thickness
t_c = rotor.t ./ result.c;
a = t_c(find((rotor.r/rotor.R) > r_R,1));
result.c((rotor.r/rotor.R) > r_R) = rotor.t((rotor.r/rotor.R) > r_R) ./ a;

% smooth chord transitions ------------->>> SOS 
splineX = [rotor.r(24) rotor.r(93) rotor.r(195)];
%splineY = linspace(result.c(4),result.c(185),19)
splineY = [result.c(24) result.c(93)-0.245 result.c(195)];                            
xq = rotor.r(24:195);
yy = spline(splineX,splineY,xq);
result.c(24:195) = yy;
if max(result.c) > c_max % Check max chord doesn't exceed c_max
    fprintf('Warning: Chord exceeds c_max, adjust splines')
end

%% Plot geometry
figure
subplot(3,1,1)
plot(rotor.r/rotor.R,result.c); hold on
plot(DTU.r/DTU.R,fnval(DTU.c,DTU.r),'x'); hold off
ylabel('Chord [m]');
legend('Redesign','DTU 10MW RWT')
grid on
subplot(3,1,2)
plot(rotor.r/rotor.R,rad2deg(result.beta)); hold on
plot(DTU.r/DTU.R,fnval(DTU.beta,DTU.r),'x'); hold off
ylabel('Twist, \beta [deg]');
legend('Redesign','DTU 10MW RWT')
grid on
subplot(3,1,3)
plot(rotor.r/rotor.R,(rotor.t./result.c)*100); hold on
plot(DTU.r/DTU.R,fnval(DTU.that,DTU.r),'x'); hold off
ylabel('t/c [%]'); xlabel('Non-dimensional radius [-]')
legend('Redesign','DTU 10MW RWT')
grid on

if length(tsr) ~= 1
    figure
    plot(tsr,result.CP)
    xlabel('TSR [-]'); ylabel('C_P [-]')
    grid on
    % Optimal TSR
    [CPmax,CPmax_idx] = max(result.CP);
    tsr_opt   = tsr(CPmax_idx);
end

%% Rerun residuals with correct axial induction factor


%% Plot geometry
figure
subplot(3,1,1)
plot(rotor.r/rotor.R,result.c); hold on
plot(DTU.r/DTU.R,fnval(DTU.c,DTU.r),'x'); hold off
ylabel('Chord [m]');
legend('Redesign','DTU 10MW RWT')
grid on
subplot(3,1,2)
plot(rotor.r/rotor.R,rad2deg(result.beta)); hold on
plot(DTU.r/DTU.R,fnval(DTU.beta,DTU.r),'x'); hold off
ylabel('Twist, \beta [deg]');
legend('Redesign','DTU 10MW RWT')
grid on
subplot(3,1,3)
plot(rotor.r/rotor.R,(rotor.t./result.c)*100); hold on
plot(DTU.r/DTU.R,fnval(DTU.that,DTU.r),'x'); hold off
ylabel('t/c [%]'); xlabel('Non-dimensional radius [-]')
legend('Redesign','DTU 10MW RWT')
grid on

if length(tsr) ~= 1
    figure
    plot(tsr,result.CP)
    xlabel('TSR [-]'); ylabel('C_P [-]')
    grid on
    % Optimal TSR
    [CPmax,CPmax_idx] = max(result.CP);
    tsr_opt   = tsr(CPmax_idx);
end

%% Interpolate geometry for HAWC input
HAWC_in.n = 27;

% change from rotor radius to blade length
HAWC_in.r = linspace(0,rotor.R-rotor.r_hub,HAWC_in.n);
% chord, twist, thickness
HAWC_in.c = interp1(rotor.r,result.c,(HAWC_in.r+rotor.r_hub));
HAWC_in.beta = interp1(rotor.r,result.beta,(HAWC_in.r+rotor.r_hub));
HAWC_in.t = interp1(rotor.r,rotor.t,(HAWC_in.r+rotor.r_hub));

figure
subplot(3,1,1)
plot(rotor.r,result.c); hold on
plot(HAWC_in.r+rotor.r_hub,HAWC_in.c); hold off
ylabel('Chord [m]');
legend('Original','HAWC')
grid on
subplot(3,1,2)
plot(rotor.r,result.beta); hold on
plot(HAWC_in.r+rotor.r_hub,HAWC_in.beta); hold off
ylabel('Twist, \beta [deg]');
legend('Original','HAWC')
grid on
subplot(3,1,3)
plot(rotor.r,(rotor.t./result.c)*100); hold on
plot(HAWC_in.r+rotor.r_hub,(HAWC_in.t./HAWC_in.c)*100); hold off
ylabel('t/c [%]'); xlabel('NRadius [-]')
legend('Original','HAWC')
grid on
sgtitle('HAWC geometry');

% Make column vectors for ease of use
HAWC_in.r = HAWC_in.r';
HAWC_in.c = HAWC_in.c';
HAWC_in.beta = HAWC_in.beta';
HAWC_in.that = (HAWC_in.t'./HAWC_in.c);

%% Save variables for post-processing
save('aero_design','aerofoil','DTU','HAWC_in','result','rotor','p','p1','p2','t_max');