% 46320 LAC Course
% Assignment 1
% Authors: Nikolaos Stamatopoulos, Nils Joseph Gaukroger, Stefanos Masalas,
% Eduardo Fredrich de Miranda
% 18th September 2021

close all; clear variables; clc

%% Add functions folder to path
addpath('../../functions');

%% DTU rotor parameters
DTU.R           = 89.1660;           % rotor radius [m]
DTU.r_hub       = 2.8;               % hub radius [m]
DTU.bladeLength = DTU.R - DTU.r_hub; % blade length [m]
DTU.V_rated     = 11.4;              % rated wind speed [m/s]
DTU.TI          = 0.16;              % IEC turbulence intensity [-]

%% Design polynomials from DTU 10MW reports
dis1 = 101; % number of points in DTU.r discretisation
[DTU.c,DTU.that,DTU.beta,DTU.t,DTU.r] = DTU10MW_des(1,dis1);

%% New rotor radius
redesign.TI     = 0.14;              % IEC turbulence intensity [-]
[redesign.R,redesign.V_rated,~] = rotorScaling(DTU.V_rated,DTU.TI,DTU.R,redesign.TI);

% % adjustment for individual assignment
% rotorScale = 8; % manual rotor radius scaling [%]
% redesign.R = (1+rotorScale/100) * DTU.R;
% redesign.V_rated = ((DTU.V_rated^3 * DTU.R^2) / (redesign.R^2))^(1/3);

%% Aerofoil data
aerofoil = polars();

%% New design polynomials
x = [0.3, 0.35, 0.47, 0.5]; % cl,max - cl,des
% x = [0.3, 0.4, 0.35, 0.5]; % Edu's values
[p1,p2,cl_max,~,x_desi] = desPolys(aerofoil,x,4);

%% cl vs alpha
n = 4; % ignore 60.0% and cylinder
limits = [0, 90; 0, 2]; % [xlim, ylim]
plot_polars(aerofoil,limits,n,1,x_desi);

%% cl vs cd
limits = [0, 90]; % limits for alpha
plot_polars(aerofoil,limits,n,2,x_desi);

%% New absolute thickness
p = newThickness(DTU.t,DTU.r,DTU.R,redesign.R,1); % output as coeffs of fitted polynomial

%% Residual
redesign.B     = 3;   % number of blades [-]
redesign.a     = 1/3; % axial induction [-]
% tsr     = 5:0.1:10; % tsr(s) [-] (NB: minimum 5)
tsr_opt     = 7.1; % optimal tsr [-]
tsr         = tsr_opt;

% Constraints
t_max    = 5.38; % maximum absolute thickness [m]
c_max    = (redesign.R/DTU.R) * 6.17447367; % scaled max chord of DTU 10 MW [m]
beta_max = (redesign.R/DTU.R) * 20;   % maximum twist [deg]
% in future, maybe limit twist up to same non-dimensional radius as DTU
% 10MW

% Radial discretisation
dis2        = 200; % discretisation for spanwise discretisation [m]
redesign.r_hub = 2.8; % hub radius [m]
redesign.r = linspace(redesign.r_hub,redesign.R,dis2); % blade span [m]
redesign.bladeLength = redesign.R - redesign.r_hub;

% Preallocation
[redesign.t, result.c, result.phi, result.alpha, result.beta,...
    result.cl, result.cd, result.ap, result.cp, result.ct, that_tsr]...
    = deal(NaN(length(tsr),length(redesign.r))); % spanwise values
[result.CP, result.CT] = deal(NaN(length(tsr),1)); % global values

% Root transition start and initial guess for root transition end
lroot=redesign.R*0.03; %change value
that=redesign.t./result.c;
crtstart = -1;
for i = 1:length(redesign.r)
    if redesign.r(i)>lroot+redesign.r(1) && crtstart == -1
        crtstart=i; % root transition start
    end
    if redesign.r(i)>10*lroot+redesign.r(1)
        crtendstart=i; % lowest root transition end
        break
    end
end

% Least-squares parameters
lb = [0, 0]; % lower bounds
ub = [inf, 1]; % upper bounds
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
    for i = 1:length(redesign.r)
        x0 = lsqnonlin(@(x)residuals(x,redesign,tsr(j),i,p,p1,p2,t_max),x0,lb,ub,opts);
        [~,values] = residuals(x0,redesign,tsr(j),i,p,p1,p2,t_max);
        redesign.t(j,i)     = values(1);
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
    result.CP(j,1) = (2/redesign.R^2) * trapz(redesign.r(2:end),redesign.r(2:end).*result.cp(j,2:end));
%     result.CT(j,1) % can't remember the equation for this right now

    % Fixing geometry
    lroot = redesign.R*0.03; %change value
    that_tsr(j) = redesign.t(j)./result.c(j);
    for i = 1:length(redesign.r)
        if redesign.r(i) > lroot + redesign.r(1)
            crtstart = i;
            break;
        end
        that_tsr(j,i) = 1;
        result.c(j,i) = redesign.t(j,1);
    end
    
    % Smooth chord transition from cylinder
    crtend = crtendstart;
    limitchord = true;
    while limitchord
        cslope = (result.c(j,crtend+1)-result.c(j,crtend))/(redesign.r(crtend+1)-redesign.r(crtend));
        result.c(j,crtstart:crtend) = spline([redesign.r(crtstart) redesign.r(crtend)], [0 [result.c(j,crtstart-1) result.c(j,crtend)] cslope], redesign.r(crtstart:crtend));
        if max(result.c(j,crtstart:crtend)) < c_max
            limitchord = false;
            crtend = crtendstart;
        else
            crtend = crtend + 10;
            if crtend > length(result.c(j,:))
                crtend = length(result.c(j,:));
                break;
            end
            %crtend = crtendlast;
        end
    end
    for i = crtstart:crtend
        that_tsr(j,i) = redesign.t(j,i) / result.c(j,i);
    end
end

%% Apply general constraints

% cap twist at beta_max
result.beta(result.beta > deg2rad(beta_max)) = deg2rad(beta_max);


% set lower limit on relative thickness of 24.1%
result.c(result.c > (redesign.t / 0.241)) =...
    redesign.t(result.c > (redesign.t / 0.241)) / 0.241;

% remove flick after r/R > 0.98 for twist
r_R = 0.98;
% result.beta = flattenTip(result.beta,result.beta,rotor,r_R);

% remove flick after r/R > 0.98 for relative thickness
for j = 1:length(tsr)
    t_c = redesign.t(j,:) ./ result.c(j,:);
    a = t_c(find((redesign.r/redesign.R) > r_R,1));
    result.c(j,(redesign.r/redesign.R) > r_R) = redesign.t(j,(redesign.r/redesign.R) > r_R) ./ a;
end

%% Residuals for a and a'
% Least-squares parameters
lb = [0, 0];     % lower bounds
ub = [1, 1];   % upper bounds
opts = optimset('display','off'); % suppress lsqnonlin messages

if length(tsr) == 1
    fprintf('Rerunning residuals for TSR = %.2f\n',tsr)
else
    fprintf('Rerunning residuals for %d TSRs between %.2f and %.2f\n',length(tsr),tsr(1),tsr(end))
end
for j = 1:length(tsr)
    if length(tsr) ~= 1
        fprintf('TSR = %.1f\n',tsr(j))
    end
    x0 = [1/3, 0.001]; % initial guess
    for i = 1:length(redesign.r)
        x0 = lsqnonlin(@(x)residuals_induction(x,redesign,tsr(j),i,p,p1,p2,t_max,redesign.R,result.c(j,i)),x0,lb,ub,opts);
        [~,values] = residuals_induction(x0,redesign,tsr(j),i,p,p1,p2,t_max,redesign.R,result.c(j,i));
        redesign.t(j,i)     = values(1);
        result2.a(j,i)     = values(2);
        result2.phi(j,i)   = values(3);
        result2.alpha(j,i) = values(4);
        result2.beta(j,i)  = values(5);
        result2.cl(j,i)    = values(6);
        result2.cd(j,i)    = values(7);
        result2.ap(j,i)    = values(8);
        result2.cp(j,i)    = values(9);
        result2.ct(j,i)    = values(10);
    end
    result2.CP(j,1) = (2/redesign.R^2) * trapz(redesign.r(2:end),redesign.r(2:end).*result2.cp(j,2:end));
%     result.CT(j,1) % can't remember the equation for this right now
end

%% Reapply constraint for twist
result2.beta(result2.beta > deg2rad(beta_max)) = deg2rad(beta_max);

%% Spline for beta
new_beta = betaSpline(redesign,result2,[0.2, 0.265, 0.35],[0.8, 1.7]);

figure
subplot(3,1,2)
plot(redesign.r/redesign.R,result2.beta); hold on
plot(redesign.r/redesign.R,new_beta); hold off
ylabel('Twist, \beta [deg]'); xlabel('Non-dimensional radius [-]')
legend('Original','with spline')
grid on
box on

result2.beta = new_beta;

%% Plot geometry
figure
subplot(3,1,1)
plot(redesign.r/redesign.R,result.c); hold on
plot(DTU.r/DTU.R,fnval(DTU.c,DTU.r)); hold off
ylabel('Chord [m]');
legend('Redesign','DTU 10MW RWT')
grid on
subplot(3,1,2)
plot(redesign.r/redesign.R,rad2deg(result2.beta)); hold on
plot(DTU.r/DTU.R,fnval(DTU.beta,DTU.r)); hold off
ylabel('Twist, \beta [deg]');
legend('Redesign','DTU 10MW RWT')
grid on
subplot(3,1,3)
plot(redesign.r/redesign.R,(redesign.t./result.c)*100); hold on
plot(DTU.r/DTU.R,fnval(DTU.that,DTU.r)); hold off
ylabel('t/c [%]'); xlabel('Non-dimensional radius [-]')
legend('Redesign','DTU 10MW RWT')
grid on

if length(tsr) ~= 1
    figure
    plot(tsr,result2.CP)
    xlabel('TSR [-]'); ylabel('C_P [-]')
    grid on
    % Optimal TSR
    [CPmax,CPmax_idx] = max(result2.CP);
    tsr_opt = tsr(CPmax_idx);
end

%% Interpolate geometry for HAWC input
HAWC_in.n = 27;

% change from rotor radius to blade length
HAWC_in.r = linspace(0,redesign.R-redesign.r_hub,HAWC_in.n);
% chord, twist, thickness
HAWC_in.c = interp1(redesign.r,result.c,(HAWC_in.r+redesign.r_hub));
HAWC_in.beta = interp1(redesign.r,result2.beta,(HAWC_in.r+redesign.r_hub));
HAWC_in.t = interp1(redesign.r,redesign.t,(HAWC_in.r+redesign.r_hub));

figure
subplot(3,1,1)
plot(redesign.r,result.c); hold on
plot(HAWC_in.r+redesign.r_hub,HAWC_in.c); hold off
ylabel('Chord [m]');
legend('Original','HAWC')
grid on
subplot(3,1,2)
plot(redesign.r,rad2deg(result2.beta)); hold on
plot(HAWC_in.r+redesign.r_hub,rad2deg(HAWC_in.beta)); hold off
ylabel('Twist, \beta [deg]');
legend('Original','HAWC')
grid on
subplot(3,1,3)
plot(redesign.r,(redesign.t./result.c)*100); hold on
plot(HAWC_in.r+redesign.r_hub,(HAWC_in.t./HAWC_in.c)*100); hold off
ylabel('t/c [%]'); xlabel('Non-dimensional radius [-]')
legend('Original','HAWC')
grid on
sgtitle('HAWC geometry');

% Make column vectors for ease of use
HAWC_in.r = HAWC_in.r';
HAWC_in.c = HAWC_in.c';
HAWC_in.beta = HAWC_in.beta';
HAWC_in.that = (HAWC_in.t'./HAWC_in.c);
HAWC_in.tsr_opt = tsr_opt;

redesign.tsr_opt = tsr_opt;

%% Save variables for post-processing
save('../../mat/aero_design','aerofoil','DTU','HAWC_in','result','result2','redesign','p','p1','p2','t_max');