% 46320 LAC Course
% Assignment 1
% Nils Joseph Gaukroger
% 18th September 2021

close all; clear variables; clc

%% Design polynomials from DTU 10MW reports
dis1 = 101; % number of points in DTU.r discretisation
[DTU.c,DTU.that,DTU.beta,DTU.t,DTU.r] = DTU10MW_des(1,dis1);

%% New rotor radius
DTU.R = 89.1660;
[rotor.R,rotor.V_rated,~] = rotorScaling(11.4,0.16,DTU.R,0.14);

%% Aerofoil data
aerofoil = polars();

%% New design polynomials
x = [0.3, 0.35, 0.47, 0.5]; % cl,max - cl,des
% x = [0.3, 0.4, 0.35, 0.5]; % Edu's values
[p1,p2,cl_max,x_desi] = desPolys(aerofoil,x,4);

%% cl vs alpha
n = 4; % ignore 60.0% and cylinder
limits = [0, 90; 0, 2]; % [xlim, ylim]
plot_polars(aerofoil,limits,n,1,x_desi);

%% cl vs cd
limits = [0, 90]; % limits for alpha
plot_polars(aerofoil,limits,n,2,x_desi);

%% Absolute thickness
figure
plot(DTU.r,DTU.t)
hold on
grid on

%% New absolute thickness
p = newThickness(DTU.t,DTU.r,DTU.R,rotor.R,1); % output as coeffs of fitted polynomial

%% Residual
rotor.B     = 3;   % number of blades [-]
rotor.a     = 1/3; % axial induction [-]
tsr     = 5:0.1:10; % tsr(s) [-] (NB: minimum 5)
% tsr     = 6.61; % optimal tsr [-]

% Constraints
t_max    = 5.38; % maximum absolute thickness [m]
c_max    = (rotor.R/DTU.R) * 6.17447367; % scaled max chord of DTU 10 MW [m]
beta_max = (rotor.R/DTU.R) * 20;   % maximum twist [deg]
% in future, maybe limit twist up to same non-dimensional radius as DTU
% 10MW

% Radial discretisation
dis2        = 200; % discretisation for spanwise discretisation [m]
rotor.r_hub = 2.8; % hub radius [m]
rotor.r = linspace(rotor.r_hub,rotor.R,dis2); % blade span [m]
rotor.bladeLength = rotor.R - rotor.r_hub;

% Preallocation
[rotor.t, result.c, result.phi, result.alpha, result.beta,...
    result.cl, result.cd, result.ap, result.cp, result.ct, that_tsr]...
    = deal(NaN(length(tsr),length(rotor.r))); % spanwise values
[result.CP, result.CT] = deal(NaN(length(tsr),1)); % global values

% Root transition start and initial guess for root transition end
lroot=rotor.R*0.03; %change value
that=rotor.t./result.c;
crtstart = -1;
for i = 1:length(rotor.r)
    if rotor.r(i)>lroot+rotor.r(1) && crtstart == -1
        crtstart=i; % root transition start
    end
    if rotor.r(i)>10*lroot+rotor.r(1)
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
%     result.CT(j,1) % can't remember the equation for this right now

    % Fixing geometry
    lroot = rotor.R*0.03; %change value
    that_tsr(j) = rotor.t(j)./result.c(j);
    for i = 1:length(rotor.r)
        if rotor.r(i) > lroot + rotor.r(1)
            crtstart = i;
            break;
        end
        that_tsr(j,i) = 1;
        result.c(j,i) = rotor.t(j,1);
    end
    
    % Smooth chord transition from cylinder
    crtend = crtendstart;
    limitchord = true;
    while limitchord
        cslope = (result.c(j,crtend+1)-result.c(j,crtend))/(rotor.r(crtend+1)-rotor.r(crtend));
        result.c(j,crtstart:crtend) = spline([rotor.r(crtstart) rotor.r(crtend)], [0 [result.c(j,crtstart-1) result.c(j,crtend)] cslope], rotor.r(crtstart:crtend));
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
        that_tsr(j,i) = rotor.t(j,i) / result.c(j,i);
    end
end

%% Apply general constraints

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

% % smooth chord transitions ------------->>> SOS 
% splineX = [rotor.r(24) rotor.r(93) rotor.r(195)];
% %splineY = linspace(result.c(4),result.c(185),19)
% splineY = [result.c(24) result.c(93)-0.245 result.c(195)];                            
% xq = rotor.r(24:195);
% yy = spline(splineX,splineY,xq);
% result.c(24:195) = yy;
% if max(result.c) > c_max % Check max chord doesn't exceed c_max
%     fprintf('Warning: Chord exceeds c_max, adjust splines')
% end

%% Residuals for a and a'
% Least-squares parameters
lb = [0, 0];     % lower bounds
ub = [1, 1];   % upper bounds
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
    x0 = [1/3, 0.001]; % initial guess
    for i = 1:length(rotor.r)
        x0 = lsqnonlin(@(x)residuals_induction(x,rotor,tsr(j),i,p,p1,p2,t_max,rotor.R,result.c(j,i)),x0,lb,ub,opts);
        [~,values] = residuals_induction(x0,rotor,tsr(j),i,p,p1,p2,t_max,rotor.R,result.c(j,i));
        rotor.t(j,i)     = values(1);
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
    result2.CP(j,1) = (2/rotor.R^2) * trapz(rotor.r(2:end),rotor.r(2:end).*result2.cp(j,2:end));
%     result.CT(j,1) % can't remember the equation for this right now
end

%% Reapply constraint for twist
result2.beta(result2.beta > deg2rad(beta_max)) = deg2rad(beta_max);

%% Plot geometry
figure
subplot(3,1,1)
plot(rotor.r,result.c)
ylabel('Chord [m]');
grid on
subplot(3,1,2)
plot(rotor.r,rad2deg(result2.beta))
ylabel('Twist, \beta [deg]');
grid on
subplot(3,1,3)
plot(rotor.r,(rotor.t./result.c)*100)
ylabel('t/c [%]'); xlabel('Radius [m]')
grid on

figure
plot(tsr,result2.CP)
xlabel('TSR [-]'); ylabel('C_P [-]')
grid on

%% Plot geometry
figure
subplot(3,1,1)
plot(rotor.r/rotor.R,result.c); hold on
plot(DTU.r/DTU.R,fnval(DTU.c,DTU.r),'x'); hold off
ylabel('Chord [m]');
legend('Redesign','DTU 10MW RWT')
grid on
subplot(3,1,2)
plot(rotor.r/rotor.R,rad2deg(result2.beta)); hold on
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
HAWC_in.r = linspace(0,rotor.R-rotor.r_hub,HAWC_in.n);
% chord, twist, thickness
HAWC_in.c = interp1(rotor.r,result.c,(HAWC_in.r+rotor.r_hub));
HAWC_in.beta = interp1(rotor.r,result2.beta,(HAWC_in.r+rotor.r_hub));
HAWC_in.t = interp1(rotor.r,rotor.t,(HAWC_in.r+rotor.r_hub));

figure
subplot(3,1,1)
plot(rotor.r,result.c); hold on
plot(HAWC_in.r+rotor.r_hub,HAWC_in.c); hold off
ylabel('Chord [m]');
legend('Original','HAWC')
grid on
subplot(3,1,2)
plot(rotor.r,result2.beta); hold on
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

%% Residual function for a and a'
function [out,varargout] = residuals_induction(x,rotor,tsr,idx,p,p1,p2,t_max,Rnew,c)
% Calculate the residuals for chord and tangential induction factor

% unpack inputs
R = Rnew; B = rotor.B;
% spanwise position
r = rotor.r(idx);

% unpack x
%c  = x(1);
a  = x(1);
ap = x(2);

% calculate that, cl and cl/cd from design polynomials
t    = thickness(r,p,t_max,rotor.R);
that = t ./ c;
cl   = x_des(that*100,p1(1,:),p2(1,:));
alpha = x_des(that*100,p1(2,:),p2(2,:));
clcd = x_des(that*100,p1(3,:),p2(3,:));

% calculate intermediate variables
phi   = atan(((1-a)/(1+ap)) * (R/(r*tsr)));
cd    = cl/clcd;
cy    = cl * cos(phi) + cd * sin(phi);
cx    = cl * sin(phi) - cd * cos(phi);
f     = (B/2) * ((R-r)/(r * sin(phi)));
F     = (2/pi) * acos(exp(-f));
sigma = (c*B) ./ (2*pi*r);

% calculate residuals
res_c  = 4*pi*r * sin(phi)^2 * F * ((2*a) / ...
    (cy * B * (1-a))) - c;

% res_a = 1/((4*F*sin(phi)*cos(phi))/(sigma*cx)-1) - a;
res_ap = 1 / ((4 * F * sin(phi) * cos(phi) / ...
    (sigma * cx)) - 1) - ap;

% pack output
out = [res_c, res_ap];

% optional outputs
if nargout == 2
    beta = phi - deg2rad(alpha);
    cp = ((1-a)^2 + (tsr*(r/R))^2 * (1+ap)^2) * tsr * (r/R) * sigma * cx;
    ct = ((1-a)^2 + (tsr*(r/R))^2 * (1+ap)^2) * sigma * cy;
%     result.a = 1 / (((4*F*np.sin(phi)**2) / (sigma * cy)) + 1)
    varargout{1} = [t,a,phi,alpha,beta,cl,cd,ap,cp,ct];
%     varargout{1} = [t,c,phi,alpha,beta,cl,cd,a,ap,cp,ct];
end
end