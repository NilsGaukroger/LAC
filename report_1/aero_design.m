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
[rotor.R,~] = rotorScaling(11.4,0.16,DTU.R,0.14);

%% Aerofoil data
aerofoil = polars();

%% cl vs alpha
n = 4; % ignore 60.0% and cylinder
limits = [0, 90; 0, 2]; % [xlim, ylim]
plot_polars(aerofoil,limits,n,1);

%% cl vs cd
limits = [0, 90]; % limits for alpha
plot_polars(aerofoil,limits,n,2);

%% New design polynomials
x = [0.3, 0.35, 0.47, 0.5]; % cl,max - cl,des
[p1,p2] = desPolys(aerofoil,x,4);

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
r_hub   = 2.8; % hub radius [m]
rotor.r = (r_hub:spacing:rotor.R-(spacing*2)); % blade span [m]
% rotor.r = (0:spacing:rotor.R-r_hub); % blade span [m]

% Preallocation
[rotor.t, result.c, result.phi, result.alpha, result.beta,...
    result.cl, result.cd, result.ap, result.cp, result.ct]...
    = deal(NaN(length(tsr),length(rotor.r))); % spanwise values
[result.CP, result.CT] = deal(NaN(length(tsr),1)); % global values

% Least-squares parameters
x0 = [3, 0.001]; % initial guess
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
        x0 = lsqnonlin(@(x)residuals(x,rotor,tsr(j),i,p,p1,p2,t_max,rotor.R),x0,lb,ub,opts);
        [~,values] = residuals(x0,rotor,tsr(j),i,p,p1,p2,t_max,rotor.R);
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

%% smooth chord transitions
[maxi,idx] = max(fnval(DTU.c,DTU.r));
r_R_max = DTU.r(idx)/DTU.R;
result.c((rotor.r/rotor.R) <= 0.1) = t_max;
point1 = rotor.r(find((rotor.r/rotor.R) < 0.1,1,'last'));
point2 = rotor.r(find((rotor.r/rotor.R) < ((r_R_max - 0.1)/2),1,'last'));
point3 = rotor.r(find((rotor.r/rotor.R) > r_R_max,1));

section1 = (rotor.r > point1) & (rotor.r < point3);
x = [point1, point2, point3];
y(1) = result.c(rotor.r == x(1));
y(2) = result.c(rotor.r == x(2));
y(3) = result.c(rotor.r == x(3));
xq = rotor.r(section1);
s = pchip(x,y,xq);

figure
plot(xq,s)

%% Edu's changes
% Fixing geometry 
%lroot=Rnew*0.03; %change value
%that=result.t./result.c;
%for i = 1:length(rotor.radii)
%    if rotor.radii(i)>lroot+rotor.radii(1)
%        crtstart=i;
%        break;
%     end
%     that(i)=1;
%     result.c(i)=result.t(1);
% end

% %Smooth chord transition from cylinder:

% for i = 1:length(rotor.radii)
%     if rotor.radii(i)>10*lroot+rotor.radii(1)
%         crtend=i;
%         break
%     end
% end
% cslope=(result.c(crtend+1)-result.c(crtend))/(rotor.radii(crtend+1)-rotor.radii(crtend));
% result.c(crtstart:crtend) = spline([rotor.radii(crtstart) rotor.radii(crtend)], [0 [result.c(crtstart-1) result.c(crtend)] cslope], rotor.radii(crtstart:crtend));

% for i = crtstart:crtend
%     that(i) = result.t(i) / result.c(i);
% end

%% Apply tweaks
% remove kink in twist


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

%% Design polynomial splines for the DTU 10MW
function [c,that,beta,t,x] = DTU10MW_des(plt,dis)
% Chord
c_breaks = [2.8000 8.1960 19.9548 28.0124 38.2220 55.0271 70.0576 78.1586 85.0000 86.2521 88.6595 88.9861 89.1660];
c_coeffs = [
    0.00000000E+00 0.00000000E+00 0.00000000E+00 5.38000000E+00;
    -4.61301443E-04 1.04228497E-02 0.00000000E+00 5.38000000E+00;
    9.54064544E-05 -5.85020456E-03 5.37688193E-02 6.07113924E+00;
    7.94309215E-05 -3.54396342E-03 -2.19256287E-02 6.17447367E+00;
    2.42467032E-05 -1.11108961E-03 -6.94518582E-02 5.66574440E+00;
    8.59761392E-06 1.11315203E-04 -8.62531671E-02 4.29988830E+00;
	4.71207167E-06 4.98994511E-04 -7.70799070E-02 3.05780226E+00;
	-4.39083118E-04 6.13511989E-04 -6.80674918E-02 2.46863016E+00;
    -1.32866948E-03 -8.39831774E-03 -1.21326462E-01 1.89106967E+00;
    -1.47600180E-02 -1.33891816E-02 -1.48606495E-01 1.72338280E+00;
	-6.79818589E+00 -1.19989176E-01 -4.69702133E-01 1.08209177E+00;
	1.04359293E+01 -6.78085172E+00 -2.72351677E+00 6.79055450E-01];
c = ppmak(c_breaks,c_coeffs,1);

% Relative thickness
that_breaks = [2.8 4.8 18.8310 27.1510 37.424 63.615 89.166];
that_coeffs = 100*[
    0.00000000E+00 0.00000000E+00 0.00000000E+00 1.00000000E+00
    1.63195479E-04 -4.64026965E-03 0.00000000E+00 1.00000000E+00
    -6.94374457E-05 2.22911763E-03 -3.38308740E-02 5.37264647E-01
    -1.19564154E-05 4.95958990E-04 -1.11582365E-02 3.70105515E-01
    -9.34760029E-07 1.27474224E-04 -4.75370706E-03 2.94855128E-01
    0.00000000E+00 0.00000000E+00 0.00000000E+00 2.41000000E-01
    ];
that = ppmak(that_breaks,that_coeffs,1);

% Twist
beta_breaks = [2.8 9.119 15.5543 21.3482 30.1888 43.4352 57.5373 89.1660];
beta_coeffs = [
	0.00000000E+00 0.00000000E+00 0.00000000E+00 1.45000000E+01
	-1.94790140E-03 -1.76633017E-02 0.00000000E+00 1.45000000E+01
	 5.39990815E-03 -5.52692913E-02 -4.69343116E-01 1.32493815E+01
	 -1.35567415E-03 3.85902922E-02 -5.65979568E-01 9.72497036E+00
	 -1.20591403E-04 2.63537363E-03 -2.01519947E-01 6.80074227E+00
	 9.58895526E-05 -2.15683227E-03 -1.95180996E-01 4.31345834E+00
	 -1.48680703E-05 1.89989991E-03 -1.98804282E-01 1.40098858E+00
	];
beta = ppmak(beta_breaks,beta_coeffs,1);

% Plot polynomials
if plt == 1
    poly = [c,that,beta];
    x = linspace(2.8,89.166,dis);
    labels = ["Chord [m]", "Relative thickness [%]", "Twist [deg]"];
    figure
    for i = 1:length(poly)
        subplot(3,1,i)
        plot(x,fnval(poly(i),x))
        ylabel(labels(i));
        grid on
    end
    xlabel('Radius [m]')
    sgtitle('DTU 10MW Geometry')
else
end

% Absolute thickness
t = fnval(c,x) .* (fnval(that,x)/100);
end

%% Aerofoil data input
function aerofoil = polars()
S = dir('polars\*.txt');
S = S(1:6,1); % remove 'hints.txt'
aerofoil = cell(2,length(S));
aerofoil(1,:) = {S.name};
vars = {'alpha','cl','cd','fsst','clinv','clfs','cm'};

for i = 1:size(aerofoil,2)
    aerofoil{2,i} = readtable("polars\" + aerofoil{1,i});
    aerofoil{2,i}.Properties.VariableNames = vars;
    aerofoil{1,i} = erase(aerofoil{1,i},"_ds.txt");
end
end

%% Polar plotting
function plot_polars(aerofoil,limits,n,mode)
figure
for i = 1:n
    subplot(ceil(n/2),2,i)
    if mode == 1
        plot(aerofoil{2,i}.alpha,aerofoil{2,i}.cl)
        xlabel('\alpha [deg]'); ylabel('c_l [-]');
        ylim(limits(2,:))
        xlim(limits(1,:));
    elseif mode == 2
        mask = (aerofoil{2,1}.alpha >= 0 & aerofoil{2,1}.alpha <= 90);
        plot(aerofoil{2,i}.cd(mask),aerofoil{2,i}.cl(mask))
        xlabel('c_d [-]'); ylabel('c_l [-]');
    end
    title(aerofoil{1,i});
    grid on
end
figure
markers = ['o','x','^','s'];
for i = 1:n
    if mode == 1
        plot(aerofoil{2,i}.alpha,aerofoil{2,i}.cl,'Marker',markers(i))
        hold on
    elseif mode == 2
        plot(aerofoil{2,i}.cd(mask),aerofoil{2,i}.cl(mask),'Marker',markers(i));
        hold on
    end
end
if mode == 1
    xlabel('\alpha [deg]'); ylabel('c_l [-]');
    xlim(limits(1,:)); ylim(limits(2,:));
elseif mode == 2
    xlabel('c_d [-]'); ylabel('c_l [-]');
end
legend(aerofoil{1,1:n});
grid on
hold off
end

%% Creating design polynomials
function [p1,p2] = desPolys(aerofoil,x1,n)
% Find values
mask = (aerofoil{2,1}.alpha >= 0 & aerofoil{2,1}.alpha <= 90);
idx0 = find(aerofoil{2,1}.alpha == 0);

t_c = [24.1, 30.1, 36.0, 48.0, 100];

mask2 = false(size(aerofoil{2,1},1),4);

cl_des = NaN(1,5);
alpha_des = NaN(1,5);
cd_des = NaN(1,5);
clcd_des = NaN(1,5);
for i = 1:n
    clmax_idx = find(diff(aerofoil{2,i}.cl(mask)) <= 0, 1, 'First');
    clmax = aerofoil{2,i}.cl(clmax_idx + idx0 - 1);
    alpha_max = aerofoil{2,i}.alpha(clmax_idx + idx0 - 1);
    cl_des(i) = clmax - x1(i);
    mask2(:,i) = aerofoil{2,i}.alpha >= 0 & aerofoil{2,i}.alpha < alpha_max;
    alpha_des(i) = interp1(aerofoil{2,i}.cl(mask2(:,i)),aerofoil{2,i}.alpha(mask2(:,i)),cl_des(i));
    cd_des(i) = interp1(aerofoil{2,i}.cl(mask2(:,i)),aerofoil{2,i}.cd(mask2(:,i)),cl_des(i));
    clcd_des(i) = cl_des(i) / cd_des(i);
end
cl_des(1,5) = 1e-6;
alpha_des(1,5) = 1e-6;
% cd_des(1,5) = 1e-6;
clcd_des(1,5) = 1e-6;

des = [cl_des; alpha_des; clcd_des];

% Fit cubics
p1 = NaN(size(des,1),4); % cubic poly coeffs
p2 = NaN(size(des,1),2); % straight line coeffs
N = 101;
x1 = linspace(t_c(1),t_c(4),N);
x2 = linspace(t_c(4),t_c(5),N);
y1 = NaN(size(des,1),N);
y2 = NaN(size(des,1),N);
for i = 1:size(des,1)
    p1(i,:) = polyfit(t_c(1:4),des(i,1:4),3);
    y1(i,:) = polyval(p1(i,:),x1);
    p2(i,:) = polyfit(t_c(4:5),des(i,4:5),1);
    y2(i,:) = polyval(p2(i,:),x2);
end

% Plot polynomials
figure
labels = ["c_{l,des} [-]", "\alpha_{des} [deg]", "{(c_l/c_d)}_{des} [-]"];

for i = 1:3
    subplot(3,1,i)
    plot(t_c,des(i,:),'x')
    hold on
    plot(x1,y1(i,:),'r');
    plot(x2,y2(i,:),'r');
    hold off
    ylabel(labels(i))
    axis tight
    grid on
end
xlabel('Relative thickness [%]')
sgtitle('New design polynomials')
end

%% Design polynomials
function p = newThickness(t_old,r_old,R_old,R_new,plt)
% Plot original absolute thickness distribution
if plt == 1
    figure
    plot(r_old,t_old)
    hold on
    grid on
end

% Create new absolute thickness distribution
scale   = R_new / R_old; % scaling factor
r_new   = r_old * scale; % new radial positions
t_new   = t_old * scale; % new thicknesses

% Constraints
mask1 = t_new < max(t_old); % cap new thickness distribution with old maximum
t_new(t_new > max(t_old)) = max(t_old);

% Create splines for new distribution and plot
if plt == 1
    plot(r_new,t_new)
    splines = csapi(r_new,t_new);
    plot(r_new,fnval(splines,r_new),'marker','x')
end

% Fit polynomial to new thickness distribution
x = r_new(mask1);
y = t_new(mask1);
p = polyfit(x,y,6); % 1MW also uses 6th order

% remove any fitted values that now exceed maximum of old distribution
t_new_fitted = polyval(p,r_new(mask1));
mask2 = t_new_fitted < max(t_old);
t_new_fitted = t_new_fitted(mask2);
r_new_fitted = x(mask2);

% Remove last point for fit (forced to zero in original design)
t_new_fitted = t_new_fitted(1:end-1);
r_new_fitted = r_new_fitted(1:end-1);

% Refit polynomial
p = polyfit(r_new_fitted,t_new_fitted,6);

% Plot fitted polynomial
if plt == 1
    plot(r_new_fitted,t_new_fitted,'marker','s')
    legend('DTU 10MW','Redesign','Redesign (splines)','Redesign (6th order polyfit)');
    xlabel('Radius [m]'); ylabel('Absolute thickness [m]');
end
end

function t = thickness(r,p,t_max,Rnew)
% Absolute thickness [m] as a function of radius [m] for redesigned blade
t = p(1)*r.^6 + p(2)*r.^5 + p(3)*r.^4 + p(4)*r.^3 + p(5)*r.^2 + p(6)*r + p(7);
t(t >= t_max) = t_max;
t(r == Rnew) = 1e-2;
end

function x = x_des(that,p1,p2)
% Design cl / clcd / alpha [ - / deg / - ] as a function of t/c [%]
x = NaN(length(that));
that(that < 24.1) = 24.1;
that(that > 100)  = 100;
for i = 1:length(that)
    if that(i) <= 48
        x(i) = p1(1)*that(i)^3 + p1(2)*that(i)^2 + p1(3)*that(i) + p1(4);
    else
        x(i) = p2(1)*that(i) + p2(2);
    end
end
end

%% Residual function for c and a'
function [out,varargout] = residuals(x,rotor,tsr,idx,p,p1,p2,t_max,Rnew)
% Calculate the residuals for chord and tangential induction factor

% unpack inputs
R = rotor.R; B = rotor.B; a = rotor.a;

% spanwise position
r = rotor.r(idx);

% unpack x
c  = x(1);
ap = x(2);

% calculate that, cl and cl/cd from design polynomials
t    = thickness(r,p,t_max,Rnew);
that = t / c;
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
sigma = (c*B) / (2*pi*r);

% calculate residuals
res_c  = 4*pi*r * sin(phi)^2 * F * ((2*a) / ...
    (cy * B * (1-a))) - c;
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
    varargout{1} = [t,c,phi,alpha,beta,cl,cd,ap,cp,ct];
%     varargout{1} = [t,c,phi,alpha,beta,cl,cd,a,ap,cp,ct];
end
end

%% Constraints
function [out,a] = flattenTip(var1,var2,rotor,r_R)
a = var1( find( (rotor.r/rotor.R) > r_R, 1 ) );
var1((rotor.r/rotor.R) > r_R) = ones(1,length(var2((rotor.r/rotor.R) > r_R))) .* a;
out = var1;
end