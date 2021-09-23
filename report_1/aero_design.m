% 46320 LAC Course
% Assignment 1
% Nils Joseph Gaukroger
% 18th September 2021

close all; clear variables; clc

%% Design polynomials from DTU 10MW reports
[DTU.c,DTU.that,DTU.beta,DTU.t,DTU.r] = DTU10MW_des(1);

%% New rotor radius
DTU.R = 89.1660;
[Rnew,~] = rotorScaling(11.4,0.16,DTU.R,0.14);

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
x = [0.3, 0.4, 0.35, 0.5]; % cl,max - cl,des
[p1,p2] = desPolys(aerofoil,x,4);

%% Absolute thickness
figure
plot(DTU.r,DTU.t)
hold on
grid on

% new absolute thickness
scale = Rnew / DTU.R;
new.r = DTU.r * scale;
new.t = DTU.t * scale;
mask = new.t < max(DTU.t);
new.t(new.t > max(DTU.t)) = max(DTU.t);
plot(new.r,new.t)
ff = csapi(new.r,new.t);
plot(new.r,fnval(ff,new.r),'marker','x')

x = new.r(mask);
y = new.t(mask);
p = polyfit(x(2:end-1),y(2:end-1),6);
r_plot = new.r(mask);
plot(r_plot(2:end-1),polyval(p,r_plot(2:end-1)),'marker','s')
legend('DTU 10MW','Redesign','Splines','polyfit');
xlabel('Radius [m]'); ylabel('Absolute thickness [m]');

%%
t_max = 5.38;
figure
plot(new.r,thickness(new.r,p,t_max,Rnew))

%% Residual
rotor.R     = Rnew;  % length of blade [m]
rotor.B     = 3;   % number of blades [-]
rotor.a     = 1/3; % axial induction [-]
rotor.tsr   = 9; % default tsr

spacing     = 0.2; % increment for spanwise discretisation [m]
rotor.radii = (2.8:spacing:rotor.R-(spacing*2));   % blade span [m]
[result.t, result.c, result.phi, result.alpha, result.beta, result.cl, result.cd, ...
    result.ap, result.cp, result.ct] = deal(NaN(1,length(rotor.radii)));

x0 = [3, 0.001]; % initial guess
lb = [0, 0]; % lower bounds
ub = [inf, 1]; % upper bounds

for i = 1:length(rotor.radii)
    x0 = lsqnonlin(@(x)residuals(x,rotor,rotor.tsr,i,p,p1,p2,t_max,Rnew),x0,lb,ub);
    [~,values] = residuals(x0,rotor,rotor.tsr,i,p,p1,p2,t_max,Rnew);
    result.t(i)     = values(1);
    result.c(i)     = values(2);
    result.phi(i)   = values(3);
    result.alpha(i) = values(4);
    result.beta(i)  = values(5);
    result.cl(i)    = values(6);
    result.cd(i)    = values(7);
    result.ap(i)    = values(8);
    result.cp(i)    = values(9);
    result.ct(i)    = values(10);
end
result.CP = (2/rotor.R^2) * trapz(rotor.radii,rotor.radii.*result.cp);

%% Geometry
figure
subplot(3,1,1)
plot(rotor.radii,result.c)
ylabel('Chord [m]')
grid on
subplot(3,1,2)
plot(rotor.radii,rad2deg(result.beta))
ylabel('Twist [deg]');
grid on
subplot(3,1,3)
plot(rotor.radii,(result.t./result.c)*100)
xlim([2.8 inf])
grid on
ylabel('Relative thickness [%]');
xlabel('Radius [m]');

figure
plot(rotor.radii,result.cp)
xlabel('Radius [m]'); ylabel('c_{p,loc} [-]')
grid on

figure
plot(rotor.radii,result.ct)
xlabel('Radius [m]'); ylabel('c_{t,loc} [-]')
grid on

%% Add constraints on geometry
% Absolute thickness
% Chord
% Relative thickness
% Twist
 
%% Examining designs for varying TSR
tsr_lst = 5:1:12;
[result_tsr.t, result_tsr.c, result_tsr.phi, result_tsr.alpha, result_tsr.beta,...
    result_tsr.cl, result_tsr.cd, result_tsr.ap, result_tsr.cp, result_tsr.ct]...
    = deal(NaN(length(tsr_lst),length(rotor.radii)));
result_tsr.CP = NaN(length(tsr_lst),1);
[chord, twist, that] = deal(NaN(length(tsr_lst),length(rotor.radii)));
for j = 1:length(tsr_lst)
    x0 = [3, 0.001]; % initial guess
    for i = 1:length(rotor.radii)
        x0 = lsqnonlin(@(x)residuals(x,rotor,tsr_lst(j),i,p,p1,p2,t_max,Rnew),x0,lb,ub);
        [~,values_tsr] = residuals(x0,rotor,tsr_lst(j),i,p,p1,p2,t_max,Rnew);
        result_tsr.t(j,i)     = values_tsr(1);
        result_tsr.c(j,i)     = values_tsr(2);
        result_tsr.phi(j,i)   = values_tsr(3);
        result_tsr.alpha(j,i) = values_tsr(4);
        result_tsr.beta(j,i)  = values_tsr(5);
        result_tsr.cl(j,i)    = values_tsr(6);
        result_tsr.cd(j,i)    = values_tsr(7);
        result_tsr.ap(j,i)    = values_tsr(8);
        result_tsr.cp(j,i)    = values_tsr(9);
        result_tsr.ct(j,i)    = values_tsr(10);
    end
    result_tsr.CP(j,1) = (2/rotor.R^2) * trapz(rotor.radii(2:end),rotor.radii(2:end).*result_tsr.cp(j,2:end));
end

figure
subplot(3,1,1)
plot(rotor.radii,result_tsr.c)
ylabel('Chord [m]');
grid on
subplot(3,1,2)
plot(rotor.radii,rad2deg(result_tsr.beta))
ylabel('Twist, \beta [deg]');
grid on
subplot(3,1,3)
plot(rotor.radii,(result_tsr.t./result_tsr.c)*100)
ylabel('t/c [%]'); xlabel('Radius [m]')
grid on

figure
plot(tsr_lst,result_tsr.CP)
xlabel('TSR [-]'); ylabel('C_P [-]')
grid on

%% Design polynomial splines for the DTU 10MW
function [c,that,beta,t,x] = DTU10MW_des(plt)
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
    x = linspace(2.8,89.166,101);
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
    clmax = aerofoil{2,i}.cl(clmax_idx + idx0);
    alpha_max = aerofoil{2,i}.alpha(clmax_idx + idx0 - 1);
    cl_des(i) = clmax - x1(i);
    mask2(:,i) = aerofoil{2,i}.alpha >= 0 & aerofoil{2,i}.alpha < alpha_max;
    alpha_des(i) = interp1(aerofoil{2,i}.cl(mask2(:,i)),aerofoil{2,i}.alpha(mask2(:,i)),cl_des(i));
    cd_des(i) = interp1(aerofoil{2,i}.cl(mask2(:,i)),aerofoil{2,i}.cd(mask2(:,i)),cl_des(i));
    clcd_des(i) = cl_des(i) / cd_des(i);
end
cl_des(1,5) = 1e-6;
alpha_des(1,5) = 1e-6;
cd_des(1,5) = 1e-6;
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
% function new_t = new_t(t,r,R)
% 
% 
% 
% figure
% plot(DTU.r,DTU.t)
% hold on
% grid on
% 
% % new absolute thickness
% scale = Rnew / DTU.R;
% new.r = DTU.r * scale;
% new.t = DTU.t * scale;
% new.t(new.t > max(DTU.t)) = max(DTU.t);
% plot(new.r,new.t)
% legend('DTU 10MW','Redesign');
% xlabel('Radius [m]'); ylabel('Absolute thickness [m]');
% ff = csapi(new.r,new.t);
% plot(new.r,fnval(ff,new.r));
% end

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
r = rotor.radii(idx);

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