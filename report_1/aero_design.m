% 46320 LAC Course
% Assignment 1
% Nils Joseph Gaukroger
% 18th September 2021

close all; clear variables; clc

%% Design polynomials from DTU 10MW reports
[DTU.c,DTU.that,DTU.beta] = DTU10MW_des();
fn = fieldnames(DTU);
x = linspace(2.8,89.166,101);
figure
for i = 1:length(fn)
    subplot(3,1,i)
    plot(x,fnval(DTU.(fn{i}),x))
end

%% New rotor radius
[Rnew,~] = rotorScaling(11.4,0.16,178.3/2,0.14);

%% Aerofoil data
S = dir('polars\*.txt');
aerofoil = cell(2,length(S));
aerofoil(1,:) = {S.name};
vars = {'alpha','cl','cd','fsst','clinv','clfs','cm'};

for i = 1:size(aerofoil,2)
    aerofoil{2,i} = readtable("polars\" + aerofoil{1,i});
    aerofoil{2,i}.Properties.VariableNames = vars;
end

%% cl vs alpha
figure
xlimits = [0 90];
for i = 1:4 % ignoring 60.0% and cylinder
    subplot(2,2,i)
    plot(aerofoil{2,i}.alpha,aerofoil{2,i}.cl)
    xlabel('\alpha [deg]'); ylabel('c_l [-]');
%     xlim([0 90]);
    title(aerofoil{1,i});
    grid on
end

%% cl vs cd
figure
mask = (aerofoil{2,1}.alpha >= 0 & aerofoil{2,1}.alpha <= 90);
for i = 1:4 % ignoring 60.0% and cylinder
    subplot(2,2,i)
    plot(aerofoil{2,i}.cd(mask),aerofoil{2,i}.cl(mask))
    xlabel('c_d [-]'); ylabel('c_l [-]');
    title(aerofoil{1,i});
    grid on
end

%% the clever bit
idx0 = find(aerofoil{2,1}.alpha == 0);

t_c = [24.1, 30.1, 36.0, 48.0, 100];

x = [0.4, 0.4, 0.4, 0.4];
mask2 = false(size(aerofoil{2,1},1),4);
cl_des = NaN(1,4);
alpha_des = NaN(1,4);
cd_des = NaN(1,4);
clcd_des = NaN(1,4);
for i = 1:4
    clmax_idx = find(diff(aerofoil{2,i}.cl(mask)) <= 0, 1, 'First');
    clmax = aerofoil{2,i}.cl(clmax_idx + idx0);
    alpha_max = aerofoil{2,i}.alpha(clmax_idx + idx0 - 1);
    cl_des(i) = clmax - x(i);
    mask2(:,i) = aerofoil{2,i}.alpha >= 0 & aerofoil{2,i}.alpha < alpha_max;
    alpha_des(i) = interp1(aerofoil{2,i}.cl(mask2(:,i)),aerofoil{2,i}.alpha(mask2(:,i)),cl_des(i));
    cd_des(i) = interp1(aerofoil{2,i}.cl(mask2(:,i)),aerofoil{2,i}.cd(mask2(:,i)),cl_des(i));
    clcd_des(i) = cl_des(i) / cd_des(i);
end

figure
des = [cl_des; alpha_des; clcd_des];

for i = 1:3
    subplot(3,1,i)
    plot(t_c,[des(i,:) 0],'x')
    axis tight
end
xlabel('Relative thickness [%]')

%% 



%% Define rotor
% rotor.R     = 35;  % length of blade [m]
% rotor.tsr   = 9.0; % TSR [-]
% rotor.B     = 3;   % number of blades [-]
% rotor.a     = 1/3; % axial induction [-]
% 
% spacing     = 0.5; % increment for spanwise discretisation [m]
% rotor.radii = (5:spacing:rotor.R-(spacing*2));   % blade span [m]
% [result.t, result.c, result.phi, result.alpha, result.beta, result.cl, result.cd, ...
%     result.ap, result.cp, result.ct] = deal(NaN(1,length(rotor.radii)));
% 
% x0 = [3, 0.001]; % initial guess
% lb = [0, 0]; % lower bounds
% ub = [inf, 1]; % upper bounds
% 
% for i = 1:length(rotor.radii)
%     x0 = lsqnonlin(@(x)residuals(x,rotor,i),x0,lb,ub);
%     [~,values] = residuals(x0,rotor,i);
%     result.t(i)     = values(1);
%     result.c(i)     = values(2);
%     result.phi(i)   = values(3);
%     result.alpha(i) = values(4);
%     result.beta(i)  = values(5);
%     result.cl(i)    = values(6);
%     result.cd(i)    = values(7);
%     result.ap(i)    = values(8);
%     result.cp(i)    = values(9);
%     result.ct(i)    = values(10);
% end
% result.CP = (2/rotor.R^2) * trapz(rotor.radii.*result.cp);
% 
% %% Examining designs for varying TSR
% tsr_lst = 2:0.5:12;
% [result_tsr.t, result_tsr.c, result_tsr.phi, result_tsr.alpha, result_tsr.beta,...
%     result_tsr.cl, result_tsr.cd, result_tsr.ap, result_tsr.cp, result_tsr.ct]...
%     = deal(NaN(length(tsr_lst),length(rotor.radii)));
% result_tsr.CP = NaN(length(tsr_lst),1);
% [chord, twist, rel_t] = deal(NaN(length(tsr_lst),length(rotor.radii)));
% for j = 1:length(tsr_lst)
%     x0 = [3, 0.001]; % initial guess
%     for i = 1:length(rotor.radii)
%         x0 = lsqnonlin(@(x)residuals(x,rotor,i,tsr_lst(j)),x0,lb,ub);
%         [~,values_tsr] = residuals(x0,rotor,i,tsr_lst(j));
%         result_tsr.t(j,i)     = values_tsr(1);
%         result_tsr.c(j,i)     = values_tsr(2);
%         result_tsr.phi(j,i)   = values_tsr(3);
%         result_tsr.alpha(j,i) = values_tsr(4);
%         result_tsr.beta(j,i)  = values_tsr(5);
%         result_tsr.cl(j,i)    = values_tsr(6);
%         result_tsr.cd(j,i)    = values_tsr(7);
%         result_tsr.ap(j,i)    = values_tsr(8);
%         result_tsr.cp(j,i)    = values_tsr(9);
%         result_tsr.ct(j,i)    = values_tsr(10);
%     end
%     result_tsr.CP(j,1) = (2/rotor.R^2) * trapz(rotor.radii.*result_tsr.cp(j,:));
% end
% 
% figure
% subplot(3,1,1)
% plot(rotor.radii,result_tsr.c)
% ylabel('Chord [m]');
% grid on
% subplot(3,1,2)
% plot(rotor.radii,rad2deg(result_tsr.beta))
% ylabel('Twist, \beta [deg]');
% grid on
% subplot(3,1,3)
% plot(rotor.radii,(result_tsr.t./result_tsr.c)*100)
% ylabel('t/c [%]'); xlabel('Radius [m]')
% grid on
% 
% figure
% plot(tsr_lst,result_tsr.CP)
% xlabel('TSR [-]'); ylabel('C_P [-]')
% grid on

%% Functions

function [out,varargout] = residuals(x,rotor,idx,varargin)
% Calculate the residuals for chord and tangential induction factor

% unpack inputs
R = rotor.R; tsr = rotor.tsr;
B = rotor.B; a = rotor.a;

% spanwise position
r = rotor.radii(idx);

% optional arguments
if nargin > 3
    tsr = varargin{1};
    rotor.tsr = tsr;
end

% unpack x
c  = x(1);
ap = x(2);

% calculate that, cl and cl/cd from design polynomials
t    = thickness(r);
that = t / c;
cl   = cl_des(that*100);
clcd = clcd_des(that*100);

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
    alpha = deg2rad(alpha_des(that*100));
    beta = phi - deg2rad(alpha_des(that*100));
    cp = ((1-a)^2 + (tsr*(r/R)) * (1+ap)^2) * tsr * (r/R) * sigma * cx;
    ct = ((1-a)^2 + (tsr*(r/R)) * (1+ap)^2) * sigma * cy;
%     result.a = 1 / (((4*F*np.sin(phi)**2) / (sigma * cy)) + 1)
    varargout{1} = [t,c,phi,alpha,beta,cl,cd,ap,cp,ct];
%     varargout{1} = [t,c,phi,alpha,beta,cl,cd,a,ap,cp,ct];
end
end

%% Design polynomial splines for the DTU 10MW
function [c,that,beta] = DTU10MW_des()
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
that_coeffs = [
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
end