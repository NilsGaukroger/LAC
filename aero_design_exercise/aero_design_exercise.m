% 46320 LAC Course
% Code for aerodynamic design
% Nils Joseph Gaukroger
% 11th September 2021

close all; clear; clc;

%% Part 1: Define a residual function
inputs.R     = 35;  % length of blade [m]
inputs.tsr   = 9.0; % TSR [-]
inputs.B     = 3;   % number of blades [-]
inputs.a     = 1/3; % axial induction [-]
inputs.r     = 8;   % blade span [m]
inputs.t     = 1.5; % absolute thickness...
             % at blade span [m]
inputs.cl    = 0.6; % design lift coefficient [-]
inputs.clcd  = 120; % design cl/cd [-]
inputs.alpha = 4;   % design AoA [deg]

% test the values
x = [3, 0.001]; % [c (m), a']
out = part1_residuals(x,inputs);

disp(out)

%% Part 2: Solve BEM equations at one span with fixed aerodynamic values
x0 = [3, 0.001]; % initial guess [c (m), a']
lb = [0, 0]; % lower bounds
ub = [inf, 1]; % upper bounds

res = lsqnonlin(@(x)part1_residuals(x,inputs),x0,lb,ub);
c = res(1); ap = res(2);

fprintf('c = %.4f m\n', c);
fprintf('a'' = %.4f\n', ap);

%% Part 3: Thickness and design polynomials

r    = linspace(0,35,201);
that = linspace(15,100,201);

figure

% absolute thickness vs. blade length
subplot(2,2,1)
plot(r,thickness(r))
xlabel('Blade length [m]'); ylabel('Thickness [m]')
grid on

% cl,des versus relative thickness
subplot(2,2,2)
plot(that,cl_des(that))
xlabel('Relative thickness [%]'); ylabel('c_{l,des} [-]')
grid on

% clcd,des versus relative thickness
subplot(2,2,3)
plot(that,clcd_des(that))
xlabel('Relative thickness [%]'); ylabel('(c_l/c_d)_{des} [-]')
grid on

% alpha,des versus relative thickness
subplot(2,2,4)
plot(that,alpha_des(that))
xlabel('Relative thickness [%]'); ylabel('\alpha_{des} [deg]')
grid on

%% Part 4: Solve BEM equations at one span with design polynomials
x0 = [3, 0.001];
residuals(x0,inputs)

x0 = [3, 0.001]; % initial guess [c (m), a']
lb = [0, 0]; % lower bounds
ub = [inf, 1]; % upper bounds

result = lsqnonlin(@(x)residuals(x,inputs),x0,lb,ub);
c = result(1); ap = result(2);

fprintf('c = %.4f m\n', c);
fprintf('a'' = %.4f\n', ap);

%% Part 5: Solve BEM equations at multiple spans with design polynomials
radii = 5:0.5:34.5;
out = NaN(length(radii), 2);
for i = 1:length(radii)
    res = lsqnonlin(@(x)residuals(x,inputs,radii(i)),x0,lb,ub);
    out(i,1) = res(1); out(i,2) = res(2);
    x0 = out(i,:);
end

out(1:5,:)

%% Part 6: Calculate other parameters from output
c = out(:,1); ap = out(:,2);

t = NaN(length(radii),1);
that = NaN(length(radii),1);
alpha = NaN(length(radii),1);
for i = 1:length(radii)
    t(i) = thickness(radii(i));
    that(i) = t(i) / c(i);
    alpha(i) = alpha_des(that(i)*100); % [deg]
end

phi   = atan(((1-inputs.a)./(1+ap)) .* (inputs.R./(radii'.*inputs.tsr)));
theta = rad2deg(phi) - alpha;

figure
ax1 = subplot(3,1,1);
plot(radii,c)
ylabel('Chord [m]')
grid on

ax2 = subplot(3,1,2);
plot(radii,theta)
ylabel('Twist [deg]')
grid on

ax3 = subplot(3,1,3);
plot(radii,that)
xlabel('Radius [m]'); ylabel('Rel. thickness [%]')
grid on

linkaxes([ax1,ax2,ax3],'x')

%% Functions
% Part 1
function out = part1_residuals(x,inputs)
% Calculate the residuals for chord and tangential induction factor

% unpack inputs
R = inputs.R; tsr = inputs.tsr;
B = inputs.B; a = inputs.a;
r = inputs.r; cl = inputs.cl;
clcd = inputs.clcd;

% unpack x
c  = x(1);
ap = x(2);

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
end

% Part 2

% Part 3
function t = thickness(r)
% Absolute thickness [m] as a function of radius [m] for 35m blade
r(r < 5) = 5;
t = 9.35996E-08*r.^6 - 1.2911E-05*r.^5 + 7.15038E-04*r.^4 - 2.03735E-02*r.^3 + 3.17726E-01*r.^2 - 2.65357E+00*r + 1.02616E+01;
end

function cl = cl_des(that)
% Design cl [-] as a function of t/c [%]
cl = NaN(length(that));
for i = 1:length(that)
    if that(i) < 36
        cl(i) = -7.34862E-07*that(i)^5 + 1.10229E-04*that(i)^4 - 6.40432E-03*that(i)^3 + 1.79563E-01*that(i)^2 - 2.43397E+00*that(i) + 1.36000E+01;
    else
        cl(i) = -0.0094*that(i) + 0.9375;
    end
end
end

function clcd = clcd_des(that)
% Design cl [-] as a function of t/c [%]
clcd = NaN(length(that));
for i = 1:length(that)
    if that(i) < 36
        clcd(i) = -8.10212E-03*that(i)^4 + 8.73313E-01*that(i)^3 - 3.41180E+01*that(i)^2 + 5.66297E+02*that(i) - 3.24932E+03;
    else
        clcd(i) = -0.8906*that(i) + 89.063;
    end
end
end

function alpha = alpha_des(that)
% Design AoA [deg] as a function of t/c [%]
alpha = NaN(length(that));
for i = 1:length(that)
    if that(i) < 36
        alpha(i) = 2.32706E-05*that(i)^5 - 2.87331E-03*that(i)^4 + 1.36343E-01*that(i)^3 - 3.10470E+00*that(i)^2 + 3.38460E+01*that(i) - 1.36500E+02;
    else
        alpha(i) = -0.0078*that(i) + 0.7813;
    end
end
end

% Part 4
function out = residuals(x,inputs,varargin)
% Calculate the residuals for chord and tangential induction factor

% unpack inputs
R = inputs.R; tsr = inputs.tsr;
B = inputs.B; a = inputs.a;

% spanwise position
if nargin == 2
    r = inputs.r;
elseif nargin > 2
    r = varargin{1};
end

% unpack x
c  = x(1);
ap = x(2);

% calculate that, cl and cl/cd from design polynomials
that = thickness(r) / c;
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
end