% 46320 LAC Course
% Streamlined code for aerodynamic design
% Nils Joseph Gaukroger
% 18th September 2021

close all; clear; clc

%% Define rotor
rotor.R     = 35;  % length of blade [m]
rotor.tsr   = 9.0; % TSR [-]
rotor.B     = 3;   % number of blades [-]
rotor.a     = 1/3; % axial induction [-]


spacing     = 0.5; % increment for spanwise discretisation [m]
rotor.r_lst = 5:spacing:rotor.R-spacing;   % blade span [m]




rotor.t     = 1.5; % absolute thickness...
             % at blade span [m]
rotor.cl    = 0.6; % design lift coefficient [-]
rotor.clcd  = 120; % design cl/cd [-]
rotor.alpha = 4;   % design AoA [deg]

%% Functions

function [out,varargout] = residuals(x,rotor,varargin)
% Calculate the residuals for chord and tangential induction factor

% unpack inputs
R = rotor.R; tsr = rotor.tsr;
B = rotor.B; a = rotor.a;

% spanwise position
if nargin == 2
    r = rotor.r;
elseif nargin > 2
    r = rotor{1};
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
    rotor.t = t;
    rotor.c = c;
    rotor.phi = phi;
    rotor.alpha = alpha_des(that);
    rotor.beta = phi - alpha_des(that);
    rotor.cl = cl;
    rotor.cd = cd;
    rotor.ap = ap;
    rotor.cp = ((1-a)^2 + (tsr*(r/R)) * (1+ap)^2) * tsr * (r/R) * sigma * cx;
    rotor.ct = ((1-a)^2 + (tsr*(r/R)) * (1+ap)^2) * sigma * cy;
    % rotor.a     = 1 / (((4*F*np.sin(phi)**2) / (sigma * cy)) + 1)
    varargout{1} = rotor;
end
end

%% Design polynomials for the 35-m blade

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