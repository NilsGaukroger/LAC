function [R2,V2,n] = rotorScaling(V1,I1,R1,I2,varargin)

% 46320 LAC Course
% Rotor scaling algorithm
% Nils Joseph Gaukroger
% 11th September 2021

% DTU 10MW rated speed
% V1 = 11.4 m/s
% D1 = 178.3 m

% IEC turbulence classes
% A: I1 = 0.16
% B: I2 = 0.14

V1_max = V1 * (1 + 2*I1);

if nargin == 4
    eps      = 1e-6;
    R2_guess = (I1/I2) * R1;
elseif nargin == 5
    eps      = varargin{1};
    R2_guess = (I1/I2) * R1;
elseif nargin == 6
    eps      = varargin{1};
    R2_guess = varargin{2};
end

conv = 0.1; n = 0;

while conv > eps
    V2 = ((V1^3 * R1^2) / (R2_guess^2))^(1/3);
    R2 = (((V1_max)^2 / (V2*(1+2*I2))^2) * R1^2)^(1/2);
    
    % Check convergence
    conv = abs(R2 - R2_guess);
    R2_guess = R2;
    n = n + 1;
end

end