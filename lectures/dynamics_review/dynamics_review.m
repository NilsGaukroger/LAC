% 46320 LAC Course
% A review of modal analysis for linear, time-invariant systems.
% Author: Nils Joseph Gaukroger
% 5th October 2021

close all; clear variables; clc

%% System parameters (mimicking DTU 10MW)
m1 = 551556;  % (nacelle mass + hub mass) [kg]
m2 = 122442;  % (rotor mass - hub mass) [kg]
c1 = 6866;    % equivalent tower damping [kgs]
c2 = 3599;    % equivalent tower damping [kgs]
k1 = 1745883; % equivalent tower stiffness [N/m]
k2 = 1503050; % equivalent rotor stiffness [N/m]

%% Assemble the mass, stiffness and damping matrices
M = [m1 0;
     0 m2];
C = [c1+c2, -c2;
      -c2,   c2];
K = [k1+k2, -k2;
      -k2,   k2];

%% Assemble the state space matrices
I = eye(2); O = zeros(2);
A = [O, I;
    -K,-C];
B = [I, O;
     O, M];

%% Q1 Solving the undamped system (no state space required)

% solve generalised eigenvalue problem
[V,D] = eig(-K,M);

% convert D from diagonal matrix to column vector (only possible for uncoupled systems)
D = diag(D);

% sort eigenvalues and vectors in increasing order of frequency
[~,idx] = sort(abs(D));
D       = D(idx);
V       = V(:,idx);

% normalise eigenvectors (not required in Python)
for i = 1:size(V,2)
    V(:,i) = V(:,i)/norm(V(:,i));
end

% calculate frequencies
omegas = abs(sqrt(D));    % [rad/s]
freqs  = omegas / (2*pi); % [Hz]

fprintf('UNDAMPED\n')
fprintf('    Mode 1    Mode 2\n')
fprintf('fk  %.4f    %.4f\n',freqs(1),freqs(2))
fprintf('x1 %.4f   %.4f\n',V(1,1),V(1,2))
fprintf('x2 %.4f    %.4f\n',V(2,1),V(2,2))

% plot modes
plot_modes(omegas,V);

%% Q2 Solving the damped system using state space

% solve the generalised eigenvalue problem
[P,L] = eig(A,B);

% convert L from diagonal matrix to column vector (only possible for uncoupled systems)
L = diag(L);

% isolate the states that correspond to displacements
L = L(1:2:end); % take every other eigenvalue, starting at 1
P = P(1:2:3);   % take every other eigenvector, ending at 2

% sort eigenvalues and vectors in increasing order of frequency
[~,idx] = sort(abs(L));
L       = L(idx);
P       = P(:,idx);

% normalise eigenvectors (not required in Python)
for i = 1:size(P,2)
    P(:,i) = P(:,i)/norm(P(:,i));
end

% calculate the natural frequencies and damping
omegas = abs(L);            % [rad/s]
freqs  = omegas / (2*pi);   % [Hz]
zetas  = cos(angle(L))*100; % If you want to convert zeta to percent critical, multiply by 100!

% intermediates: magnitude and phase of modes
mags = abs(P);
pha  = angle(P);

fprintf('DAMPED\n')
print('    Mode 1    Mode 2')
% print('fk  {freqs[0]:.2f}   {freqs[1]:.2f}')
% print('zk  {zetas[0]:.2f}   {zetas[1]:.2f}')
% print('x1  {mags[0, 0]:.2f}<{pha[0, 0]:.2f}   {mags[0, 1]:.2f}<{pha[0, 1]:.2f}')
% print('x2  {mags[1, 0]:.2f}<{pha[1, 0]:.2f}   {mags[1, 1]:.2f}<{pha[1, 1]:.2f}')

% plot mode oscillation vs. time
plot_modes(omegas,P)
sgtitle('Damped')