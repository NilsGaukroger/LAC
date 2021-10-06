function plot_modes(omegas,V,varargin)
% Plot the oscillations of the two masses for the two mode shapes.
if nargin == 2
    t = linspace(0,4,201);
elseif nargin == 3
    t = varargin{1};
elseif nargin > 3
    error('Too many input arguments')
end

% Plot the oscillating undamped mode
figure
for i = 1:2 % loop over modes
    subplot(2,1,i)
    mode = real(V(:,i) * exp(1i * omegas(i) * t));
    plot(t,mode(1,:)); hold on
    plot(t,mode(2,:)); hold off
    grid on
    tit = ['Mode ', string(i)];
    title(tit)
    xlim([0 max(t)])
    xlabel('Time [s]')
    if i == 1
        legend('x_1','x_2')
    end
end
end