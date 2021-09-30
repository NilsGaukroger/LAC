function p = newThickness(t_old,r_old,R_old,R_new,plt)
% Function for creating new absolut thickness distribution.

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