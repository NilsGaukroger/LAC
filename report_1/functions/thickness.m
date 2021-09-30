function t = thickness(r,p,t_max,Rnew)
% Absolute thickness [m] as a function of radius [m] for redesigned blade
t = p(1)*r.^6 + p(2)*r.^5 + p(3)*r.^4 + p(4)*r.^3 + p(5)*r.^2 + p(6)*r + p(7);
t(t >= t_max) = t_max;
t(r == Rnew) = 1e-2;
end