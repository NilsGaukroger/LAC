function [p1,p2,cl_max,x_des] = desPolys(aerofoil,x1,n)
% Function to create new design polynomials from aerofoil polars, first
% finds cl_max, then cl_des based on x1, then corresponding alpha_des and
% cd_des (and (cl/cd)_des). Then fits splines and plots the result.

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
    cl_max_idx = find(diff(aerofoil{2,i}.cl(mask)) <= 0, 1, 'First');
    cl_max = aerofoil{2,i}.cl(cl_max_idx + idx0 - 1);
    alpha_max = aerofoil{2,i}.alpha(cl_max_idx + idx0 - 1);
    cl_des(i) = cl_max - x1(i);
    mask2(:,i) = aerofoil{2,i}.alpha >= 0 & aerofoil{2,i}.alpha < alpha_max;
    alpha_des(i) = interp1(aerofoil{2,i}.cl(mask2(:,i)),aerofoil{2,i}.alpha(mask2(:,i)),cl_des(i));
    cd_des(i) = interp1(aerofoil{2,i}.cl(mask2(:,i)),aerofoil{2,i}.cd(mask2(:,i)),cl_des(i));
    clcd_des(i) = cl_des(i) / cd_des(i);
end
cl_des(1,5) = 1e-6;
alpha_des(1,5) = 1e-6;
% cd_des(1,5) = 1e-6;
clcd_des(1,5) = 1e-6;

x_des = [cl_des; alpha_des; clcd_des];

% Fit cubics
p1 = NaN(size(x_des,1),4); % cubic poly coeffs
p2 = NaN(size(x_des,1),2); % straight line coeffs
N = 101;
x1 = linspace(t_c(1),t_c(4),N);
x2 = linspace(t_c(4),t_c(5),N);
y1 = NaN(size(x_des,1),N);
y2 = NaN(size(x_des,1),N);
for i = 1:size(x_des,1)
    p1(i,:) = polyfit(t_c(1:4),x_des(i,1:4),3);
    y1(i,:) = polyval(p1(i,:),x1);
    p2(i,:) = polyfit(t_c(4:5),x_des(i,4:5),1);
    y2(i,:) = polyval(p2(i,:),x2);
end

% Plot polynomials
figure
labels = ["c_{l,des} [-]", "\alpha_{des} [deg]", "{(c_l/c_d)}_{des} [-]"];

for i = 1:3
    subplot(3,1,i)
    plot(t_c,x_des(i,:),'x','color','k')
    hold on
    plot(x1,y1(i,:),'r');
    plot(x2,y2(i,:),'r');
    hold off
    ylabel(labels(i))
    axis tight
    grid on
end
xlabel('Relative thickness [%]')
% sgtitle('New design polynomials')
end