function plot_polars(aerofoil,limits,n,mode,x_des)
% function for plotting either cl vs. alpha (mode == 1)
% or cl vs. cd (mode == 2) with specified axis limits

figure
for i = 1:n
    subplot(ceil(n/2),2,i)
    if mode == 1
        plot(aerofoil{2,i}.alpha,aerofoil{2,i}.cl); hold on
        plot(x_des(2,i),x_des(1,i),'x','color','k','MarkerSize',10); hold off
        xlabel('\alpha [deg]'); ylabel('c_l [-]');
        ylim(limits(2,:))
        xlim(limits(1,:))
    elseif mode == 2
        mask = (aerofoil{2,1}.alpha >= 0 & aerofoil{2,1}.alpha <= 90);
        plot(aerofoil{2,i}.cd(mask),aerofoil{2,i}.cl(mask)); hold on
        plot(x_des(1,i)/x_des(3,i),x_des(1,i),'x','color','k','MarkerSize',10); hold off
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