function plotGeo(orig,comp)
    figure
    subplot(3,1,1)
    plot(comp.r/comp.R,comp.c); hold on
    plot(orig.r/orig.R,fnval(orig.c,orig.r),'x'); hold off
    ylabel('Chord [m]');
    legend('Redesign','DTU 10MW RWT')
    grid on
    subplot(3,1,2)
    plot(rotor.r/rotor.R,rad2deg(result.beta)); hold on
    plot(DTU.r/DTU.R,fnval(DTU.beta,DTU.r),'x'); hold off
    ylabel('Twist, \beta [deg]');
    legend('Redesign','DTU 10MW RWT')
    grid on
    subplot(3,1,3)
    plot(rotor.r/rotor.R,(rotor.t./result.c)*100); hold on
    plot(DTU.r/DTU.R,fnval(DTU.that,DTU.r),'x'); hold off
    ylabel('t/c [%]'); xlabel('Non-dimensional radius [-]')
    legend('Redesign','DTU 10MW RWT')
    grid on

    if length(tsr) ~= 1
        figure
        plot(tsr,result.CP)
        xlabel('TSR [-]'); ylabel('C_P [-]')
        grid on
        % Optimal TSR
        [CPmax,CPmax_idx] = max(result.CP);
        tsr_opt   = tsr(CPmax_idx);
    end
end