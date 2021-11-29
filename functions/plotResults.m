function plotResults(filename)
data = readtable(filename);
idx  = [15,3,10,100]; scale = [1,1,1,1e-3];
ylabels = ["Wind speed [m/s]", "Pitch angle [deg]", "Rotational speed [rad/s]", "Electrical power [kW]"];

for i = 1:4
    subplot(2,2,i)
    plot(data{:,1},data{:,idx(i)}.*scale(i))
    xlabel('Time [s]'); ylabel(ylabels(i))
    box on
    grid minor
end
end