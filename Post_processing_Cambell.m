clear all
clc
close all

%% Outputs from HAWC2Stab2.exe
Aeroelastic=readcell('Aeroelastic_frequencies_and_Damping_vs_V.txt'); 
Structural=readcell('Modal_frequencies_and_Damping_vs_V.txt');

%% Post processing
for i=2:23
    V(i-1)=Aeroelastic{i,1};
    for j=3:13
        Naturalfreq.aeroelastic(i-1,j-2)=Aeroelastic{i,j};
        Naturalfreq.structural(i-1,j-2)=Structural{i,j};
    end
end

for i=2:23
    for j=15:25
        Damping.aeroelastic(i-1,j-14)=Aeroelastic{i,j};
        Damping.structural(i-1,j-14)=Structural{i,j};
    end
end

vars={'-o','-+','-*','-.','-x','-_','-|','-s','-d','-^','-v'};
figure
for i=1:length(Damping.structural(1,:))
    plot(V,Naturalfreq.aeroelastic(:,i),vars{i})
    
    hold on
end
grid on
grid minor
xlabel('V(m/s)');
ylabel('Nat. frequency(Hz)');
title('Cambell diagram for aeroelastic mode shapes');
legend('1st Twr FA','1st Twr SS','1st BW flap','1st SYM flap','1st FW flap','1st BW edge','1st FW edge','2nd BW flap','2nd FW flap','2nd SYM flap','1st COL edge','Location','northwest','NumColumns',1);

figure
for i=1:length(Damping.structural(1,:))
    plot(V,Naturalfreq.structural(:,i),vars{i})
    
    hold on
end
grid on
grid minor
xlabel('V(m/s)');
ylabel('Nat. frequency(Hz)');
title('Cambell diagram for structural mode shapes');
legend('1st Twr FA','1st Twr SS','1st BW flap','1st SYM flap','1st FW flap','1st BW edge','1st FW edge','2nd BW flap','2nd FW flap','2nd SYM flap','1st COL edge','Location','northwest','NumColumns',1);

figure
for i=1:length(Damping.structural(1,:))
    plot(V,Damping.structural(:,i),vars{i})
    
    hold on
end
grid on
grid minor
xlabel('V(m/s)');
ylabel('Damping[% critical]');
title('Damping values for structural mode shapes');
legend('1st Twr FA','1st Twr SS','1st BW flap','1st SYM flap','1st FW flap','1st BW edge','1st FW edge','2nd BW flap','2nd FW flap','2nd SYM flap','1st COL edge','Location','northwest','NumColumns',1);

figure
for i=1:length(Damping.structural(1,:))
    plot(V,Damping.aeroelastic(:,i),vars{i})
    
    hold on
end
grid on
grid minor
xlabel('V(m/s)');
ylabel('Damping[% critical]');
title('Damping values for aeroelastic mode shapes');
legend('1st Twr FA','1st Twr SS','1st BW flap','1st SYM flap','1st FW flap','1st BW edge','1st FW edge','2nd BW flap','2nd FW flap','2nd SYM flap','1st COL edge','Location','northwest','NumColumns',1);