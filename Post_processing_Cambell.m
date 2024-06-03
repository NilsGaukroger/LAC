clear all
clc
close all

%% Outputs from HAWC2Stab2.exe for new rotor
Aeroelastic=readcell('Aeroelastic_frequencies_and_Damping_vs_V.txt'); 
Structural=readcell('Modal_frequencies_and_Damping_vs_V.txt');

%% Outputs from pwr for DTU 10 MW and new rotor
turbine = {'DTU_10MW_flexible_hawc2s','DTU_10MW_flexible_hawc2s_newrotor'};
%elas    = {'flexible','rigid'};
%file    = {'/DTU_10MW_struct_','/redesign_struct_'};
pwr     = cell(1,2); % array for storing .pwr files, rows = turbine, cols = elasticity

for i = 1:length(turbine)
  %  for j = 1:size(pwr,2)
        pwr{i} = import_pwr(strcat(turbine{i},".pwr"));
  % end
end

Freq1P=(pwr{1,2}.('Speed').*pi./30)./(2.*pi); % Rotational frequency(Hz)
DTU10MW.Freq1P=(pwr{1,1}.('Speed').*pi./30)./(2.*pi); % Rotational frequency(Hz)

Freq3P=3.*Freq1P;
Freq6P=6.*Freq1P;
DTU10MW.Freq3P=3.*DTU10MW.Freq1P;
DTU10MW.Freq6P=6.*DTU10MW.Freq1P;


DTU10MW.Structural=readcell('Modal_frequencies_and_Damping_vs_V_DTU10MW.txt');
DTU10MW.Aeroelastic=readcell('Aeroelastic_frequencies_and_Damping_vs_V_DTU10MW.txt');

%% Post processing New rotor
for i=2:23
    V(i-1)=Aeroelastic{i,1};
    for j=3:14
        Naturalfreq.aeroelastic(i-1,j-2)=Aeroelastic{i,j};
        Naturalfreq.structural(i-1,j-2)=Structural{i,j};
    end
end

for i=2:23
    for j=16:27
        Damping.aeroelastic(i-1,j-15)=Aeroelastic{i,j};
        Damping.structural(i-1,j-15)=Structural{i,j};
    end
end

vars={'-o','-+','-*','-.','-x','-_','-|','-s','-d','-^','-v','->','-<','-p','-h'};
figure
for i=1:length(Damping.structural(1,:))
    plot(V,Naturalfreq.aeroelastic(:,i),vars{i})
    
    hold on
end
plot(V,Freq1P)
hold on
plot(V,Freq3P)
hold on
plot(V,Freq6P)
grid on
grid minor
xlabel('V(m/s)');
ylabel('Nat. frequency(Hz)');
title('Campbell diagram for aeroelastic mode shapes(new rotor)');
legend('1st Twr FA','1st Twr SS','1st BW flap','1st SYM flap','1st FW flap','1st BW edge','1st FW edge','2nd BW flap','2nd FW flap','2nd SYM flap','1st COL edge','3rd BW flap','1P','3P','6P','Location','northwest','NumColumns',1);

figure
for i=1:length(Damping.structural(1,:))
    plot(V,Naturalfreq.structural(:,i),vars{i})
    
    hold on
end
grid on
grid minor
xlabel('V(m/s)');
ylabel('Nat. frequency(Hz)');
title('Campbell diagram for structural mode shapes(new rotor)');
legend('1st Twr FA','1st Twr SS','1st BW flap','1st SYM flap','1st FW flap','1st BW edge','1st FW edge','2nd BW flap','2nd FW flap','2nd SYM flap','1st COL edge','3rd BW flap','Location','northwest','NumColumns',1);

figure
for i=1:length(Damping.structural(1,:))
    plot(V,Damping.structural(:,i),vars{i})
    
    hold on
end
grid on
grid minor
xlabel('V(m/s)');
ylabel('Damping[% critical]');
title('Damping values for structural mode shapes(new rotor)');
legend('1st Twr FA','1st Twr SS','1st BW flap','1st SYM flap','1st FW flap','1st BW edge','1st FW edge','2nd BW flap','2nd FW flap','2nd SYM flap','1st COL edge','3rd BW flap','Location','northwest','NumColumns',1);

figure
for i=1:length(Damping.structural(1,:))
    plot(V,Damping.aeroelastic(:,i),vars{i})
    
    hold on
end
grid on
grid minor
xlabel('V(m/s)');
ylabel('Damping[% critical]');
title('Damping values for aeroelastic mode shapes(new rotor)');
legend('1st Twr FA','1st Twr SS','1st BW flap','1st SYM flap','1st FW flap','1st BW edge','1st FW edge','2nd BW flap','2nd FW flap','2nd SYM flap','1st COL edge','3rd BW flap','Location','northwest','NumColumns',1);

%% Post processing DTU 10 MW

for i=2:23
    DTU10MW.V(i-1)=DTU10MW.Aeroelastic{i,1};
    for j=3:14
        DTU10MW.Naturalfreq.aeroelastic(i-1,j-2)=DTU10MW.Aeroelastic{i,j};
        DTU10MW.Naturalfreq.structural(i-1,j-2)=DTU10MW.Structural{i,j};
    end
end

for i=2:23
    for j=16:27
        DTU10MW.Damping.aeroelastic(i-1,j-15)=DTU10MW.Aeroelastic{i,j};
        DTU10MW.Damping.structural(i-1,j-15)=DTU10MW.Structural{i,j};
    end
end

vars={'-o','-+','-*','-.','-x','-_','-|','-s','-d','-^','-v','->','-<','-p','-h'};

figure
for i=1:length(DTU10MW.Damping.structural(1,:))
    plot(DTU10MW.V,DTU10MW.Naturalfreq.aeroelastic(:,i),vars{i})
    
    hold on
end
plot(V,DTU10MW.Freq1P)
hold on
plot(V,DTU10MW.Freq3P)
hold on
plot(V,DTU10MW.Freq6P)
grid on
grid minor
xlabel('V(m/s)');
ylabel('Nat. frequency(Hz)');
title('Campbell diagram for aeroelastic mode shapes(DTU 10MW)');
legend('1st Twr FA','1st Twr SS','1st BW flap','1st SYM flap','1st FW flap','1st BW edge','1st FW edge','2nd BW flap','2nd FW flap','2nd SYM flap','1st COL edge','3rd BW flap','1P','3P','6P','Location','northwest','NumColumns',1);

figure
for i=1:length(DTU10MW.Damping.structural(1,:))
    plot(DTU10MW.V,DTU10MW.Naturalfreq.structural(:,i),vars{i})
    
    hold on
end
grid on
grid minor
xlabel('V(m/s)');
ylabel('Nat. frequency(Hz)');
title('Campbell diagram for structural mode shapes(DTU 10MW)');
legend('1st Twr FA','1st Twr SS','1st BW flap','1st SYM flap','1st FW flap','1st BW edge','1st FW edge','2nd BW flap','2nd FW flap','2nd SYM flap','1st COL edge','3rd BW flap','Location','northwest','NumColumns',1);

figure
for i=1:length(DTU10MW.Damping.structural(1,:))
    plot(DTU10MW.V,DTU10MW.Damping.structural(:,i),vars{i})
    
    hold on
end
grid on
grid minor
xlabel('V(m/s)');
ylabel('Damping[% critical]');
title('Damping values for structural mode shapes(DTU 10MW)');
legend('1st Twr FA','1st Twr SS','1st BW flap','1st SYM flap','1st FW flap','1st BW edge','1st FW edge','2nd BW flap','2nd FW flap','2nd SYM flap','1st COL edge','3rd BW flap','Location','northwest','NumColumns',1);

figure
for i=1:length(DTU10MW.Damping.structural(1,:))
    plot(DTU10MW.V,DTU10MW.Damping.aeroelastic(:,i),vars{i})
    
    hold on
end
grid on
grid minor
xlabel('V(m/s)');
ylabel('Damping[% critical]');
title('Damping values for aeroelastic mode shapes(DTU 10MW)');
legend('1st Twr FA','1st Twr SS','1st BW flap','1st SYM flap','1st FW flap','1st BW edge','1st FW edge','2nd BW flap','2nd FW flap','2nd SYM flap','1st COL edge','3rd BW flap','Location','northwest','NumColumns',1);

%% Table creation for the new rotor
[row,col]=find(V==14);

Tower_signal=Naturalfreq.aeroelastic(col,1:7);
Tower_signal(end+1)=Naturalfreq.aeroelastic(col,11);
N=7.4322;
omega=pi.*N./30;
frequency=omega./(2*pi);

%Blade_signal=[Naturalfreq.aeroelastic(col,1)+frequency Naturalfreq.aeroelastic(col,1)-frequency 
    
% for i=1:(2*length(Tower_signal)-2)
% 
% if i==1 
%     Blade_signal(i)=Tower_signal(i)+frequency; % 1st Twr FA
%     Blade_signal(i+1)=Tower_signal(i)-frequency;
% elseif i==3
%     Blade_signal(i)=Tower_signal(i-1)+frequency; % 1st Twr SS
%     Blade_signal(i+1)=Tower_signal(i-1)-frequency;
% 
% elseif i==5
%     Blade_signal(i)=Tower_signal(i-2)+frequency; % 1st BW flap
%     Blade_signal(i+1)=Tower_signal(i-2)-frequency;
% elseif i==7
%     Blade_signal(i)=Tower_signal(i-3); % 1st SYM flap
% 
% elseif i==8
%     Blade_signal(i)=Tower_signal(i-3)+frequency; %1st FW flap
%     Blade_signal(i+1)=Tower_signal(i-3)-frequency;
%     
% elseif i==10
%     Blade_signal(i)=Tower_signal(i-4)+frequency; %1st BW edge
%     Blade_signal(i+1)=Tower_signal(i-4)-frequency;
%     
% elseif i==12
%     Blade_signal(i)=Tower_signal(i-5)+frequency; %1st FW edge
%     Blade_signal(i+1)=Tower_signal(i-5)-frequency;
%     
% elseif i==14
%     Blade_signal(i)=Tower_signal(i-6); %1st COL edge
% end
% end

for i=1:(length(Tower_signal)-1)
    Blade_signal(2*(i-1)+1)=Tower_signal(i)+frequency;
    Blade_signal(2*(i-1)+2)=Tower_signal(i)-frequency;
end
Blade_signal(7)=Tower_signal(4);
Blade_signal(14)=Tower_signal(8);


Blade_signal=Blade_signal';
Tower_signal=Tower_signal';


    
