clear all
clc
close all

%% Save plots parameter
save_figs=0; %save_figs=0=0: Save plots with title, save_figs=0=1: Save plots without title

%% Outputs from HAWC2Stab2.exe for new rotor

Aeroelastic=readtable('redesign_ael.cmb','FileType','text'); 
Structural=readtable('redesign_st.cmb','FileType','text');
Aeroelastic.Var2=[];
Natfreq.aeroelastic=Aeroelastic(:,1:13);
Structural.Var2=[];
Natfreq.structural=Structural(:,1:13);

Damp.aeroelastic(:,1)=Aeroelastic(:,1);
Damp.aeroelastic(:,2:13)=Aeroelastic(:,15:26);
Damp.structural(:,1)=Structural(:,1);
Damp.structural(:,2:13)=Structural(:,15:26);

Natfreq.aeroelastic.Properties.VariableNames = {'V [m/s]' '1st Twr FA' '1st Twr SS' '1st BW flap' '1st SYM flap' '1st FW flap' '1st BW edge' '1st FW edge' '2nd BW flap' '2nd FW flap' '2nd SYM flap' '1st COL edge' '3rd BW flap'};
Natfreq.structural.Properties.VariableNames = {'V [m/s]' '1st Twr FA' '1st Twr SS' '1st BW flap' '1st SYM flap' '1st FW flap' '1st BW edge' '1st FW edge' '2nd BW flap' '2nd FW flap' '2nd SYM flap' '1st COL edge' '3rd BW flap'};
Damp.aeroelastic.Properties.VariableNames = {'V [m/s]' '1st Twr FA' '1st Twr SS' '1st BW flap' '1st SYM flap' '1st FW flap' '1st BW edge' '1st FW edge' '2nd BW flap' '2nd FW flap' '2nd SYM flap' '1st COL edge' '3rd BW flap'};
Damp.structural.Properties.VariableNames = {'V [m/s]' '1st Twr FA' '1st Twr SS' '1st BW flap' '1st SYM flap' '1st FW flap' '1st BW edge' '1st FW edge' '2nd BW flap' '2nd FW flap' '2nd SYM flap' '1st COL edge' '3rd BW flap'};

var1 = table2array(Natfreq.aeroelastic(1,2:13));
var2 = table2array(Natfreq.structural(1,2:13));



[~,idx1] = sort(var1);
[~,idx2] = sort(var2);


idx1 = [1, idx1+1, 15:size(Natfreq.aeroelastic,2)];
idx2 = [1, idx2+1, 15:size(Natfreq.structural,2)];

Natfreq.aeroelastic_sorted = Natfreq.aeroelastic(:,idx1);
Natfreq.structural_sorted = Natfreq.structural(:,idx2);

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


DTU10MW.Structural=readtable('DTU_10MW_st.cmb','FileType','text');
DTU10MW.Aeroelastic=readtable('DTU_10MW_ael.cmb','FileType','text');
DTU10MW.Aeroelastic.Var2=[];
DTU10MW.Natfreq.aeroelastic=DTU10MW.Aeroelastic(:,1:13);
DTU10MW.Structural.Var2=[];
DTU10MW.Natfreq.structural=DTU10MW.Structural(:,1:13);

DTU10MW.Damp.aeroelastic(:,1)=DTU10MW.Aeroelastic(:,1);
DTU10MW.Damp.aeroelastic(:,2:13)=DTU10MW.Aeroelastic(:,15:26);
DTU10MW.Damp.structural(:,1)=DTU10MW.Structural(:,1);
DTU10MW.Damp.structural(:,2:13)=DTU10MW.Structural(:,15:26);

DTU10MW.Natfreq.aeroelastic.Properties.VariableNames = {'V [m/s]' '1st Twr FA' '1st Twr SS' '1st BW flap' '1st SYM flap' '1st FW flap' '1st BW edge' '1st FW edge' '2nd BW flap' '2nd FW flap' '1st COL edge' '2nd SYM flap' '3rd BW flap'};
DTU10MW.Natfreq.structural.Properties.VariableNames = {'V [m/s]' '1st Twr FA' '1st Twr SS' '1st BW flap' '1st SYM flap' '1st FW flap' '1st BW edge' '1st FW edge' '2nd BW flap' '2nd FW flap' '2nd SYM flap' '1st COL edge' '3rd BW flap'};
DTU10MW.Damp.aeroelastic.Properties.VariableNames = {'V [m/s]' '1st Twr FA' '1st Twr SS' '1st BW flap' '1st SYM flap' '1st FW flap' '1st BW edge' '1st FW edge' '2nd BW flap' '2nd FW flap' '2nd SYM flap' '1st COL edge' '3rd BW flap'};
DTU10MW.Damp.structural.Properties.VariableNames = {'V [m/s]' '1st Twr FA' '1st Twr SS' '1st BW flap' '1st SYM flap' '1st FW flap' '1st BW edge' '1st FW edge' '2nd BW flap' '2nd FW flap' '2nd SYM flap' '1st COL edge' '3rd BW flap'};

DTU10MW.Damp.structural(2,1)={5};
DTU10MW.Damp.structural(7,1)={10};
DTU10MW.Damp.structural(8,1)={11};

%% DTU 10MW reordering



var3 = table2array(DTU10MW.Natfreq.aeroelastic(1,2:13));
var4=table2array(DTU10MW.Natfreq.structural(1,2:13));



[~,idx3] = sort(var3);
[~,idx4] = sort(var4);


idx3 = [1, idx3+1, 15:size(DTU10MW.Natfreq.aeroelastic,2)];
idx4 = [1, idx4+1, 15:size(DTU10MW.Natfreq.structural,2)];



DTU10MW.Natfreq.aeroelastic_sorted = DTU10MW.Natfreq.aeroelastic(:,idx3);
DTU10MW.Natfreq.structural_sorted = DTU10MW.Natfreq.structural(:,idx4);

DTU10MW.Natfreq.structural_sorted(2,1)={5};
DTU10MW.Natfreq.structural_sorted(7,1)={10};
DTU10MW.Natfreq.structural_sorted(8,1)={11};


%% Swaping for 1st COL edge and 2nd SYM flap
[row1,col1]=find(DTU10MW.Natfreq.structural_sorted.('V [m/s]')==13);

DTU10MW.save_flap=DTU10MW.Natfreq.structural_sorted.('2nd SYM flap')(row1:end);
DTU10MW.save_edge=DTU10MW.Natfreq.structural_sorted.('1st COL edge')(row1:end);

DTU10MW.Natfreq.structural_sorted.('2nd SYM flap')(row1:end)=DTU10MW.save_edge;
DTU10MW.Natfreq.structural_sorted.('1st COL edge')(row1:end)=DTU10MW.save_flap;

[row2,col2]=find(Natfreq.structural_sorted.('V [m/s]')==13);

save_flap=Natfreq.structural_sorted.('2nd SYM flap')(row2:end);
save_edge=Natfreq.structural_sorted.('1st COL edge')(row2:end);

Natfreq.structural_sorted.('2nd SYM flap')(row2:end)=save_edge;
Natfreq.structural_sorted.('1st COL edge')(row2:end)=save_flap;

%% DTU 10 MW plots

if save_figs==0


vars={'-o','-+','-*','-s','-x','-o','-|','-s','-d','-^','-v','->','-<','-p','-h'};

figure
for i=2:width(DTU10MW.Natfreq.aeroelastic(1,:))
   plot(DTU10MW.Natfreq.aeroelastic_sorted.('V [m/s]'),DTU10MW.Natfreq.aeroelastic_sorted.(DTU10MW.Natfreq.aeroelastic_sorted.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end
plot(DTU10MW.Natfreq.aeroelastic_sorted.('V [m/s]'),DTU10MW.Freq1P,'k--','LineWidth',2.0)
hold on
plot(DTU10MW.Natfreq.aeroelastic_sorted.('V [m/s]'),DTU10MW.Freq3P,'k-.','LineWidth',2.0)
hold on
plot(DTU10MW.Natfreq.aeroelastic_sorted.('V [m/s]'),DTU10MW.Freq6P,'k:','LineWidth',2.0)
set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Natural frequency [Hz]');
title('Campbell diagram for aeroelastic mode shapes(DTU 10MW)');
%legend('1st Twr FA','1st Twr SS','1st BW flap','1st SYM flap','1st FW flap','1st BW edge','1st FW edge','2nd BW flap','2nd FW flap','1st COL edge','2nd SYM flap','3rd BW flap','1P','3P','6P','Location','eastoutside','NumColumns',1);
legend(DTU10MW.Natfreq.aeroelastic_sorted.Properties.VariableNames{2:end},'1P','3P','6P','Location','eastoutside','NumColumns',1)



figure
for i=2:width(DTU10MW.Natfreq.structural(1,:))
   plot(DTU10MW.Natfreq.structural_sorted.('V [m/s]'),DTU10MW.Natfreq.structural_sorted.(DTU10MW.Natfreq.structural_sorted.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end
set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Natural frequency [Hz]');
title('Campbell diagram for structural mode shapes(DTU 10MW)');
legend(DTU10MW.Natfreq.structural_sorted.Properties.VariableNames{2:end},'Location','eastoutside','NumColumns',1)

figure
for i=2:width(DTU10MW.Damp.structural(1,:))
   plot(DTU10MW.Damp.structural.('V [m/s]'),DTU10MW.Damp.structural.(DTU10MW.Damp.structural.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end
set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Damping[% critical]');
title('Damping values for structural mode shapes(DTU 10MW)');
legend(DTU10MW.Damp.structural.Properties.VariableNames{2:end},'Location','eastoutside','NumColumns',1)


figure
for i=2:width(DTU10MW.Damp.aeroelastic(1,:))
   plot(DTU10MW.Damp.aeroelastic.('V [m/s]'),DTU10MW.Damp.aeroelastic.(DTU10MW.Damp.aeroelastic.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end
set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Damping[% critical]');
title('Damping values for aeroelastic mode shapes(DTU 10MW)');
legend(DTU10MW.Damp.aeroelastic.Properties.VariableNames{2:end},'Location','eastoutside','NumColumns',1)

%% New rotor plots
figure
for i=2:width(Natfreq.aeroelastic_sorted(1,:))
   plot(Natfreq.aeroelastic_sorted.('V [m/s]'),Natfreq.aeroelastic_sorted.(Natfreq.aeroelastic_sorted.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end
plot(Natfreq.aeroelastic_sorted.('V [m/s]'),Freq1P,'k--','LineWidth',2.0)
hold on
plot(Natfreq.aeroelastic_sorted.('V [m/s]'),Freq3P,'k-.','LineWidth',2.0)
hold on
plot(Natfreq.aeroelastic_sorted.('V [m/s]'),Freq6P,'k:','LineWidth',2.0)
set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Natural frequency [Hz]');
title('Campbell diagram for aeroelastic mode shapes(Redesign)');
legend(Natfreq.aeroelastic_sorted.Properties.VariableNames{2:end},'1P','3P','6P','Location','eastoutside','NumColumns',1)


figure
for i=2:width(Natfreq.structural_sorted(1,:))
   plot(Natfreq.structural_sorted.('V [m/s]'),Natfreq.structural_sorted.(Natfreq.structural_sorted.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end

set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Natural frequency [Hz]');
title('Campbell diagram for structural mode shapes(Redesign)');
legend(Natfreq.structural_sorted.Properties.VariableNames{2:end},'Location','eastoutside','NumColumns',1)



figure
for i=2:width(Damp.aeroelastic(1,:))
   plot(Damp.aeroelastic.('V [m/s]'),Damp.aeroelastic.(Damp.aeroelastic.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end

set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Damping[% critical]');
title('Damping values for aeroelastic mode shapes(Redesign)');
legend(Damp.aeroelastic.Properties.VariableNames{2:end},'Location','eastoutside','NumColumns',1)



figure
for i=2:width(Damp.structural(1,:))
   plot(Damp.structural.('V [m/s]'),Damp.structural.(Damp.structural.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end

set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Damping[% critical]');
title('Damping values for structural mode shapes(Redesign)');
legend(Damp.structural.Properties.VariableNames{2:end},'Location','eastoutside','NumColumns',1)

elseif save_figs==1
    
    vars={'-o','-+','-*','-s','-x','-o','-|','-s','-d','-^','-v','->','-<','-p','-h'};
    
    %% DTU 10 MW plots

figure
for i=2:width(DTU10MW.Natfreq.aeroelastic(1,:))
   plot(DTU10MW.Natfreq.aeroelastic_sorted.('V [m/s]'),DTU10MW.Natfreq.aeroelastic_sorted.(DTU10MW.Natfreq.aeroelastic_sorted.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end
plot(DTU10MW.Natfreq.aeroelastic_sorted.('V [m/s]'),DTU10MW.Freq1P,'k--','LineWidth',2.0)
hold on
plot(DTU10MW.Natfreq.aeroelastic_sorted.('V [m/s]'),DTU10MW.Freq3P,'k-.','LineWidth',2.0)
hold on
plot(DTU10MW.Natfreq.aeroelastic_sorted.('V [m/s]'),DTU10MW.Freq6P,'k:','LineWidth',2.0)
set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Natural frequency [Hz]');
%title('Campbell diagram for aeroelastic mode shapes(DTU 10MW)');
%legend('1st Twr FA','1st Twr SS','1st BW flap','1st SYM flap','1st FW flap','1st BW edge','1st FW edge','2nd BW flap','2nd FW flap','1st COL edge','2nd SYM flap','3rd BW flap','1P','3P','6P','Location','eastoutside','NumColumns',1);
legend(DTU10MW.Natfreq.aeroelastic_sorted.Properties.VariableNames{2:end},'1P','3P','6P','Location','eastoutside','NumColumns',1)



figure
for i=2:width(DTU10MW.Natfreq.structural(1,:))
   plot(DTU10MW.Natfreq.structural_sorted.('V [m/s]'),DTU10MW.Natfreq.structural_sorted.(DTU10MW.Natfreq.structural_sorted.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end
set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Natural frequency [Hz]');
%title('Campbell diagram for structural mode shapes(DTU 10MW)');
legend(DTU10MW.Natfreq.structural_sorted.Properties.VariableNames{2:end},'Location','eastoutside','NumColumns',1)

figure
for i=2:width(DTU10MW.Damp.structural(1,:))
   plot(DTU10MW.Damp.structural.('V [m/s]'),DTU10MW.Damp.structural.(DTU10MW.Damp.structural.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end
set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Damping[% critical]');
%title('Damping values for structural mode shapes(DTU 10MW)');
legend(DTU10MW.Damp.structural.Properties.VariableNames{2:end},'Location','eastoutside','NumColumns',1)


figure
for i=2:width(DTU10MW.Damp.aeroelastic(1,:))
   plot(DTU10MW.Damp.aeroelastic.('V [m/s]'),DTU10MW.Damp.aeroelastic.(DTU10MW.Damp.aeroelastic.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end
set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Damping[% critical]');
%title('Damping values for aeroelastic mode shapes(DTU 10MW)');
legend(DTU10MW.Damp.aeroelastic.Properties.VariableNames{2:end},'Location','eastoutside','NumColumns',1)

%% New rotor plots
figure
for i=2:width(Natfreq.aeroelastic_sorted(1,:))
   plot(Natfreq.aeroelastic_sorted.('V [m/s]'),Natfreq.aeroelastic_sorted.(Natfreq.aeroelastic_sorted.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end
plot(Natfreq.aeroelastic_sorted.('V [m/s]'),Freq1P,'k--','LineWidth',2.0)
hold on
plot(Natfreq.aeroelastic_sorted.('V [m/s]'),Freq3P,'k-.','LineWidth',2.0)
hold on
plot(Natfreq.aeroelastic_sorted.('V [m/s]'),Freq6P,'k:','LineWidth',2.0)
set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Natural frequency [Hz]');
%title('Campbell diagram for aeroelastic mode shapes(Redesign)');
legend(Natfreq.aeroelastic_sorted.Properties.VariableNames{2:end},'1P','3P','6P','Location','eastoutside','NumColumns',1)


figure
for i=2:width(Natfreq.structural_sorted(1,:))
   plot(Natfreq.structural_sorted.('V [m/s]'),Natfreq.structural_sorted.(Natfreq.structural_sorted.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end

set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Natural frequency [Hz]');
%title('Campbell diagram for structural mode shapes(Redesign)');
legend(Natfreq.structural_sorted.Properties.VariableNames{2:end},'Location','eastoutside','NumColumns',1)



figure
for i=2:width(Damp.aeroelastic(1,:))
   plot(Damp.aeroelastic.('V [m/s]'),Damp.aeroelastic.(Damp.aeroelastic.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end

set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Damping[% critical]');
%title('Damping values for aeroelastic mode shapes(Redesign)');
legend(Damp.aeroelastic.Properties.VariableNames{2:end},'Location','eastoutside','NumColumns',1)



figure
for i=2:width(Damp.structural(1,:))
   plot(Damp.structural.('V [m/s]'),Damp.structural.(Damp.structural.Properties.VariableNames{i}),vars{i},'LineWidth',2.0)
    
    hold on
end

set(gca,'FontSize',20)
xlim([4 25]);
grid on
grid minor
xlabel('Wind speed [m/s]');
ylabel('Damping[% critical]');
%title('Damping values for structural mode shapes(Redesign)');
legend(Damp.structural.Properties.VariableNames{2:end},'Location','eastoutside','NumColumns',1)

end

%% Table creation for the new rotor
[row,col]=find(Natfreq.aeroelastic_sorted.('V [m/s]')'==14);

Natfreq_aeroelastic_sort=table2array(Natfreq.aeroelastic_sorted); 
Tower_signal=Natfreq_aeroelastic_sort(col,2:8);
Tower_signal(end+1)=Natfreq_aeroelastic_sort(col,12);
N=pwr{1,2}.('Speed')(end);
omega=pi.*N./30;
frequency=omega./(2*pi);

for i=1:(length(Tower_signal)-1)
    Blade_signal(2*(i-1)+1)=Tower_signal(1,i)+frequency;
    Blade_signal(2*(i-1)+2)=Tower_signal(1,i)-frequency;
end
Blade_signal(7)=Tower_signal(4);
Blade_signal(14)=Tower_signal(8);


Blade_signal=Blade_signal';
Tower_signal=Tower_signal';