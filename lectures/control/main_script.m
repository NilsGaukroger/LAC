%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%             From here on you don't need to change anything!
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% -------------------------------------------------------------------------
%                   Wind Speed File replacing!
% -------------------------------------------------------------------------

pathCell    =   regexp(path, pathsep, 'split'); 
onPath      =   any(strcmpi([cd, '/lib'], pathCell));
if ~onPath;     addpath([cd, '/lib']);  end
addpath(genpath(pwd));

SimParams.Ts                = 0.1;          % Controller sampling time in seconds (Normally you don't need to change this!)

MakeWSP;
InitWTModel;
PIControllerInit;

% -------------------------------------------------------------------------
%                       Initialize variables
% -------------------------------------------------------------------------
% keyboard;
OmegaInit   = 1;
x           = zeros(size(A,1),NSim);
y           = zeros(size(C,1),NSim);
x(1:2)      = OmegaInit;%[OmegaInit; OmegaInit; 0; 0; 0];
u           = zeros(2,NSim);


for k=1:NSim
    [x(:,k+1),y(:,k)]   = WT_Nonlinear(x(:,k),u(:,k), wsp(k),WT);
    u(:,k+1)            = PI_Controller(y(1,k));    
    counterfunc(t(k));
end


%%

StartTime = 0;
nStart  = find(t>StartTime,1);
Time    = t(nStart:NSim)-StartTime;

Vars.OmegaR   = x(1,nStart:NSim);

switch WT.Model
    case 'WT0'
        Vars.OmegaG   = 50*Vars.OmegaR;
        Vars.Vt       = zeros(size(Vars.OmegaR));
    case 'WT1'
        Vars.OmegaG   = WT.Params.GEARBOX_RATIO*x(2,nStart:NSim);
        Vars.Vt       = zeros(size(Vars.OmegaR));
    case 'WT2'
        Vars.OmegaG   = WT.Params.GEARBOX_RATIO*x(2,nStart:NSim);
        Vars.Vt       = x(5,nStart:NSim);
    otherwise
        error('Wrong model type!');
end

Vars.BldPitch = u(1,nStart:NSim);
Vars.GenTq    = u(2,nStart:NSim)*1e-3;

Vars.WSP      = wsp(1,nStart:NSim);

Vars.Pe       = Vars.OmegaG.*Vars.GenTq*1e-3/WT.Params.GEARBOX_RATIO;
Vars.Cp       = (Vars.Pe*1e6/WT.Params.Efficiency)./(1/2*WT.Params.Air_density*WT.Params.Rotor_Swept_Area*Vars.WSP.^3);

VarNames    = fieldnames(Vars);
titles      = {'Rotor speed [rad/s]','Generator speed transfered to LSS [rad/s]','Tower fore-aft velocity [m/s]',...
    'Blade pitch [degrees]','Generator torque [KNM]','Wind Speed [m/s] ','Generated power [MW]','Cp [-]'};

%%
color_vector = rand(1,3)*0.8;

for i=1:length(VarNames)
    figure(i);
    hold on;
    a = plot(Time,Vars.(VarNames{i}),'-');
    set(a,'color',color_vector,'LineWidth',2);
    grid on; box on;
    title(titles{i});
    xlabel('Time [s]')
end

figure(1);
v = axis;
plot(v(1:2),[1.005 1.005],'r--','LineWidth',2)
%%
plot_on_cp = 1;
if plot_on_cp == 1
    Lambda = WT.Params.Rotor_Radius*Vars.OmegaR./Vars.WSP;
    figure(11);
    hold on;
    contour(WT.CPCTTable.Pitch,WT.CPCTTable.Lambda,max(0,WT.CPCTTable.CP),30);
    a = plot(Vars.BldPitch,Lambda,'.');
    set(a,'color',color_vector);
    v = axis;
    axis([v(1:2) 2 18])
end

xlabel('Pitch angle [degrees]');
ylabel('Tip speed ratio [-]');

