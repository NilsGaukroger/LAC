function [camp,damp,realPart] = import_cmb(n_modes,filename,VariableNames,swap)
%% Import file
cmb = readtable(filename,'FileType','text');

%% Correct wind speed if required
if any(cmb.(1)) % if any wind speed == 0
    cmb(:,1) = num2cell(linspace(cmb.(1)(1),cmb.(1)(end),length(cmb.(1)))');
end

%% Make separate tables for frequencies, damping and real parts
% Extract frequencies for Campbell diagram (ignoring rigid body mode)
out.camp = cmb(:,[1, 2+1:2+(n_modes)]);
% Extract damping values
out.damp = cmb(:,[1, 4+n_modes:3+(n_modes*2)]);
% Extract real parts of aeroelastic eigenvalues
if size(cmb,2) == 1+((n_modes+1)*3) % if this is the aeroelastic analysis
    out.realPart = cmb(:,[1, 5+(n_modes*2):4+(n_modes*3)]);
end
fn = fieldnames(out);
for i = 1:length(fn)
    % Apply column headings
    out.(fn{i}).Properties.VariableNames = [ {'Wind Speed'}, VariableNames];
end

%% Organise tables in order of increasing frequency
% Find order of frequencies
[~,idx] = sort(table2array(out.camp(1,2:end)));
idx = [1, idx+1]; % don't sort the wind speed column
for i = 1:length(fn)
    out.(fn{i}) = out.(fn{i})(:,idx);
end

%% Mode 10 & 11 correction
if swap
    fprintf('Correction\n')
end

%% Create outputs
camp = out.camp;
damp = out.damp;
if size(cmb,2) == 1+((n_modes+1)*3) % if this is the aeroelastic analysis
    realPart = out.realPart;
end
end