function [u] = import_u(filepath,file,idx)
names  = {'_u','_d','_f'};
output = cell(length(names),1);

for i = 1:3
    output{i} = dir(strcat(filepath, file, names{i}, '*'));
end
filename = {output{1}.name};

% Parse results as tables
u_vars   = ["s","A","AP","PHI0","ALPHA0","U0","FX0","FY0","M0","UX0","UY0","UZ0",...
    "Twist","X_AC0","Y_AC0","Z_AC0","CL0","CD0","CM0","CLp0","CDp0","CMp0",...
    "F0","F'","CL_FS0","CLFS'","V_a","V_t","Tors.","vx","vy","chord",...
    "CT","CP","angle","v_1","v_2","v_3"];

u = readtable(strcat(filepath, filename{idx}),'FileType','text');
u.Properties.VariableNames = u_vars;
end