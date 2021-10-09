close all; clear variables; clc

table = cell(2,2);
table{1,1} = readtable('your_model\results_dtu10mw\stab\DTU_10MW_st.cmb','FileType','text');
table{1,2} = readtable('your_model\results_dtu10mw\stab\DTU_10MW_ael.cmb','FileType','text');
table{2,1} = readtable('your_model\results_redesign\stab\redesign_st.cmb','FileType','text');
table{2,2} = readtable('your_model\results_redesign\stab\redesign_ael.cmb','FileType','text');

var = table2array(table{1,2}(1,2:14))

[~,idx] = sort(var)

idx = [1, idx+1, 15:size(table{1,2},2)]

sorted = table{1,2}(:,idx)
% var(idx)