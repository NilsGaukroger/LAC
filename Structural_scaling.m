% import numpy as np
% #
% # Output file name
%f_new_st = 'st_new.dat'
% # ============================================================================
%# INPUTS
%# ============================================================================
clear
clc
close all
f_original_st = importdata('scaling_inputs.txt');



initial.r=f_original_st(:,1);
initial.m=f_original_st(:,2);
initial.x_cg=f_original_st(:,3);
initial.y_cg=f_original_st(:,4);
initial.ri_x=f_original_st(:,5);
initial.ri_y=f_original_st(:,6);
initial.x_sh=f_original_st(:,7);
initial.y_sh=f_original_st(:,8);
initial.E=f_original_st(:,9);
initial.G=f_original_st(:,10);
initial.I_x=f_original_st(:,11);
initial.I_y=f_original_st(:,12);
initial.I_p=f_original_st(:,13);
initial.k_x=f_original_st(:,14);
initial.k_y=f_original_st(:,15);
initial.A=f_original_st(:,16);
initial.pitch=f_original_st(:,17);
initial.x_e=f_original_st(:,18);
initial.y_e=f_original_st(:,19);

r_tip_10MW=89.1660;
r_tip_newrotor=97.7893;
r_root=2.8;

Blade_length_10MW=r_tip_10MW-r_root;
Blade_length_newrotor=r_tip_newrotor-r_root;
s=Blade_length_newrotor./Blade_length_10MW;

%% s^1

new.r=s.*initial.r;
new.x_cg=s.*initial.x_cg;
new.y_cg=s.*initial.y_cg;
new.ri_x=s.*initial.ri_x;
new.ri_y=s.*initial.ri_y;
new.x_sh=s.*initial.x_sh;
new.y_sh=s.*initial.y_sh;
new.x_e=s.*initial.x_e;
new.y_e=s.*initial.y_e;


%% s^2

new.m=(s.^2).*initial.m;
new.A=(s.^2).*initial.A;

%% s^4

new.I_x=(s.^4).*initial.I_x;
new.I_y=(s.^4).*initial.I_y;
new.I_p=(s.^4).*initial.I_p;

%% st_new.dat update

f_new_st=[new.r new.m new.x_cg new.y_cg new.ri_x new.ri_y new.x_sh new.y_sh initial.E initial.G new.I_x new.I_y new.I_p initial.k_x initial.k_y new.A initial.pitch new.x_e new.y_e];

C={1, 'number of sets,', 'Nset', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',' ',' ',' ';'-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------',' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',' ',' ',' ',' ',' '; '#1', 'user_mads', 'generated blade',' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',' ',' ',' '; 'r', 'm', 'x_cg', 'y_cg', 'ri_x', 'ri_y', 'x_sh', 'y_sh', 'E', 'G', 'I_x', 'I_y', 'I_p', 'k_x', 'k_y', 'A', 'pitch', 'x_e', 'y_e'; '$1', 51, ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',' ',' ',' ',' '};
for i=1:length(f_new_st(:,1))
    for j=1:length(f_new_st(1,:))
          C{i+5,j}=f_new_st(i,j);
    end
end
writecell(C,'st_new.dat','Delimiter','tab')

%% c2new.dat update

f_original_c2 = importdata('c2_inputs.txt');
HAWCinp=load('aero_design.mat');

initial.Nsec=f_original_c2(:,1);
initial.x_pos=f_original_c2(:,2);
initial.y_pos=f_original_c2(:,3);
initial.z_pos=f_original_c2(:,4);
initial.twist=f_original_c2(:,5);

new.x_pos=s.*initial.x_pos;
new.y_pos=s.*initial.y_pos;
new.z_pos=s.*initial.z_pos;
new.twist=-rad2deg(HAWCinp.HAWC_in.beta);



c2new=[initial.Nsec new.x_pos new.y_pos new.z_pos new.twist];

for i=1:length(c2new(:,1))
    D{i,1}='sec';
    D{i,2}=i;
    for j=3:(length(c2new(1,:))+1)
      D{i,j}=c2new(i,j-1);
    end
end

writecell(D,'c2_new.dat','Delimiter','tab');

        
