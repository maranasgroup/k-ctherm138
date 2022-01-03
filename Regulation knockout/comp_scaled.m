%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

clear, clc

%res={};
load optimal_reinit1_06_21_19_CF_83.mat


model=modcompile(strcat(pwd,'/Cth_model_cofac_new_WT.xlsx'),strcat(pwd,'/ctherm_mechanism_20_jan_125.xlsx'),strcat(pwd,'/Ctherm_data_all_v3.xlsx'));
load xbest_45.mat
%model.d.vpert(24,3)=1;
x = res.reinit_data;

res = compileresult(x,model)


save('optimal_corrected_83.mat','model','res')
%[k,~,~,~] = calc_k_act( model,x );
%updated_rxn = 35;
