%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
tic
res={};
model=modcompile(strcat(pwd,'/Cth_model_cofac_new_WT.xlsx'),strcat(pwd,'/ctherm_mechanism_20_jan_new.xlsx'),strcat(pwd,'/Ctherm_data_all_v3.xlsx'));
model.d.vpert(19,2)=0.5;
model.d.vpert(39,2)=0.5;
model.d.vpert(19,3)=0.2;
model.d.vpert(39,3)=0.2;
model.d.vpert(26,3)=1;
model.d.vpert(12,3)=5;
model.d.vpert(19,4)=0.5;
model.d.vpert(39,4)=0.5;
model.d.vpert(19,5)=1;
model.d.vpert(43,5)=0.1;
model.d.vpert(43,9)=0.5;
%model.options.multistarts=3;
res=kineticestimate_v4(model,res);
tout = toc;

save('optimal_1_4_19_wt_all_3.mat','res','model','tout');

