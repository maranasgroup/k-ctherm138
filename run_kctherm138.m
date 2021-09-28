tic
res={};
model=modcompile(strcat(pwd,'/kctherm138_model.xlsx'),strcat(pwd,'/kctherm138_mechanism.xlsx'),strcat(pwd,'/kctherm138_data.xlsx'));
model.d.vpert(19,2)=0.5;
model.d.vpert(39,2)=0.5;
model.d.vpert(19,3)=0.2;
model.d.vpert(39,3)=0.2;
model.d.vpert(24,3)=1;
model.d.vpert(12,3)=5;
model.d.vpert(19,4)=0.5;
model.d.vpert(39,4)=0.5;
model.d.vpert(19,5)=1;
model.d.vpert(43,5)=0.1;
model.d.vpert(43,9)=0.5;

load x_kctherm138.mat
model.fstore = fbest;
res=kineticestimate(model,res);
tout = toc;

