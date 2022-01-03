clear, clc

res={};
load xbest.mat
%load optimal_reinit1_06_21_19_CF_62.mat
model=modcompile(strcat(pwd,'/Cth_model_v8_juice.xlsx'),strcat(pwd,'/ctherm_mechanism_v8_juice.xlsx'),strcat(pwd,'/Ctherm_data_all_v8_juice.xlsx'));
%x = res.reinit_data;
x = xbest;
%for i = 3:11
%    model.d.vpert(:,3) = [];
%end
model.d.vpert(:,2) = model.d.vpert(:,1);


%x = res.reinit_data;
[k,~,~,~] = calc_k_act( model,x );
vstore = cell(50,1);
cstore = cell(50,1);
[~,~,nums] = xlsread('parameters.xlsx',1,'A2:D5036');
for i = 1:size(nums,1)
    for j = 1:size(nums,2)
        if isnan(nums{i,j}) == 1
            nums{i,j} = [];
        end
    end
end
bbb = 1;
for i = 1:19
    k2 = k;
    numcat2 = [];
    for j = 1:size(nums,2)
        numcat1 = nums{i,j};
        if ~isa(numcat1,'double')
            numcat1 = strsplit(numcat1,',');
            numcat1 = str2double(numcat1);
        end
        numcat2 = [numcat2 numcat1];
    end
    numcat2 = sort(numcat2,'descend');
    for ll = 1:length(numcat2)
        k2(numcat2(ll)) = [];
    end
    model2=modcompile(strcat(pwd,'/Cth_model_v8_juice.xlsx'),strcat(pwd,'/ctherm_mechanism_v8_juice',num2str(i),'.xlsx'),strcat(pwd,'/Ctherm_data_all_v8_juice.xlsx'));
    [x,~,~,~] = initialize_v4(model2);
    [r,W,J,vop,cs] = rescalc_reg(x,model2,k2);
    vstore{bbb} = vop(:,1);
    cstore{bbb} = cs(:,1);
    bbb = bbb + 1;
end
for i = 1:length(vstore)
    vv = vstore{i};
    cc = cstore{i};
    etoh_v(i,:) = vv(7,:);
    etoh_c(i,:) = cc(29,:);
    vv_store(:,i) = vv;
    cc_store(:,i) = cc;
end
% for i = 1:length(res.predictions)
%     etoh_modv(1,i) = res.predictions(i).fluxes(7).val;
%     etoh_modc(1,i) = res.predictions(i).concentrations(29).val;
% end
% etoh_v = [etoh_modv;etoh_v];
% etoh_c = [etoh_modc;etoh_c];

save('ethanol_full_1.mat','vv_store','cc_store');