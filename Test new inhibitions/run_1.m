%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

clear, clc

res={};
load optimal.mat
model=modcompile(strcat(pwd,'/Cth_model_v8_juice.xlsx'),strcat(pwd,'/ctherm_mechanism_v8_juice.xlsx'),strcat(pwd,'/Ctherm_data_all_v8_juice.xlsx'));
load xbest.mat
x = xbest;
res.reinit_data = x;
[k,~,~,~] = calc_k_act( model,x );
updated_rxn = 39;
kI = 100;
done = false;

for i = 1:length(model.rid)
    nact(i,1) = length(model.kinetic(i).act);
end

while ~done
    for i = 1:length(model.ensemble.ne)
        ne(i,1) = model.ensemble.ne{i} - 1;
    end
    for i = 1:length(model.ensemble.nei)
        nei(i,1) = model.ensemble.nei{i};
    end
    for i = 1:length(model.ensemble.nvr)
        nvr(i,1) = model.ensemble.nvr{i};
    end

    Vnet = model.d.flx{1}(updated_rxn);
    kblocks = model.p.kblocks';
    k_pull(:,1) = k(kblocks(updated_rxn)+1:kblocks(updated_rxn+1));
    k_use(1:2*length(model.kinetic(updated_rxn).subs)+2*length(model.kinetic(updated_rxn).pdts)+2*length(model.kinetic(updated_rxn).act)+2,1) = k_pull(1:2*length(model.kinetic(updated_rxn).subs)+2*length(model.kinetic(updated_rxn).pdts)+2*length(model.kinetic(updated_rxn).act)+2);
    k_use(2*length(model.kinetic(updated_rxn).subs)+2*length(model.kinetic(updated_rxn).pdts)+2*length(model.kinetic(updated_rxn).act)+3,1) = kI;
    k_use(2*length(model.kinetic(updated_rxn).subs)+2*length(model.kinetic(updated_rxn).pdts)+2*length(model.kinetic(updated_rxn).act)+4:2*length(model.kinetic(updated_rxn).subs)+2*length(model.kinetic(updated_rxn).pdts)+2*length(model.kinetic(updated_rxn).act)+nei(updated_rxn)+3,1) = k_pull(2*length(model.kinetic(updated_rxn).subs)+2*length(model.kinetic(updated_rxn).pdts)+2*length(model.kinetic(updated_rxn).act)+3:end);


    
    cnt = 1;
    if length(model.kinetic(updated_rxn).act) == 0
        Aeq = zeros(ne(updated_rxn)+nei(updated_rxn)+3,ne(updated_rxn)+nei(updated_rxn)+2);
        for i = 1:nvr(updated_rxn)-1
            Aeq(i,i) = k_use(cnt);
            Aeq(i,i+1) = -k_use(cnt+1);
            cnt = cnt + 2;
        end
        Aeq(nvr(updated_rxn),nvr(updated_rxn)) = k_use(cnt);
        Aeq(nvr(updated_rxn),1) = -k_use(cnt+1);
        Aeq(nvr(updated_rxn)+1,nvr(updated_rxn)+1) = 1;
        Aeq(nvr(updated_rxn)+1,1) = -kI;
    else
        Aeq = zeros(ne(updated_rxn)+nei(updated_rxn)+3,ne(updated_rxn)+nei(updated_rxn)+2);
        for i = 1:nvr(updated_rxn)-2
            Aeq(i,i) = k_use(cnt);
            Aeq(i,i+1) = -k_use(cnt+1);
            cnt = cnt + 2;
        end
        Aeq(nvr(updated_rxn)-1,nvr(updated_rxn)-1) = k_use(cnt);
        Aeq(nvr(updated_rxn)-1,1) = -k_use(cnt+1);
        cnt = cnt + 2;
        for i = 1:length(model.kinetic(updated_rxn).act)
            Aeq(nvr(updated_rxn)-1+i,nvr(updated_rxn)-1+i) = k_use(cnt);
            Aeq(nvr(updated_rxn)-1+i,1) = -k_use(cnt+1);
            cnt = cnt+2;
        end
        Aeq(nvr(updated_rxn)+i,nvr(updated_rxn)+i) = 1;
        Aeq(nvr(updated_rxn)+i,1) = -kI;
    end
    
    

    kI_cnt = nvr(updated_rxn) + 1;
    if nei(updated_rxn) ~= 0
        for i = 1:nei(updated_rxn)+1
            Aeq(kI_cnt,1) = -k_use(2*length(model.kinetic(updated_rxn).subs)+2*length(model.kinetic(updated_rxn).pdts)+2*nact(updated_rxn)+2+i);
            Aeq(kI_cnt,kI_cnt) = 1;
            kI_cnt = kI_cnt + 1;
        end
    end

    for i = 1:size(Aeq,2)
        Aeq(end,i) = 1;
    end
    
    if length(model.kinetic(updated_rxn).act) == 0
        for i = 1:nvr(updated_rxn)
            beq(i,1) = Vnet;
        end
        for i = nvr(updated_rxn)+1:nvr(updated_rxn)+nei(updated_rxn)+1
            beq(i,1) = 0;
        end
        beq(size(Aeq,1),1) = 1;
    else
        for i = 1:nvr(updated_rxn)-1
            beq(i,1) = Vnet;
        end
        beq(nvr(updated_rxn),1) = 0;
        for i = nvr(updated_rxn)+1:nvr(updated_rxn)+nei(updated_rxn)+1
            beq(i,1) = 0;
        end
        beq(size(Aeq,1),1) = 1;
    end
        
    
    %pull updated variable vector up to reaction to be updated 
    x_update(1:sum(ne(1:updated_rxn -1))+sum(nei(1:updated_rxn - 1)) + sum(nvr(1:updated_rxn-1)),1) = res.reinit_data(1:sum(ne(1:updated_rxn -1))+sum(nei(1:updated_rxn - 1)) + sum(nvr(1:updated_rxn-1))) ;

    %pull values to use as initial enzyme fractions in lp problem
    x_pull(:,1) = res.reinit_data(length(x_update)+1:length(x_update)+ne(updated_rxn));
    x_pull(length(x_pull)+1,1) = 0;
    x_pull = [1-sum(x_pull); x_pull];

    %finish updated variable vector with regulation
    %x_update(length(x_update)+1:length(x_update)+ne(updated_rxn),1) = res.reinit_data(length(x_update)+1:length(x_update)+ne(updated_rxn));
    %x_update(length(x_update)+1:length(res.reinit_data)+1,1) = res.reinit_data(length(x_update):end);

    lb = zeros(1,ne(updated_rxn)+nei(updated_rxn)+2);
    ub = ones(1,ne(updated_rxn)+nei(updated_rxn)+2);
    %x0 = x_pull';
    x0 = rand(ne(updated_rxn)+nei(updated_rxn)+2,1);
    x0 = x0./sum(x0);
    x0 = x0';

    fun = @(x) x(nvr(updated_rxn)+1);
    options = optimoptions('fmincon','StepTolerance',1e-30,'Algorithm','active-set','MaxFunctionEvaluations',100000);
    [x,fval,exitflag] = fmincon(fun,x0,[],[],Aeq,beq,lb,ub);
    if exitflag ~= 1
        kI = kI/10;
    elseif exitflag == 1
        done = true;
        x_store = x';
    end
    clear x_pull
end

x_update(length(x_update)+1:length(x_update)+length(x_store)-1,1) = x_store(2:end);
for i = 1:nvr(updated_rxn)
    if i<nvr(updated_rxn)
        x_update(length(x_update)+1) = k_use(i*2)*x_store(i+1);
    else
        x_update(length(x_update)+1) = k_use(i*2)*x_store(1);
    end
end
x_update(length(x_update)+1:length(res.reinit_data)+1,1) = res.reinit_data(length(x_update):end);

clear model

model=modcompile(strcat(pwd,'/Cth_model_v8_juice.xlsx'),strcat(pwd,'/ctherm_mechanism_v8_juice_1.xlsx'),strcat(pwd,'/Ctherm_data_all_v8_juice.xlsx'));
res.reinit_data = x_update;

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

model.options.reinit =0;
model.fstore = res.fmin;
model.fstore = fbest;
model.starter = 1;
res=kineticestimate_v4_1(model,res);
save('optimal_1.mat','res','model');