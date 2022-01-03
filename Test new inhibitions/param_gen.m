function [MMparam] = param_gen(model)


substrate = {'A','B','C','D','E','F','G','H','I','J','L','M'};
product = {'P','Q','R','S','T','U','V','W','X','Y','Z','N'};
inhibitor = {'a','b','c','d','e','f','g','h','i','j','l','m'};

%[ k,dKdp,Revs,dRevsdp ] = calc_k( model,x0 );
for i = 1:size(model.ensemble.ne,1)
    dkdc{i,1} = model.kinetic(i).dkdc;
end
p = 1;
for i = 1:size(dkdc,1)
    currentrxn = dkdc{i};
    for j = 1:size(currentrxn,2)
        a = find(currentrxn(:,j));
        if ~isempty(a)
            kinetic_name{p,1} = model.metprop(a).metid;
            p = p + 1;
        else
            kinetic_name{p,1} = model.kinetic(i).id;
            p = p + 1;
        end
    end   
end
p = 1;
for i = 1:size(dkdc,1)
    for j = 1:model.kinetic(i).nk
        RXN_NAME{p,1} = model.kinetic(i).id;
        p = p + 1;
    end
end



MMparam = struct('id',{''});
for i = 1:length(model.kinetic)
    MMparam(i).id = model.kinetic(i).id;
end
%MMparam.id = model.kinetic.id;

%[MMparam.id(1:length(model.kinetic))] = deal(struct('param',{''}));
for i = 1:length(model.kinetic)
    location = find(strcmp(RXN_NAME,model.kinetic(i).id));
    l = 1;
    k = 1;
    param = cell(2*(length(model.kinetic(i).subs)+length(model.kinetic(i).pdts))+2,1);
    MMparam(i).id = model.kinetic(i).id;
    for j = 1:length(model.kinetic(i).subs)
        param{k} = strcat('k',num2str(k),'*',substrate{j});
        participant{l} = substrate{j};
        l = l + 1;
        param{k+1} = strcat('k',num2str(k+1));
        k = k + 2;
    end
    for j = 1:2
        param{k} = strcat('k',num2str(k));
        k = k + 1;
    end
    for j = 1:length(model.kinetic(i).pdts)
        param{k} = strcat('k',num2str(k));
        param{k+1} = strcat('k',num2str(k+1),'*',product{length(model.kinetic(i).pdts)-(j-1)});
        participant{l} = product{j};
        l = l + 1;
        k = k + 2;
    end
    kk = 1;
    s = 1;
    if ~isempty(model.kinetic(i).c_in)
        for j = 1:length(model.kinetic(i).c_in)
            inhib_param{kk} = strcat('k',num2str(k),'*',inhibitor{s});
            inhib_param{kk+1} = strcat('k',num2str(k+1));
            participant{l} = inhibitor{s};
            k = k + 2;
            kk = kk + 2;
            s = s + 1;
            l = l + 1;
        end
    end
    if ~isempty(model.kinetic(i).uc_in)
        for j = 1:length(model.kinetic(i).uc_in)
            inhib_param{kk} = strcat('k',num2str(k),'*',inhibitor{s});
            inhib_param{kk+1} = strcat('k',num2str(k+1));
            participant{l} = inhibitor{s};
            k = k + 2;
            kk = kk + 2;
            s = s + 1;
            l = l + 1;
        end
    end
    if ~isempty(model.kinetic(i).nc_in)
        for j = 1:length(model.kinetic(i).nc_in)
            inhib_param{kk} = strcat('k',num2str(k),'*',inhibitor{s});
            inhib_param{kk+1} = strcat('k',num2str(k+1));
            inhib_param{kk+2} = strcat('k',num2str(k+2),'*',inhibitor{s});
            inhib_param{kk+3} = strcat('k',num2str(k+3));
            participant{l} = inhibitor{s};
            k = k + 4;
            kk = kk + 4;
            s = s + 1;
            l = l + 1;
        end
    end
        
        
    participant = participant';
    MMparam(i).param = param;
    if exist('inhib_param','var') == 1
        MMparam(i).inhib_param = inhib_param';
    else
        MMparam(i).inhib_param = [];
    end
    MMparam(i).ne = (length(param))/2;
    MMparam(i).np = MMparam(i).ne-1;
    MMparam(i).participant = participant;
    param = [];
    inhib_param = [];
    participant = [];
end
end