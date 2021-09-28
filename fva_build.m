clear, clc

[~,rxn,~] = xlsread('core_model_4_16_WT.xlsx','bounds','A2:A74');
[~,metab,~] = xlsread('core_model_4_16_WT.xlsx','Metabolites','A2:A149');
[,~,stoic,~] = xlsread('core_model_4_16_WT.xlsx','bounds','C2:C74');
[LB,~,~] = xlsread('core_model_4_16_WT.xlsx','bounds','D2:D74');
[UB,~,~] = xlsread('core_model_4_16_WT.xlsx','bounds','E2:E74');

fileID = fopen('sij.txt','w');
for i = 1:length(stoic)
    for_use = strsplit(stoic{i},' ');
    for j = 1:length(for_use)
        a = strfind(for_use{j},'>');
        if ~isempty(a)
            st = j;
        end
    end
    for j = 1:length(for_use)
        if ~isempty(strfind(for_use{j},')')) && j < st
            fprintf(fileID,'''%s''.''%s'' -%s\n',for_use{j+1},rxn{i},strrep(strrep(for_use{j},'(',''),')',''));
        elseif ~isempty(strfind(for_use{j},')')) && j > st
            fprintf(fileID,'''%s''.''%s'' %s\n',for_use{j+1},rxn{i},strrep(strrep(for_use{j},'(',''),')',''));
        end
    end
end
fclose(fileID);
fileID = fopen('lower_bound.txt','w');
for i = 1:length(rxn)
    %fprintf(fileID,'''%s'' -1000\n',rxn{i});
    fprintf(fileID,'''%s'' %f\n',rxn{i},LB(i));
end
fclose(fileID);
fileID = fopen('upper_bound.txt','w');
for i = 1:length(rxn)
    %fprintf(fileID,'''%s'' 1000\n',rxn{i});
    fprintf(fileID,'''%s'' %f\n',rxn{i},UB(i));
end
fclose(fileID);
fileID = fopen('reactions.txt','w');
for i = 1:length(rxn)
    fprintf(fileID,'''%s''\n',rxn{i});
end
fclose(fileID);
fileID = fopen('metabolites.txt','w');
for i = 1:length(metab)
    fprintf(fileID,'''%s''\n',metab{i});
end
fclose(fileID);

    

