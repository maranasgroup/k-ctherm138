[~,~,reg_info] = xlsread('reg_id_031021.xlsx','Sheet1','A2:F20');
[~,~,mechanism] = xlsread('ctherm_mechanism_v8_juice.xlsx','A1:J139');
%parameter_names = ['one','two','three','four','five','six','seven','eight','nine','ten','eleven','twelve','thirteen','fourteen','fifteen','sixteen','seventeen','eighteen','nineteen'];
reg_ind = cell2mat(reg_info(:,5));
aa = cell(19,1);
for i = (1:19)
    aa{i,1} = combnk([1:19],i);
end
counter = 1;
 for i = 1:4%length(aa)
     aa_pull = aa{i};
     for j = 1:size(aa_pull,1)
         Sheet1 = mechanism;
         for k = 1:size(aa_pull,2)
            ind = find(reg_ind==aa_pull(j,k));
            effector = reg_info(ind,3);
            effector = effector{1};
            type = reg_info(ind,4);
            type = type{1};
            RXN = reg_info(ind,2);
            RXN = RXN{1};
            k_ID = reg_info(ind,6);
            k_ID = k_ID{1};
            rxn_ind = find(strcmp(RXN,Sheet1(:,1)));
            parameters{counter,k} = k_ID;
            if strcmp(type,'i') == 1
                Sheet1(rxn_ind,5) = strrep(Sheet1(rxn_ind,5),effector,'');
                Sheet1(rxn_ind,5) = strrep(Sheet1(rxn_ind,5),';','');
            else
                Sheet1(rxn_ind,8) = strrep(Sheet1(rxn_ind,8),effector,'');
            end
         end
         for ii = 1:size(Sheet1,1)
             for jj = 1:size(Sheet1,2)
                 if isnan(Sheet1{ii,jj}) == 1
                     Sheet1{ii,jj} = [];
                 end
             end
         end
         
         Sheet1 = cell2table(Sheet1(2:end,:),'VariableNames',Sheet1(1,:));
         writetable(Sheet1,strcat('ctherm_mechanism_v8_juice',num2str(counter),'.xlsx'))
         %csvwrite(strcat('ctherm_mechanism_v8_juice',num2str(counter),'.xlsx'),Sheet1);
         counter = counter + 1;
     end
 end

parameters = cell2table(parameters(1:end,:));
writetable(parameters,'parameters.xlsx')
 