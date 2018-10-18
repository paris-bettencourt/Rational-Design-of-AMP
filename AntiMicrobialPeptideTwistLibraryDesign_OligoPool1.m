%% settings
clear; close all; clc
Library = {};
while numel(Library) > 12000 | numel(Library) <11950 % generate a large enough library
% extract the wild-type antimicrobial peptide amino acid sequences
AMPFile = [pwd,'/AMPs.xlsx'];
[~,AMP_list,~] = xlsread(AMPFile); % second output of the xlsread gives the text format
AMP_wt = cell2table(AMP_list(2:end,:));
AMP_wt.Properties.VariableNames = AMP_list(1,:);
AMP_seq = string(table2array(AMP_wt(:,2)));
% define 20 commonly-used amino acids, and special groups:
AA = 'ARNDCQEGHILKMFSTWYV';
Hydrophobic = 'AILMFVW';
Positive = 'RKH';
Negative = 'ED';
Polar = 'QNHSTYC';
%% --- locating position of positive aa in AMP --- %%
Index_Pos = cell(numel(AMP_seq),numel(Positive));
for i = 1:numel(Positive)
    Index_Pos(:,i) = strfind(AMP_seq,Positive(i));
end
% Merge the positive charge positions
% Note: only consider the lysine (K) and arginine (R) as "real" positive
Index_RK = cell(numel(AMP_seq),1);
% statistic: how many positive charge amino acids exists in the AMPs
Count_RK = zeros(numel(AMP_seq),1);
for i = 1:numel(AMP_seq)
   Index_RK{i} = [cell2mat(Index_Pos(i,1)),cell2mat(Index_Pos(i,2))];
   Count_RK(i) = length(Index_RK{i});
end

%% replacing positive aa by similar size neutral aa
% generate the full list of all the possible combinations:
AMP_Neut_xy = {};AMP_Neut = {};
for i = 1:numel(AMP_seq)
    for j = 1:Count_RK(i) % j define how many positive aa will be replaced
        C = combnk(1:Count_RK(i),j); % j positions chosen from all positions
        for k = 1:(numel(C)/j)
            AMP_temp = AMP_seq{i};
            AMP_temp(intersect(Index_RK{i}(C(k,:)),Index_Pos{i,1})) = 'x';% replace R by x
            AMP_temp(intersect(Index_RK{i}(C(k,:)),Index_Pos{i,2})) = 'y';% replace K by y
            AMP_Neut_xy = [AMP_Neut_xy;AMP_temp];
        end
    end
end
% apply different replace rules
X = 'QL';Y = 'NV';
for i = 1:2
    AMP_x = strrep(AMP_Neut_xy,'x',X(i));
    for j = 1:2
        AMP_y = strrep(AMP_x,'y',Y(j));
        AMP_Neut = [AMP_Neut;AMP_y];
    end
end
%% mutate in vicinity of positive charges
% If the positive charged amino acid is lysine (K) or arginine (R), the
% adjacent amino acids will be changed.
AMP_charge_x = {};AMP_charge = {};
for i = 1:numel(AMP_seq)
    for j = 1:Count_RK % j define how many positive aa will be chosen as the center of the mutation
        for k = 1:258 % replicates
                distance = round(2*randn(j,1));% randomly choose the distance between the mutation and positive aa in wt
                Index_charge = randsample(Index_RK{i},j) + distance'; % for each chosen positive aa, mutate one aa close to it.
                Index_charge(find(Index_charge <1)) = 1;
                Index_charge(find(Index_charge > length(AMP_seq{i})))= length(AMP_seq{i}); % cut off
                AMP_temp = AMP_seq{i};
                AMP_temp(Index_charge) = 'x'; % replace aa by x
                AMP_charge_x = [AMP_charge_x;AMP_temp];
        end
    end
end
AMP_charge_x = unique(AMP_charge_x);
% apply different replace rules
for i = 1:2
    AMP_charge = [AMP_charge;strrep(AMP_charge_x,'x',Positive(i))];
end
for i = 1:2
    AMP_charge = [AMP_charge;strrep(AMP_charge_x,'x',Negative(i))];
end
AMP_charge = unique(AMP_charge);
%     % statistic: how many charge amino acids exists in the AMPs
%     Index_Posnew = cell(numel(AMP_charge),numel(Positive));
%     for i = 1:numel(Positive)
%         Index_Posnew(:,i) = strfind(AMP_charge,Positive(i));
%     end
%     Index_Negnew = cell(numel(AMP_charge),numel(Negative));
%     for i = 1:numel(Negative)
%         Index_Negnew(:,i) = strfind(AMP_charge,Negative(i));
%     end
%     Index_RKnew = cell(numel(AMP_charge),1);
%     Index_EDnew = cell(numel(AMP_charge),1);
%     Count_EDnew = zeros(numel(AMP_charge),1);
%     Count_RKnew = zeros(numel(AMP_charge),1);
%     Density = zeros(numel(AMP_charge),1);
%     for i = 1:numel(AMP_charge)
%        Index_RKnew{i} = [cell2mat(Index_Posnew(i,1)),cell2mat(Index_Posnew(i,2))];
%        Count_RKnew(i) = length(Index_RKnew{i});
%        Index_EDnew{i} = [cell2mat(Index_Negnew(i,1)),cell2mat(Index_Negnew(i,2))];
%        Count_EDnew(i) = length(Index_EDnew{i});
%        Density(i) = (Count_RKnew(i) - Count_EDnew(i))./length(AMP_charge{i});
%     end
AMP_hydrophobic = {};
AMP_hydrophobic = [AMP_hydrophobic;strrep(AMP_charge_x,'x','L')];
AMP_hydrophobic = unique(AMP_hydrophobic);
%% Find the unique list of all the mutants along with wild-type.
AMP_all = [cellstr(AMP_seq);AMP_Neut;AMP_charge;AMP_hydrophobic];
AMP_all = unique(AMP_all);
%% generate the DNA library
SeqNT = cell(size(AMP_all));
for i = 1:numel(AMP_all)
    SeqNT(i) = {aa2nt(AMP_all{i})};
end

% quality control of the DNA sequences (cut sites)
BbsI = 'GAAGAC'; BsaI = 'GGTCTC';BsmBI = 'CGTCTC';
BbsIrc = seqrcomplement(BbsI);BsaIrc = seqrcomplement(BsaI);BsmBIrc = seqrcomplement(BsmBI);
Forbidden = {BbsI, BsaI, BsmBI, BbsIrc, BsaIrc, BsmBIrc};
isSafe = zeros(numel(SeqNT),6);% store if the sequence is safe to use, 1 == no cut site and safe to use

while prod(prod(isSafe)) == 0 % finish the loop when all the sequences are safe to use
idx = {};
for i = 1:6
    idx = [idx,strfind(SeqNT,Forbidden{i})];
end
for i = 1:numel(idx(:,1))
    for j = 1:6
    isSafe(i,j) = isempty(idx{i,j});
    end
end
index_Breakingrule = find(prod(isSafe') == 0);
for i = index_Breakingrule % for those sequences containing the forbidden cut sites
    SeqNT(i) = {aa2nt(AMP_all{i})}; % return another random mapping of nucleotide
end
end





% define the overhang DNA sequences
Upstream = 'CATTCGTTGAAGACGGAATG';
Downstream = 'GGAGTGAGACGAAATCAAGAGA';
Library = [char(ones(size(SeqNT))*Upstream),char(SeqNT),char(ones(size(SeqNT)))*Downstream]; % add the up and downstream overhang
Library = strrep(cellstr(Library),' ','');% remove the spaces
numel(Library)
end

save AMPseq.mat
xlswrite('TwistOligoPool1.xls',Library)