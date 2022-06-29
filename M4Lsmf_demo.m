clear all;
clc;
funspath=[pwd,filesep,'data',filesep];
addpath(funspath);
fprintf('load LncRNA dataset\n');
%% load data to construct the relation matrix R
load('lncRNAGene')
R12 = LGasso;
load('lncRNAMiA');
R13 = lncMI;
load('LncDOs')
load('lncDisease2')
LncDO = lncDisease + LncCancer;
LncDO(LncDO>1)=1;
R14 = LncDO;
load('GeneDisease');
R24 = GDasso;
load('miRNAGene');
R32 = MGasso;
load('MiDOs');
R34 = miDOs;

%% filter the no instanse miRNA.
fun_stat_miRNA = sum(R13,1);
sel_miRNA_idx = find(fun_stat_miRNA>0);
R13 = R13(:,sel_miRNA_idx);
R32 = R32(sel_miRNA_idx,:);
R34 = R34(sel_miRNA_idx,:);

Rcell = {R12,R13,R14,R24,R32,R34};

load('ThetaCell');
fprintf('ThetaCell is load complecated\n');
ThetaCell{3,1} = ThetaCell{3,1}(sel_miRNA_idx,sel_miRNA_idx);
ThetaCell{3,2} = ThetaCell{3,2}(sel_miRNA_idx,sel_miRNA_idx);
ThetaCell{3,3} = ThetaCell{3,3}(sel_miRNA_idx,sel_miRNA_idx);
ThetaCell{3,4} = ThetaCell{3,4}(sel_miRNA_idx,sel_miRNA_idx);
ThetaCell{3,5} = ThetaCell{3,5}(sel_miRNA_idx,sel_miRNA_idx);
ThetaCell{3,6} = ThetaCell{3,6}(sel_miRNA_idx,sel_miRNA_idx);

R1 = [R12 R14];
R2 = R24;
R3 = [R32 R34];
R4 = [R14' R24' R34'];

rank_k=30;

fprintf('Start SVD....\n');

[G1,~] = NNDSVD(abs(R1),rank_k,0);
[G2,~] = NNDSVD(abs(R2),rank_k,0);
[G3,~] = NNDSVD(abs(R3),rank_k,0);
[G4,~] = NNDSVD(abs(R4),rank_k,0);


fprintf('End SVD....\n');
Gcell = {G1,G2,G3,G4};
%%
instanseIdx = {2,3,4,8,10,12};   
nTypes=4;
max_iter = 30;

fprintf('Start M4Lsmf....\n');
M4Lsmf(  nTypes,instanseIdx,Gcell,Rcell,ThetaCell,max_iter  );