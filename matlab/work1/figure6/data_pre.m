% obtain absolute variations of thermosphere density to SB
load SBtype
load density
load daily_GSM_Bz
load daily_SW_v
SBEQNSCIR=intersect(SBEQNS,SByesCIR,'rows');
SBEQNSNOCIR=intersect(SBEQNS,SBnoCIR,'rows');
SBEQSNCIR=intersect(SBEQSN,SByesCIR,'rows');
SBEQSNNOCIR=intersect(SBEQSN,SBnoCIR,'rows');

fdenreCIR=cell(5,2,2);% paramters,CIR,NS
type1='actual';
type2='perc';
%% GSM,Bz
fdenreCIR{1,1,1}=resp2SB(SBEQNSCIR,daily_GSM_Bz,type1);
fdenreCIR{1,2,1}=resp2SB(SBEQNSNOCIR,daily_GSM_Bz,type1);
fdenreCIR{1,1,2}=resp2SB(SBEQSNCIR,daily_GSM_Bz,type1);
fdenreCIR{1,2,2}=resp2SB(SBEQSNNOCIR,daily_GSM_Bz,type1);
%% solar wind velocity
fdenreCIR{2,1,1}=resp2SB(SBEQNSCIR,daily_SW_v,type1);
fdenreCIR{2,2,1}=resp2SB(SBEQNSNOCIR,daily_SW_v,type1);
fdenreCIR{2,1,2}=resp2SB(SBEQSNCIR,daily_SW_v,type1);
fdenreCIR{2,2,2}=resp2SB(SBEQSNNOCIR,daily_SW_v,type1);
%% density 250 km
fdenreCIR{3,1,1}=resp2SB(SBEQNSCIR,obser250_actual,type2);
fdenreCIR{3,2,1}=resp2SB(SBEQNSNOCIR,obser250_actual,type2);
fdenreCIR{3,1,2}=resp2SB(SBEQSNCIR,obser250_actual,type2);
fdenreCIR{3,2,2}=resp2SB(SBEQSNNOCIR,obser250_actual,type2);
%% density 400km
fdenreCIR{4,1,1}=resp2SB(SBEQNSCIR,obser400_actual,type2);
fdenreCIR{4,2,1}=resp2SB(SBEQNSNOCIR,obser400_actual,type2);
fdenreCIR{4,1,2}=resp2SB(SBEQSNCIR,obser400_actual,type2);
fdenreCIR{4,2,2}=resp2SB(SBEQSNNOCIR,obser400_actual,type2);
%% density 550km
fdenreCIR{5,1,1}=resp2SB(SBEQNSCIR,obser550_actual,type2);
fdenreCIR{5,2,1}=resp2SB(SBEQNSNOCIR,obser550_actual,type2);
fdenreCIR{5,1,2}=resp2SB(SBEQSNCIR,obser550_actual,type2);
fdenreCIR{5,2,2}=resp2SB(SBEQSNNOCIR,obser550_actual,type2);
%% save
dir='F:\mywork\matlabworkspace\ther_dens_sect\figure6\';
nam='fdenreCIR.mat';
save([dir,nam],'fdenreCIR')
