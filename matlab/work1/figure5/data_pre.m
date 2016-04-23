% obtain absolute variations of thermosphere density to SB
load SBtype
load density
SBEQNSH=intersect(SBEQNS,SBF107H,'rows');
SBEQNSM=intersect(SBEQNS,SBF107M,'rows');
SBEQNSL=intersect(SBEQNS,SBF107L,'rows');
SBEQSNH=intersect(SBEQSN,SBF107H,'rows');
SBEQSNM=intersect(SBEQSN,SBF107M,'rows');
SBEQSNL=intersect(SBEQSN,SBF107L,'rows');
% absolute change of mass density
fdenres=cell(3,3,2);% altitude,solar activity,NS
% 250 km
type='perc';
fdenres{1,1,1}=resp2SB(SBEQNSL,obser250_actual,type);
fdenres{1,2,1}=resp2SB(SBEQNSM,obser250_actual,type);
fdenres{1,3,1}=resp2SB(SBEQNSH,obser250_actual,type);
fdenres{1,1,2}=resp2SB(SBEQSNL,obser250_actual,type);
fdenres{1,2,2}=resp2SB(SBEQSNM,obser250_actual,type);
fdenres{1,3,2}=resp2SB(SBEQSNH,obser250_actual,type);
% 400km
fdenres{2,1,1}=resp2SB(SBEQNSL,obser400_actual,type);
fdenres{2,2,1}=resp2SB(SBEQNSM,obser400_actual,type);
fdenres{2,3,1}=resp2SB(SBEQNSH,obser400_actual,type);
fdenres{2,1,2}=resp2SB(SBEQSNL,obser400_actual,type);
fdenres{2,2,2}=resp2SB(SBEQSNM,obser400_actual,type);
fdenres{2,3,2}=resp2SB(SBEQSNH,obser400_actual,type);
% 550km
fdenres{3,1,1}=resp2SB(SBEQNSL,obser550_actual,type);
fdenres{3,2,1}=resp2SB(SBEQNSM,obser550_actual,type);
fdenres{3,3,1}=resp2SB(SBEQNSH,obser550_actual,type);
fdenres{3,1,2}=resp2SB(SBEQSNL,obser550_actual,type);
fdenres{3,2,2}=resp2SB(SBEQSNM,obser550_actual,type);
fdenres{3,3,2}=resp2SB(SBEQSNH,obser550_actual,type);
dir='F:\mywork\matlabworkspace\ther_dens_sect\figure5\';
nam='fdenres.mat';
save([dir,nam],'fdenres')
