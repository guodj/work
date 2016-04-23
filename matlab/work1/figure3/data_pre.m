% obtain absolute variations of thermosphere density to SB
load SBtype
load density
% absolute change of mass density
fdenab=cell(3,4,2);% altitude,season,orien
% 250 km
type='residual';
fdenab{1,1,1}=resp2SB(SBATME,obser250_actual,type);
fdenab{1,2,1}=resp2SB(SBATSE,obser250_actual,type);
fdenab{1,3,1}=resp2SB(SBATJS,obser250_actual,type);
fdenab{1,4,1}=resp2SB(SBATDS,obser250_actual,type);
fdenab{1,1,2}=resp2SB(SBTAME,obser250_actual,type);
fdenab{1,2,2}=resp2SB(SBTASE,obser250_actual,type);
fdenab{1,3,2}=resp2SB(SBTAJS,obser250_actual,type);
fdenab{1,4,2}=resp2SB(SBTADS,obser250_actual,type);
% 400km
fdenab{2,1,1}=resp2SB(SBATME,obser400_actual,type);
fdenab{2,2,1}=resp2SB(SBATSE,obser400_actual,type);
fdenab{2,3,1}=resp2SB(SBATJS,obser400_actual,type);
fdenab{2,4,1}=resp2SB(SBATDS,obser400_actual,type);
fdenab{2,1,2}=resp2SB(SBTAME,obser400_actual,type);
fdenab{2,2,2}=resp2SB(SBTASE,obser400_actual,type);
fdenab{2,3,2}=resp2SB(SBTAJS,obser400_actual,type);
fdenab{2,4,2}=resp2SB(SBTADS,obser400_actual,type);
% 550km
fdenab{3,1,1}=resp2SB(SBATME,obser550_actual,type);
fdenab{3,2,1}=resp2SB(SBATSE,obser550_actual,type);
fdenab{3,3,1}=resp2SB(SBATJS,obser550_actual,type);
fdenab{3,4,1}=resp2SB(SBATDS,obser550_actual,type);
fdenab{3,1,2}=resp2SB(SBTAME,obser550_actual,type);
fdenab{3,2,2}=resp2SB(SBTASE,obser550_actual,type);
fdenab{3,3,2}=resp2SB(SBTAJS,obser550_actual,type);
fdenab{3,4,2}=resp2SB(SBTADS,obser550_actual,type);
dir='F:\mywork\matlabworkspace\ther_dens_sect\figure3\';
nam='fdenab.mat';
save([dir,nam],'fdenab')
