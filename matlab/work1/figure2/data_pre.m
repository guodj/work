% solar parameters and ap variation.
load daily_GSE_Bx.mat
load daily_GSE_By.mat
load daily_GSE_Bz.mat
load daily_GSM_Bz.mat
load daily_SW_v.mat
load daily_Ap.mat
load daily_F107.mat
load SBtype.mat
% fpic
fpic=cell(7,4,2);% solar parameters, season, orien
SBtype='actual';
%% variation of GSE Bx
% away_toward
fpic{1,1,1}=resp2SB(SBATME,daily_GSE_Bx,SBtype);
fpic{1,2,1}=resp2SB(SBATSE,daily_GSE_Bx,SBtype);
fpic{1,3,1}=resp2SB(SBATJS,daily_GSE_Bx,SBtype);
fpic{1,4,1}=resp2SB(SBATDS,daily_GSE_Bx,SBtype);
% toward_away
fpic{1,1,2}=resp2SB(SBTAME,daily_GSE_Bx,SBtype);
fpic{1,2,2}=resp2SB(SBTASE,daily_GSE_Bx,SBtype);
fpic{1,3,2}=resp2SB(SBTAJS,daily_GSE_Bx,SBtype);
fpic{1,4,2}=resp2SB(SBTADS,daily_GSE_Bx,SBtype);
%% variation of GSE By
% away_toward
fpic{2,1,1}=resp2SB(SBATME,daily_GSE_By,SBtype);
fpic{2,2,1}=resp2SB(SBATSE,daily_GSE_By,SBtype);
fpic{2,3,1}=resp2SB(SBATJS,daily_GSE_By,SBtype);
fpic{2,4,1}=resp2SB(SBATDS,daily_GSE_By,SBtype);
% toward_away
fpic{2,1,2}=resp2SB(SBTAME,daily_GSE_By,SBtype);
fpic{2,2,2}=resp2SB(SBTASE,daily_GSE_By,SBtype);
fpic{2,3,2}=resp2SB(SBTAJS,daily_GSE_By,SBtype);
fpic{2,4,2}=resp2SB(SBTADS,daily_GSE_By,SBtype);
%% variation of GSE Bz
% away_toward
fpic{3,1,1}=resp2SB(SBATME,daily_GSE_Bz,SBtype);
fpic{3,2,1}=resp2SB(SBATSE,daily_GSE_Bz,SBtype);
fpic{3,3,1}=resp2SB(SBATJS,daily_GSE_Bz,SBtype);
fpic{3,4,1}=resp2SB(SBATDS,daily_GSE_Bz,SBtype);
% toward_away
fpic{3,1,2}=resp2SB(SBTAME,daily_GSE_Bz,SBtype);
fpic{3,2,2}=resp2SB(SBTASE,daily_GSE_Bz,SBtype);
fpic{3,3,2}=resp2SB(SBTAJS,daily_GSE_Bz,SBtype);
fpic{3,4,2}=resp2SB(SBTADS,daily_GSE_Bz,SBtype);
%% variation of GSM Bz
% away_toward
fpic{4,1,1}=resp2SB(SBATME,daily_GSM_Bz,SBtype);
fpic{4,2,1}=resp2SB(SBATSE,daily_GSM_Bz,SBtype);
fpic{4,3,1}=resp2SB(SBATJS,daily_GSM_Bz,SBtype);
fpic{4,4,1}=resp2SB(SBATDS,daily_GSM_Bz,SBtype);
% toward_away
fpic{4,1,2}=resp2SB(SBTAME,daily_GSM_Bz,SBtype);
fpic{4,2,2}=resp2SB(SBTASE,daily_GSM_Bz,SBtype);
fpic{4,3,2}=resp2SB(SBTAJS,daily_GSM_Bz,SBtype);
fpic{4,4,2}=resp2SB(SBTADS,daily_GSM_Bz,SBtype);
%% variation of solar wind velocity
% away_toward
fpic{5,1,1}=resp2SB(SBATME,daily_SW_v,SBtype);
fpic{5,2,1}=resp2SB(SBATSE,daily_SW_v,SBtype);
fpic{5,3,1}=resp2SB(SBATJS,daily_SW_v,SBtype);
fpic{5,4,1}=resp2SB(SBATDS,daily_SW_v,SBtype);
% toward_away
fpic{5,1,2}=resp2SB(SBTAME,daily_SW_v,SBtype);
fpic{5,2,2}=resp2SB(SBTASE,daily_SW_v,SBtype);
fpic{5,3,2}=resp2SB(SBTAJS,daily_SW_v,SBtype);
fpic{5,4,2}=resp2SB(SBTADS,daily_SW_v,SBtype);
%% variation of Ap
% away_toward
fpic{6,1,1}=resp2SB(SBATME,daily_Ap,SBtype);
fpic{6,2,1}=resp2SB(SBATSE,daily_Ap,SBtype);
fpic{6,3,1}=resp2SB(SBATJS,daily_Ap,SBtype);
fpic{6,4,1}=resp2SB(SBATDS,daily_Ap,SBtype);
% toward_away
fpic{6,1,2}=resp2SB(SBTAME,daily_Ap,SBtype);
fpic{6,2,2}=resp2SB(SBTASE,daily_Ap,SBtype);
fpic{6,3,2}=resp2SB(SBTAJS,daily_Ap,SBtype);
fpic{6,4,2}=resp2SB(SBTADS,daily_Ap,SBtype);
%% variation of F107
% away_toward
fpic{7,1,1}=resp2SB(SBATME,daily_F107,SBtype);
fpic{7,2,1}=resp2SB(SBATSE,daily_F107,SBtype);
fpic{7,3,1}=resp2SB(SBATJS,daily_F107,SBtype);
fpic{7,4,1}=resp2SB(SBATDS,daily_F107,SBtype);
% toward_away
fpic{7,1,2}=resp2SB(SBTAME,daily_F107,SBtype);
fpic{7,2,2}=resp2SB(SBTASE,daily_F107,SBtype);
fpic{7,3,2}=resp2SB(SBTAJS,daily_F107,SBtype);
fpic{7,4,2}=resp2SB(SBTADS,daily_F107,SBtype);
%% save
dir='F:\mywork\matlabworkspace\ther_dens_sect\figure2\';
nam='fig2dat.mat';
save([dir,nam],'fpic')
