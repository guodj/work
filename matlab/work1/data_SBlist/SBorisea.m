SBorien
SBseason
SBATME=intersect(SBAT,SBME,'rows');
SBATSE=intersect(SBAT,SBSE,'rows');
SBATJS=intersect(SBAT,SBJS,'rows');
SBATDS=intersect(SBAT,SBDS,'rows');
SBTAME=intersect(SBTA,SBME,'rows');
SBTASE=intersect(SBTA,SBSE,'rows');
SBTAJS=intersect(SBTA,SBJS,'rows');
SBTADS=intersect(SBTA,SBDS,'rows');
dir='F:\mywork\matlabworkspace\ther_dens_sect\data_SBlist\';
nam='SBtype.mat';
save([dir,nam],'SBATME','SBATSE','SBATJS','SBATDS',...
    'SBTAME','SBTASE','SBTAJS','SBTADS','-append')