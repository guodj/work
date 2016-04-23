% Find SBs in different seasons.
load 'SBlist67_07'
SBME=findSBseason(SBlist(:,1:2),'me');
SBSE=findSBseason(SBlist(:,1:2),'se');
SBJS=findSBseason(SBlist(:,1:2),'js');
SBDS=findSBseason(SBlist(:,1:2),'ds');
save('F:\mywork\matlabworkspace\ther_dens_sect\data_SBlist\SBtype.mat',...
    'SBSE','SBME','SBJS','SBDS','-append')