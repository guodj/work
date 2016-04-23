% find A-T SB and T-A SB
% -1 means away to toward. 1 means toward to away
load 'SBlist67_07.mat'
SBAT=SBlist(SBlist(:,3)==-1,1:2);
SBTA=SBlist(SBlist(:,3)==1,1:2);
save('F:\mywork\matlabworkspace\ther_dens_sect\data_SBlist\SBtype.mat',...
    'SBAT','SBTA','-append')