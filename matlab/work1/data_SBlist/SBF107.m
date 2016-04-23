load SBlist67_07.mat
% obtain the daily average F10.7 (DailyF107) during 19650101~20101231
diurnal_avarages_of_f107;
% smooth f10.7
DailyF107_smooth=DailyF107;
DailyF107_smooth(:,3)=smooth(DailyF107_smooth(:,3),11);
SBF107_i=find(ismember(DailyF107_smooth(:,1:2),SBlist(:,1:2),'rows'));
SBF10767_07=DailyF107_smooth(SBF107_i,:);
SBF107H=SBF10767_07(SBF10767_07(:,3)>160,1:2);
SBF107L=SBF10767_07(SBF10767_07(:,3)<100,1:2);
SBF107M=SBF10767_07(SBF10767_07(:,3)>100 & SBF10767_07(:,3)<160,1:2);
save('F:\mywork\matlabworkspace\ther_dens_sect\data_SBlist\SBtype.mat',...
    'SBF107H','SBF107M','SBF107L','-append')