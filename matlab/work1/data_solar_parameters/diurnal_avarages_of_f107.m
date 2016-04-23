% purpose: 
% calculate daily average F10.7 at the Earth-Sun distance,save them in
% DailyF107
%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 'indexnew.mat'
DailyF107=zeros(size(indexnew,1),3);
DailyF107(:,1)=indexnew(:,1);
DailyF107(:,2)=date2doy(indexnew(:,1:3));
DailyF107(:,3)=indexnew(:,14);
save('F:\mywork\matlabworkspace\ther_dens_sect\data_solar_parameters\DailyF107.mat',...
    'DailyF107')