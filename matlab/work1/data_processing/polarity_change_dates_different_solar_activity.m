% purpose:
% classify 'away_toward_date' and 'toward_away_date' according to F10.7 (solar
% activity)
%
% define variables:
% intput: away_toward_date, toward_away_date
% outputs: away_toward_date_meh ...
%
% other instruction:
% solar maximum~medium:107=160; solar medium~minimum: F107=100
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine the SB crossing dates
sector_polarity_change_dates_SBlist;
% obtain the daily average F10.7 (DailyF107) during 19650101~20101231
diurnal_avarages_of_f107;
% smooth f10.7
DailyF107_smooth=DailyF107;
DailyF107_smooth(:,3)=smooth(DailyF107_smooth(:,3),11);

isSBat=ismember(DailyF107_smooth(:,1:2),away_toward_date,'rows');
isSBta=ismember(DailyF107_smooth(:,1:2),toward_away_date,'rows');
F107_at=DailyF107_smooth(isSBat,:);
F107_ta=DailyF107_smooth(isSBta,:);

hm=160;
ml=100;

% at
ishigh_at=F107_at(:,3)>hm;
ismedium_at=F107_at(:,3)<=hm & F107_at(:,3)>=ml;
islow_at=F107_at(:,3)<ml;
away_toward_dateh=F107_at(ishigh_at,1:2);
away_toward_datem=F107_at(ismedium_at,1:2);
away_toward_datel=F107_at(islow_at,1:2);

% ta
ishigh_ta=F107_ta(:,3)>hm;
ismedium_ta=F107_ta(:,3)<=hm & F107_ta(:,3)>=ml;
islow_ta=F107_ta(:,3)<ml;
toward_away_dateh=F107_ta(ishigh_ta,1:2);
toward_away_datem=F107_ta(ismedium_ta,1:2);
toward_away_datel=F107_ta(islow_ta,1:2);
% solar maximum
% clasify the SB lists according to seasons.
away_toward_date_meh=findSBseason(away_toward_dateh,'me');
away_toward_date_seh=findSBseason(away_toward_dateh,'se');
away_toward_date_jsh=findSBseason(away_toward_dateh,'js');
away_toward_date_dsh=findSBseason(away_toward_dateh,'ds');

toward_away_date_meh=findSBseason(toward_away_dateh,'me');
toward_away_date_seh=findSBseason(toward_away_dateh,'se');
toward_away_date_jsh=findSBseason(toward_away_dateh,'js');
toward_away_date_dsh=findSBseason(toward_away_dateh,'ds');

% solar medium
% clasify the SB lists according to seasons.
away_toward_date_mem=findSBseason(away_toward_datem,'me');
away_toward_date_sem=findSBseason(away_toward_datem,'se');
away_toward_date_jsm=findSBseason(away_toward_datem,'js');
away_toward_date_dsm=findSBseason(away_toward_datem,'ds');

toward_away_date_mem=findSBseason(toward_away_datem,'me');
toward_away_date_sem=findSBseason(toward_away_datem,'se');
toward_away_date_jsm=findSBseason(toward_away_datem,'js');
toward_away_date_dsm=findSBseason(toward_away_datem,'ds');

% solar minimum
% clasify the SB lists according to seasons.
away_toward_date_mel=findSBseason(away_toward_datel,'me');
away_toward_date_sel=findSBseason(away_toward_datel,'se');
away_toward_date_jsl=findSBseason(away_toward_datel,'js');
away_toward_date_dsl=findSBseason(away_toward_datel,'ds');

toward_away_date_mel=findSBseason(toward_away_datel,'me');
toward_away_date_sel=findSBseason(toward_away_datel,'se');
toward_away_date_jsl=findSBseason(toward_away_datel,'js');
toward_away_date_dsl=findSBseason(toward_away_datel,'ds');
