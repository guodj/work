% test F10.7 response to solar activity for different F10.7
sector_polarity_change_dates_SBlist
diurnal_avarages_of_f107
diurnal_averages_of_Ap
% smooth F107 with span 11
DailyF107_smooth=DailyF107;
DailyF107_smooth(:,3)=smooth(DailyF107_smooth(:,3),11);

isSBat=ismember(DailyF107_smooth(:,1:2),away_toward_date,'rows');
isSBta=ismember(DailyF107_smooth(:,1:2),toward_away_date,'rows');
DailyF107_at=DailyF107_smooth(isSBat,:);
DailyF107_ta=DailyF107_smooth(isSBta,:);

hm=160;
ml=100;

ishigh_at=DailyF107_at(:,3)>hm;
ismedium_at=DailyF107_at(:,3)<=hm & DailyF107_at(:,3)>=ml;
islow_at=DailyF107_at(:,3)<ml;
away_toward_date_h=DailyF107_at(ishigh_at,1:2);
away_toward_date_m=DailyF107_at(ismedium_at,1:2);
away_toward_date_l=DailyF107_at(islow_at,1:2);

ishigh_ta=DailyF107_ta(:,3)>hm;
ismedium_ta=DailyF107_ta(:,3)<=hm & DailyF107_ta(:,3)>=ml;
islow_ta=DailyF107_ta(:,3)<ml;
toward_away_date_h=DailyF107_ta(ishigh_ta,1:2);
toward_away_date_m=DailyF107_ta(ismedium_ta,1:2);
toward_away_date_l=DailyF107_ta(islow_ta,1:2);

[F107_ath,F107_ath_num]=resp2SB_mean(away_toward_date_h,DailyF107);
[F107_atm,F107_atm_num]=resp2SB_mean(away_toward_date_m,DailyF107);
[F107_atl,F107_atl_num]=resp2SB_mean(away_toward_date_l,DailyF107);

[F107_tah,F107_tah_num]=resp2SB_mean(toward_away_date_h,DailyF107);
[F107_tam,F107_tam_num]=resp2SB_mean(toward_away_date_m,DailyF107);
[F107_tal,F107_tal_num]=resp2SB_mean(toward_away_date_l,DailyF107);

[Ap_ath,Ap_ath_num]=resp2SB_mean(away_toward_date_h,DailyAp);
[Ap_atm,Ap_atm_num]=resp2SB_mean(away_toward_date_m,DailyAp);
[Ap_atl,Ap_atl_num]=resp2SB_mean(away_toward_date_l,DailyAp);

[Ap_tah,Ap_tah_num]=resp2SB_mean(toward_away_date_h,DailyAp);
[Ap_tam,Ap_tam_num]=resp2SB_mean(toward_away_date_m,DailyAp);
[Ap_tal,Ap_tal_num]=resp2SB_mean(toward_away_date_l,DailyAp);

figure(1)
subplot(1,3,1)
plot(-5:5,F107_ath(2,:),-5:5,F107_tah(2,:))
title('solar maximum')
ylim([-10,10])

subplot(1,3,2)
plot(-5:5,F107_atm(2,:),-5:5,F107_tam(2,:))
title('solar medium')
ylim([-10,10])

subplot(1,3,3)
plot(-5:5,F107_atl(2,:),-5:5,F107_tal(2,:))
title('solar minimum')
ylim([-10,10])

figure(2)
subplot(1,3,1)
plot(-5:5,Ap_ath(2,:),-5:5,Ap_tah(2,:))
title('solar maximum')
ylim([-5,5])

subplot(1,3,2)
plot(-5:5,Ap_atm(2,:),-5:5,Ap_tam(2,:))
title('solar medium')
ylim([-5,5])

subplot(1,3,3)
plot(-5:5,Ap_atl(2,:),-5:5,Ap_tal(2,:))
title('solar minimum')
ylim([-5,5])






