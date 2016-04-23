% test the variations of F10.7 during CIRs
clear
clc
load 'final.mat'
load 'streaminterface.mat'
% 1, obtain daily average F10.7:DailyF107
diurnal_avarages_of_f107
% smooth f107
DailyF107(:,3)=smooth(DailyF107(:,3),5);
% 2, obtain datenum of DailyF107, CIRs and SBs.
dn_F107=datenum(DailyF107(:,1),1,DailyF107(:,2));
dn_CIR=datenum(streaminterface1);
dn_SB=datenum(SBlist(:,2),1,SBlist(:,3));
% 3,find dates of CIRs in f107
find_dnCIRinf107=find(ismember(dn_F107,dn_CIR));
interval=15
find_dnCIRinf107(:,2)=find_dnCIRinf107(:,1)-interval;
find_dnCIRinf107(:,3)=find_dnCIRinf107(:,1)+interval;
dele_row=find_dnCIRinf107(:,2)<=0 | find_dnCIRinf107(:,3)>length(dn_F107);
find_dnCIRinf107(dele_row,:)=[];
% 4,F107 response to CIR
f1072CIR=[];
for irow=1:length(find_dnCIRinf107)
    addcol=DailyF107(find_dnCIRinf107(irow,2):find_dnCIRinf107(irow,3),3);
    f1072CIR=[f1072CIR,addcol];
end
plot(-interval:interval,nanmean(f1072CIR,2))
% ylim([0,130])
% conclusion: F10.7 does not show strong response to CIRs or SBs.
