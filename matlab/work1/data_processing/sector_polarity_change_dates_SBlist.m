% use the alrady found SB lists to determine the SB dates.
% output:away_toward_date_me ... SB_ie_eq...
SBtype=menu('choose a SB type','1):all','2):ynCIR',...
    '3):noCIR','4):yesCIR');
switch SBtype
    case 2
        find_SBlist_CIR
        mySB=SBlist_all;
    case 3
        find_SBlist_CIR
        mySB=SBlistnoCIR;
    case 4
        find_SBlist_CIR
        mySB=SBlistyesCIR;
    case 1
        mySB=SBlist;
    otherwise
        error('wrong input')
end
away_toward_date=mySB(mySB(:,1)==-1,2:3);
toward_away_date=mySB(mySB(:,1)==1,2:3);
% % % % % % % % away_toward_date=indep_SBlist(indep_SBlist(:,3)==-1,1:2);
% % % % % % % % toward_away_date=indep_SBlist(indep_SBlist(:,3)==1,1:2);
% clasify the SBs according to seasons.
away_toward_date_me=findSBseason(away_toward_date,'me');
away_toward_date_se=findSBseason(away_toward_date,'se');
away_toward_date_js=findSBseason(away_toward_date,'js');
away_toward_date_ds=findSBseason(away_toward_date,'ds');

toward_away_date_me=findSBseason(toward_away_date,'me');
toward_away_date_se=findSBseason(toward_away_date,'se');
toward_away_date_js=findSBseason(toward_away_date,'js');
toward_away_date_ds=findSBseason(toward_away_date,'ds');
% clasify SBs due to geoeffectiveness.
SB_ie_eq=[away_toward_date_me;toward_away_date_se];
SB_ie_so=[away_toward_date_js;toward_away_date_ds];
SB_ei_eq=[toward_away_date_me;away_toward_date_se];
SB_ei_so=[toward_away_date_js;away_toward_date_ds];
