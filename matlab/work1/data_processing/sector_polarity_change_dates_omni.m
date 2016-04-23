% 4 kinds of polarities may exist in the IMF polarity table: 'away' marked 
% as 1, 'toward' marked as -1, neither of the two polarities marked as 0, 
% no data marked as NaN or 9.
%
% input: ipt, 29 columns: year, doy, 27-day polarities. 
% output: away_toward_date, toward_away_date, 2 columns: year,doy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ipt1=ipt(:,3:29);                                           %remove dates,leave behind sector polarities.
ipt2=zeros(numel(ipt1),3);
% running average of ipt2
ipt2(:,3)=reshape(ipt1',numel(ipt1),1);
ipt2(ipt2(:,3)==9,3)=0;
ipt2(:,3)=smooth(ipt2(:,3),5,'moving');
ipt2(ipt2(:,3)>0,3)=1;ipt2(ipt2(:,3)<0,3)=-1;     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set dates for ipt2. 
year=ipt(1,1); doy=ipt(1,2);        %begining date
for p=1:numel(ipt1)
    ipt2(p,1)=year;ipt2(p,2)=doy;
    doy=doy+1;
    if doy==days(year)+1;
        doy=1;year=year+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
away_toward_date=[];toward_away_date=[];  % output
p=1;     % p points to ipt2
while p<=length(ipt2)-9       
    if sum(ipt2(p:p+4,3)==1)>=4 && sum(ipt2(p+5:p+9,3)==-1)>=4 ...
            && ipt2(p+4,3)~=-1 && ipt2(p+5,3)~=1
        away_toward_date=[away_toward_date;ipt2(p+5,1:2)];%away-toward date
        p=p+4;
    elseif sum(ipt2(p:p+4,3)==-1)>=4 && sum(ipt2(p+5:p+9,3)==1)>=4 ...
            && ipt2(p+4,3)~=1 && ipt2(p+5,3)~=-1
        toward_away_date=[toward_away_date;ipt2(p+5,1:2)];%toward-away date
        p=p+4;
    else
        p=p+1;
    end
end
% clasify the SB lists according to seasons.
away_toward_date_me=findSBseason(away_toward_date,'me');
away_toward_date_se=findSBseason(away_toward_date,'se');
away_toward_date_js=findSBseason(away_toward_date,'js');
away_toward_date_ds=findSBseason(away_toward_date,'ds');

toward_away_date_me=findSBseason(toward_away_date,'me');
toward_away_date_se=findSBseason(toward_away_date,'se');
toward_away_date_js=findSBseason(toward_away_date,'js');
toward_away_date_ds=findSBseason(toward_away_date,'ds');



