% use iptrussia to determine the SB crossing dates
% 3 kinds of polarities may exist in the IMF polarity table: 'away' marked 
% as 1, 'toward' marked as -1, neither of the two polarities marked as 0, 
%
% input: iptrussia, 56 columns: year, doy, 27-day polarities. 
% output: away_toward_date, toward_away_date, 2 columns: year,doy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ipt1=iptrussia(:,3:56);                              %remove dates,leave behind sector polarities.
len=size(ipt1,1)*size(ipt1,2);
ipt2=zeros(len,3);
ipt2(:,3)=reshape(ipt1',len,1);
away_toward_date=[];toward_away_date=[];                  % output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set dates for ipt2. 
year=iptrussia(1,1); doy=iptrussia(1,2);             %begining date
for m=1:len
    if mod(m,2)==1
        ipt2(m,1)=year;ipt2(m,2)=doy;
        doy=doy+1;
        if doy==days(year)+1;
            doy=1;year=year+1;
        end
    elseif mod(m,2)==0
        ipt2(m,1)=ipt2(m-1,1);ipt2(m,2)=ipt2(m-1,2);
    end
end
% until this step,ipt2 is the same as iptrussia except format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  count SB crossing dates
p=1;          % p points to ipt2
while p<=length(ipt2)-19
    if sum(ipt2(p:p+9,3)==1)>=8 && sum(ipt2(p+10:p+19,3)==-1)>=8 ...
            && sum(ipt2(p+8:p+9,3)==-1)==0 && sum(ipt2(p+10:p+11,3)==1)==0
        away_toward_date=cat(1,away_toward_date,ipt2(p+10,1:2));
        p=p+8;
    elseif sum(ipt2(p:p+9,3)==-1)>=8 && sum(ipt2(p+10:p+19,3)==1)>=8 ...
            && sum(ipt2(p+8:p+9,3)==1)==0 && sum(ipt2(p+10:p+11,3)==-1)==0
        toward_away_date=cat(1,toward_away_date,ipt2(p+10,1:2));
        p=p+8;
    else
        p=p+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clasify the SB lists according to seasons.
away_toward_date_me=findSBseason(away_toward_date,'me');
away_toward_date_se=findSBseason(away_toward_date,'se');
away_toward_date_js=findSBseason(away_toward_date,'js');
away_toward_date_ds=findSBseason(away_toward_date,'ds');

toward_away_date_me=findSBseason(toward_away_date,'me');
toward_away_date_se=findSBseason(toward_away_date,'se');
toward_away_date_js=findSBseason(toward_away_date,'js');
toward_away_date_ds=findSBseason(toward_away_date,'ds');