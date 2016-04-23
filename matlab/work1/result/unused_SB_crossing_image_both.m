% smooth ipt and iptrussia.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smooth ipt
ipt_a=ipt(:,3:29);ipt_a(ipt_a==9)=0;
% '0' and '9' cases are the same in this process.
ipt_a=reshape(ipt_a',numel(ipt_a),1);
% matlab arranges matrix according to columns, so ipt_a' is used.
ipt_a=smooth(ipt_a,5,'moving');
ipt_a(ipt_a>0)=1;ipt_a(ipt_a<0)=-1;
% find SB crossing positions in ipt_smooth
pta=[];pat=[];
p=1;     % p points to ipt_a
while p<=length(ipt_a)-9       
    if sum(ipt_a(p:p+4)==1)>=4 && sum(ipt_a(p+5:p+9)==-1)>=4 ...
            && ipt_a(p+4)~=-1 && ipt_a(p+5)~=1
        pat=cat(1,pat,p+5);%away-toward date
        p=p+4;
    elseif sum(ipt_a(p:p+4)==-1)>=4 && sum(ipt_a(p+5:p+9)==1)>=4 ...
            && ipt_a(p+4)~=1 && ipt_a(p+5)~=-1
        pta=cat(1,pta,p+5);%toward-away date
        p=p+4;
    else
        p=p+1;
    end
end
%%%%%%%%%%%%%%%%%%%%
zz1=ipt_a;zz1(:,:)=0;zz1(pat)=1;zz1(pta)=-1;
zz1=reshape(zz1,27,numel(zz1)/27)';
[y1,x1]=find(zz1==1);[y2,x2]=find(zz1==-1);
%%%%%%%%%%%%%%%%%%%%%%%%
ipt_a=reshape(ipt_a,27,numel(ipt_a)/27)';
ipt_smooth=[ipt(:,1:2),ipt_a];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% smooth iptrussia
iptrussia_a=iptrussia(:,3:56);
iptrussia_a=reshape(iptrussia_a',numel(iptrussia_a),1);
% matlab arranges matrix according to columns, so iptrussia_a' is used.
iptrussia_a=smooth(iptrussia_a,10,'moving');
iptrussia_a(iptrussia_a>0)=1;iptrussia_a(iptrussia_a<0)=-1;
% find SB crossing positions in iptrussia_smooth
patrussia=[];ptarussia=[];
p=1;          % p points to iptrussia_a
while p<=length(iptrussia_a)-19
    if sum(iptrussia_a(p:p+9)==1)>=8 && sum(iptrussia_a(p+10:p+19)==-1)>=8 ...
            && sum(iptrussia_a(p+8:p+9)==-1)==0 && sum(iptrussia_a(p+10:p+11)==1)==0
        patrussia=[patrussia;p+10];
        p=p+8;
    elseif sum(iptrussia_a(p:p+9)==-1)>=8 && sum(iptrussia_a(p+10:p+19)==1)>=8 ...
            && sum(iptrussia_a(p+8:p+9)==1)==0 && sum(iptrussia_a(p+10:p+11)==-1)==0
        ptarussia=[ptarussia;p+10];
        p=p+8;
    else
        p=p+1;
    end
end
%%%%%%%%%%%%%%%%%%%%
zz2=iptrussia_a;zz2(:,:)=0;zz2(patrussia)=1;zz2(ptarussia)=-1;
zz2=reshape(zz2,54,numel(zz2)/54)';
[y11,x11]=find(zz2==1);[y22,x22]=find(zz2==-1);
%%%%%%%%%%%%%%%%%%%%%%%%
iptrussia_a=reshape(iptrussia_a,54,numel(iptrussia_a)/54)';
iptrussia_smooth=[iptrussia(:,1:2),iptrussia_a];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin to draw pictures.
figure(1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iptomni
% color index of image() function begins from 1
h1=subplot(1,2,1);
ipt_smooth(:,3:29)=ipt_smooth(:,3:29)+2;
image(ipt_smooth(:,3:29))

colormap([0 0 0;.5 .5 .5;1 1 1]);
shading flat

set(gca,'xtick',5:5:25,'xticklabel',5:5:25)

title('Satellite Data')
ylabel('Year'),xlabel('Time(Day)')

hold on
plot(x1,y1,'b.')
plot(x2,y2,'r.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% iptrussia
h2=subplot(1,2,2);
iptrussia_smooth(:,3:56)=iptrussia_smooth(:,3:56)+2;
image(iptrussia_smooth(:,3:56))

colormap([0 0 0;.5 .5 .5;1 1 1])
shading flat

set(gca,'xtick',9:10:54,'xticklabel',5:5:25)

title('Ground Data')
xlabel('Time(Day)')

hold on
plot(x11,y11,'b.')
plot(x22,y22,'r.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set yaxis limit
year1=1975;year2=1980;
p1b=find(ipt_smooth(:,1)==year1,1);p1e=find(ipt_smooth(:,1)==year2,1,'last');
p2b=find(iptrussia_smooth(:,1)==year1,1);p2e=find(iptrussia_smooth(:,1)==year2,1,'last');
for p=year1:year2
    iptytick(p-year1+1)=find(ipt(:,1)==p,1);
    iptrussiaytick(p-year1+1)=find(iptrussia(:,1)==p,1);
end
ylim(h1,[p1b,p1e])
set(h1,'ytick',iptytick,'yticklabel',year1:year2)
ylim(h2,[p2b,p2e])
set(h2,'ytick',iptrussiaytick,'yticklabel',year1:year2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveas(gcf,'/home/gdj/study/matlabworkspace/graduation_project/picture/crossing_dates_both_smooth.eps','psc2');

