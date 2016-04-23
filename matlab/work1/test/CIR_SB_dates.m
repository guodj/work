% test the relationship between CIR dates and SB dates
load streaminterface.mat
load SBlist67_07.mat
stream_datenum1=datenum(streaminterface1,'dd-mmm-yyyy');
% stream_datenum2=datenum(streaminterface2,'dd-mmm-yyyy');

densityenddate=datenum(2007,1,279);
% limdate=find(stream_datenum2>max(stream_datenum1) &...
%     stream_datenum2<densityenddate);
% stream_datenum2=stream_datenum2(limdate);
% CIR_datenum=[stream_datenum1;stream_datenum2];
CIR_datenum=stream_datenum1;

SB_datenum=datenum(SBlist(:,1),1,SBlist(:,2));
% CIR_datenum=unique(CIR_datenum);
myday=5;
diffday=-myday:myday;
for i=-myday:myday
    iscoin=ismember(SB_datenum,CIR_datenum+i);
    sumcoin=sum(iscoin);
    diffday(myday+1-i)=sumcoin;
end
plot(-myday:myday,diffday,'k*:','linewidth',1.5)
set(gca,'fontsize',12,'xtick',-5:5)
xlabel('Epoch Time (Days)')
ylabel('Number of CIRs')
text(-4,80,['Total CIRs: ',num2str(size(CIR_datenum,1),'%3d')])
set(gcf,'paperpositionmode','auto')
figdir=cd;
print(gcf,'-depsc2',[figdir,'\CIR_dis.eps'])
