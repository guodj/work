% SB distribution due to different seasons, solar 10.7cm flux and
% geoeffectiveness.
%---------------------------------------------------------------%
load SBtype.mat
load DailyF107.mat
% smooth f10.7
DailyF107_smooth=DailyF107;
DailyF107_smooth(:,3)=smooth(DailyF107_smooth(:,3),11);
% 11-day average F10.7 at SBs
isSBat=ismember(DailyF107_smooth(:,1:2),SBAT,'rows');
isSBta=ismember(DailyF107_smooth(:,1:2),SBTA,'rows');
SBATF107=DailyF107_smooth(isSBat,:);
SBTAF107=DailyF107_smooth(isSBta,:);

dtnSBATF107=datenum(SBATF107(:,1),1,SBATF107(:,2));
dtnSBTAF107=datenum(SBTAF107(:,1),1,SBTAF107(:,2));
% nums of SBs for A-T and T-A events
AT_sum=num2str(length(SBAT));
TA_sum=num2str(length(SBTA));
% nums of SBs at different solar activity
strlow_sum=num2str(length(SBF107L));
strmed_sum=num2str(length(SBF107M));
stract_sum=num2str(length(SBF107H));
% nums of SBs at different seasons
me_num=num2str(length(SBME));
se_num=num2str(length(SBSE));
js_num=num2str(length(SBJS));
ds_num=num2str(length(SBDS));
% % nums of SBs associated with CIR
% all_num=num2str(length(SBall));
% CIR_num=num2str(length(SByesCIR));
% rat=num2str(100*(length(SByesCIR)/length(SBall)),'%.0f%%');
% % time interval of CIR
% SBall_l=datenum(SBall(1,1),1,SBall(1,2));
% lenall=length(SBall);
% SBall_r=datenum(SBall(lenall,1),1,SBall(lenall,2));
%% plot,xylim
myplot(dtnSBATF107,SBATF107(:,3),'b.')
hold on
myplot(dtnSBTAF107,SBTAF107(:,3),'r.')
xran=datenum([1967,2008],1,1);
yran=[50,300];
set(gca,'xlim',xran,'ylim',yran)
%% title,label,tick,legend
% title('Distribution of SBs','fontweight','bold')
xlabel('Year')
ylabel('F107')
xticlab=1967:5:2007;
xtic=datenum(xticlab,1,1);
set(gca,'xtick',xtic,'xticklabel',xticlab)
hl=legend(['\color[rgb]{0 0 1}Away-Toward \color[rgb]{0 0 0}SBs ( '...
    ,AT_sum,' )'],['\color[rgb]{1 0 0}Toward-Away \color[rgb]{0 0 0}SBs ( '...
    ,TA_sum,' )']);
set(hl,'box','off',...
    'fontweight','bold','position',[0.18 0.8 0.2 0.1])
%% add reflines
xl=get(gca,'xlim');
yl=get(gca,'ylim');
a=myplot(xl,[100,100],'k--');
set(a,'color',[0.5 0.5 0.5])
a=myplot(xl,[160,160],'k--');
set(a,'color',[0.5 0.5 0.5])
% a=myplot([SBall_l,SBall_l],yl,'--');
% set(a,'color',[0.5 0.5 0.5])
% a=myplot([SBall_r,SBall_r],yl,'--');
% set(a,'color',[0.5 0.5 0.5])
%% text
plot_text(0.7,0.95,['Solar Max:',stract_sum],10)
plot_text(0.7,0.90,['Solar Med:',strmed_sum],10)
plot_text(0.7,0.85,['Solar Min:',strlow_sum],10)
% plot_text(0.76,0.88,['JS: ',js_num,',  ','DS: ',ds_num],10)
% plot_text(0.76,0.93,['ME: ',me_num,', ','SE: ',se_num],10)
% plot_text(0.75,0.8,['SB_{CIR}:',CIR_num,'(',rat,')'],10)
% plot_text(0.75,0.9,['SB_{total}:',all_num],10)
%% minorxy
minorxy([5,5])
%% save
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[226   198   704   581])
dir='D:\mywork\ther_dens_sect\matlabworkspace\figure1\';
nam='figure1.eps';
print(gcf,'-depsc2',[dir,nam])
