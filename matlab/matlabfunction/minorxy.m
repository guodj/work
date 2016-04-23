function minorxy(num,lengthf)
% MINORXY add minortick to the axes
%% Input:
% num is the number of minortick 1x2;
% Lengthf is the factor of parent ticklabel for minortick;
% Example:
%     figure; plot(1:10);
%     minorxy([10 5]);
% writed by Lei J. H., Oct 2005.
if nargin==1,lengthf=0.6; end
if length(num)==1, num=[num num]; end

%
h1=gca;set(h1,'position',get(h1,'position'));
xk=get(h1,'xtick');yk=get(h1,'ytick');
x1=xlim; y1=ylim;
fp=xk<x1(1) | xk>x1 (end); xk(fp)=[];
fp=yk<y1(1) | yk>y1 (end); yk(fp)=[];

xb=xk(1);xe=xk(end);
yb=yk(1);ye=yk(end);
%%
if xb>x1(1), aa=xk(2):-(xk(2)-xk(1))/num(1):x1(1); xb=aa(end); end
if xe<x1(end), aa=xk(1):(xk(2)-xk(1))/num(1):x1(end); xe=aa(end); end
if yb>y1(1), aa=yk(2):-(yk(2)-yk(1))/num(2):y1(1); yb=aa(end); end
if ye<y1(end), aa=yk(1):(yk(2)-yk(1))/num(2):y1(end); ye=aa(end); end
if num(1)>1,xk=xb:(xk(2)-xk(1))/num(1):xe;else xk=[]; end
if num(end)>1,yk=yb:(yk(2)-yk(1))/num(end):ye;else yk=[]; end
len=get(gca,'Ticklength');length_manul=len*lengthf;
%% add axes;
h2=axes('position',get(h1,'position'),'color','none');
set(h2,'xlim',get(h1,'xlim'),'ylim',get(h1,'ylim'),'xtick',xk,'ytick',yk,'xticklabel',[],'yticklabel',[],'ycolor',get(h1,'ycolor'),'xcolor',get(h1,'xcolor'),'ticklength',length_manul,'box',get(h1,'box'),'Clipping',get(h1,'Clipping'));
set(h1,'layer','top')
set(gcf,'currentaxes',h1)
%% use line to add minor xy.

%%
%h2=axes('position',get(h1,'position'),'color','none');
%set(h2,'xlim',get(h1,'xlim'),'ylim',get(h1,'ylim'),'xtick',xk,'ytick',yk,'xticklabel',[],'yticklabel',[],'ticklength',length_manul)
%set(gca,'XAxisLocation','top','YAxisLocation','right');
