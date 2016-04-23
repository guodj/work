function hvl=vline(x0)
% hvl=vline(x0) draws a vertial line: x=x0
yl=get(gca,'ylim');
h=myplot(gca,[x0,x0],[-1e16,1e16],'--');
set(gca,'ylim',yl)
if nargout , hvl=h; end
end