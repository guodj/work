% add legend to the current axis(gca)
function varargout = mylegend(varargin )
hleg=legend(varargin{1:end-2});
if nargout==1, varargout{1}=hleg; end
sizeleg=get(hleg,'position');
pleg=get(gca,'position');
xf=varargin{end-1};
yf=varargin{end};
set(hleg,'position',[pleg(1)+xf*pleg(3),pleg(2)+yf*pleg(4),sizeleg(3),sizeleg(4)])
end
