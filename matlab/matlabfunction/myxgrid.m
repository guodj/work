function varargout = myxgrid(varargin )
%XGRID turns on xgrid with lines.
%   Make grids as lines, then they can be set
hold on,
xt=get(gca,'xtick');
yl=ylim;
for ii=1:length(xt)
    hg=plot([xt(ii),xt(ii)],yl,varargin{:});
end
if nargout,varargout{1}=hg; end
end

