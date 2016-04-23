function varargout = myplot(varargin )
%------------------------------------------------------------------------%
% MYPLOT is the same as 'plot' but with linewidth 1.5
% It is suggested one plot contains one line.
%------------------------------------------------------------------------%
    pic=plot(varargin{:},'linewidth',1.5);
    if nargout==1,varargout{1}=pic; end
end

