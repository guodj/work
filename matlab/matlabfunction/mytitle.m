function tit = mytitle( varargin)
%MYTITLE is the same as title,but change some defaults
% title position is a ratio(>1) of ylim
    for ii=1:length(varargin)
        if strcmp(varargin{ii},'position') || strcmp(varargin{ii},'Position') 
            tit_pos=varargin{ii+1};
            if length(tit_pos)==1, tit_pos=[0.5,tit_pos]; end
            xl=get(gca,'xlim');
            yl=get(gca,'ylim');
            tit_pos(1)=xl(1)+tit_pos(1)*(xl(2)-xl(1));
            tit_pos(2)=yl(1)+tit_pos(2)*(yl(2)-yl(1));
            varargin{ii+1}=[tit_pos(1),tit_pos(2)];
            break
        end
    end
    tit=title(varargin{:});
    set(tit,'fontweight','bold','fontsize',11)
end

