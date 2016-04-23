function h=get_contourf_rho_mdays(lrdays,satellite,updown)
    % Draw a contour plot of density in the current axis
    % with ut and latitude as x and y axes
    switch satellite
        case {'champ', 'grace'}
            [rhoup,rhodown]=get_champ_grace_mdays(lrdays,satellite);
        case 'goce'
            [rhoup,rhodown]=get_goce_mdays(lrdays);
        otherwise
            error('wrong input of satellite!');
    end
    switch updown
        case 'up'
            rho=rhoup;
        case 'down'
            rho=rhodown;
        otherwise
            error('wrong input of updown!')
    end
    if isempty(rho)
        h=[];
        return
    end
    switch satellite
        case {'champ', 'grace'}
            x0=rho(:,3);
            y0=rho(:,4);
            z0=rho(:,13);
        case 'goce'
            x0=rho(:,1);
            y0=rho(:,4);
            z0=rho(:,6);
        otherwise
            error('wrong input of satellite!');
    end

    x=0:0.5/24:ceil(max(x0));
    y=-90:3:90;
    [xx,yy]=meshgrid(x,y);
    zz=griddata(x0,y0,z0,xx,yy);
    dx=abs(x(2)-x(1));
    dy=abs(y(2)-y(1));

    for k=1:length(x)
        fp=abs(x(k)-x0)<(2*dx);
        if sum(fp)<1
            zz(:,k)=nan;
        end
    end

    for k=1:length(y)
        fp=abs(y(k)-y0)<(2*dy);
        if sum(fp)<1
            zz(k,:)=nan;
        end
    end
    h=contourf(xx,yy,zz,11,...
        'linecolor','none');
    caxis(prctile(zz(:),[1,99]));

    lt=nanmedian(rho(:,8));
    lt=round(lt*10)/10;
    title([upper(satellite),' LT=',num2str(lt)],'color','red')
end
