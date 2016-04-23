function [ha,h1,h2,h3,h4,h5,h6]=get_fig_1234(epoch,lrdays,satellite,updown)
    epoch=get_dates_satellite(epoch,satellite);
    if isempty(epoch)
        return
    end
    lrdates=epoch2lrdates(epoch,lrdays);

    lyear=lrdates(1);
    ldoy=lrdates(2);
    ryear=lrdates(3);
    rdoy=lrdates(4);

    total_days=datenum(ryear,1,rdoy)-datenum(lyear,1,ldoy)+1;
    yeardoy=ones(total_days+1,2);
    yeardoy(:,1)=lyear;
    yeardoy(:,2)=ldoy : (ldoy+total_days);
    days_lyear=yeardays(lyear);
    index=yeardoy(:,2)>days_lyear;
    yeardoy(index,1)=lyear+1;
    yeardoy(index,2)=yeardoy(index,2)-days_lyear;
    [month,dom]=doy2date(yeardoy);
    xtickstr=cell(total_days+1,1);
    for k=1:total_days+1
        xtickstr{k,1}=[num2str(month(k)),'/', num2str(dom(k))];
    end
    figure
    set(gcf,'position',[64   340   706   533])
    ha=sub_creat_axes(3,1,[0.05,0.15,0.1,0.1],[0.05,0.1]);
    %% Bx and By variations
    set(gcf,'currentaxes',ha(1))
    hold on
    box on
    imf=get_imf_mdays(lrdates);
    xut=imf(:,1);
    bx=imf(:,2);
    bym=imf(:,5);
    bzm=imf(:,6);
    bt=sqrt(bx.^2+bym.^2+bzm.^2);
    % to feature the daily variation, the imf data are smoothed.
    h1=myplot(xut,bx,'r');
    h2=myplot(xut,bym,'b');
    h3=myplot(xut,bt,'k');
    set(gca,'xlim',[0,lrdays(1)+lrdays(2)+1],'xtick',0:(lrdays(1)+lrdays(2)+1),'xticklabel',[])
    ylabel('GSM {\color{red}Bx}, {\color{blue}By}, {\color{black}Bt}')
    %------------------------Bz variation, added on Oct. 19--------------------------------------
    %% Bz variation
    set(gcf,'currentaxes',ha(2))
    hold on
    box on
    imf=get_imf_mdays(lrdates);
    xut=imf(:,1);
    bx=imf(:,2);
    bye=imf(:,3);
    bze=imf(:,4);
    bym=imf(:,5);
    bzm=imf(:,6);
    % to feature the daily variation, the imf data are smoothed.
    h4=myplot(xut,bzm,'r');
    h5=myplot(xut,bze,'b');
    set(gca,'xlim',[0,lrdays(1)+lrdays(2)+1],'xtick',0:(lrdays(1)+lrdays(2)+1),'xticklabel',[])
    ylabel('{\color{red}GSM}, {\color{blue}GSE} Bz')
    %--------------------------------------------------------------------------
    %% rho variation
    set(gcf,'currentaxes',ha(3));
    h6=get_contourf_rho_mdays(lrdates,satellite,updown);
    set(gca,'xlim',[0,lrdays(1)+lrdays(2)+1],'xtick',0:(lrdays(1)+lrdays(2)+1),'xticklabel',xtickstr)
    xlabel('Date of 2003')
    ylabel('Latitude (deg)')
end
