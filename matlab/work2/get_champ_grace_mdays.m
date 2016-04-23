function [uporbit,downorbit]=get_champ_grace_mdays(lrdays,satellite)
    % [uporbit,downorbit]=get_rho_mdays(lrdays,satellite) returns champ or grace data during specified days
    % input:
    % lrdays=[lyear,ldoy,ryear,rdoy]
    % satellite='champ' or 'grace'
    % output:
    % uporbit=[1year,2doy,3sec/24/3600 ( from 0 ),...
    %     4lat3,5lat,6lon,7height,8lt,9mlat,10mlon,11mlt,...
    %     12rho,13rho400,14rho410,15rhomsis];
    % note: the date range should be small so that the up and down orbits can be in the same lt.
    if isempty(lrdays)
        uporbit=[];
        downorbit=[];
        return
    end

    lyear=lrdays(1);
    ldoy=lrdays(2);
    ryear=lrdays(3);
    rdoy=lrdays(4);

    total_days=datenum(ryear,1,rdoy)-datenum(lyear,1,ldoy)+1;
    yeardoy=ones(total_days,2);
    yeardoy(:,1)=lyear;
    yeardoy(:,2)=ldoy : (ldoy+total_days-1);
    days_lyear=yeardays(lyear);
    index=yeardoy(:,2)>days_lyear;
    yeardoy(index,1)=lyear+1;
    yeardoy(index,2)=yeardoy(index,2)-days_lyear;

    uporbit=[];
    downorbit=[];
    for k=1:total_days
        str_year=num2str(yeardoy(k,1));
        str_doy=num2str(yeardoy(k,2));
        if yeardoy(k,2)<10
            str_doy=['00',str_doy];
        elseif yeardoy(k,2)<100
            str_doy=['0',str_doy];
        end
        switch satellite
            case 'champ'
            filename=['/data/CHAMP23/',str_year,...
                '/ascii/Density_3deg_' str_year(3:4) '_' str_doy '.ascii'];
            case 'grace'
            filename=['/data/Grace23/',str_year,...
                '/ascii/Density_graceA_3deg_' str_year(3:4) '_' str_doy '.ascii'];
            otherwise
                error('wrong input of satellite!');
        end
        if ~exist(filename,'file')
            continue
        end

        fid=fopen(filename);
        %1,Two-digit Year (years);
        %2,Day of the Year (days);
        %3,Second of the Day (GPS time,sec);
        %4,Center Latitude of 3-degree Bin (deg);
        %5,Satellite Geodetic Latitude (deg);
        %6,Satellite Longitude (deg);
        %7,Satellite Height (km);
        %8,Satellite Local Time (hours);
        %9,Satellite Quasi-Dipole Latitude (deg);
        %10,Satellite Magnetic Longitude (deg);
        %11,Satellite Magnetic Local Time (hours)
        %12,Neutral Density (kg/m^3);
        %13,Neutral Density Normalized to 400km using NRLMSISe00;
        %14,Neutral Density Normalized to 410km using NRLMSISe00;
        %15,NRLMSISe00 Neutral Density at Satellite Height;
        %16,Uncertainty in Neutral Density (kg/m^3);
        %17,Number of Data Points in Current Averaging Bin;
        %18,Number of Points in Current Averaging Bin that Required Interpolation;
        %19,Average Coefficient of Drag Used in Current Averaging Bin
        fdata=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d %d %f',...
            'headerlines',2);
        fc=fclose(fid);
        [fyear,fdoy,fsec,...
            flat3,flat,flon,fheight,flt,fmlat,fmlon,fmlt,...
            frho,frho400,frho410,frhomsis,~,~,~,~]...
            =fdata{:};
        fyear=fyear+2000;
        [fsec,p]=unique(fsec);
        fyear=fyear(p);
        fdoy=fdoy(p);
        flat3=flat3(p);
        flat=flat(p);
        flon=flon(p);
        fheight=fheight(p);
        flt=flt(p);
        fmlat=fmlat(p);
        fmlon=fmlon(p);
        fmlt=fmlt(p);
        frho=frho(p);
        frho400=frho400(p);
        frho410=frho410(p);
        frhomsis=frhomsis(p);
        fsec=fsec+(k-1)*24*3600;

        diff_lat3=diff(flat3);
        diff_lat3(end+1)=diff_lat3(end);
        u0=diff_lat3>0;
        d0=diff_lat3<0;
        e0=(diff_lat3==0);
        switch satellite
            case 'champ'
                u0(e0 & flat3==87)=true;
                d0(e0 & flat3==-87)=true;
            case 'grace'
                u0(e0 & flat3==90)=true;
                d0(e0 & flat3==-90)=true;
            otherwise
                error('wrong input of satellite!');
        end

        % u0 and d0 may not cover all the data
        uporbit=[uporbit;[fyear(u0),fdoy(u0),fsec(u0)/24/3600,...
            flat3(u0),flat(u0),flon(u0),fheight(u0),flt(u0),fmlat(u0),fmlon(u0),fmlt(u0),...
            frho(u0),frho400(u0),frho410(u0),frhomsis(u0)]];
        downorbit=[downorbit;[fyear(d0),fdoy(d0),fsec(d0)/24/3600,...
            flat3(d0),flat(d0),flon(d0),fheight(d0),flt(d0),fmlat(d0),fmlon(d0),fmlt(d0),...
            frho(d0),frho400(d0),frho410(d0),frhomsis(d0)]];
    end
end
