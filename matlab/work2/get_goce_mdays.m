function [uporbit,downorbit]=get_goce_mdays(lrdays)
    % [uporbit,downorbit]=get_rho_mdays_goce(lrdays) returns goce data during specified continuous days
    % input:
    % lrdays=[lyear,ldoy,ryear,rdoy]
    % output:
    % uporbit=[1epochday(from 0),2faltitude,3flon,4flat,5flt,6frho,7fwinde,8fwindn,9fwindu]
    % [] if no data exists
    if isempty(lrdays)
        uporbit=[];
        downorbit=[];
        return
    end

    lyear=lrdays(1);
    ldoy=lrdays(2);
    ryear=lrdays(3);
    rdoy=lrdays(4);
    [lmonth, ldom]=doy2date(lyear,ldoy);
    [rmonth, rdom]=doy2date(ryear,rdoy);

    tmonths=(ryear-lyear)*12+rmonth-lmonth+1;
    yearmonth=ones(tmonths,2);
    yearmonth(:,1)=lyear;
    yearmonth(:,2)=lmonth : (lmonth+tmonths-1);
    index=yearmonth(:,2)>12;
    while sum(index)>=1;
        yearmonth(index,1)=yearmonth(index,1)+1;
        yearmonth(index,2)=yearmonth(index,2)-12;
        index=yearmonth(:,2)>12;
    end

    uporbit=[];
    downorbit=[];
    for k=1:tmonths
        str_year=num2str(yearmonth(k,1));
        str_month=num2str(yearmonth(k,2));
        if yearmonth(k,2)<10
            str_month=['0',str_month];
        end
        filename=['/data/GOCE/data/goce_denswind_v1_3_',...
            str_year,'-',str_month,'.txt'];
        if ~exist(filename,'file')
            continue
        end
        fid=fopen(filename);
        % Column  1:         Date (yyyy-mm-dd)
        % Column  2:         Time (hh:mm:ss.sss)
        % Column  3:         Time system
        % Column  4:  f10.3  Altitude (m)
        % Column  5:   f8.3  Geodetic longitude (deg)
        % Column  6:   f7.3  Geodetic latitude (deg)
        % Column  7:   f6.3  Local solar time (h)
        % Column  8:   f7.3  Argument of latitude (deg)
        % Column  9:  e12.6  Density (kg/m3)
        % Column 10:   f8.2  Eastward component of crosswind speed (m/s)
        % Column 11:   f8.2  Northward component of crosswind speed (m/s)
        % Column 12:   f8.2  Upward component of crosswind speed (m/s)
        % Format string: (a27,1x,f10.3,1X,f8.3,1X,f7.3,1X,f6.3,1X,f7.3,1X,e12.6,1X,f8.2,1X,f8.2,1X,f8.2)
        fdata=textscan(fid,'%s %s %s %f %f %f %f %f %f %f %f %f','commentstyle','#');
        fc=fclose(fid);
        [fdate,ftime,~,...
            faltitude,flon,flat,flt,~,...
            frho,fwinde,fwindn,fwindu]=fdata{:};

        fdate=char(fdate);
        ftime=char(ftime);
        fdate(:,end+1)=' ';
        ftime_temp=[fdate,ftime];
        fmattime=datenum(ftime_temp,'yyyy-mm-dd HH:MM:SS');
        epochtime=fmattime-datenum(lyear,1,ldoy);

        [epochtime,p]=unique(epochtime);
        faltitude=faltitude(p);
        flon=flon(p);
        flat=flat(p);
        flt=flt(p);
        frho=frho(p);
        fwinde=fwinde(p);
        fwindn=fwindn(p);
        fwindu=fwindu(p);

        ceillat=ceil(flat);
        index1=ceillat-flat<=0.25;
        index2=ceillat-flat>0.25 & ceillat-flat<=0.75;
        index3=ceillat-flat>0.75;
        flat(index1)=ceillat(index1);
        flat(index2)=ceillat(index2)-0.5;
        flat(index3)=ceillat(index3)-1;

        diff_lat=diff(flat);
        diff_lat(end+1)=diff_lat(end);
        u0=diff_lat>0;
        d0=diff_lat<0;
        e0=(diff_lat==0);
        % u0 and d0 may not cover all the data

        uporbit=[uporbit;[epochtime(u0),...
            faltitude(u0),flon(u0),flat(u0),flt(u0),...
            frho(u0),fwinde(u0),fwindn(u0),fwindu(u0)]];
        downorbit=[downorbit;[epochtime(d0),...
            faltitude(d0),flon(d0),flat(d0),flt(d0),...
            frho(d0),fwinde(d0),fwindn(d0),fwindu(d0)]];
    end
    if isempty(uporbit) & isempty(downorbit)
        return% should return here 
    end

    lepochtime=0;
    repochtime=datenum(ryear,1,rdoy)-datenum(lyear,1,ldoy)+1;
    if ~isempty(uporbit)
        index=uporbit(:,1)>=lepochtime & uporbit(:,1)<=repochtime;
        uporbit=uporbit(index,:);
    end
    if ~isempty(downorbit)
        index=downorbit(:,1)>=lepochtime & downorbit(:,1)<=repochtime;
        downorbit=downorbit(index,:);
    end
end
