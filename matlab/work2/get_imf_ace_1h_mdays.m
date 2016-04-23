function imf_ace=get_imf_ace_1h_mdays(lrdates)
    % Get hourly averaged real-time IMF values from ACE-Magnetometer
    % input:
    % lrdates=[lyear,ldoy,ryear,rdoy]
    % output:
    % imf_ace=[1day(from 0),2Bx,3By,4Bz,5Bt]
    imf_ace=[];
    if isempty(lrdates)
        return
    end

    lyear=lrdates(1);
    ldoy=lrdates(2);
    ryear=lrdates(3);
    rdoy=lrdates(4);
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

    imf_ace_tmp=[];
    for k=1:tmonths
        str_year=num2str(yearmonth(k,1));
        str_month=num2str(yearmonth(k,2));
        if yearmonth(k,2)<10
            str_month=['0',str_month];
        end
        filename=['/data/ACE_mag_1h/',str_year,str_month,'_ace_mag_1h.txt'];
        if ~exist(filename,'file')
            continue
        end
        fid=fopen(filename);
        % Column  1:         Year
        % Column  2:         Month
        % Column  3:         Day of Month
        % Column  4:         Time: HHMM
        % Column  5:         Modified Julian day
        % Column  6:         Seconds of the day
        % Column  7:         status: 0=normal, 1 to 8 = bad data record, 9 = no data
        % Column  8:         GSM Bx
        % Column  9:         GSM By
        % Column 10:         GSM Bz
        % Column 11:         Bt
        % Column 12:         Lat.
        % Column 13:         Long.
        fdata=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f','headerlines',20);
        fc=fclose(fid);
        [fyear,fmonth,fdom,fhour,~,~,...
            status,fbx,fby,fbz,fbt,flat,flong]=fdata{:};
        fhour=fhour/100;
        imf_ace_tmp=[imf_ace_tmp;[fyear,fmonth,fdom,fhour,status,fbx,fby,fbz,fbt]];
    end

    lindex=find(imf_ace_tmp(:,1)==lyear & imf_ace_tmp(:,2)==lmonth &...
        imf_ace_tmp(:,3)==ldom,1,'first');
    rindex=find(imf_ace_tmp(:,1)==ryear & imf_ace_tmp(:,2)==rmonth &...
        imf_ace_tmp(:,3)==rdom,1,'last');
    imf_ace_tmp=imf_ace_tmp(lindex:rindex,:);

    imf_ace_day=datenum(imf_ace_tmp(:,1:3))-datenum(lyear,lmonth,ldom)+imf_ace_tmp(:,4)/24;
    imf_ace=[imf_ace_day,imf_ace_tmp(:,6:9)];
    index=imf_ace_tmp(:,5)~=0;
    imf_ace(index,2:5)=nan;
end
