function solar_wind_parameters=get_solar_wind_parameters_mdays(lrdates)
    % solar_wind_parameters=get_solar_wind_mdays(lrdates) returns solar wind parameter values during specified days
    % input:
    % lrdates=[lyear,ldoy,ryear,rdoy]
    % output:
    % solar_wind_parameters=[1day, 2temperature, 3density, 4speed, 5longitude, 6latitude, 7pressure];

    lyear=lrdates(1);
    ldoy=lrdates(2);
    ryear=lrdates(3);
    rdoy=lrdates(4);

    year=[];
    doy=[];
    hour=[];
    temperature=[];
    density=[];
    speed=[];
    longitude=[];
    latitude=[];
    pressure=[];

    for k=lyear:ryear
        str_year=num2str(k);
        fname=['/data/SW_IMF_1h/SW_l1/plasma', str_year, '.txt'];
        if ~exist(fname,'file')
            continue
        end
        fid = fopen(fname);
        fdata = textscan(fid,'%f %f %f %f %f %f %f %f %f');
        fc = fclose(fid);
        [fyear,fdoy,fhour,ftemperature,fdensity,fspeed,flongitude,flatitude,fpressure] = fdata{:};
        ftemperature(ftemperature==9999999.) = nan;
        fdensity(fdensity==999.9) = nan;
        fspeed(fspeed==9999.)=nan;
        flongitude(flongitude==999.9) = nan;
        flatitude(flatitude==999.9) = nan;
        fpressure(fpressure==99.99) = nan;

        year=[year;fyear];
        doy=[doy;fdoy];
        hour=[hour;fhour];
        temperature=[temperature;ftemperature];
        density=[density; fdensity];
        speed=[speed;fspeed];
        longitude=[longitude; flongitude];
        latitude=[latitude; flatitude];
        pressure=[pressure; fpressure];
    end

    if isempty(year)
        solar_wind_parameters=[];
        return
    end

    lp=find( (year==lyear & doy==ldoy), 1, 'first');
    rp=find( (year==ryear & doy==rdoy), 1, 'last');
    year=year(lp:rp);
    doy=doy(lp:rp);
    hour=hour(lp:rp);
    temperature=temperature(lp:rp);
    density=density(lp:rp);
    speed=speed(lp:rp);
    longitude=longitude(lp:rp);
    latitude=latitude(lp:rp);
    pressure=pressure(lp:rp);

    day=datenum(year,1,doy)-datenum(lyear,1,ldoy)+hour/24;
    solar_wind_parameters=[day,temperature,density,speed,longitude,latitude,pressure];
end
