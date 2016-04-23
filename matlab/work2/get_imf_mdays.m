function imf=get_imf_mdays(lrdates)
    % imf=get_imf_mdays(lrdates) returns imf values during specified days
    % input:
    % lrdates=[lyear,ldoy,ryear,rdoy]
    % output:
    % imf=[1day(from 0),2bx,3bye,4bze,5bym,6bzm]

    lyear=lrdates(1);
    ldoy=lrdates(2);
    ryear=lrdates(3);
    rdoy=lrdates(4);

    year=[];
    doy=[];
    hour=[];
    Bx=[];
    By_GSE=[];
    Bz_GSE=[];
    By_GSM=[];
    Bz_GSM=[];

    for k=lyear:ryear
        str_year=num2str(k);
        fname=['/data/SW_IMF_1h/IMF/imf', str_year, '.txt'];
        if ~exist(fname,'file')
            continue
        end
        fid = fopen(fname);
        fdata = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f');
        fc = fclose(fid);
        [fyear,fdoy,fhour,~,~,~,fBx,fBy_GSE,fBz_GSE,fBy_GSM,fBz_GSM] = fdata{:};
        fBx(fBx==999.9) = nan;
        fBy_GSE(fBy_GSE==999.9) = nan;
        fBz_GSE(fBz_GSE==999.9) = nan;
        fBy_GSM(fBy_GSM==999.9) = nan;
        fBz_GSM(fBz_GSM==999.9) = nan;

        year=[year;fyear];
        doy=[doy;fdoy];
        hour=[hour;fhour];
        Bx=[Bx;fBx];
        By_GSE=[By_GSE; fBy_GSE];
        Bz_GSE=[Bz_GSE; fBz_GSE];
        By_GSM=[By_GSM; fBy_GSM];
        Bz_GSM=[Bz_GSM; fBz_GSM];
    end

    if isempty(year)
        imf=[];
        return
    end

    lp=find( (year==lyear & doy==ldoy), 1, 'first');
    rp=find( (year==ryear & doy==rdoy), 1, 'last');
    year=year(lp:rp);
    doy=doy(lp:rp);
    hour=hour(lp:rp);
    Bx=Bx(lp:rp);
    By_GSE=By_GSE(lp:rp);
    Bz_GSE=Bz_GSE(lp:rp);
    By_GSM=By_GSM(lp:rp);
    Bz_GSM=Bz_GSM(lp:rp);

    day=datenum(year,1,doy)-datenum(lyear,1,ldoy)+hour/24;
    imf=[day,Bx,By_GSE,Bz_GSE,By_GSM,Bz_GSM];
end
