function dates=get_dates_satellite(dates_list,satellite)
    % return dates within the range of satellite data
    % input:
    % dates_list=[year,doy,...];
    % satellite='champ', 'grace' or 'goce'
    % output:
    % dates=selected dates_list that are within the date range of the satellite
    lchamp=[2001, 135];
    rchamp=[2010, 182];
    lgrace=[2002, 213];
    rgrace=[2010, 177];
    lgoce=[2009, 305];
    rgoce=[2012, 152];

    matlchamp=datenum(2001,1, 135);
    matrchamp=datenum(2010,1, 182);
    matlgrace=datenum(2002,1, 213);
    matrgrace=datenum(2010,1, 177);
    matlgoce=datenum(2009,1, 305);
    matrgoce=datenum(2012,1, 152);

    matdates=datenum(dates_list(:,1),1,dates_list(:,2));
    switch satellite
        case 'champ'
            index = (matdates>matlchamp & matdates<matrchamp);
        case 'grace'
            index = (matdates>matlgrace & matdates<matrgrace);
        case 'goce'
            index = (matdates>matlgoce & matdates<matrgoce);
        otherwise
            error('wrong input of satellite!');
    end

    dates=dates_list(index,:);
end
