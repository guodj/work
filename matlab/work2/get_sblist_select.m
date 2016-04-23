function sblist_select=get_sblist_select(sblist,sbtype,season,yncir)
    switch sbtype
        case 'at'
            index1=is_away_toward(sblist);
        case 'ta'
            index1=is_toward_away(sblist);
        otherwise
            error('wrong input of sbtype!');
    end
    switch season
        case 'me'
            index2=isME(sblist);
        case 'se'
            index2=isSE(sblist);
        case 'js'
            index2=isJS(sblist);
        case 'ds'
            index2=isDS(sblist);
        otherwise
            error('wrong input of season!');
    end
    switch yncir
        case 'all'
            index3=~isnan(sblist(:,1));
        case 'ncir'
            [~,index3]=yncir(sblist);
        case 'ycir'
            [index3,~]=yncir(sblist);
        otherwise
            error('wrong input of cirtype!');
    end
    sblist_select=sblist(index1 & index2 & index3,:);
end
