function [ycir,ncir] = yncir( sblist )
    % find sblist accompanied with and without cir.
    fname='/data/CIRlist/streaminterfacelist.txt';
    fid=fopen(fname);
    fdata=textscan(fid,'%s %s %f %f %f %f %f %f %f %f','commentstyle','%');
    fc=fclose(fid);

    [~,~,mattime_cirlist,~,~,~,~,~,~,~]=fdata{:};  
    mattime_cirlist=fix(mattime_cirlist);
    mattime_cirlist=sort(mattime_cirlist);
    mattime_cirlist_extend=[mattime_cirlist; mattime_cirlist+1; mattime_cirlist-1];

    mattime_sblist=datenum(sblist(:,1),1,sblist(:,2));

    ycir=ismember(mattime_sblist,mattime_cirlist_extend);

    ncir=~ycir & (mattime_sblist>=mattime_cirlist(1)) & (mattime_sblist<=mattime_cirlist(end));
end
