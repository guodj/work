function yeardoy=numdate(matday)
    % this function transforms matlat day to year and doy.
    matday=reshape(matday,[length(matday),1]);

    yeardoy=ones(length(matday),2)*nan;
    yeardoy(:,1)=2000;
    yeardoy(:,2)=matday-datenum(2000,1,1)+1;

    p=yeardoy(:,2)>yeardays(2000) | yeardoy(:,2)<1;

    while sum(p)>=1
        yeardoy(p,1)=yeardoy(p,1)+floor(yeardoy(p,2)/365);
        yeardoy(p,2)=matday(p)-datenum(yeardoy(p,1),1,1)+1;

        p1 = (yeardoy(:,2)==0);
        if sum(p1)>=1
            yeardoy(p1,1)=yeardoy(p1,1)-1;
            yeardoy(p1,2)=matday(p1)-datenum(yeardoy(p1,1),1,1)+1;
        end


        % pp=yeardoy(p,2)>yeardays(yeardoy(p,1)) | yeardoy(p,2)<1;
        % p=p(pp);
        p=yeardoy(:,2)>yeardays(yeardoy(:,1)) | yeardoy(:,2)<1;
    end
end
