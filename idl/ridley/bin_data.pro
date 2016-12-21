
pro bin_data, datax, datay, nPts, datas

    xRange = max(datax) - min(datax)
    dx = xRange/nPts
    newx = findgen(nPts)
    newy = fltarr(nPts)
    news = fltarr(nPts)

    n = 0
    for i=0,nPts-2 do begin
        lowx = i*dx + min(datax)
        highx = lowx + dx
        loc = where(datax ge lowx and datax lt highx, count)
        if (count gt 0) then begin
            newx(n) = median(datax(loc))
            newy(n) = median(datay(loc))
            if (count gt 1) then news(n) = stddev(datay(loc)) else news(n)=0.0
            n = n + 1
        endif
    endfor

    datax = newx(0:n-1)
    datay = newy(0:n-1)
    datas = news(0:n-1)
    
end
