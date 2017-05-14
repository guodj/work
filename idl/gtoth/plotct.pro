;*****************************************************************************

pro plotct, pos, maxmin

;******************************************************************************

    !p.title = ' '
    !y.tickname=strarr(60)
    !y.title = ' '
    !x.title = ' '
    xrange=!x.range & yrange=!y.range & !x.range=0 & !y.range=0

    maxi = max(maxmin)
    mini = min(maxmin)

    array = findgen(10,256)
    for i=0,9 do array(i,*) = findgen(256)/(256-1)*(maxi-mini) + mini

    levels=(findgen(60)-1)/(58-1)*(maxi-mini)+mini

    contour, array, /noerase, /cell_fill, xstyle = 5, ystyle = 5, $
        levels = levels, pos=pos

    plot, maxmin, /noerase, pos = pos, xstyle=1,ystyle=1, /nodata,$
          xtickname = [' ',' '], xticks = 1, xminor=1  , $
          ytickname = strarr(60) + ' ', yticklen = 0.25
    axis, 1, ystyle=1, /nodata, yax=1, charsize=0.9*(!p.charsize > 1.)

    !x.range=xrange & !y.range=yrange

  return

end

