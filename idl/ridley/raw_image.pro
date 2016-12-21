pro read_raw, filelist, images, times

  flist = findfile(filelist)

  n = n_elements(flist)

  if n gt 0 then begin

    images = fltarr(n,200,228)
    times = strarr(n)

    for i=0,n-1 do begin

      openr,1,flist(i)

      x=0
      y=0
      readf,1,x,y
      dtime = ''
      readf,1,dtime

      im = fltarr(x,y)
      readu,1,im
      readu,1,im
      readu,1,im

      images(i,0:x-1,0:y-1) = im(0:x-1,0:y-1)
      times(i) = dtime

      close,1

    endfor

  endif

  return

end

print, 'Enter file list :'
filelist = ''
read,filelist

print, 'Enter ps file name (return for screen) :'
psfile = ''
read,psfile

print, 'Enter color table name :'
ctname = ''
read,ctname

print, 'Enter maximum photons/(cm2s) (500-1000) :'
maxct = 0.0
read,maxct

print, 'Enter number of plots per page:'
ppp = 0
read,ppp

readct, ncolors,ctname

if strlen(psfile) gt 0 then setdevice, psfile, 'l', 4, 0.95

space = 0.04
pos_space,ppp,space,sizes
fac = 200.0/228.0

read_raw, filelist, images, times

images = float(ncolors)*images/maxct
n = n_elements(times)

for i=0,n-1 do begin

  pn = i mod ppp
  if pn eq 0 then begin
    plotdumb
    pos = [0.97,0.3,1.0,0.7]
    title = 'Photons/(cm!E2!N/s)'
    plotct, ncolors, pos, maxct, title, /right
  endif

  get_position,ppp,space,sizes,pn,pos
  xpos = pos(0)
  ypos = pos(1)
  xs = (pos(2)-pos(0))*fac
  ys = (pos(3)-pos(1))

  tv,images(i,*,*),xpos,ypos,xsize = xs,ysize=ys,/norm

  xyouts,xpos,ypos+ys+0.005,times(i), /norm

endfor

if !d.name eq 'PS' then begin
  device, /close
  set_plot,'X'
endif

end 

