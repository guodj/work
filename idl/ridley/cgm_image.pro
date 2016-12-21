
print, 'Enter file list which you would like plotted : '
flist = ''
read,flist

min_lat = 50.0

print, 'Enter number of times to skip in between each image : '
print, '(0 for no skipping, 1 for every other image, ...)'
nskip = 0
read, nskip
nskip = nskip + 1

print, 'Enter color table which you would like :'
ct = ''
read,ct

print, 'Enter maximum reyleighs for color table (200-600) :'
maxct = 200
read,maxct

print, 'Enter number of plots per page : '
ppp = 4
read, ppp

print, 'Enter ps file name (return for screen) :'
psfile = ''
read,psfile

if strlen(psfile) gt 0 then begin
  print, 'Would you like a different ps file for each page (y/n) ?'
  que=''
  read,que
  if (strmid(que,0,1) eq 'y') then begin
    moreppp = 1
    strl = strpos(psfile,'.ps')
    if strl gt 0 then psfile = strmid(psfile,0,strl)
    outfile = psfile+'_01'+'.ps'
  endif else begin
    moreppp = 0
    outfile = psfile
  endelse
  setdevice,outfile,'land',4,0.95
endif else moreppp = 0

readct, ncolors, ct

read_images, flist, nskip, uttimes, inttimes, image_save, min_lat

nimage = n_elements(image_save(*,0,0))

space = 0.02

pos_space, ppp, space, sizes
get_position, ppp, space, sizes, 0, pos  
xoff = pos(0)

maxran = 90.0-min_lat

for i= 0, nimage-1 do begin

  pn = i mod ppp
  pagen = i/ppp + 1

  if (pagen gt 1) and (pn eq 0) then begin
    if (moreppp) then begin
      device, /close
      outfile = psfile+'_'+chopr('0'+tostr(pagen),2)+'.ps'
      setdevice,outfile,'land',4,0.95
    endif
    if (!d.name eq 'X') then prompt_for_next 
  endif

  if (pn eq 0) then begin
    plotdumb
    itime = intarr(6)
    itime(0) = strmid(uttimes(i),0,2)
    itime(1) = strmid(uttimes(i),2,2)
    itime(2) = strmid(uttimes(i),4,2)
    c_a_to_s, itime, stime
    date = strmid(stime,0,9)
    xyouts, 0.5, 1.0, date, /norm, alignment = 0.5
    ctpos = [0.97,0.3,1.0,0.7]
    plotct, ncolors, ctpos, maxct, 'photons/(cm2 s)',/right
  endif

  get_position, ppp, space, sizes, pn, pos  

  pos([0,2]) = pos([0,2]) - xoff
  xs = pos(2) - pos(0)
  ys = pos(3) - pos(1)

  plot, [-maxran,maxran],[-maxran,maxran],		$
	pos = pos, xstyle = 5, ystyle = 5,		$
	/nodata, /noerase

  image = float(ncolors)*image_save(i,*,*)/float(maxct)

  tv, image, pos(0), pos(1), xsize = xs, ysize = ys, /norm

  hms = strmid(uttimes(i),7,4)+':'+strmid(uttimes(i),11,2)+' UT'
  xyouts, -maxran, 0.9*maxran, hms, color = 255, charsize=0.8

  tp = sizes.nby*sizes.nbx

  if ((nimage-1)-i lt sizes.nbx) or 				$
     (pn+1 gt tp-sizes.nbx) then no00 = 0 else no00 = 1
  if (pn mod sizes.nbx eq sizes.nbx-1) or 			$
     (i eq nimage-1) then no06 = 0 else no06 = 1
  if (pn lt sizes.nbx) then no12 = 0 else no12 = 1
  if (pn mod sizes.nbx eq 0) or 				$
     (i eq 0) then no18 = 0 else no18 = 1

  plotmlt, 90.0-min_lat, /white, 			$
	no00 = no00, no06 = no06, no12 = no12, no18 = no18

endfor

end
