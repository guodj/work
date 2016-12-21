
filein = ask('data file from T-GCM run','NO_eq.bin')
bora = ''
read,'Is this file a binary or ascii file (b/a) [a] : ',bora

psfile = ask('PS file name (return for screen)','')
if strlen(psfile) eq 0 then window else setdevice, psfile, 'p', 4, 0.9

if (strmid(bora,0,1) eq 'b') or (strmid(bora,0,1) eq 'B') then 		$
  rd_tgcm_bin, filein, in_data, vars, utdes, axdes, xax, yax		$
else									$
  rd_tgcm_ascii, filein, in_data, vars, utdes, axdes, xax, yax, 200, run

print, '1. Plot normal data'
print, '2. Plot log(data)'
print, '3. Plot 10^data'
plotdata = ask('Plotting option (1,2,3)','1')

if (plotdata eq '2') then in_data = alog10(in_data)
if (plotdata eq '3') then in_data = 10^(in_data)

nva = n_elements(in_data(*,0,0))
nla = n_elements(in_data(0,*,0))
nal = n_elements(in_data(0,0,*))

xpos  = fltarr(nla,nal)
ypos  = fltarr(nla,nal)
datai = fltarr(nla,nal)

for ial = 0, nal-1 do xpos(*,ial) = xax
for ila = 0, nla-1 do ypos(ila,*) = yax

h_string = 'Geographic Coordinates'

; read in color table

ct_dir = getenv('IDL_EXTRAS')
;ctname = ct_dir+'joule4.ct'
ctname = ct_dir+'bw.ct'
readct,ncolors,ctname

dx = 1.25*float(!d.y_ch_size)/float(!d.y_size)
ppp = 1
space = dx*2.0
pos_space,ppp,space,sizes
get_position, ppp, space, sizes, 0, pos
pos([0,2]) = pos([0,2]) - pos(0)/2.0

spawn,'date',date

for j=0,nva-1 do begin

  mini = min(in_data(j,*,*)) 
  maxi = max(in_data(j,*,*))

  plotdumb

  xyouts, 0.0, 0.0, strcompress(utdes(j)), /norm
  xyouts, 0.0, -dx*1.5, h_string, /norm
  title = strcompress(vars(j))

  datai(*,*) = in_data(j,*,*)

  image = ncolors*(datai-mini)/(maxi-mini)

  levels = float(ncolors)*findgen(30)/29

  contour, image, xpos, ypos, nlevels=30, 		$
           /noerase, pos = pos, xstyle = 1, 		$
           ystyle = 1, /fill,levels=levels, 		$
           c_colors=levels, 				$
           xtitle =strcompress(axdes(0,j)),		$
           ytitle =strcompress(axdes(1,j))

  contour, datai, xpos, ypos, /follow, nlevels=15, 		$
           /noerase, pos = pos, xstyle = 5, 		$
           ystyle = 5

  xyouts, mean([pos(0),pos(2)]), pos(3)+0.01, title, /norm, alignment =0.5

  xyouts, -0.02, -0.02,filein+'  '+date,				$
	orientation=90,charsize=0.5,/norm

  posct = [pos(2)+0.01,pos(1)+0.1,pos(2)+0.04,pos(3)-0.1]
  plotct, ncolors, posct, [mini,maxi], strcompress(vars(j)), /right,color=0

endfor

closedevice

end