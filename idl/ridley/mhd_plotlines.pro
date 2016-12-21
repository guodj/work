
nVars = 14
nPts  = 100

Files = findfile("*_p0.dat")
nFiles = n_elements(Files)

data = fltarr(nVars,nFiles,nPts)
Vars = strarr(nVars)

line = ''
tmp  = fltarr(nVars)

for iFile = 0, nFiles-1 do begin

   openr,1,Files(iFile)
   for i=0,nVars -1 do begin
      readf,1,line
      s = strpos(line,'"')
      Vars(i) = strmid(line,s+1,strlen(line)-s-2)
   endfor
   readf,1,line
   readf,1,line

   for i=0,nPts-1 do begin
      readf,1,tmp
      data(*,iFile,i) = tmp
      data(1,iFile,i) = float(iFile)*30.0/60.0
   endfor

   close,1

endfor

setdevice, 'loc0.ps','p',5

ppp = 4
space = 0.01
pos_space, ppp, space, sizes, ny = ppp

t = reform(data(1,*,0))
bz = reform(data(9,*,0))

rho_si = data(3,*,0)*1.0e6 * cmp_
v_si = abs(data(4,*,0)*1000.0)
b_si = abs(bz)*1.0e-9
beq_si = 31100.0*1.0e-9
cmu0_ = 4*!pi*1.0e-7

dp = rho_si*v_si^2

rnorm = (sqrt(2*cmu0_*dp)/(2*beq_si))^(-1.0/3.0)
rnormwb = (sqrt(2*cmu0_*dp + b_si^2)/(2*beq_si))^(-1.0/3.0)

shue98 = (dp*1.0e9)^(-1.0/6.6) * (10.22 + 1.29*tanh(0.84*(bz+8.14)))
shue97 = (dp*1.0e9)^(-1.0/6.6) * (11.40 + 0.14*bz)

get_position, ppp, space, sizes, 0, pos, /rect
pos(0) = pos(0) + 0.1

plot, t, bz, pos = pos, /noerase, xstyle = 1, yrange = [-20.0,20.0], $
      xtickname = strarr(10)+' ', ytitle = 'Bz at 32 Re (nT)'

get_position, ppp, space, sizes, 1, pos1, /rect
get_position, ppp, space, sizes, 2, pos2, /rect
pos = [pos1(0) + 0.1, pos2(1), pos1(2), pos1(3)]

makect,'mid'
contour, data(12,*,*), data(1,*,*), data(0,*,*), yrange = [4.0,12.0], $
         /fill, nlevels = 31, ytitle = 'X (Re)', $
         pos = pos, xstyle = 1, /noerase, xtickname = strarr(10)+' '

mp = fltarr(nFiles)

for i=0,nFiles-1 do begin

   jy = reform(data(12,i,*))
   l  = where(jy eq max(jy))
   mp(i) = data(0,i,l(0))

endfor

get_position, ppp, space, sizes, 3, pos, /rect
pos(0) = pos(0) + 0.1

oplot, t, mp

plot, t, mp, yrange = [4.0, 12.0], ytitle = 'X (Re)', xtitle = 'Time (hours)', $
         pos = pos, /noerase, xstyle = 1

oplot, t, rnorm, linestyle = 1
oplot, t, rnormwb, linestyle = 2
oplot, t, shue98, linestyle = 3
oplot, t, shue97, linestyle = 4

closedevice

end
