

if (n_elements(amiefile) eq 0) then begin
    amiefile = findfile('b*_eq')
    amiefile = amiefile(0)
endif

amiefile = ask('file to plot',amiefile)

read_amie_binary, amiefile, data, lats, mlts, time, fields, 	$
  imf, ae, dst, hp, cpcp, version

nTimes = n_elements(time)
nLats = n_elements(lats)
nLons = n_elements(mlts)

efielde = reform(data(*,0,1:nLons-1,*)) - reform(data(*,0,0:nLons-2,*))

for iTime=0,nTimes-1 do for iLon=0,nLons-2 do $
  efielde(iTime,iLon,*) = efielde(iTime,iLon,*)/(15*111.0*cos(lats*!dtor)+0.1)

plotdumb

ppp = 8
space = 0.01
pos_space, ppp, space, sizes, ny = ppp

for i=0,ppp-1 do begin

    get_position, ppp, space, sizes, i, pos, /rect
    pos(0) = pos(0) + 0.1

    plot, efielde(*,0,(i+1)*5-1)*1000.0, pos = pos, /noerase

endfor

end
