pro get_temp_ps, psfile, broadps

  if !version.os eq 'vms' then begin

    spawn, 'show process',list

    path = 'zasu$dkb500:[software.plotting.idl.temp]'
    name = strcompress(strmid(list(1),32,10),/remove_all)
    date = strcompress(strmid(list(1),0,2)+				$
		       strmid(list(1),3,3)+				$
		       strmid(list(1),9,2)+				$
		       strmid(list(1),12,2)+				$
		       strmid(list(1),15,2)+				$
		       strmid(list(1),18,2), /remove_all)
     ext = '.ps'

     psfile = path+name+date+ext
     broadps = path+name+'*'+ext+';*'

   endif else begin

     psfile = 'temp.ps'
     broadps = 'temp.ps'

   endelse

   return

end    
