
pro download_dst, itime, dstfile

  syyyy = tostr(itime(0))
  smm = chopr('0'+tostr(itime(1)),2)

  found = 0

  dstfile = findfile('dstfinal'+syyyy+smm+'.html')

  if (strlen(dstfile) eq 0) then begin

     spawn, 'wget http://wdc.kugi.kyoto-u.ac.jp/dst_final/'+syyyy+smm+'/index.html'
     spawn, 'mv index.html dstfinal'+syyyy+smm+'.html'
     dstfile = findfile('dstfinal'+syyyy+smm+'.html')
     if (strlen(dstfile) eq 0) then begin
        ; probably provisional dst file
        spawn,'wget http://wdc.kugi.kyoto-u.ac.jp/dst_provisional/'+$
              syyyy+smm+'/index.html'
        test = findfile('index.html')
        if (strlen(test) eq 0) then begin
           spawn,'wget http://wdc.kugi.kyoto-u.ac.jp/dst_realtime/'+$
                 syyyy+smm+'/index.html'
           test = findfile('index.html')
           if (strlen(test) gt 0) then found = 1
        endif else found = 1
        if (found) then begin
           dstfile = 'dstfinal'+syyyy+smm+'.html'
           spawn, 'mv index.html '+dstfile
        endif else begin
           print, 'cant download dst file. Stopping'
           stop
        endelse
     endif
  endif

end
