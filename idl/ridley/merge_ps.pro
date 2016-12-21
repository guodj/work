function params, nplots

  ny = fix(sqrt(nplots))
  nx = nplots/ny

  if nx*ny lt nplots then nx = nx + 1

  if nx gt ny then rotation = 90 else rotation = 0

  return, [nx,ny,rotation]

end


pro merge_ps, filelist=filelist

  if n_elements(filelist) eq 0 then begin

    print, 'Enter number of files : '
    nps = 0
    read, nps

    filelist = strarr(nps)

    for i=0,nps-1 do begin

      dum$=''
      print, 'Enter file number ',strcompress(string(i+1),/remove),' :'
      read, dum$
      filelist(i) = dum$

    endfor

  endif

  nps = n_elements(filelist)

  openw,1, 'TEMP_FILE.ps'
  line = ''
  nplots = 0

  for i=0,nps-1 do begin

    openr, 2, filelist(i)

    while not eof(2) do begin

      readf,2,line
      if strpos(line,'showpage') ne -1 then nplots = nplots + 1
      printf,1,line

    endwhile

    close, 2

  endfor

  close, 1

  p = params(nplots)
  dx = fix(25400/p(0))
  dy = fix(17780/p(1))
  if p(2) eq 90 then scale = 0.75/p(1) else scale = 1.0/p(1)

  rot = strcompress(string(p(2)),/remove_all)+' rotate'

  print, 'Enter new ps file name : '
  newps = ''
  read, newps

  openw,4,newps
  openr,3,'TEMP_FILE.ps'

  flag = 0
  plotn = 0
  l = ''
  while not eof(3) do begin
    readf,3,l
    if flag eq 1 and strpos(l,'showpage') ne -1 then flag=2
    if flag ne 2 then printf,4,l
    if flag eq 0 and strpos(l,'%%EndPageSetup') ne -1 then begin
      printf,4,'gsave                     %% save graphics state'
      x = dx*(float(plotn mod p(0)) + 0.75)
      y = dy*(plotn/p(0))
      t = strcompress(string(x)+string(y)+' translate')
      print, plotn mod p(0), plotn/p(0)
      printf,4,t
      printf,4,rot
      printf,4,strcompress(string(scale),/remove_all)+' dup scale'
      plotn = plotn+1
      flag=1
    endif
    if flag eq 2 and strpos(l,'%%EndPageSetup') ne -1 then begin
      printf,4,'grestore                  %% return graphics'
      printf,4,'gsave                     %% save graphics state'
      x = dx*(float(plotn mod p(0)) + 0.75)
      y = dy*(plotn/p(0))
      print, plotn mod p(0), plotn/p(0)
      t = strcompress(string(x)+string(y)+' translate')
      printf,4,t
      printf,4,rot
      printf,4,strcompress(string(scale),/remove_all)+' dup scale'
      plotn = plotn + 1
      if plotn eq nplots then flag=3 else flag = 1
    endif
  endwhile
  close,3,4

  spawn, '/bin/rm TEMP_FILE.ps'

  return

end
