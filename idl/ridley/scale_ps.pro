pro scale_ps,file,scale=scale

;	SCALE [,FILE] [,SCALE=SCALE]

;	This scales an existing postscript file


;    *** Choose a file ***

	if n_elements(file) ne 1 then begin
	  file=pickfile(filter='*.ps')
	  file=file(0)
	endif


;    *** Choose a scaling ***

	if n_elements(scale) ne 1 then begin
	  scale = 0.0
	  read,'What scale? (1.0 is full size.) ',scale
	endif

	fileout = strmid(file,0,strlen(file)-5)
	fileout = fileout+tostr(100.0*scale)+'.ps'


;    *** Open files ***

	openr,1,file
	openw,2,fileout


;    *** Read through file,writing to fileout ***

	l=''
	scaled = 0
	while not eof(1) do begin

	  readf,1,l
	  printf,2,l


;    *** Look for 'dup scale' and add additional scaling after ***

	  if (strpos(l,'dup scale',0) ne -1) and not scaled then begin
	    scaled = 1
	    printf,2,strcompress(scale,/rem)+' dup scale'
	  endif

	endwhile


;    *** close file and fileout ***

	close,1,2

	print,'New scaled PostScript file created:'
	print,'	'+fileout

	return
end

