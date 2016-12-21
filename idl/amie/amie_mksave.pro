
;--------------------------------------------------------------
; Get Inputs from the user
;--------------------------------------------------------------

initial_guess = 'b?????????'
amie_file = ask('AMIE binary file name',initial_guess)

filelist = findfile(amie_file)

nfiles = n_elements(filelist)

for i=0,nfiles-1 do begin

  amie_file = filelist(i)
  savefile = amie_file+'.save'

  print, "Reading File : ", amie_file

  read_amie_binary, amie_file, data, lats, mlts, time, fields, 		$
                  imf, ae, dst, hp, cpcp, version

  print, "Saving File : ", savefile

  save, amie_file, data, lats, mlts, time, fields, 		$
                  imf, ae, dst, hp, cpcp, version, filename = savefile

endfor

end


