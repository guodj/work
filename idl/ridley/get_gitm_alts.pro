
pro get_gitm_alts, gitmalts

  filelist = findfile('3DALL*.bin')
  file = filelist(0)

  AddedVars = [2]
  gitm_read_bin_1var, file, gitmcoords, gitmtime1, nVars, Vars, ver,$
                      VarsToGet = AddedVars
  gitmAlts = reform(gitmcoords(0,0,0,*))/1000.0

end
