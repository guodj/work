
pro read_single_file, infile, VarToRead

  openr,1,infile

  i = 0
  j = 0
  nPts = 0
  mlt_set = 0.0
  mlat_set = 0.0
  B0 = 0.0
  Bsme = 0.0
  Bt1 = 0.0
  Bt2 = 0.0
  Rsq_e = 0.0

  while not eof(1) do begin

     readf,1,i,j,npnts, mlt_set, mlat_set, B0
     readf,1,Bsme,Bt1,Bt2,Rsq_e
  
     if (npnts lt 160) then begin
        B0 = 0.0
        Bsme = 0.0
        Bt1 = 0.0
        Bt2 = 0.0
     endif

     VarToRead(0,i,j) = B0 
     VarToRead(1,i,j) = Bsme
     VarToRead(2,i,j) = Bt1
     VarToRead(3,i,j) = Bt2
     VarToRead(4,i,j) = npnts

  endwhile

  close,1

end


pro load_ovation_sme, OvationSmeInternals

  nMlts = 96
  nMlats = 80

  beta_je_diff = fltarr(5, nMlts+1, nMLats)
  beta_je_mono = fltarr(5, nMlts+1, nMLats)
  beta_je_wave = fltarr(5, nMlts+1, nMLats)
  beta_je_iono = fltarr(5, nMlts+1, nMLats)
  beta_jn_diff = fltarr(5, nMlts+1, nMLats)
  beta_jn_mono = fltarr(5, nMlts+1, nMLats)
  beta_jn_wave = fltarr(5, nMlts+1, nMLats)
  beta_jn_iono = fltarr(5, nMlts+1, nMLats)

  dir = '/Users/ridley/Gitm/GITM2/srcData/Aurora/'

  read_single_file, dir+'Diff_je_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt', beta_je_diff
  beta_je_diff(*,nMlts,*) = beta_je_diff(*,0,*)
  read_single_file, dir+'Mono_je_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt', beta_je_mono
  beta_je_mono(*,nMlts,*) = beta_je_mono(*,0,*)
  read_single_file, dir+'Wave_je_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt', beta_je_wave
  beta_je_wave(*,nMlts,*) = beta_je_wave(*,0,*)
  read_single_file, dir+'Iono_je_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt', beta_je_iono
  beta_je_iono(*,nMlts,*) = beta_je_iono(*,0,*)

  read_single_file, dir+'Diff_jn_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt', beta_jn_diff
  beta_jn_diff(*,nMlts,*) = beta_jn_diff(*,0,*)
  read_single_file, dir+'Mono_jn_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt', beta_jn_mono
  beta_jn_mono(*,nMlts,*) = beta_jn_mono(*,0,*)
  read_single_file, dir+'Wave_jn_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt', beta_jn_wave
  beta_jn_wave(*,nMlts,*) = beta_jn_wave(*,0,*)
  read_single_file, dir+'Iono_jn_B0_Bsme_Bt1_Bt2_Rsq_by_mlt_mlat_ch_180_sqsme.txt', beta_jn_iono
  beta_jn_iono(*,nMlts,*) = beta_jn_iono(*,0,*)

  OvationSmeInternals = {nMlts: nMlts, $
                         nMLats: nMLats, $
                         beta_je_diff: beta_je_diff, $
                         beta_je_mono: beta_je_mono, $
                         beta_je_wave: beta_je_wave, $
                         beta_je_iono: beta_je_iono, $
                         beta_jn_diff: beta_jn_diff, $
                         beta_jn_mono: beta_jn_mono, $
                         beta_jn_wave: beta_jn_wave, $
                         beta_jn_iono: beta_jn_iono, $
                         minlat: 50.0}

end
