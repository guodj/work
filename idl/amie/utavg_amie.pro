pro utavg_amie, amie_data, n_amie_times, amie_date, amie_time,nmlts,nlats

ta = intarr(n_amie_times)

for i=0, n_amie_times-1 do begin

  stime = strmid(amie_date(i),4,2)+'-'+strmid(amie_date(i),0,3)+'-'+$
          strmid(amie_date(i),10,2)
  stime = stime + ' ' + strmid(amie_time(i),0,5)
  ;Parse String Time into a string array
  c_s_to_a,itime,stime
  ta(i) = itime(3)

endfor

new_amie_time = strarr(24)
new_amie_date = strarr(24)
new_amie_data = fltarr(24, nmlts,nlats)


for cnt=0, 23 do begin
	ndx = where (ta EQ cnt)
	new_amie_data(cnt,*,*) = TOTAL (amie_data(ndx,*,*),1)/size (ndx, /n_elements)
	; Now we need to get the right time stamp on the data.
	new_amie_date(cnt) = amie_date(ndx(0))
	new_amie_time(cnt) = String( cnt, FORMAT='(I2)') + ':00'
endfor

amie_date = new_amie_date
amie_data = new_amie_data
amie_time = new_amie_time

end