
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; Procedure remove_mean
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PRO remove_mean, m_array, mean

  mean = 0.0
  nele = n_elements(m_array)

  for n=0,nele-1 do mean = mean + m_array(n)

  mean = mean/nele

  m_array = m_array - mean

END

