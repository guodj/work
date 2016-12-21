;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
; Procedure remove_spike
;
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

PRO remove_spike, m_array, max_val, removed, nrem 

  removed = where(abs(m_array) gt max_val, nrem)

  if nrem gt 0 then m_array(removed) = 0.0

END
