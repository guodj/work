;===========================================================================
function curl_rz,u,v,r,z
;
; Take curl of "u,v" vector with respect to "r,z" on a structured 2D grid
; with axial symmetry in the ignored direction.
; Curl is the contour integral of edge_vector.(u,v)*R_edge_averaged
; divided by R*cell_volume for each cells - v/R.
; See also comments for the div_rz function on edge average and edge cells.
;
;===========================================================================

if n_elements(u) eq 0 or n_elements(v) eq 0 $
   or n_elements(r) eq 0 or n_elements(z) eq 0 then begin
   print,'Missing arguments in function curl_rz'
   retall
endif

gengrid,'curl',r,z,rc,zc,vol2,u,v

return,filledge(-((intedge_rz(v,rc,zc)+intedge(u,rc^2))/vol2 - v)/2/r)

end

