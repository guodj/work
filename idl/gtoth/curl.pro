;===========================================================================
function curl,u,v,x,y
;
; Take curl of "u,v" vector with respect to "x,y" on a structured 2D grid.
; Curl is the contour integral of edge_vector.(u,v)_edge_averaged
; divided by cell_volume for each cells. See also comments for div function.
;
;===========================================================================

if n_elements(u) eq 0 or n_elements(v) eq 0 $
   or n_elements(x) eq 0 or n_elements(y) eq 0 then begin
   print,'Missing arguments in function curl'
   retall
endif

gengrid,'curl',x,y,xc,yc,vol2,u,v

return,filledge((intedge(u,xc)+intedge(v,yc))/vol2)

end

