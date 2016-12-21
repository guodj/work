;===========================================================================
function div,u,v,x,y
;
; Take divergence of "u,v" vector with respect to "x,y" on a structured 2D grid
; Divergence is the contour integral of edge_normal.(u,v)_edge_averaged
; divided by cell_volume for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center cancels.
; Divergence can be calculated for inner points only, edge values are
; copies of inner neighbors.
;===========================================================================

if n_elements(u) eq 0 or n_elements(v) eq 0 $
   or n_elements(x) eq 0 or n_elements(y) eq 0 then begin
   print,'Missing arguments in function div'
   retall
endif

gengrid,'div',x,y,xc,yc,vol2,u,v

return,filledge((intedge(u,yc)-intedge(v,xc))/vol2)

end

