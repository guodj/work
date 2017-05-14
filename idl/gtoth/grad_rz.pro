;===========================================================================
function grad_rz,idir,f,r,z
;
; Take gradient of "f" in direction "idir" on the "r,z" structured 2D grid
; assuming axial symmetry in the ignored direction.
; Gradient is the contour integral of edge_normal_idir*f*R_edge_averaged
; divided by R*cell_volume - f/r for each cells. The cell corners are at the
; averaged coordinates of the four neighboring cell centers.
; However there is no need for edge averaging since the contribution of
; the value in the cell center cancels for idir=2, or equals +f/2R for idir=1.
; Gradient can be calculated for inner points only, edge values are
; copies of inner neighbors.
;===========================================================================

if n_elements(ndir) eq 0 or n_elements(f) eq 0 $
   or n_elements(r) eq 0 or n_elements(z) eq 0 then begin
   print,'Missing arguments in function grad_rz'
   retall
endif

gengrid,'grad_rz',r,z,rc,zc,vol2,f

if idir eq 1 then return,filledge( (intedge_rz(f,rc,zc)/vol2 - f)/2/r ) $
else              return,filledge( -intedge(f,rc^2)/vol2/2/r)

end

