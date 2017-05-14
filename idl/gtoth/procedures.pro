;^CFG COPYRIGHT VAC_UM
; Written by G. Toth for the Versatile Advection Code and BATSRUS
; Some improvements by Aaron Ridley.
;
; Procedures for
;
; reading ascii and binary data produced by VAC and VACINI:
;    openfile,gettype,gethead,getpict,getpict_asc,getpict_bin
; reading numbers and strings from input:
;    asknum, askstr, str2arr, readplotpar, readlimits
; transforming initial data:
;    regulargrid, polargrid, unpolargrid, spheregrid, getaxes
; calculating functions of the data
;    getfunc, getlimits
; plotting
;    plotfunc, plotgrid
; calculating cell corners and cell volumes for general 2D grids
;    gengrid
; comparing two w or x arrays for relative differences
;    compare
; checking index ranges for functions quadruplet and triplet
;    checkdim
; procedure "quit" as an alias for "exit"
;    quit
;
; Functions for
;
; calculating derivatives in 2D for Cartesian grids to 2nd,3rd,4th order
;    diff2,diff3,diff4
; calculate minmod limited slope
;    minmod
; calculating symmetric differences with respect to some mirror plane
;    symmdiff
; calculating derivatives in 2D for general grids
;    grad,div,curl,grad_rz,div_rz,curl_rz, filledge,intedge,intedge_rz
; taking a part of an array or coarsen an array
;    triplet, quadruplet, coarse
; eliminating degenerate dimensions from an array
;    reform2


