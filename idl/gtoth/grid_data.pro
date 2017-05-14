;==============================================================================
function grid_data,x,y,data,nxreg,xreglimits,triangles,wregpad

return,griddata(x,y,data,$
        dimension=nxreg,$
	start=xreglimits(0:1),$
	delta=[(xreglimits(2)-xreglimits(0))/(nxreg(0)-1),  $
               (xreglimits(3)-xreglimits(1))/(nxreg(1)-1)] ,$
        triangles=triangles,$
;        method='NearestNeighbor',$
;        method='Linear',$
        method='InverseDistance',$
       smoothing=0.5,$
       max_per_sector=4, $
        missing=wregpad $
        )
end

