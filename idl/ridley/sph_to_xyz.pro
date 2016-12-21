
; Convert spherical coordinate to Cartesian

function sph_to_xyz, r, theta, phi

Xyz_D = fltarr(3)

Xyz_D(0) = r * sin(theta * !pi/180.) * cos(phi * !pi/180.) 
Xyz_D(1) = r * sin(theta * !pi/180.) * sin(phi * !pi/180.)
Xyz_D(2) = r * cos(theta * !pi/180.)

return, Xyz_D

end
