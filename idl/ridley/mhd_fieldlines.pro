file = "trace_1.dat"
file = ask("file name to process",file)

spawn, "wc "+file, result

print, result

n = long(result) - 8
n = n(0)
xyz = dblarr(3,n)

dummy = ""
openr,1,file
for i=1,8 do readf,1,dummy
readf,1,xyz
close,1

nPtsInCircle = 7

e = n-nPtsInCircle
p = ulonarr(5,e)
p(0,*) = 4
for i=0L,e-1 do begin
  p(1,i) = i
  p(2,i) = i+nPtsInCircle
  if (i mod nPtsInCircle eq nPtsInCircle-1) then begin
      p(3,i) = i+1
      p(4,i) = i-(nPtsInCircle-1)
  endif else begin
      p(3,i) = i+(nPtsInCircle+1)
      p(4,i) = i+1
  endelse
endfor

t = e*4+e
edges = ulonarr(t)

edges = p

opoly = obj_new('IDLgrPolygon',xyz,polygon=edges)
omodel = obj_new('IDLgrModel')
oview = obj_new('IDLgrView')
omodel->Add,opoly
oview->Add,omodel

ovrml = obj_new('IDLgrVRML',filename=file+'.wrl',quality=2)
ovrml->Draw,oview

obj_destroy, opoly
obj_destroy, omodel
obj_destroy, oview
obj_destroy, ovrml

end
