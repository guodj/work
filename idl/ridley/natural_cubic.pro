
PRO natural_cubic, a, cubic_a

  if n_elements(a) lt 3 then a=randomu(s,10)
  x=findgen(n_elements(a))

  n     = n_elements(x) -1
  h     = fltarr(n+1)
  alpha = fltarr(n)
  l     = fltarr(n+1)
  mu    = fltarr(n+1)
  z     = fltarr(n+1)
  b     = fltarr(n+1)
  c     = fltarr(n+1)
  d     = fltarr(n+1)

  for i=0,n-1 do h(i) = x(i+1)-x(i)
  for i=1,n-1 do alpha(i) = 3.0*(a(i+1)*h(i-1) -			      $
		 a(i)*(x(i+1)-x(i-1)) + a(i-1)*h(i))/(h(i-1)*h(i))

  l(0) = 1.0
  mu(0) = 0.0
  z(0) = 0.0

  for i=1, n-1 do begin

    l(i) = 2.0*(x(i+1)-x(i-1))-h(i-1)*mu(i-1)
    mu(i) = h(i)/l(i)
    z(i) = (alpha(i)-h(i-1)*z(i-1))/l(i)

  endfor

  l(n) = 1.0
  z(n) = 0.0
  c(n) = 0.0

  for j=n-1,0,-1 do begin

    c(j) = z(j) - mu(j)*c(j+1)
    b(j) = (a(j+1) - a(j))/h(j) - h(j)*(c(j+1) + 2.0*c(j))/3.0
    d(j) = (c(j+1) - c(j))/(3.0*h(j))

  endfor

  xx = findgen(10*(n+1))/(10.0*(n+1)-1.0)*(max(x)-min(x)) + min(x)
  yy = fltarr(10*(n+1))
  xjpo = x(findgen(n)+1)
  for i=0,10*(n+1)-1 do begin
    j = where(x le xx(i) and xx(i) le xjpo)
    j = j(0)
    xmxj = xx(i) - x(j)
    yy(i) = a(j) + b(j)*xmxj + c(j)*xmxj^2.0 + d(j)*xmxj^3.0
  endfor

  cubic_a = yy

END
