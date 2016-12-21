
if n_elements(speed) eq 0 then speed = 400.0
if n_elements(density) eq 0 then density = 5.0
if n_elements(b) eq 0 then b = 5.0
if n_elements(t) eq 0 then t = 150000.0

speed   = float(ask('speed',string(speed)))
density = float(ask('density',string(density)))
b       = float(ask('b',string(b)))
t       = float(ask('t',string(t)))

gamma = 5.0/3.0

mp = 1.6726e-27
k  = 1.3807e-23
mu0 = 4*!pi*1e-7

n_u = density * 100.0^3
rho_u = n_u * mp
b_u = b * 1.0e-9
v_u = speed*1000.0
t_u = t

p_u = n_u * k * T_u

Cs = sqrt(gamma * p_u / rho_u)
Ca = B_u / sqrt(mu0 * rho_u)

Ms  = V_u / Cs
Ma  = V_u / Ca

print, "Alfven Speed : ", Ca/1000.0, " km/s"
print, "Sound Speed  : ", Cs/1000.0, " km/s"
print, "Alfven Mach Number : ", Ma
print, "Sound Mach Number  : ", Ms


end


