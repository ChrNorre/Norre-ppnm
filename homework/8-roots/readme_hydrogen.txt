I change the shrodinger equation from a second order diff eq to two first order diff eqs
-(1/2)f'' -(1/r)f= E*f 
becomes f=y1, f'=y1'=y2, f''=y2'
y1' = y2
y2' = -2(E + 1/r)*y1

I then use my ode routines to solve from small radius (0.01) to large radius (10)
The function i give to my rootfinder is the entire function which initiates the ode_solver
The 'x' value (the value to change in order to reach 0) is 'E' in the diff eq
The 'fx' value (which should become 0) is the function value at large radius (effectively the result of ode_solver)




proving f(r->0) = r-r^2

f'' = -2(E+1/r)*f
insert proposed solution to check that it works
(0-2)  =  (r-r^2)''     =      -2(E+1/r)*(r-r^2)  =  -2Er+2Er^2 - 2+2r  =  -2 + (-2E + 2Er + 2)*r
0 = (-2E + 2Er + 2)*r
This is true when r->0