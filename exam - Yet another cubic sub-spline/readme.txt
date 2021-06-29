Name: Christian Norre
Student id: 201806854

54 mod 22 = 10

My project is 10. Yet another cubic sub-spline


COMMENTS:

From the example plots I would say that my cubic sub-spline is somewhere between regular cubic and Akima, 
	in the sense that my sub-spline doesn't eliminate wiggles imediately, like Akima, but one point over. 
This is no surprise, as for Akima, the Si'(xi) is calculated as a weighted linear combination of 
	nearby slopes, whereas my Si'(xi) is specified to just be my pi's. (second condition on assignment)

The good news is that I find pi's with a temporary quadratic spline, which means that I only have two
	parameters when I do cubic sub-spline, which allows me to easily solve the systems without considering 
	the Matrix method for regular cubic (see Homework "1-Interpolation").

I would rate my work as [6 + 3 + 0] = 9 out of 10 (on the same scale used for homework)
Considering how long we got for this exam, and how easy it would be, 
	I might as well have made my own Akima sub-spline, though, 
	I see that that would be another exam question.






WORKING OUT:

Quadratic spline:
Spline at point xi depends on both left and right point (according to assigment), which is different from our Homework.

I propose continuity in Si to left and right (previously we did continuity in Si and Si' to the right only)
Si(xi+1) = yi+1 = yi + bi(xi+1 - xi) + ci(xi+1 - xi)**2
Si(xi-1) = yi-1 = yi + bi(xi-1 - xi) + ci(xi-1 - xi)**2

I can solve for bi and ci

yi-1 - yi = ym,      yi+1 - yi = yp
xi-i - xi = xm,		 xi+1 - xi = xp
xp,yp for step i becomes -xm,-ym for step i+1

bi =  (yp*xm^2 - ym*xp^2) / (xm*xp(xm - xp))		Solved symbolically with MATLAB
ci = -(yp*xm   - ym*xp  ) / (xm*xp(xm - xp))

pi's can easily be found:
Si(x) = yi + bi(x-xi) +ci(x-xi)**2
Si'(x) = bi + 2*ci(x-xi)
pi = Si'(xi) = bi
p1 = S1'(x1) = b2 + 2*c2 * (x1-x2)
pn = Sn'(xn) = b_(n-1) + 2*c_(n-1) * (xn-x_(n-1))





Cubic sub-spline:

Si(x)= yi + bi(x-xi) + ci(x-xi)**2 + di(x-xi)**3

Conditions:
(1) Si(x(i+1)) = y(i+1)		Continuity to the right
(2) Si'(xi) = pi 			Derivative matches what is already known
(3) Si'(x(i+1)) = p(i+1)	Continuity of derivative to the right

(2) ->   Si'(xi) = bi = pi,	    pi is known, so bi is known

xp = x(i+1) - xi
(1) ->     yi + pi*xp +   ci*xp**2 +   di*xp**3  =  y(i+1)
(3) ->          pi    + 2*ci*xp    + 3*di*xp**2  =  p(i+1)

Two equations with two unknowns (ci and di)
yp = y(i+1) - yi

ci = (3*yp - (p(i+1) + 2*pi)*xp) / xp**2			Solved symbolically with MATLAB
di = ((p(i+1) + pi)*xp - 2*yp) / xp**3
There might be a way to calculate ci and di from c(i-1) and d(i-1), thus reducing calculations,
but i cannot see it

Each spline is between two points, i and i+1, so i need n-1 coefficients





Differentiation at point z:
find i, so xi < z < x(i+1)
Si'(z) = pi + 2*ci*(z-xi) + 3*di*(z-xi)**2



Integration from x0 to z:

for each interval xi to x(i+1), the integral gives
yi*xp + 1/2 * pi*xp**2 + 1/3 *ci*xp**3 + 1/4 *di*xp**4			(given xp = x(i+1)-xi)
which factors into:
1/12 * xp * [12*yi + 6*pi*xp + 4*ci*xp**2 + 3*di*xp**3]		(eq1)
given condition (1):
1/12 * xp * [3*y(i+1) + 9*yi  + 3*pi*xp +   ci*xp**2]		(eq2)

find i, so xi < z < x(i+1)
integral(Si(z)) = 
for j<i
	xp = x(j+1) - x(j)
	+ 1/12 * xp * [3*y(j+1) + 9*yj  + 3*pj*xp +   cj*xp**2]

xp = z - xi
+ 1/12 * xp * [3*Si(z) + 9*yi  + 3*pi*xp +   ci*xp**2]

unsurprisingly, eq1 and eq2 gives the same result









