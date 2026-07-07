#import "template.typ": *
#import "macros.typ": *
#import "@preview/xarrow:0.3.1": xarrow
#import "@preview/fletcher:0.4.5" as fletcher: diagram, node, edge

#show: project.with(
  //title: "Dirac-Bergmann and Hamiltonian Field Theory",
  title: "Things I've derived before",
  authors: (
    (name: "Niels Slotboom", email: "slotboom.n@gmail.com"),
  ),
)

#outline(indent:auto)

= Finite Difference Stencils
== First Derivatives
Two-point forward difference stencil with error:
$
  f'(x) &= (f(x+epsilon)-f(x))/epsilon \ & quad- epsilon/2 f''(x) + cal(O)(epsilon²)
$
Three-point central difference stencil with error: 
$
  f'(x) &= (f(x+epsilon) - f(x-epsilon))/(2epsilon) \
  & quad - epsilon^2/6 f^((3))(x) + cal(O)(epsilon^3)
$
Five-point central difference stencil with error:
$
  f'(x) &= (-f(x+2epsilon) + 8f(x+epsilon) - 8f(x-epsilon) +f(x-2epsilon))/(12 epsilon)\
  &quad + epsilon^4/30 f^((5))(x) + cal(O)(epsilon^6)
$
== Second Derivatives
Three-point stencil with error:
$
  f''(x) &= (f(x+epsilon) - 2f(x) + f(x-epsilon))/epsilon^2\ &quad - epsilon^2/12 f^((4))(x) + cal(O)(epsilon^4)
$
Five-point stencil with error:
$
  f''(x) &= (-f(x+2epsilon) + 16 f(x+epsilon) - 30 f(x) + 16 f(x-
  epsilon) - f(x-2epsilon))/(12 epsilon^2)\ &quad- epsilon^4/90 f^((6))(x) + cal(O)(epsilon^6)
$

== Laplacian
27-point stencil: 
$
  Delta f = 1/epsilon^2 ((-8 lambda -6) dot "(center)" + (4 lambda+1) dot "(faces)" - 2lambda dot "(edges)" + lambda dot "(corners)")
$
where, using $f_(+0-) = f(x+epsilon,y,z-epsilon)$ etc., 
$
  "(center)" &= f_(000)\
  "(faces)" &= f_(+00) + f_(-00) + f_(0+0) + f_(0-0) +  f_(0 0 +) + f_(0 0 -)\
  "(edges)" &= f_(++0) + f_(+-0) + f_(-+0) + f_(--0)\
  &quad f_(+0+) + f_(+0-) + f_(-0+) + f_(-0-)\
  &quad f_(0++) + f_(0+-) + f_(0-+) + f_(0--)\
  "(corners)" &= f_(+++) + f_(++-) + f_(+-+) + f_(-++)\
  &quad f_(+--) + f_(-+-) + f_(--+) + f_(---).
$
The error is $cal(O)(epsilon^4)$, with leading contribution given by
$
  epsilon^2/12 ((diff^4 f)/(diff x^4) + (diff^4 f)/(diff y^4) + (diff^4 f)/(diff z^4)).
$
The choice of $lambda$ doesn't change this leading error order, but it changes the anisotropy of the error. Useful choices are $lambda = 1/22$ or $lambda = 1/26$. The latter minimises the anisotropy for Fourier modes, and its full expression reads
$
  Delta f = (-164 dot "(center)" + 30 dot "(faces)" - 2 dot "(edges)" + 1 dot "(corners)")/(26 epsilon^2) + cal(O)(epsilon^2)_"iso" + cal(O)(epsilon^4)_"aniso". wide
$