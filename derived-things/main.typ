#import "template.typ": *
#import "macros.typ": *
#import "@preview/xarrow:0.3.1": xarrow
#import "@preview/fletcher:0.4.5" as fletcher: diagram, node, edge

#show: project.with(
  //title: "Dirac-Bergmann and Hamiltonian Field Theory",
  title: "Things I've derived and can trust",
  authors: (
    (name: "Niels Slotboom", email: "slotboom.n@gmail.com"),
  ),
)

#outline(indent:auto)

= Finite Difference Stencils
== Second Derivatives
Three-point stencil with error:
$
  f''(x) &= (f(x+epsilon) - 2f(x) + f(x-epsilon))/epsilon^2\ &quad - epsilon^2/12 f^((4))(x) + cal(O)(epsilon^4)
$
Five-point stencil with error:
$
  f''(x) &= (-f(x+2epsilon) + 16 f(x+epsilon) - 30 f(x) + 16 f(x-
  epsilon) - f(x-2epsilon))/(12 epsilon^2)\ &quad+ epsilon^4/90 f^((6))(x) + cal(O)(epsilon^6)
$