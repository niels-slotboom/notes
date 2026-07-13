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
The error is $cal(O)(epsilon^2)$, with leading contribution given by
$
  epsilon^2/12 ((diff^4 f)/(diff x^4) + (diff^4 f)/(diff y^4) + (diff^4 f)/(diff z^4)).
$
The choice of $lambda$ doesn't change this leading error order, but it changes the anisotropy of the error. Useful choices are $lambda = 1/22$ or $lambda = 1/26$. The latter minimises the anisotropy for Fourier modes, and its full expression reads
$
  Delta f = (-164 dot "(center)" + 30 dot "(faces)" - 2 dot "(edges)" + 1 dot "(corners)")/(26 epsilon^2) + cal(O)(epsilon^2)_"iso" + cal(O)(epsilon^4)_"aniso". wide
$
= The $3+1$-Formalism
== Foliations and Projectors
On a Lorentzian 3-manifold $(fM,g)$, given a function $t:fM->RR$ such that $g(dt,dt) <= 0$ everywhere, we call the covering of $fM$ by the sets 
$
  Sigma_t_0 = {p in fM | t(p) = t_0}
$
a _spacelike foliation_ of $(fM,g)$. If $(fM,g)$ is time-orientable and $dt$ future-oriented, then the foliation ${Sigma_t}$ defines a unique timelike unit normal vector field $n$, perpendicular to the foliation everywhere, by
$
  n = -alpha dt^sharp,
$
where $alpha>0$ is a function fixed by the normalisation condition
$
  -1 = g(n,n) = alpha^2 g^(t t) quad <=> quad g^(t t) = -1/alpha^2.
$
In components, we have
$
  n_mu = -alpha delta^t_mu, quad n^mu = -alpha g^(mu t).
$
The vector $n$ can be used to define an induced metric $gamma$ / a projection operator $P$,
$
  gamma = g + n^flat otimes n^flat quad <=> quad P = delta + n otimes n^flat
$
or equivalently in components,
$
  gamma_(mu nu) = g_(mu nu) + n_mu n_nu quad <=> quad tensor(P,+mu,-nu) = tensor(delta,+mu,-nu) + n^mu n_nu.
$
Note: although related by metric-induced isomorphism, the fully contra- and covariant object is denoted by $gamma_(mu nu)$ and $gamma^(mu nu)$, respectively, wherease the (1,1)-tensor is denoted $P$. This is because in its fully contra- and covariant forms, it is best interpreted as an induced metric, whereas the (1,1)-form is better interpreted as a projector.

This projector has the following properties:
+ $tensor(P,+mu,-nu) n^nu = 0$;

+ If $g(X,n) = 0$ ($<=> X in Gamma(T Sigma)$) then $tensor(P,+mu,-nu)X^nu =X^mu$;

+ $tensor(P,+mu,-nu) tensor(P,+nu,-lambda) = tensor(P,+mu,-lambda)$;

+ $gamma_(mu nu) = gamma_(nu mu)$;

+ $tensor(P,+mu,-mu) = 3$.

+ For $X,Y in Gamma(T Sigma)$, it holds that $g(X,Y) = gamma(X,Y)$.

These properties establish that pointwise, $P:T_p fM -> T_p Sigma_t(p)$ is a surjective orthogonal projection. In particular, $P|_(T Sigma) = id_(T Sigma)$, and $P|_(N Sigma) = 0$. The last property in the above shows that $gamma$ is indeed the induced metric on $Sigma$ when restricting arguments to $T Sigma$.

At this point, we introduce some additional notation. Given a tensor $T in T^((r,s))_p fM$, its projection $P T in T_p^((r,s)) Sigma$ is defined by
$
  tensor((P T),+mu...,-nu...) = tensor(P,+mu,-lambda) ... tensor(P,+rho,-nu)...tensor(T,+lambda...,-rho...),
$
that is, $P T$ has as components those of $T$ contracted with a projector on each index. In this notation we have, for example,
$
  P n = 0 quad "and" quad P g = gamma.
$
== Extrinsic Curvature
In the codimension 1 case, the extrinsic curvature tensor is defined by
$
  K(X,Y) = g(n,nabla_(P X) P Y) quad <=> quad K_(mu nu) X^mu Y^nu = n_mu (P X)^nu nabla_nu (P Y)^mu
$
for $X,Y in Gamma(T fM)$ and $nabla$ the Levi-Civita connection on $(fM,g)$. For $X,Y in Gamma(T Sigma)$, this reduces to
$
  K(X,Y) = g(n,nabla_X Y) = -g(nabla_X n, Y),
$
which hence measures the failure of parallel transport of $Y$ along $X$ to remain tangent. The second expression shows that equivalently, it measures the rate of change of the normal vector $n$ along $X$, as projected onto $Y$. 

Defining the acceleration $a$ as
$
  a = nabla_n n quad <=> quad a^mu = n^nu nabla_nu n^mu,
$
we have the following properties:
+ $K_(mu nu) = -tensor(P,+lambda,-mu) nabla_lambda n_nu = -nabla_mu n_nu - n_mu a_nu$;

+ $K_(mu nu) = K_(nu mu)$ as a consequence of vanishing torsion;

+ $n^mu K_(mu nu) = K_(mu nu) n^nu = 0$.

These follow straightforwardly either from the definition or the alternative characterisation of the extrinsic curvature above.

We denote the trace of the extrinsic curvature by $K$, and note that
$
  K = g^(mu nu) K_(mu nu) = gamma^(mu nu) K_(mu nu)
$
due to the transverse property 3.

Finally, we can relate $K_(mu nu)$ to the Lie-derivative of $gamma$ along $n$ as follows:
$
  K_(mu nu) = -1/2 \(fL_n gamma\)_(mu nu).
$
This follows by expanding the right-hand side and making use of 
$
  fL_n n^flat = a^flat.
$
Note that the Lie-derivative of the _vector_ $n$ would be zero, $fL_n n = [n,n] = 0$---in the identity above, we take the Lie-derivative of the _covector_ $n^flat$.

== Induced Connection and Intrinsic Curvature