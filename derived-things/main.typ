#import "template.typ": *
#import "macros.typ": *
#import "@preview/xarrow:0.3.1": xarrow
#import "@preview/fletcher:0.4.5" as fletcher: diagram, node, edge

#show: project.with(
  //title: "Dirac-Bergmann and Hamiltonian Field Theory",
  title: "Stuff I've derived before",
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
= Runge-Kutta
== General Structure
Runge-Kutta is the name of a family of integrators for first-order differential equations of the form
$
  y'(t) = F(t,y(t)).
$
Generally speaking, given a timestep $epsilon$, one evaluates the right-hand side of the ODE $s$ times, at different times $t_i = epsilon c_i  in (t_0,t_0+epsilon)$ and for different guesses for $y$. Concretely, starting from $y(t_0) = y_0$, one evaluates
$
  k_1 &= F(t_0,y_0),\
  k_2 &= F(t_0 + epsilon c_2, y_0 + epsilon a_(2 1) k_1),\
  k_3 &= F(t_0 + epsilon c_3, y_0 + epsilon( a_(3 2) k_2 + a_(3 1) k_1)),\
  &#h(0.29em) dots.v\
  k_s &= F(t_0 + epsilon c_s, y_0 + epsilon(a_(s,s-1) k_(s-1) + ... + a_(s 1) k_1).
$
These are then used to approximate $y(t+epsilon)$ as
$
  y(t+epsilon) approx y_0 + epsilon (b_1 k_1 + ... + b_s k_s).
$
An RK integrator is hence defined by picking a value of $s$, quadrature points $c_i in [0,1]$, $i = 2,...,s$, as well as the coefficients $b_i$, $i=1,...,s$ and $a_(i j)$, $i > j$. These coefficients are not chosen at random, but rather picked such that the error of $y(t+epsilon)$ is as small as possible, i.e. of as high of an order $epsilon$ as can be arranged.

The coefficients can be thought of as a matrix $vA in RR^(s times s)$ and a pair of vectors $vb,vc in RR^s$, with the coefficients arranged such that
$
  vA = mat(0,,,dots.h,0;a_(21),0,,dots.h,0;a_31, a_32, 0,dots.h,0;dots.v,dots.v,dots.down,dots.down,dots.v;a_(s 1),a_(s 2),dots.h,a_(s,s-1),0), quad vb = mat(b_1;dots.v;b_s), quad vc = mat(0;c_2;dots.v;c_s).
$
Further writing $vk = (k_1,...,k_s)^top$, and denoting row evalution of a vector matrix by $[ dot ]_i$, we can write the above as follows:
$
  k_i = F(t_0 + epsilon[vc]_i, y_0 + epsilon [vA#h(0em) vk]_i), wide
  y(t+epsilon) approx y_0 + epsilon vb^top #h(0em) vk.
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
due to the transverse property (iii).

Finally, we can relate $K_(mu nu)$ to the Lie derivative of $gamma$ along $n$ as follows:
$
  K_(mu nu) = -1/2 \(fL_n gamma\)_(mu nu).
$
This follows by expanding the right-hand side and making use of 
$
  fL_n n^flat = a^flat.
$
Note that the Lie derivative of the _vector_ $n$ would be zero, $fL_n n = [n,n] = 0$---in the identity above, we take the Lie derivative of the _covector_ $n^flat$.

== Induced Connection and Intrinsic Curvature
=== Induced Connection
Given the Levi-Civita connection $nabla$ on $(fM,g)$, we define the _induced/three-dimensional/spatial covariant derivative_ $mnabla$ on the foliation $Sigma = {Sigma_t}$ as 
$
  mnabla T = P (nabla T) quad <=>quad tensor((mnabla_lambda T),+mu...,-nu...) = tensor(P,+alpha,-lambda) tensor(P,+mu,-beta) ... tensor(P,+gamma,-nu)... nabla_alpha tensor(T,+beta...,-gamma...) 
$
for any $T in Gamma(T^((r,s)) Sigma)$. Note that $T$ must be tangent to $Sigma$---i.e. have no normal components---for this to define a connection on $T Sigma$. The direction $X$ of the derivative does not necessarily have to be in $T Sigma$, but if it is, it holds that
$
  mnabla_X T = P (nabla_X T), quad "if" X in Gamma(T Sigma).
$
Otherwise, for general $X in Gamma(TM))$, we have
$
  mnabla_X T = P(nabla_(P X)T).
$
The fact that $mnabla$ defines a connection on $T Sigma$ is clear from the fact that $nabla$ is a connection, and that $P$ is $C^infty (fM)$-linear. In particular, we have the Leibniz rule
$
  mnabla_X (T otimes S) &= P (nabla_(P X) (T otimes S)) = P (nabla_(P X) T) otimes P S + P T otimes P(nabla_(P X) S)\
  &= (mnabla_X T) otimes S + T otimes (mnabla_X S).
$
Hence, for example, on a $1$-form $eta in Gamma(T^* Sigma)$,
$
  (mnabla_X eta)(Y) = X(eta(Y)) - eta(mnabla_X Y).
$

Further, $mnabla$ is torsion-free and metric-compatible with $gamma$:
- Torsion: For any $X,Y in Gamma(T Sigma)$, we have
  $
    macron(T)(X,Y) &= mnabla_X Y - mnabla_Y X - underbrace([X,Y],in Gamma(T Sigma))= P(nabla_X Y - nabla_Y X - [X,Y]) \ 
    &= P T(X,Y) = 0
  $

- Metricity: For any $X,Y,Z in Gamma(T Sigma)$,
  $
    X(gamma(Y,Z)) &= X(g(Y,Z)) = g(nabla_X Y,Z) + g(Y,nabla_X Z) \
    &= g(P nabla_X Y, Z) + g(Y,P nabla_X Z) = gamma(mnabla_X Y,Z) + gamma(Y,mnabla_X Z).
  $
=== Intrinsic Curvature
Since $mnabla$ defines a connection on $T Sigma$, and by that on each of the leaves of the foliation, it comes with an associated Riemannian curvature tensor, defined by
$
  macron(R)(X,Y)Z = mnabla_X mnabla_Y Z - mnabla_Y mnabla_X Z - mnabla_[X,Y]Z
$
In components, this reads
$
  (macron(R)(X,Y)Z)^mu = tensor(macron(R),+mu,-nu rho sigma) Z^nu X^rho Y^sigma.
$
We further have an associated Ricci tensor and scalar,
$
  macron(R)_(mu nu) = tensor(macron(R),+lambda,-mu lambda nu), quad macron(R) = g^(mu nu) macron(R)_(mu nu) = gamma^(mu nu) macron(R)_(mu nu),
$
where the last equality holds because $tensor(macron(R),+mu,-nu rho sigma)$ vanishes when contracted with the foliation normal $n$ on any index. 

Since $mnabla$ is the Levi-Civita connection with respect to $gamma$, $macron(R)_(mu nu rho sigma)$ enjoys the same symmetries as the ambient Riemann tensor. Further, the Ricci identity
$
  (mnabla_mu mnabla_nu thin - thin mnabla_nu mnabla_mu)V^lambda = tensor(macron(R),+lambda,-rho mu nu) V^rho
$
holds.

== Gauss-Codazzi-Ricci Equations
#theorem(name: "Gauss Equation")[
  With the definitions given in the above, the _Gauss equation_
  $
    tensor((P R),+lambda,-rho mu nu) = tensor(macron(R),+lambda,-rho mu nu) + tensor(K,+lambda,-mu) K_(rho nu) - tensor(K,+lambda,-nu) K_(rho mu)
  $
  holds. Contracting over $lambda mu$ leads to (after renaming indices)
  $
    (P R)_(mu nu) + R_(lambda mu rho nu) n^lambda n^rho = macron(R)_(mu nu) + K K_(mu nu) - tensor(K,+lambda,-mu) K_(nu lambda).
  $
  Tracing with respect to $gamma^(mu nu)$ (or equivalently, $g^(mu nu)$) leads to
  $
    R + 2 R_(mu nu) n^mu n^nu = macron(R) + K^2 - K_(mu nu) K^(mu nu).
  $
]
#proof[(Hint) To show this, first compute $mnabla_mu mnabla_nu X^sigma$ for $X in Gamma(T Sigma)$. Use $K_(mu nu) = -tensor(P,+lambda,-mu) nabla_lambda n_nu$ to trade derivatives of $n$ for the extrinsic curvature.]

#theorem(name: "Codazzi Equation")[
  Under the same assumptions, the _Codazzi equation_
  $
    P(R_(rho lambda mu nu ) n^rho) = mnabla_mu K_(nu lambda) - mnabla_nu K_(mu lambda)
  $
  holds. Tracing over $lambda mu$ with respect to $gamma^(lambda mu)$ (or equivalently, $g^(lambda mu)$) leads to
  $
    tensor(P,+mu,-lambda) n^nu R_(mu nu) = mnabla_lambda K - mnabla_mu tensor(K,+mu,-lambda) 
  $
]

#lemma[
  With the setup from the above, we have:
  + The following Lie derivatives of the projector/induced metric:
    $
      fL_n gamma^(mu nu) &= 2 K^(mu nu) + n^mu a^nu + a^mu n^nu,\
      fL_n tensor(P,+mu,-nu) &= n^mu a_nu,\
      fL_n gamma_(mu nu) &= -2 K_(mu nu)
    $
  + If $T in Gamma(T^((0,s)) Sigma)$, that is, if $T$ is a foliation-tangent covariant tensor, then so is $fL_n T$. That is,
    $
      fL_n T = P (fL_n T), quad T in Gamma(T^((0,s))Sigma).
    $
  + The acceleration covector $a_mu$ and the lapse $alpha$ satisfy the identity
    $
      a_mu = tensor(P,+lambda,-mu) nabla_lambda log alpha = mnabla_mu log alpha.
    $
]

#theorem(name: "Ricci Equation")[
  Again with the same setup, it holds the _Ricci equation_
  $
    tensor(P,+alpha,-mu) n^rho tensor(P,+beta,-nu) n^sigma R_(alpha rho beta sigma) = fL_n K_(mu nu) + K_(mu lambda) tensor(K,+lambda,-nu) + 1/alpha mnabla_mu mnabla_nu alpha
  $
  Note that on the left-hand side, the projections of the $alpha$ and $beta$ indices are superfluous, since the $n otimes n^flat$-part of $P$ vanishes when contracted due to the symmetries of the Riemann tensor. That is, we may write the Ricci equation more compactly as
  $
    R_(mu rho nu sigma) n^rho n^sigma = fL_n K_(mu nu) + K_(mu lambda) tensor(K,+lambda,-nu) + 1/alpha mnabla_mu mnabla_nu alpha.
  $
  Since both sides yield zero when contracted with $n$, the traces with respect to $g^(mu nu)$ and $gamma^(mu nu)$ are identical and read
  $
    R_(mu nu) n^mu n^nu = fL_n K - K_(mu nu) K^(mu nu) + 1/alpha mnabla^mu mnabla_mu alpha.
  $
]

With the triad of the Gauss-Codazzi-Ricci equations, we can now isolate the final projection of the Ricci tensor, the purely spatial one, and find an expression for the full Ricci scalar in terms of the intrinsic and extrinsic quantities.
#theorem[
  The above results combine into 
  $
    (P R)_(mu nu) &= - fL_n K_(mu nu) + K K_(mu nu) - 2 K_(mu lambda) tensor(K,+lambda,-nu) + macron(R)_(mu nu) - 1/alpha mnabla_mu mnabla_nu alpha,\
    R &= -2 fL_n K + K^2 + K_(mu nu) K^(mu nu) + macron(R) - 2/alpha mnabla^mu mnabla_mu alpha
  $<eqProjectionsRicci>
]
== The Einstein Equations in the $3+1$-Formalism
We work with the Einstein equations in the convention where
$
  R_(mu nu) - 1/2 R g_(mu nu) + Lambda g_(mu nu) = 8 pi T_(mu nu).
$
Taking the trace of this equation with respect to $g^(mu nu)$, we obtain
$
  1/2 R = 2 Lambda - 1/2 8 pi T ,
$
with $T = g^(mu nu) T_(mu nu)$ the trace of the energy-momentum tensor. Inserting this into the Einstein equations yields their trace-reversed counterpart,
$
  R_(mu nu) = 8pi (T_(mu nu) - 1/2 T g_(mu nu)) + Lambda g_(mu nu).
$
Note that a second trace-reversal brings us back to the original form, implying that the Einstein equations and their trace-reverse are equivalent.

We define the energy density $rho$, the energy current $j^mu$, and the stress tensor $S_(mu nu)$ as measured by Eulerian observers traveling along the integral lines of $n$ by
$
  rho &= n^mu n^nu T_(mu nu),\
  j_mu &= - tensor(P,+alpha,-mu) n^nu T_(alpha nu)\
  S_(mu nu) &= tensor(P,+alpha,-mu) tensor(P,+beta,-nu) T_(alpha beta)
$
These are simply the time-time, time-space and space-space projections of the energy-momentum tensor, and allow us to decompose it as
$
  T_(mu nu) = rho n_mu n_nu + j_mu n_nu + j_nu n_mu + S_(mu nu).
$
Its trace $T$ is given by
$
  T = S - rho
$
where $S = g^(mu nu) S_(mu nu) = gamma^(mu nu) S_(mu nu)$.

Together with the results from the previous section, these definitions now allow us to cast the Einstein equations into a first-order system for a second-order in time evolution of the spatial metric $gamma_(mu nu)$, with the extrinsic curvature $K_(mu nu)$ serving as an auxiliary variable:
#theorem[
  The Einstein equations written in terms of $alpha, n, gamma, K, rho, j$ and $S$ read
  #bottom-number[$
    cal(H) &:= macron(R) + K^2 - K_(mu nu) K^(mu nu) - 2Lambda - 16 pi rho = 0,\ \
    cal(M)_mu &:= mnabla_mu K - mnabla_nu tensor(K,+nu,-mu) + 8 pi j_mu = 0,\ \
    fL_n gamma_(mu nu) &= -2 K_(mu nu),\ \
    fL_n K_(mu nu) &= K K_(mu nu) - 2 K_(mu lambda) tensor(K,+lambda,-nu) + macron(R)_(mu nu) - 1/alpha mnabla_mu mnabla_nu alpha\ &wide - Lambda gamma_(mu nu) - 8pi ((S_(mu nu) - 1/2 S gamma_(mu nu)) + 1/2 rho gamma_(mu nu)).
  $<eqEvSys>]
]
The first two equations are derived by taking the normal-normal and normal-tangential projections of the original Einstein equations and subsequently identifying the definitions of $rho$ and $j$ as well as applying the Gauss and Codazzi equations, respectively. These two equations involve no time derivatives, and are hence a set of four constraints. 

The remaining two equations are first order in time (by the presence of the $fL_n$ normal derivatives), and hence dynamical. The first of the two is simply a consequence of the definition of the extrinsic curvature. The latter is the result of projecting the trace-reversed Einstein equations onto $T Sigma$ by contracting with $P$ on both indices, and subsequently applying the first equation in @eqProjectionsRicci[] to re-express the projected Ricci tensor $(P R)_(mu nu)$.

== Adapted Coordinates
=== The Metric and Bases
We can amend the time function $t:fM->RR$ (for which $dt != 0$) which defines the foliation $Sigma_t$ to a coordinate system _adapted to $Sigma_t$_ by introducing three additional functions $x^i : fM -> RR$ such that $dx^i != 0$ and $x^i|_Sigma_t$ label points uniquely on any leaf $Sigma_t$. This defines a basis ${diff_t, diff_i}$ of $TM$ and an associated dual basis ${dt,dx^i}$ of $T^*fM$. 

A natural question is to ask how the two time directions that we now have are related---we have the (coordinate-independent) timelike normal vector $n$, and additionally, the coordinate time direction $diff_t$. To this end, we introduce the shift vector $beta$, defined by
$
  beta := diff_t - alpha n quad <=> quad diff_t = alpha n + beta quad<=> quad n = 1/alpha (diff_t - beta).
$
This vector is tangent to the foliation, since
$
  dt(beta) = dt(diff_t) thin underbrace(- thin alpha dt,=n^flat)(n) = 1 - g(n,n) = 0.
$
Hence, it quantifies the failure of $diff_t$ to be normal to the foliation; that is, how much the coordinate lines of $t$ _shift_ with respect to the geometric normal direction specified by $n$.

With the lapse $alpha$, the shift $beta$, and the induced metric $gamma$, we now have all the objects we need for the 3+1-decomposition of the metric tensor $g$. Its components in the basis ${diff_t,diff_i}$  read
$
  g_(i j) &= g(diff_i,diff_j) = gamma(diff_i,diff_j) =: gamma_(i j)\
g_(t i) &= g(diff_t,diff_i) =  gamma_(i j) beta^j =: beta_i, \
  g_(t t) &= g(diff_t,diff_t) = -alpha^2 + beta_i beta^i.\
$
In matrix form, this reads
$
  [g_(mu nu)]_(mu,nu in {t,i}) = mat(-alpha^2 + beta_k beta^k, beta_j;beta_i,gamma_(i j)),
$
and tensorially, we may write
$
  g &= g_(mu nu) dx^mu otimes dx^nu\ &= -alpha^2 dt otimes dt + gamma_(i j) (dx^i + beta^i dt)(dx^j + beta^j dt).
$
The corresponding inverse metric has the components
$
  g^(t t) = -1/alpha^2, quad g^(t i) = beta^i/alpha^2, quad g^(i j) = gamma^(i j) - (beta^i beta^j)/alpha^2,
$
which takes the matrix form
$
  [g^(mu nu)]_(mu,nu in {t,i}) = mat(-1\/alpha^2, beta^j\/alpha^2;beta^i\/alpha^2, gamma^(i j) - beta^i beta^j\/alpha^2).
$
Here, $gamma^(i j)$ denotes the matrix inverse of $gamma_(i j)$, characterised uniquely by the condition
$
  gamma^(i j) gamma_(j k) = tensor(delta,+i,-k).
$
Tensorially, the inverse metric can be expressed compactly as 
$
  g^(-1) = -n otimes n + gamma^(i j) diff_i otimes diff_j.
$
As we have seen in the previous sections, we often care about normal projections of tensors, not their coordinate time components. Because of this, another useful basis is ${n,diff_i}$, alongside its dual ${-n^flat, dx^i} = {alpha dt, dx^i}$, with indices $perp, i$. Generically, this is _not_ a coordinate basis. However, it does give the metric a nice block-diagonal structure:
$
  g_(perp perp) = g(n,n) = -1, quad g_(perp i) = g(n,diff_i) = 0, quad g_(i j) = g(diff_i,diff_j) = gamma_(i j),
$
and
$
  g^(perp perp) = -1, quad g^(perp i) = 0, quad g^(i j) = gamma^(i j).
$
For example, in terms of this basis, the Hamiltonian and momentum constraints are nothing but the Einstein equations
$
  fH = G_(perp perp) - 8pi T_(perp perp) = 0, quad fM_i = G_(perp i) - 8pi T_(perp i).
$
In this basis, foliation-tangent tensors, i.e. objects $T in Gamma(T^((r,s)) Sigma)$, only have components where all indices are spatial. That is, a component vanishes if any of its indices is $perp$. 

Back in the coordinate basis ${diff_t,diff_i}$, the same is true for upstairs indices; for example, for a vector $X in Gamma(T Sigma)$, we have
$
  X^t = dt(X) = -1/alpha n^flat (X) = 0,
$
since $X perp n$. For downstairs indices, this is generally not true, as, for instance,
$
  X_t = g_(t mu) X^mu = g_(t i) X^i = beta_i X^i.
$
However, the $X_t$-component is not independent of the other three, $X_i = gamma_(i j) X^j$. Because of this, it is irrelevant, and hence, the full information about any foliation-tangent tensor $T in Gamma(T^((r,s)) Sigma)$ (such as $gamma_(mu nu)$, $K_(mu nu)$, $macron(R)_(mu nu rho sigma)$, $S_(mu nu)$ etc.) is stored in its purely spatial components. We also dont need the downstairs-, $t$ components for contractions; it will always be contracted with an upstairs-$t$ component, which is zero. Concretely, this implies that for example,
$
  K_(mu nu) K^(mu nu) = K_(i j) K^(i j) quad "and" quad tensor(macron(R),-i j) = tensor(macron(R),+mu,-i mu j) = tensor(macron(R),+k,-i k j).
$

=== Extrinsic and Intrinsic Curvature
#lemma[
  For any vector field $N$ that is normal to $T Sigma$ (not necessarily the unit-normal $n$), and a contravariant tangential tensor $T in Gamma(T^((0,s))Sigma)$, we have
  $
    fL_(f N) T = f fL_N T.
  $
  Note that this is not necessarily true if $T$ is not foliation-tangential or has upstairs indices.
]
This lemma might seem kind of useless, but it allows us to write down an explicit expression for the extrinsic curvature in an adapted coordinate basis: 
$
  K_(mu nu) = -1/2 fL_n gamma_(mu nu) = -1/2 fL_(1/alpha (diff_t-beta)) gamma_(mu nu) = -1/(2 alpha) (diff_t gamma_(mu nu) - fL_beta gamma_(mu nu))
$
Since the extrinsic curvature itself is also a foliation-tangent tensor, an analogous result holds. Writing directly in the adapted basis (and dropping any $t$-components), we obtain
$
  K_(i j) &= -1/(2alpha) (diff_t gamma_(i j) - fL_beta gamma_(i j)),\
  fL_n K_(i j) &= 1/alpha (diff_t K_(i j) - fL_beta K_(i j))
$
Writing out the Lie derivatives along $beta$, we obtain the following relationships:
$
  diff_t gamma_(i j) &= beta^k diff_k gamma_(i j) + 2 gamma_(k \(i) diff_(j\)) beta^k -2 alpha K_(i j), \
  diff_t K_(i j) &= beta^k diff_k K_(i j) + 2 K_(k \(i) diff_(j\)) beta^k +  alpha fL_n K_(i j).
$
These will be useful to recast the 3+1 Einstein equations we derived in the previous section into adapted coordinates, in a form where they clearly describe the first-order time evolution of $(gamma_(i j),K_(i j))$.

In the 3+1 Einstein equations, beyond the extrinsic curvature, the intrinsic connection $mnabla$ as well as its associated curavture $macron(R)$ appears. We therefore also need to write these in terms of the adapted coordinate basis. Since $mnabla$ is the Levi-Civita connection on any of the leaves with respect to the induced metric $gamma$, its components in adapted coordinates read
$
  tensor(macron(Gamma),+k,-i j) = 1/2 gamma^(k ell) (diff_i gamma_(j ell) + diff_j gamma_(i ell) - diff_ell gamma_(i j))
$
Correspondingly, the associated Riemann curvature has the components
$
  tensor(macron(R),+k,-ell i j) = diff_i tensor(macron(Gamma),+k,-ell j) - diff_j tensor(macron(Gamma),+k,-ell i) + tensor(macron(Gamma),+k,-n i) tensor(macron(Gamma),+n,-ell j) - tensor(macron(Gamma),+k,-n j) tensor(macron(Gamma),+n,-ell i)
$<eqTimeDerivGammaK>
=== The Einstein Equations in $3+1$-Adapted Coordinates
The @eqTimeDerivGammaK[equations] directly imply that the Einstein equations in adapted coordinates reads
$
  fH &= macron(R) + K^2 - K_(i j) K^(i j) - 2Lambda - 16pi rho = 0, \ \ \
  fM_i &= mnabla_i K - mnabla_j tensor(K,+j,-i) + 8pi j_i \ \ \
  diff_t gamma_(i j) &= beta^k diff_k gamma_(i j) + 2 gamma_(k \(i) diff_(j\)) beta^k - 2 alpha K_(i j), \ \ \
  diff_t K_(i j)&= beta^k diff_k K_(i j) + 2K_(k \(i) diff_(j\)) beta^k - mnabla_i mnabla_j alpha + alpha(K K_(i j) - 2 K_(i k) tensor(K,+k,-j) + macron(R)_(i j))\
  &quad - alpha Lambda gamma_(i j) - 8pi alpha (S_(i j) - 1/2 (S-rho) gamma_(i j))

$
The conservation equation $nabla_mu T^(mu nu) = 0$ further implies equations of motion for the energy density $rho$ and $j_i$. To derive them, we must first relate $rho$ and $j_i$ to components of the energy-momentum tensor in adapted coordinates. To this end, we should recall that in adapted coordinates,
$
  n^t &= 1/alpha, &quad&& n^i &= -beta^i/alpha,\
  n_t &= -alpha, &quad&& n_i &= 0.
$
Consequently, the projector $tensor(P,+mu,-nu)$ has the components
$
  tensor(P,+t,-t) &= 0, &&quad&
  tensor(P,+t,-i) &= 0,\
  tensor(P,+i,-t) &= beta^i, &&&
  tensor(P,+i,-j) &= tensor(delta,+i,-j),
$
which may be summarised as $tensor(P,+t,-mu) = 0$ and $tensor(P,+i,-mu) = tensor(delta,+i,-mu) + beta^i tensor(delta,+t,-mu)$
With this projector and the components of $n$, we can directly evaluate
$
  rho &= T^(mu nu) n_mu n_nu = alpha^2 T^(t t),\
  j^i &= -tensor(P,+i,-mu) T^(mu nu) n_nu = alpha (beta^i T^(t t) + T^(i t)),\
  S^(i j) &= tensor(P,+i,-mu) tensor(P,+j,-nu) T^(mu nu) = T^(i j) + 2 beta^(\(i) T^(j\) t) + beta^i beta^j T^(t t).
$
Here, we recall that $j^t = 0$, $S^(t mu) = 0$, and that hence, we may raise and lower indices with $gamma_(i j)$. Note that $j_t$ and $S_(t mu)$ are not necessarily zero, but by the reasoning given earlier, these components are not independent of the purely spatial ones and hence carry no additional information. 

The above puts $rho,j^i$ and $S^(i j)$ in terms of the components of $T^(mu nu)$. To evaluate projections of $nabla_mu T^(mu nu) = 0$, however, we need to invert this relationship, and express $T^(mu nu)$ in terms of $rho,j^i$ and $S^(i j)$. Carrying out this inversion leads to
$
  T^(t t) = 1/alpha^2 rho, wide T^(i t) = 1/alpha j^i - beta^i/alpha^2 rho, wide T^(i j) = S^(i j) - 2/alpha beta^(\(i) j^(j\)) + (beta^i beta^j)/alpha^2 rho.
$
Note that this is simply saying that
$
  T = rho n otimes n + j otimes n + n otimes j + S.
$
The normal projection $(nabla_mu T^(mu nu)) n_nu = 0$ then implies the equation of motion for $rho$, reading
$
  diff_t rho = beta^k diff_k rho - 2 j^k diff_k alpha + alpha (K rho + K_(i j) S^(i j) - mnabla_i j^i).
$
This is derived by first integrating by parts, $(nabla_mu T^(mu nu)) n_nu = nabla_mu (T^(mu nu) n_nu) - T^(mu nu) nabla_mu n_nu$, and then applying identities from the above to replace derivatives of $n$ with the extrinsic curvature. The standard formula $nabla_mu X^mu = 1/sqrt(g)diff_mu (sqrt(g)X^mu)$, together with $sqrt(g) = alpha sqrt(gamma)$, then allows the deduction of the result.

A similar equation of motion can be derived for $diff_t j_i$ by projecting $nabla_mu T^(mu nu) = 0$ onto $T Sigma$ using $P$, but I will not do this here. 
= Metric Derivative and Curvature Identities
== Jacobi Formula
For the variation $delta$ of a matrix $M$, we have the _Jacobi formula_
$
  delta det M = det M dot tr(M^(-1) delta M).
$
Using $g = det [g_(mu nu)]$, this implies that
$
  delta g = g g^(mu nu) delta g_(mu nu) = - g g_(mu nu) delta g^(mu nu).
$
For a derivative variation, $delta = diff_lambda$, this can be expressed in terms of Christoffel symbols as
$
  diff_lambda g = g g^(mu nu) diff_lambda g_(mu nu) = 2 g tensor(Gamma,+mu,-mu lambda).
$
This expression further implies expressions for the divergence of $X in Gamma(TM)$ and the wave operator on $f in C^infty (fM)$, reading
$
  Div X = nabla_mu X^mu = 1/sqrt(g) diff_mu (sqrt(g) X^mu)\ \ "and"\ \ Box_g f = g^(mu nu) nabla_mu nabla_nu f = 1/sqrt(g) diff_mu (sqrt(g) g^(mu nu) diff_nu f),
$
respectively.
== Connection and Curvature of Conformal Metrics
In the following, let $g_(mu nu)$ and $tilde(g)_(mu nu)$ be two metrics on a $d$-dimensional manifold $fM$, related by
$
  tilde(g)_(mu nu) = e^(4 phi.alt) g_(mu nu).
$
Trivially, their inverses and determinants are then related by
$
  tilde(g)^(mu nu) = e^(-4 phi.alt) g^(mu nu) quad "and" quad tilde(g) = e^(4d phi.alt) g.
$
Denoting their Christoffel symbols by $Gamma$ and $tilde(Gamma)$, respectively, we have the identity
$
  tensor(tilde(Gamma),+lambda,-mu nu) = tensor(Gamma,+lambda,-mu nu) + 2 (tensor(delta,+lambda,-mu) diff_nu phi.alt + tensor(delta,+lambda,-nu) diff_mu phi.alt - g_(mu nu) g^(lambda rho) diff_rho phi.alt)
$
or equivalently,
$
  tensor(Gamma,+lambda,-mu nu) = tensor(tilde(Gamma),+lambda,-mu nu) - 2 (tensor(delta,+lambda,-mu) diff_nu phi.alt + tensor(delta,+lambda,-nu) diff_mu phi.alt - g_(mu nu) g^(lambda rho) diff_rho phi.alt).
$
We note that the difference between the two connection coefficients,
$
  tensor(C,+lambda,-mu nu) := tensor(tilde(Gamma),+lambda,-mu nu) - tensor(Gamma,+lambda,-mu nu) = 2 (tensor(delta,+lambda,-mu) tnabla_nu phi.alt + tensor(delta,+lambda,-nu) tnabla_mu phi.alt - tilde(g)_(mu nu) tilde(g)^(lambda rho) tnabla_rho phi.alt)
$
is a tensor, and independent of whether one writes it in terms of $g_(mu nu)$ or $tilde(g)_(mu nu)$ in the last term. This means that the Levi-Civita connection actions on an arbitrary tensor $T in Gamma(T^((r,s))fM)$ are related by
$
  tilde(nabla)_lambda tensor(T,+mu...,-nu...) = nabla_lambda tensor(T,+mu...,-nu...) + tensor(C,+mu,-rho lambda)tensor(T,+rho...,-nu...) +... - tensor(C,+rho,-nu lambda)tensor(T,+mu...,-rho...) - ...
$<eqConformalConnection>
Concretely, using the action on vectors, one can derive that
$
  tensor(R,+lambda,-rho mu nu) = tensor(tilde(R),+lambda,-rho mu nu) - tnabla_mu tensor(C,+lambda,-rho nu) + tnabla_nu tensor(C,+lambda,-rho mu) + tensor(C,+lambda,-sigma mu) tensor(C,+sigma,-rho nu) - tensor(C,+lambda,-sigma nu) tensor(C,+sigma,-rho mu).
$
By contracting over $lambda mu$, we find a relationship between the Ricci tensors reading
$
   R_(mu nu) &= tilde(R)_(mu nu) + 2 tilde(g)_(mu nu) tilde(g)^(rho sigma) tnabla_rho tnabla_sigma phi.alt  + 2 (d-2)(tnabla_mu tnabla_nu phi.alt + 2 tnabla_mu phi.alt tnabla_nu phi.alt - 2 tilde(g)_(mu nu) tilde(g)^(rho sigma) tnabla_rho phi.alt tnabla_sigma phi.alt)
$<eqConformalRicciTensor>
For the Ricci scalars---taking note that $tilde(R) = tilde(g)^(mu nu) tilde(R)_(mu nu)$ and $R = g^(mu nu) R_(mu nu)$ are traced with respect to their respective metric---we obtain 
$
  R = e^(4phi.alt)(tilde(R) + 4(d-1)tilde(g)^(mu nu) tnabla_mu tnabla_nu phi.alt -4(d-2)(d-1) tilde(g)^(mu nu) tnabla_mu phi.alt tnabla_nu phi.alt).
$<eqConformalRicciScalar>
== Bochner Formula and the Ricci Tensor
For two scalar fields $phi$ and $psi$, the _Bochner formula_
$
  1/2 Box_g (nabla_mu phi nabla^mu psi) &= 1/2 (nabla_mu Box_g phi) nabla^mu psi + 1/2 nabla_mu phi (nabla^mu Box_g psi) \ 
  &quad+ nabla_mu nabla_nu phi nabla^mu nabla^nu psi + R^(mu nu) nabla_mu phi nabla_b psi.
$
Using the coordinate functions $phi = x^mu$, $psi = x^nu$, it follows that
$
  R^(mu nu) = 1/2 g^(lambda rho) diff_lambda diff_rho g^(mu nu) - 1/2 Gamma^lambda diff_lambda g^(mu nu) + diff^(\(mu) Gamma^(nu\)) - tensor(Gamma,+mu,-lambda rho) Gamma^(nu lambda rho),
$
where $Gamma^mu = g^(lambda rho) tensor(Gamma,+mu,-lambda rho)$.
Explicitly lowering the indices leads to
$
  R_(mu nu) &= -1/2 g^(lambda rho) diff_lambda diff_rho g_(mu nu) + g_(lambda \(mu) diff_(nu\)) Gamma^lambda\
  & quad + 1/2 Gamma^lambda diff_lambda g_(mu nu) + g^(lambda rho) g^(sigma tau) diff_lambda g_(mu sigma) diff_rho g_(nu tau) - Gamma_(mu rho sigma) tensor(Gamma,-nu,+rho sigma).
$
This expression will motivate why in BSSNOK, one promotes a version of $Gamma^mu$ to be a dynamical variable; in this case, the only second-order part in $R_(mu nu)$ is the well-behaved (in that case, elliptic) expression $g^(lambda rho) diff_lambda diff_rho g_(mu nu)$. 

For numerical implementations, we should address one issue the above has, though. This issue is that the expression we obtained for the Ricci tensor depends on both $diff g$ as well as $Gamma$. This increases register pressure in GPU-based implementations, which we would like to avoid; hence, we should make use of the identity
$
  diff_lambda g_(mu nu) = Gamma_(mu nu lambda) + Gamma_(nu mu lambda)
$
to re-express everything in terms of Christoffel symbols. Doing so yields
$
  R_(mu nu) &= -1/2 g^(lambda rho) diff_lambda diff_rho g_(mu nu) + g_(lambda\(mu) diff_(nu\)) Gamma^lambda \
  &quad+ thin Gamma^lambda Gamma_((mu nu) lambda) + Gamma_(lambda rho mu) tensor(Gamma,+lambda rho, -nu) + 2 Gamma_(lambda rho \(mu) tensor(Gamma,-nu\),+lambda rho),
$<eqRicciTensorNice>
which is now written entirely in terms of $g$, $diff^2 g$, and the two flavours of $Gamma$.
= BSSNOK
In this section, we take what we have derived so far and build atop it the BSSNOK formalism. This involves

- defining the BSSNOK variables in terms of the ADM variables $gamma_(i j)$ and $K_(i j)$, scaling out the conformal factor to isolate the conformal metric $tilde(gamma)_(i j)$ and trace-free extrinsic curvature $tilde(A)_(i j)$;

- deriving their first time derivative expressions from the known expressions for $diff_t gamma_(i j)$ and $diff_t K_(i j)$ via the ADM evolution equations, as well as re-expressing the Hamiltonian and momentum constraints in terms of the new conformal variables;

- and finally, relating the physical induced connections and curvatures to their conformal counterparts using conformal transformations and the auxiliary connection quantities.

In doing so, we will use many of the identities from the previous section, and establish a number of additional ones to handle the conformal rescaling. Unfortunately, these algebraic definitions are rather lengthy and messy, but they ultimately lead to a strongly hyperbolic formulation of general relativity well-suited to be implemented in high-performance computing (HPC) code. 

== BSSNOK Reparametrisation and Evolution Equations
#definition(name: "BSSN Variables")[The dynamical variables in the BSSN formulation of GR are

+ the _conformal factor_ $phi.alt$, defined as
  $
    phi.alt = -1/12 log gamma,
  $
+ the _conformal metric_ $tilde(gamma)_(i j)$, given by
  $
    tilde(gamma)_(i j) = e^(4 phi.alt) gamma_(i j),
  $
+ the trace of the extrinsic curvature $K$, defined by
  $
    K = gamma^(i j) K_(i j),
  $
+ the traceless conformal extrinsic curvature $tilde(A)_(i j)$, expressed as
  $
    tilde(A)_(i j) = e^(4phi.alt)(K_(i j) - 1/3 K gamma_(i j))
  $
+ and the contracted conformal connection $tilde(Gamma)^i$, defined as
  $
    tilde(Gamma)^i = tilde(gamma)^(k ell) tensor(tilde(Gamma),+i,-k ell), quad "with" quad tensor(tilde(Gamma),+i,-k ell) = 1/2 tilde(gamma)^(i j) (diff_k tilde(gamma)_(ell j) + diff_ell tilde(gamma)_(k j) - diff_j tilde(gamma)_(k ell)).
  $
]
#corollary[
  These definitions have a number of important immediate corollaries:
]
  + The determinant of the conformal metric $tilde(gamma)_(i j)$ is 1, that is,
    $
      tilde(gamma) = e^(3 dot 4 phi.alt) gamma = e^(-log gamma) gamma = 1.
    $
    This has a few important consequences; firstly, the trace of its Christoffel symbols vanishes,
    $
      tensor(tilde(Gamma),+i,- k i) = 1/2 tilde(gamma)^(i j) diff_k tilde(gamma)_(i j) = 1/sqrt(tilde(gamma)) diff_k sqrt(tilde(gamma)) = diff_k 1 = 0.
    $
    Moreover, its contracted Christoffel symbols $tilde(Gamma)^i$ may be written as
    $
      tilde(Gamma)^i = - diff_k tilde(gamma)^(k i)
    $

  + $tilde(A)_(i j)$ is traceless with respect to both $tilde(gamma)^(i j)$ and $gamma^(i j)$, that is,
    $
      tilde(gamma)^(i j) tilde(A)_(i j) = gamma^(i j) tilde(A)_(i j) = 0.
    $

  + The definitions of $phi.alt, tilde(gamma)_(i j), K$ and $tilde(A)_(i j)$ can be inverted for $gamma_(i j)$ and $K_(i j)$, with the inversions reading
    $
      gamma_(i j) = e^(-4phi.alt) tilde(gamma)_(i j) quad "and" quad K_(i j) = e^(- 4phi.alt) (tilde(A)_(i j) + 1/3 K tilde(gamma)_(i j))
    $

  + From the point of view of the dynamical evolution equations we are about to derive for $tilde(gamma)_(i j), tilde(A)_(i j)$ and $tilde(Gamma)^i$, these are two symmetric $3 times 3$-matrix- and one $3$-vector-valued functions. The above tells us that in addition to these evolution equations, we get additional constraints; algebraically, we must always have
    $
      tilde(gamma) = 1, quad tilde(gamma)^(i j) tilde(A)_(i j) = 0,
    $
    and, in addition, a differential constraint emerges from
    $
      fG^i := tilde(Gamma)^i - tilde(gamma)^(k ell) tensor(tilde(Gamma),+i,-k ell) = 0.
    $
    Although mathematically speaking, the evolution equations will preserve these, numerically, there is a slight drift away from them. The algebraic constraints are very easy to enforce numerically---during discrete evolution we can simply regularly remap
    $
      tilde(gamma)_(i j) |-> tilde(gamma)^(-1\/3) tilde(gamma)_(i j), quad tilde(A)_(i j) |-> tilde(A)_(i j) - 1/3 tilde(gamma)_(i j) tensor(tilde(A),+k,-k).
    $
    To preserve the differential constraint $fG^i = 0$ even in numerical simulations, we will artificially add it to the right-hand side of the evolution equation for $tilde(Gamma)^i$, as 
    $
      diff_t tilde(Gamma)^i = -sigma fG^i + ...
    $
    for some positive constant $sigma > 0$. This creates exponential damping/decay of $fG^i$ towards 0, so that numerical error does not diverge.
  + The trace of the square of $K_(i j)$ reads
    $
      K_(i j) K^(i j) = tilde(A)_(i j) tilde(A)^(i j) + 1/3 K^2.
    $
    Here, indices of quantities with a tilde are raised and lowered using $tilde(gamma)$, whereas for quantities without tilde, $gamma$ is used. We will continue to use this convention in the following to omit writing hundreds of instances of the induced and conformal metrics.

#proposition(name: "Preliminary BSSNOK Evolution Equations")[
  Below are the evolution equations for the BSSNOK variables $phi.alt, tilde(gamma)_(i j), K, tilde(A)_(i j)$ and $tilde(Gamma)^i$.
  #bottom-number[$
    fH &= macron(R) + 2/3 K^2 - tilde(A)_(i j) tilde(A)^(i j) - 2Lambda - 16 pi rho = 0,\ \
    fM_i &= 2/3 diff_i K - tnabla_k tensor(tilde(A),+k,-i) + 6 tensor(tilde(A),+k,-i) diff_k phi.alt + 8pi j_i = 0,\ \
    diff_t phi.alt &= beta^k diff_k phi.alt - 1/6 diff_k beta^k + 1/6 alpha K,\ \
    diff_t tilde(gamma)_(i j) &= beta^k diff_k tilde(gamma)_(i j) + 2 tilde(gamma)_(k\(i) diff_(j\)) beta^k - 2/3 tilde(gamma)_(i j) diff_k beta^k - 2alpha tilde(A)_(i j), \ \
    diff_t K&= beta^k diff_k K - gamma^(i j) mnabla_i mnabla_j alpha + alpha (tilde(A)_(i j) tilde(A)^(i j) + 1/3 K^2 - Lambda + 4pi (S+rho)),\ \
    diff_t tilde(A)_(i j) &= beta^k diff_k tilde(A)_(i j) + 2 tilde(A)_(k\(i) diff_(j\)) beta^k - 2/3 tilde(A)_(i j) diff_k beta^k + e^(4phi.alt) [alpha macron(R)_(i j) - mnabla_i mnabla_j alpha - 8 pi alpha S_(i j)]^"TF"\
    &quad  +alpha (K tilde(A)_(i j) - 2 tilde(A)_(i k) tensor(tilde(A),+k,-j))\ \ 
    diff_t tilde(Gamma)^i &= beta^k diff_k tilde(Gamma)^i - tilde(Gamma)^k diff_k beta^i + 2/3 tilde(Gamma)^i diff_m beta^m + tilde(gamma)^(j ell) diff_j diff_ell beta^i + 1/3 tilde(gamma)^(i j) diff_j (diff_k beta^k)\
    &quad- 2tilde(A)^(i j) diff_j alpha + 2alpha tensor(tilde(Gamma),+i,-j k) tilde(A)^(j k) - 4/3 alpha tilde(gamma)^(i j) diff_j K  - 12 alpha tensor(tilde(A),+i j) diff_j phi.alt - 16 pi alpha tilde(gamma)^(i k) j_k - sigma fG^i
  $<eqBSSN>
  In the above, the $#h(0em)^"TF"$ exponent refers to the trace-free part of the bracketed expression, i.e.
  $
    X_(i j)^"TF" = X_(i j) - 1/3 gamma_(i j) gamma^(k ell) X_(k ell) = X_(i j) - 1/3 tilde(gamma)_(i j) tilde(gamma)^(k ell) X_(k ell).
  $
  Noteworthily, it does not matter whether trace removal is carried out using $gamma_(i j)$ or $tilde(gamma)_(i j)$.
  ]
]
#remark[
]
+ In the derivation of the right-hand side for $diff_t K$, $alpha fH = 0$ was subtracted. This does not affect the validity of the equality on-shell, but changes the principal part of the evolution system. Nonetheless, sometimes it is useful to have the same equation without this subtraction. For such situations, we provide it below:
    $
      diff_t K = beta^k diff_k K - gamma^(i j) mnabla_i mnabla_j alpha + alpha (K^2 + macron(R)) - 3 alpha Lambda + 4 pi alpha(S-3rho)
    $
    One such situation is in the derivation of the right-hand side for $diff_t tilde(A)_(i j)$ above, where it turns out to be more convenient to use the unsubtracted version. The reason to subtract the Hamiltonian constraint from this equation is to eliminate the intrinsic Ricci scalar $macron(R)$. In the standard ADM formulation, $macron(R)$ introduces second spatial derivatives of the metric into $diff_t K$, which contributes to weak hyperbolicity and numerical instabilities. Subtracting $alpha fH$ removes all second derivatives of the metric from the right-hand side of $diff_t K$, which helps cast the BSSNOK evolution system into a strongly hyperbolic form.

+ This version of the equations is only preliminary as it contains expressions involving $macron(R)_(i j)$, $macron(R)$ and $mnabla$, which are associated with the non-conformal quantity $gamma_(i j)$ which is not available directly to the evolution code. Concretely, the offending terms in the above equations which still contain objects tied to $gamma_(i j)$ are $macron(R)$ in the Hamiltonian constraint, the Laplacian $macron(Delta)_gamma alpha = gamma^(i j) mnabla_i mnabla_j alpha$ and the intrinsic Ricci tensor $macron(R)_(i j)$. We will need to re-express these in terms of conformal quantities---this requires the results presented below.

+ Since the evolution equation for $tilde(Gamma)^i$ contains an artificially introduced damping term $-sigma fG^i$ with $sigma>0$, we present the full derivation of the equation below. We begin by expanding
  #bottom-number($
    diff_t tilde(Gamma)^i &= - diff_t diff_j tilde(gamma)^(i j) = -diff_j diff_t tilde(gamma)^(i j) = diff_j (tilde(gamma)^(i k) tilde(gamma)^(j ell) diff_t tilde(gamma)_(k ell))\
    &= diff_j (tilde(gamma)^(i k) tilde(gamma)^(j ell) (beta^m diff_m tilde(gamma)_(k ell) + tilde(gamma)_(m k) diff_ell beta^m + tilde(gamma)_(m ell) diff_k beta^m - 2/3 tilde(gamma)_(k ell) diff_m beta^m - 2alpha tilde(A)_(k ell)))\
    &= diff_j (-beta^m diff_m tilde(gamma)^(i j) + tilde(gamma)^(j ell) diff_ell beta^i + tilde(gamma)^(i ell) diff_ell beta^j - 2/3 tilde(gamma)^(i j) diff_m beta^m - 2alpha tilde(A)^(i j))\

    &= cancelr(-(diff_j beta^m) diff_m tilde(gamma)^(i j)) - beta^m diff_m underbrace(diff_j tilde(gamma)^(i j),=-tilde(Gamma)^i) + underbrace((diff_j tilde(gamma)^(j ell)),=-tilde(Gamma)^ell) diff_ell beta^ i + tilde(gamma)^(j ell) diff_j diff_ell beta^i + cancelr((diff_j tilde(gamma)^(i ell)) diff_ell beta^j) + underline(tilde(gamma)^(i ell) diff_j diff_ell beta^j)\
    &quad -2/3 underbrace((diff_j tilde(gamma)^(i j)),=-tilde(Gamma)^i) diff_m beta^m underline(-2/3 tilde(gamma)^(i j) diff_j diff_m beta^m) - 2tilde(A)^(i j) diff_j alpha - 2 alpha diff_j tilde(A)^(i j)\
    &= beta^k diff_k tilde(Gamma)^i - tilde(Gamma)^k diff_k beta^i + 2/3 tilde(Gamma)^i diff_m beta^m + tilde(gamma)^(j ell) diff_j diff_ell beta^i + 1/3 tilde(gamma)^(i j) diff_j (diff_k beta^k)\
    &quad- 2tilde(A)^(i j) diff_j alpha - 2 alpha diff_j tilde(A)^(i j)
  $)
  The last term is still in an unfortunate form, which we can improve by using the momentum constraint $fM^i = 0$. Concretely, it may be used in re-arranged form to replace $tnabla_j tilde(A)^(i j)$ as
  $
    diff_j tilde(A)^(i j) &= tnabla_j tilde(A)^(i j) - tensor( tilde(Gamma),+i,-j k) tilde(A)^(j k) - underbrace(tensor(tilde(Gamma),+j,-k j),=0) tilde(A)^(k i)\
    &= - tensor(tilde(Gamma),+i,-j k) tilde(A)^(j k) + 2/3tilde(gamma)^(i j) diff_j K  + 6 tensor(tilde(A),+i j) diff_j phi.alt + 8pi j^i
  $
  Consequently, the evolution equation for $tilde(Gamma)^i$, with the artificial damping term $-sigma fG^i$ introduced, reads
  $
    diff_t tilde(Gamma)^i &= beta^k diff_k tilde(Gamma)^i - tilde(Gamma)^k diff_k beta^i + 2/3 tilde(Gamma)^i diff_m beta^m + tilde(gamma)^(j ell) diff_j diff_ell beta^i + 1/3 tilde(gamma)^(i j) diff_j (diff_k beta^k)\
    &quad- 2tilde(A)^(i j) diff_j alpha + 2alpha tensor(tilde(Gamma),+i,-j k) tilde(A)^(j k) - 4/3 alpha tilde(gamma)^(i j) diff_j K  - 12 alpha tensor(tilde(A),+i j) diff_j phi.alt - 16pi alpha j^i - sigma fG^i
  $

+ The keen-eyed reader might have noticed that for any of the BSSN variables, its time derivative contains similarly structured spatial derivative terms of that variable involving the shift $beta$ on the right hand side. These combinations are reminiscient of Lie derivatives---and in fact, they are---with a catch. The catch comes from the terms involving the divergence $diff_i beta^i$, which in a tensorial Lie derivative do not appear. However, we should recall that $phi.alt$ is a function of the _determinant_ of the tensor $gamma_(i j)$, and that in the definitions of $tilde(gamma)_(i j)$ as well as $tilde(A)_(i j)$ include factors of 
  $
    e^(4 phi.alt) = gamma^(-1\/3).
  $
  From this, we conclude that $phi.alt$ is not a scalar but a scalar _density_ of weight $-1/6$, and that $tilde(gamma)_(i j)$ and $tilde(A)_(i j)$ are $(0,2)$-tensor _densities_ of weight $-2/3$. We recall that for a tensor density $tensor(T,+mu...,-nu...)$ of weight $w$, its Lie derivative along some vector field $X$ is given by
  $
    fL_X tensor(T,+mu...,-nu...) &= X^lambda diff_lambda tensor(T,+mu...,-nu...) - (diff_lambda X^mu) tensor(T,+lambda...,-nu...) -...\ &quad + (diff_nu X^lambda) tensor(T,+mu...,-lambda...) + ... + w (diff_lambda X^lambda) tensor(T,+mu...,-nu...)
  $<defLieDerivDensity>
  A similar expression can be obtained for connections as well. Taking $X = beta$, $T = phi.alt, tilde(gamma)_(i j), K, tilde(A)_(i j)$ and $tilde(Gamma)^i$, we identify that the @eqBSSN[BSSNOK evolution equations in] may be written in the (slightly) more compact form
  $
    diff_t phi.alt &= fL_beta phi.alt + 1/6 alpha K,\
    diff_t tilde(gamma)_(i j) &= fL_beta tilde(gamma)_(i j) - 2 alpha tilde(A)_(i j),\
    diff_t K &= fL_beta K - gamma^(i j) mnabla_i mnabla_j alpha + alpha (tilde(A)_(i j) tilde(A)^(i j) + 1/3 K^2 - Lambda + 4pi (S+rho)),\
    diff_t tilde(A)_(i j)&= fL_beta tilde(A)_(i j) + e^(4phi.alt) [alpha macron(R)_(i j) - mnabla_i mnabla_j alpha - 8pi alpha S_(i j)]^"TF" + alpha (K tilde(A)_(i j) - 2 tilde(A)_(i k) tensor(tilde(A),+k,-j)),\
    diff_t tilde(Gamma)^i &= fL_beta tilde(Gamma)^i - 2 tilde(A)^(i j) diff_j alpha + 2 alpha tensor(tilde(Gamma),+i,-j k) tilde(A)^(j k) - 4/3 alpha tilde(gamma)^(i j) diff_j K - 12 alpha tilde(A)^(i j) diff_j phi.alt - 16 pi alpha tilde(gamma)^(i k) j_k - sigma fG^i.
  $
  Note that the trace of the extrinsic curvature, $K$, is the only true scalar, so $fL_beta K = beta^i diff_i K$. Moreover, it should be noted that $fL_beta tilde(Gamma)^i$ is _not_ to be interpreted as the Lie derivative of a vector field nor vector density, but rather expanded as
  $
    fL_beta tilde(Gamma)^i = fL_beta (tensor(tilde(Gamma),+i,-j k) tilde(gamma)^(j k)) = (fL_beta tensor(tilde(Gamma),+i,-j k)) tilde(gamma)^(j k) + tensor(tilde(Gamma),+i,-j k) (fL_beta tilde(gamma)^(j k)).
  $
  The Lie derivative of the connection coefficients is then evaluated using the analogon of @defLieDerivDensity for connections, and gives rise to the terms involving second derivatives of $beta$ in the last of @eqBSSN[equations].

  Although this compactification of the BSSN evolution equations does not introduce any benefits to its numerical implementation---we still have to compute all the individual terms the Lie derivatives consist of---writing them down in this form enables a more straightforward interpretation. For any of the $T = phi.alt, tilde(gamma)_(i j), K, tilde(A)_(i j), tilde(Gamma)^i$, we have the main Lie-advection piece
  $
    diff_t T = fL_beta T.
  $
  In absence of any other terms, this implies that the BSSN variables are Lie-transported along the normal vector field $alpha n = diff_t -beta$ and thus remain invariant under its generated diffeomorphisms. However, there are other terms; these source the Lie transport by introducing curvature, matter and gauge contributions, such as $e^(4phi.alt)alpha macron(R)_(i j)$, $8pi e^(4phi.alt) alpha S_(i j)^"TF"$ and $gamma^(i j) mnabla_i mnabla_j alpha$, respectively. Put differently, this means that the failure of the BSSN variables to be Lie-transported across the leaves of the foliation is caused by global curvature, the presence of matter, and the choice of the foliation and coordinates on it.

== From Bar to Tilde
As alluded to before, the @eqBSSN[BSSN equations] are preliminary in the sense that they still contain expressions associated with the induced metric $gamma_(i j)$, which is not part of the BSSN variables. Although technically, one can always reconstruct $gamma_(i j)$ from $phi.alt$ and $tilde(gamma)_(i j)$ as $gamma_(i j) = e^(-4phi.alt) tilde(gamma)_(i j)$, this is computationally inefficient and induces more memory usage and traffic. In this section, we address the offending terms and reexpress them in terms of the BSSN variables $phi.alt, tilde(gamma)_(i j),K,tilde(A)_(i j)$ and $tilde(Gamma)^i$. Concretely, the terms that appear which we need to reformulate are
$
  macron(R)_(i j), quad macron(R), quad mnabla_i mnabla_j alpha quad "and" quad Delta_gamma alpha = gamma^(i j) mnabla_i mnabla_j alpha.
$
We begin with the last of the four, as it is the simplest. Using the standard formula for the Laplacian, we get
$
  Delta_gamma alpha &= 1/sqrt(gamma) diff_i (sqrt(gamma) gamma^(i j) diff_j alpha) = 1/sqrt(e^(-12 phi.alt)) diff_i (sqrt(e^(-12 phi.alt)) e^(4phi.alt) tilde(gamma)^(i j) diff_j alpha)\
  &= e^(6phi.alt) diff_i (e^(-2phi.alt) tilde(gamma)^(i j) diff_j alpha) = e^(4phi.alt)(diff_i (tilde(gamma)^(i j) diff_j alpha) - 2 tilde(gamma)^(i j) diff_i phi.alt diff_j alpha )\
  &= e^(4phi.alt)( Delta_tilde(gamma) alpha - 2 tilde(gamma)^(i j) diff_i phi.alt diff_j alpha).
$<eqConformalLaplacian>
In the last step, we made use of the fact that $tilde(gamma) = 1$. 

Next up, we consider the Hessian of $alpha$. For this, we recall the @eqConformalConnection[relationship], which allows us to infer
$
  mnabla_i mnabla_j alpha &= mnabla_i diff_j alpha = tnabla_i diff_j alpha + tensor(C,+k,-j i) diff_k alpha\
  &= tnabla_i tnabla_j alpha + 2 (tensor(delta,+k,-i) diff_j phi.alt + tensor(delta,+k,-j) diff_i phi.alt - tilde(gamma)_(i j) tilde(gamma)^(k ell) diff_ell phi.alt) diff_k alpha\
  &= tnabla_i tnabla_j alpha + 4 diff_(\(i) phi.alt diff_(j\)) alpha - 2 tilde(gamma)_(i j) tilde(gamma)^(k ell) diff_k phi.alt diff_ell alpha.
$
Note that since $tilde(gamma)^(i j) tilde(gamma)_(i j) =3$, this is consistent with @eqConformalLaplacian.

Now turning to $macron(R)_(i j)$, we can simply employ our @eqConformalRicciTensor[result] for the conformal Ricci tensor, which for $d=3$ implies
$
  macron(R)_(i j) = tilde(R)_(i j) + 2 (tnabla_i tnabla_j phi.alt + tilde(gamma)_(i j) tilde(gamma)^(k ell) tnabla_k tnabla_ell phi.alt) + 4 (diff_i phi.alt diff_j phi.alt - tilde(gamma)_(i j) tilde(gamma)^(k ell) diff_k phi.alt diff_ell phi.alt).
$ 
Here, according to @eqRicciTensorNice, the conformal Ricci tensor $tilde(R)_(i j)$ is given by 
$
  tilde(R)_(i j) &= -1/2 tilde(gamma)^(k ell) diff_k diff_ell tilde(gamma)_(i j) + tilde(gamma)_(k\(i) diff_(j\)) tilde(Gamma)^k \
&quad+ thin tilde(Gamma)^k tilde(Gamma)_((i j) k) + tilde(Gamma)_(k ell i) tensor( tilde(Gamma),+k ell, -j) + 2  tilde(Gamma)_(k ell \(i) tensor( tilde(Gamma),-j\),+k ell).
$
Here it becomes apparent why introducing $tilde(Gamma)^i$ as an additional dynamical variable is a good idea: since it only appears with first derivatives, the principal part of the conformal Ricci tensor is the elliptic operator
$
  tilde(gamma)^(k ell) diff_k diff_ell tilde(gamma)_(i j)
$
which is much more well-behaved numerically than the second-order terms introduced by $diff_i tilde(Gamma)^j$ would be. The introduction of $tilde(Gamma)^i$ turns exactly these terms into mere first-order contributions, so that they do not mess up the principal symbol.

Lastly, for the Ricci scalar, we may use @eqConformalRicciScalar, which for $d=3$ reads
$
  R = e^(4phi.alt)(tilde(R) + 8 tilde(gamma)^(i j) (tnabla_i tnabla_j phi.alt - diff_i phi.alt diff_j phi.alt)).
$
For quicker reference, we provide all these results again below:
$
  mnabla_i mnabla_j alpha &= tnabla_i tnabla_j alpha + 4 diff_(\(i) phi.alt diff_(j\)) alpha - 2 tilde(gamma)_(i j) tilde(gamma)^(k ell) diff_k phi.alt diff_ell alpha,\ \ \
  gamma^(i j) mnabla_i mnabla_j alpha &= e^(4phi.alt)( Delta_tilde(gamma) alpha - 2 tilde(gamma)^(i j) diff_i phi.alt diff_j alpha),\ \ \
  macron(R)_(i j) &= tilde(R)_(i j) + 2 (tnabla_i tnabla_j phi.alt + tilde(gamma)_(i j) tilde(gamma)^(k ell) tnabla_k tnabla_ell phi.alt) + 4 (diff_i phi.alt diff_j phi.alt - tilde(gamma)_(i j) tilde(gamma)^(k ell) diff_k phi.alt diff_ell phi.alt)\ \
  & quad "with"quad tilde(R)_(i j) = -1/2 tilde(gamma)^(k ell) diff_k diff_ell tilde(gamma)_(i j) + tilde(gamma)_(k\(i) diff_(j\)) tilde(Gamma)^k \ 
& #h(5.65em) quad+ thin tilde(Gamma)^k tilde(Gamma)_((i j) k) + tilde(Gamma)_(k ell i) tensor( tilde(Gamma),+k ell, -j) + 2  tilde(Gamma)_(k ell \(i) tensor( tilde(Gamma),-j\),+k ell),\ \ \
  R &= e^(4phi.alt)(tilde(R) + 8 tilde(gamma)^(i j) (tnabla_i tnabla_j phi.alt - diff_i phi.alt diff_j phi.alt)).
$

= Initial Data