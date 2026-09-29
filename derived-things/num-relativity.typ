#import "template.typ": *
#import "macros.typ": *
#import "@preview/xarrow:0.3.1": xarrow
#import "@preview/fletcher:0.4.5" as fletcher: diagram, node, edge
#import "@preview/cetz:0.5.1": canvas, draw
#import "@preview/cetz-plot:0.1.4": plot



= Numerical Relativity
== Conventions
In the following, the spacetime metric $g$ always has signature $(-+++)$. The induced metric $gamma$ and conformal metric $tilde(gamma)$ are Riemannian, that is, of signature $(+++)$. The Riemann tensor has the sign convention
$
  R(X,Y)Z= nabla_X nabla_Y Z - nabla_Y nabla_X Z - nabla_[X,Y] Z.
$
We use the Einstein equations in the form
$
  R_(mu nu)- 1/2 R g_(mu nu) + Lambda g_(mu nu) = 8 pi T_(mu nu),
$
that is, with $c= G = 1$.
== The $3+1$-Formalism
=== Foliations and Projectors
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
=== Extrinsic Curvature
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

=== Induced Connection and Intrinsic Curvature
==== Induced Connection
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
==== Intrinsic Curvature
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

=== Gauss-Codazzi-Ricci Equations
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
=== The Einstein Equations in the $3+1$-Formalism
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
    cal(H) &:= 2 cal(E)_(mu nu)n^mu n^nu = macron(R) + K^2 - K_(mu nu) K^(mu nu) - 2Lambda - 16 pi rho = 0,\ \
    cal(M)_mu &:= tensor(P,+lambda,-mu) n^nu cal(E)_(lambda nu) = mnabla_mu K - mnabla_nu tensor(K,+nu,-mu) + 8 pi j_mu = 0,\ \
    fL_n gamma_(mu nu) &= -2 K_(mu nu),\ \
    fL_n K_(mu nu) &= K K_(mu nu) - 2 K_(mu lambda) tensor(K,+lambda,-nu) + macron(R)_(mu nu) - 1/alpha mnabla_mu mnabla_nu alpha\ &wide - Lambda gamma_(mu nu) - 8pi ((S_(mu nu) - 1/2 S gamma_(mu nu)) + 1/2 rho gamma_(mu nu)).
  $<eqEvSys>]
]
The first two equations are derived by taking the normal-normal and normal-tangential projections of the original Einstein equations $cal(E) = G_(mu nu) + Lambda g_(mu nu) - 8pi T_(mu nu)=0$ and subsequently identifying the definitions of $rho$ and $j$ as well as applying the Gauss and Codazzi equations, respectively. These two equations involve no time derivatives, and are hence a set of four constraints. 

The remaining two equations are first order in time (by the presence of the $fL_n$ normal derivatives), and hence dynamical. The first of the two is simply a consequence of the definition of the extrinsic curvature. The latter is the result of projecting the trace-reversed Einstein equations onto $T Sigma$ by contracting with $P$ on both indices, and subsequently applying the first equation in @eqProjectionsRicci[] to re-express the projected Ricci tensor $(P R)_(mu nu)$.

=== Adapted Coordinates
==== The Metric and Bases
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

==== Extrinsic and Intrinsic Curvature
#lemma[
  For any vector field $N$ that is normal to $T Sigma$ (not necessarily the unit-normal $n$), and a contravariant tangential tensor $T in Gamma(T^((0,s))Sigma)$, we have
  $
    fL_(f N) T = f fL_N T.
  $<eqTensorialityNormalTangentialLieDeriv>
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
==== The Einstein Equations in $3+1$-Adapted Coordinates
The @eqTimeDerivGammaK[equations] directly imply that the Einstein equations in adapted coordinates reads
$
  fH &= macron(R) + K^2 - K_(i j) K^(i j) - 2Lambda - 16pi rho = 0, \ \ \
  fM_i &= mnabla_i K - mnabla_j tensor(K,+j,-i) + 8pi j_i =0, \ \ \
  diff_t gamma_(i j) &= beta^k diff_k gamma_(i j) + 2 gamma_(k \(i) diff_(j\)) beta^k - 2 alpha K_(i j), \ \ \
  diff_t K_(i j)&= beta^k diff_k K_(i j) + 2K_(k \(i) diff_(j\)) beta^k - mnabla_i mnabla_j alpha + alpha(K K_(i j) - 2 K_(i k) tensor(K,+k,-j) + macron(R)_(i j))\
  &quad - alpha Lambda gamma_(i j) - 8pi alpha (S_(i j) - 1/2 (S-rho) gamma_(i j))
$<eqADM>
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
== Metric Derivative and Curvature Identities
=== Jacobi Formula
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
=== Connection and Curvature of Conformal Metrics
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
$<eqConformalConnectionRelation>
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
=== Bochner Formula and the Ricci Tensor
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
== BSSNOK
In this section, we take what we have derived so far and build atop it the BSSNOK formalism. This involves

- defining the BSSNOK variables in terms of the ADM variables $gamma_(i j)$ and $K_(i j)$, scaling out the conformal factor to isolate the conformal metric $tilde(gamma)_(i j)$ and trace-free extrinsic curvature $tilde(A)_(i j)$;

- deriving their first time derivative expressions from the known expressions for $diff_t gamma_(i j)$ and $diff_t K_(i j)$ via the ADM evolution equations, as well as re-expressing the Hamiltonian and momentum constraints in terms of the new conformal variables;

- and finally, relating the physical induced connections and curvatures to their conformal counterparts using conformal transformations and the auxiliary connection quantities.

In doing so, we will use many of the identities from the previous section, and establish a number of additional ones to handle the conformal rescaling. Unfortunately, these algebraic definitions are rather lengthy and messy, but they ultimately lead to a strongly hyperbolic formulation of general relativity well-suited to be implemented in high-performance computing (HPC) code. 

=== BSSNOK Reparametrisation and Evolution Equations
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

#proposition(name: "BSSNOK Evolution Equations")[
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

=== From Bar to Tilde
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
  macron(R) &= e^(4phi.alt)(tilde(R) + 8 tilde(gamma)^(i j) (tnabla_i tnabla_j phi.alt - diff_i phi.alt diff_j phi.alt)).
$
=== #text(fill: red)[Gauge Dynamics]
== Initial Data
In this section, we discuss different approaches used to decompose and solve the constraint equations
#bottom-number[$
  fH &= macron(R) + K^2 - K_(i j) K^(i j) - 2Lambda - 16pi rho = 0, \ \ \
  fM_i &= mnabla_i K - mnabla_k tensor(K,+k,-i) + 8pi j_i =0.
$<eqConstraints>]
These are 4 equations for the total of 12 functions in $gamma_(i j)$ and $K_(i j)$, meaning that there is a significant amount of freedom in initial data. To separate the freely specifiable data from data fixed by these equations, we perform different kinds of decompositions which have been developed historically to approach the problem of computing physically meaningful initial data.

=== York-Lichnerowicz Conformal Transverse-Traceless (CTT) Decomposition
The York-Lichnerowicz conformal traceless split employs the variables
$
  gamma_(i j) &= psi^4 hat(gamma)_(i j) quad<=>quad gamma^(i j) = psi^(-4) hat(gamma)^(i j),\
  K_(i j) &= A_(i j) + 1/3 gamma_(i j) K,\
  A_(i j) &= psi^(-2) hat(A)_(i j) quad<=>quad A^(i j) = psi^(-10) hat(A)^(i j).
$
Concretely, this means that $hat(gamma)_(i j)$ is a conformal metric,
$
  hat(gamma)_(i j) = psi^(-4) gamma_(i j),
$
for which we do _not_ necessarily require $det hat(gamma) = 1$, $A_(i j)$ is the traceless extrinsic curvature
$
  A_(i j) = K_(i j) - 1/3 gamma_(i j) K,
$
and $hat(A)_(i j) = psi^2 A_(i j)$ its conformal rescaling.

We further denote by $hnabla$ the Levi-Civita connection associated with $hat(gamma)_(i j)$, whose coefficients---according to @eqConformalConnectionRelation with $psi = e^(-phi.alt) <=> phi.alt = -log psi$---are related to those of $mnabla$ by
$
  tensor(macron(Gamma),+k,-i j) = tensor(hat(Gamma),+k,-i j) + underbrace(2/psi (tensor(delta,+k,-i) diff_j psi + tensor(delta,+k,-j) diff_i psi - hat(gamma)_(i j) hat(gamma)^(k ell) diff_ell psi),=: tensor(C,+k,-i j)).
$
Besides the connection $mnabla$ needing to be re-expressed in terms of the variables $(psi,hat(gamma)_(i j), K, hat(A)_(i j))$, we also need to relate the Ricci scalar appearing in the Hamiltonian constraint to quantities associated to $hat(gamma)_(i j)$ and $psi$. In $d=3$, @eqConformalRicciScalar tells us that
$
  macron(R) = e^(4 phi.alt) (hat(R) + 8 hat(gamma)^(i j) hnabla_i hnabla_j phi.alt - 8 hat(gamma)^(i j) hnabla_i phi.alt hnabla_j phi.alt).
$
Again using $phi.alt = -log psi$, we expand 
$
  hat(gamma)^(i j) hnabla_i hnabla_j phi.alt &= -hat(gamma)^(i j)hnabla_i hnabla_j log psi = -hat(gamma)^(i j)hnabla_i (psi^(-1) hnabla_j psi)\ 
  &= -psi^(-1) hat(gamma)^(i j) hnabla_i hnabla_j psi + psi^(-2) hat(gamma)^(i j) hnabla_i psi hnabla_j psi\
  hat(gamma)^(i j)hnabla_i phi.alt hnabla_j phi.alt &= psi^(-2) hat(gamma)^(i j) hnabla_i psi hnabla_j psi,
$
whence the Ricci scalar becomes
$
  macron(R) = psi^(-4) hat(R) - 8 psi^(-5) hat(gamma)^(i j)hnabla_i hnabla_j psi.
$
With this in hand, we can write the Hamiltonian constraint as
$
  hat(fH) = 8 hat(gamma)^(i j) hnabla_i hnabla_j psi - hat(R) psi - 2/3 K^2 psi^5 + hat(A)_(i j) hat(A)^(i j) psi^(-7) + (2 Lambda + 16 pi rho) psi^5 = 0.
$
For the momentum constraint, we need to do some more work. We rewrite it slightly to become
$
  fM^j = mnabla_i (gamma^(i j) K - K^(i j)) + 8 pi j^j = 0,
$
and rewrite the first term separately as follows:
$
  mnabla_i (gamma^(i j) K - K^(i j)) = mnabla_i (gamma^(i j) K - A^(i j) - 1/3 gamma^(i j) K) = 2/3 gamma^(i j) mnabla_i K- mnabla_i A^(i j).
$
For the divergence term, we get
$
  mnabla_i A^(i j) &= mnabla_i (psi^(-10) hat(A)^(i j)) = hnabla_i (psi^(-10) hat(A)^(i j)) + psi^(-10) tensor(C,+i,-k i) hat(A)^(k j) + psi^(-10) tensor(C,+j,-k i) hat(A)^(i k)\
  &= psi^(-10) hnabla_i hat(A)^(i j) - 10 psi^(-11) hat(A)^(i j) diff_i psi + 6 psi^(-11) hat(A)^(k j) diff_k psi + 4 psi^(-11) hat(A)^(j k) diff_k psi\
  &= psi^(-10) hnabla_i hat(A)^(i j). 
$
Hence, inserting back into the momentum constraint, we get
$
  -psi^10 fM^i = hnabla_k hat(A)^(k i) -2/3 psi^6 hat(gamma)^(i k) hnabla_k K  - 8 pi psi^(10) j^i = 0.
$
Although we could leave it at this, we can make the following important observation which will allow us to split the degrees of freedom even further: the traceless conformal extrinsic curvature $hat(A)^(i j)$ appears in the momentum constraint only through its divergence $hnabla_j hat(A)^(i j)$. This means that any divergence-free components remain completely unconstrained and can be chosen arbitrarily. According to the following lemma, the remaining longitudinal component of $hat(A)^(i j)$ can be written in terms of a vector potential:

#lemma[
  Let $hat(A)^(i j)$ be a symmetric traceless tensor. Then, there exists a symmetric, traceless and transverse tensor $Q^(i j)$, i.e.
  $
    hnabla_j Q^(i j) = 0, quad tensor(Q,+i,-i) = 0,
  $
  and a vector field $X^i$ such that
  $
    hat(A)^(i j) = Q^(i j) + (LL X)^(i j) = Q^(i j) + hnabla^i X^j + hnabla^j X^i - 2/3 hat(gamma)^(i j) hnabla_k X^k.
  $
  The operator $LL$ is called the _Killing operator_, and $tensor((LL X),+i,-i) = 0$ by construction.
]

We insert this into the divergence term that appears in the momentum constraint to obtain
$
  hnabla_k hat(A)^(k i) &= hnabla_k hnabla^k X^i + hnabla_k hnabla^i X^k - 2/3 hnabla^i hnabla_k X^k \
  &= hnabla_k hnabla^k X^i +1/3 hnabla^i hnabla_k X^k +  tensor(hat(R),+i,-k) X^k.
$
Thus, the Hamiltonian and momentum constraints in this decomposition read
#bottom-number[$
  hat(fH) &= 8 hat(gamma)^(i j) hnabla_i hnabla_j psi - hat(R) psi - 2/3 K^2 psi^5 + hat(A)_(i j) hat(A)^(i j) psi^(-7) + (2 Lambda + 16 pi rho) psi^5 = 0,\ \ \
  hat(fM)^i &= hnabla_k hnabla^k X^i + 1/3 hnabla^i hnabla_k X^k + tensor(hat(R),+i,-k) X^k - 2/3 psi^6 hat(gamma)^(i k) diff_k K - 8pi psi^10 j^i = 0. 
$<eqCTTConstraints>]
The implications of this decomposition are very convenient; the data is now clearly split into eight freely specifiable functions---5 from the background conformal metric $hat(gamma)_(i j)$, one from the mean curvature $K$, and two from the symmetric transverse-traceless extrinsic curvature component $Q^(i j)$. The remaining four variables---the scalar $psi$ and the three components of $X^i$---are then fixed by the constraint equations above. Further, the momentum constraint is linear in the vector potential $X^i$, which will enable us to find analytical solutions in certain cases.
==== Brill-Lindquist Data
The most convenient way to solve an equation is always to specialise to the most simplified version of the problem. Here, this is achieved by assuming vacuum ($rho=j^i=0$), no cosmological constant ($Lambda = 0$), and time-symmetry ($K_(i j) = 0$, so $K = Q^(i j) = X^i = 0$), as well as a flat conformal background, $hat(gamma)_(i j) = delta_(i j)$. These assumptions collapse the @eqCTTConstraints[constraint equations] down to just a single scalar equation,
$
  Delta psi = 0,
$<eqBrillLindquistEqn>
where $Delta = delta^(i j) diff_i diff_j$ is the flat-space Laplacian. If we require solutions to be asymptotically flat and smooth, we must have $psi -> 1$ as $r->infty$ and $||psi||_infty < infty$, which by Liouville's theorem forces $psi equiv 1$ everywhere. However, if we allow for isolated singular points, then solutions of the form
$
  psi(vx) = 1 + C/r
$
with $vx = (x^i)$ and $r = sqrt(delta_(i j) x^i x^j)$ are allowed as well. Given that in isotropic coordinates, the Schwarzschild metric reads
$
  g_"SS" = - ((1-M/(2r))/(1+M/(2r)))^2 dt^2 + (1+M/(2r))^4 delta_(i j) dx^i dx^j,
$
we identify 
$
  psi(vx) = 1+ M/(2r),
$
whence $C = M\/2$ is related to the mass of the black hole. By linearity of  @eqBrillLindquistEqn, we can superimpose multiple such solutions to obtain multi-black hole initial data,
$
  psi(vx) = 1 + sum_(i = 1)^n M_((a))/(2|vx-vx_((a))|)
$
where $M_((a))$ is the mass and $vx_((a))$ the initial position of the $i$-th black hole. This form of initial data is known as _Brill-Lindquist data_. 

#remark[
  There are a couple of observations to be made here.
]
+ As can be seen in the Schwarzschild example, the singularities of $psi$ do not actually correspond to physical singularities inside the black hole, but rather represent the coordinate singularities at spatial infinity of the additional asymptotically flat ends that are introduced by the presence of black holes. 

+ The time-symmetry assumption $K_(i j) = 0$ implies that the mean curvature vanishes, $K=0$. This is the condition for a submanifold to be maximal---in essence, it ensures that the initial data slice has a maximal 3-volume within the four-dimensional spacetime. 

+ The masses $M_((a))$ are the _bare_ masses of the black holes, not the physical ADM masses of their associated asymptotically flat ends. Near the puncture $vx_((a))$---that is, near spatial infinity of the flat end, where $vx->vx_((a))$---the conformal factor $psi$ behaves as
  $
    psi(vx) sim 1 + M_((a))/(2|vx-vx_((a))|) + (sum_(j != i) M_((j))/(2|vx-vx_((j))|)).
  $
  Although the leading "$1+M_((a))\/2r$" ensures that $M_((a))$ appears as a term in the ADM mass of the $i$-th end, the remaining terms in parentheses provide an additional, regular contribution to the surface integral that is best interpreted as a binding energy. 
  #text(fill:red)[Explicitly, it can be computed that
  $
    M_("ADM",(i)) = M_((a))(1 + sum_(j!=i) M_((j))/(2|vx_((a))-vx_((j))|))
  $
  I did not check this explicitly yet.] Nonetheless, the total ADM mass of the "main end" at $r->infty$ is 
  $
    M_("ADM",infty) = sum_(i=1)^n M_((a)),
  $
  as is easily verified by the fact that asymptotically, $|vx-vx_((a))| sim r$ and hence,
  $
    psi(vx) sim 1 + (sum_(i=1)^n M_((a)))/(2r) quad "as" quad r->infty.
  $

+ Although Brill-Lindquist can accommodate multi-black hole initial data, it does so in a very limited way. The black holes are initially stationary (although of course, they do not remain stationary during evolution), and have no spin. In astrophysical scenarios, this is typically rather uninteresting; we usually want to simulate black hole mergers where the black holes have both linear and angular momentum. For this reason, in the next section, we will loosen our simplifying assumptions, with the goal of enabling us to place black holes with arbitrary mass, momentum and spin into our simulation domain.
==== Bowen-York Initial Data and the Puncture Method
To improve on the lack of initial linear and angular momentum of Brill-Lindquist initial data, in this section, we make our assumptions more general. Concretely, we drop time symmetry $K_(i j) = 0$, and instead only demand maximal slicing ($K=0$) and that the transverse part vanishes, $Q^(i j)=0$. To arrive at Bowen-York initial data, we keep the flat conformal background assumption $hat(gamma)_(i j) = delta_(i j)$ in place, and continue to work in vacuum with a vanishing cosmological constant. 

These assumptions turn the constraint equations into
#bottom-number[$
  hat(fH) &= 8 Delta psi + hat(A)_(i j) hat(A)^(i j) psi^(-7) = 0,\
  hat(fM)^i &= Delta X^i + 1/3 diff^i diff_k X^k = 0,
$<eqBowenYorkConstraints>]
where $Delta = delta^(i j)diff_i diff_j$ is the flat space Laplacian, indices are raised with $hat(gamma)^(i j) = delta^(i j)$, and $hat(A)^(i j)$ is computed as
$
  hat(A)^(i j) = (LL X)^(i j) = diff^i X^j + diff^j X^i - 2/3 delta^(i j) diff_k X^k.
$ 
We can see that under these assumptions, only $X^i$ appears in the momentum constraints, decoupling it fully from $psi$. Further, the momentum constraint is now a linear equation for $X^i$---which we can actually solve analytically---and which allows for solutions to be superimposed. Once $X^i$ is determined, we can evaluate $hat(A)^(i j)$, leaving us to solve the nonlinear Hamiltonian constraint for $psi$. Although this has to be approached numerically, we will consider an ansatz that splits off the analytical Brill-Lindquist terms, making it so that we only have to solve for a regular function---instead of the full $psi$ which typically contains singular points.

Let us begin by discussing the two most important types of solutions to the momentum constraint equation, starting with angular momentum.

#proposition(name: "Bowen-York Data: Angular Momentum")[
  The @eqBowenYorkConstraints[momentum constraint equation] is solved by 
  $
    X^i_"ang" = epsilon^(i j k) x_j/r^3 J_k,
  $
  which has the associated traceless conformal extrinsic curvature tensor
  $
    hat(A)^(i j)_"ang" = -6/r^5 x^(\(i)epsilon^(j\)k ell)x_k J_ell.
  $
]
#text(fill:red)[One can verify that for asymptotically flat initial data, where $psi->1$ for $r->infty$, $J_"ADM"^i = J^i$ for the AF end produced by the puncture.]
#text(fill:red)[
#proposition(name: "Bowen-York Data: Linear Momentum")[
  The @eqBowenYorkConstraints[momentum constraint equation] is also solved by 
  $
    X^i_"lin" = -1/(4r)(7P^i+(x^i x_k)/r^2 P^k)
  $
  for a constant vector $P^i$. The associated traceless conformal extrinsic curvature tensor is
  $
    hat(A)^(i j)_"lin" = 3/(2r^3) [2 P^(\(i) x^(j\)) + ((x^i x^j)/r^2 - delta^(i j))x_m P^m].
  $
]]
#text(fill:red)[For asymptotically flat initial data, where $psi->1$ for $r->infty$, $P^i_"ADM" = P^i$ for the AF end produced by the puncture, so that $X^i_"lin"$ has the interpretation of the linear momentum of a black hole.]

We now have everything we need to close the Bowen-York initial data story. Given positions $vx_((a))$, masses $M_((a))$, linear momenta $vP_((a))$ and angular momenta $vJ_((a))$, we set up the total $hat(A)^(i j)$ as
$
  hat(A)^(i j) (vx) = sum_(a = 1)^n hat(A)^(i j)_"lin" (vx-vx_((a)); vP_((a))) + sum_(a=1)^n hat(A)^(i j)_"ang" (vx-vx_((a)); vJ_((a))).
$
It now remains to solve the Hamiltonian constraint for $psi$. To assign the masses to the punctures, it is sufficient---mathematically speaking---to impose the boundary conditions that
$
  psi(vx) sim M_((a))/(2|vx-vx_((a))|) quad "for" quad vx->vx_((a)), quad psi(vx)->1 quad "for" quad r->infty
$
However, this is hard to enforce numerically, and will cause the function $psi$ which is being solved for to have singularities that are difficult to resolve on a grid. Given that we already know the Brill-Lindquist solution, we rather build on that; we take the ansatz
$
  psi(vx) = u(vx) + underbrace(1 + sum_(a=1)^n M_((a))/(2|vx-vx_((a))|),=:psi_0 (vx)),
$
and solve for the (now regular) function $u(vx)$ instead. Since away from the punctures at $vx_((a))$, the Laplacian acting on the $sim 1\/r$ terms vanishes, the equation for $u$ reads
$
  Delta u + 1/8 hat(A)_(i j) hat(A)^(i j) (u + psi_0)^(-7) = 0.
$
We have thus now reduced the problem of solving the rather complicated constraint equations to an almost fully analytical solution, with as final step solving a non-linear scalar equation for a regular function. This final reduction step is often termed the _puncture method_, due to Brandt and Brügmann. 

The solutions that result are black holes in specifiable positions, with tunable mass as well as linear and angular momentum. However, the initial data produced this way is not perfect, in the sense that it contains more than just black holes. This is evident from the following line of reasoning: It can be shown that the Kerr spacetime---i.e. the spacetime containing a singular rotating black hole---does not allow for any conformally flat spatial hypersurfaces. Our initial data, however, is able to produce a hypersurface that contains a black hole with spin, starting out from the assumption that it is a conformally flat slice ($hat(gamma)_(i j) = delta_(i j)$). This means that, besides a Kerr black hole, such a slice must contain additional gravitational radiation, as otherwise we would contradict the theorem that Kerr allows no conformally flat spatial slices. For slow spins, this radiation is typically of no concern, as it propagates away quickly and hence only contaminates the first part of the gravitational wave signal. 
=== #text(fill: red)[(Extended) Conformal Thin Sandwich ((X)CTS)] 
== An Initial Overview of Z4
=== Time Evolution of Violated Constraints in Standard GR
It is rather simple to show that mathematically speaking, if the constraints are satisfied initially, they are satisfied for the entirety of the evolution---this is a direct consequence of $nabla^mu G_(mu nu) = 0$. Unfortunately, we cannot perfectly satisfy $cal(H) = cal(M)_i = 0$ when computing initial data numerically, and machine as well as discretisation error introduce additional deviations, so that at best, $cal(H) approx cal(M)_i approx 0$. Motivated by this, let us see how nonzero values of $cal(H)$ and $cal(M)_i$ evolve under evolution.

Since we are only interested in a mathematical result for now, and not a numerically stable implementation, we can make our lives easier by starting from the Einstein equations written as $cal(E)_(mu nu) = G_(mu nu) + Lambda g_(mu nu) - 8pi T_(mu nu) = 0$, without needing to decompose or invoke adapted coordinates just yet. We recall that the constraints are defined as
$
  cal(H) = 2 cal(E)_(mu nu) n^mu n^nu , quad cal(M)_mu = tensor(P,+lambda,-mu)n^nu cal(E)_(lambda nu).
$
By the contracted Bianchi identity $nabla^mu G_(mu nu) = 0$, the metric compatibility of the connection, $nabla_lambda g_(mu nu) = 0$, and conservation of energy-momentum $nabla^mu T_(mu nu) = 0$, we have
$
  nabla^mu cal(E)_(mu nu) = 0.
$
Besides this identity, we will also need to assume that the @eqEvSys[evolution equations in] hold. We recall that these are derived from the tangential-tangential projection of the _trace-reversed_ Einstein equations, that is, from
$
  tensor(P,+lambda,-mu) tensor(P,+rho,-nu) (cal(E)_(lambda rho) - 1/2 g_(lambda rho) cal(E)) = 0.
$<eqTraceRevEEProjected>
Since in the derivations that follow, the projection $tensor(P,+lambda,-mu) tensor(P,+rho,-nu) cal(E)_(lambda rho)$ will appear, let us first examine what this evaluates to assuming that the above holds. Taking the trace with respect to $g^(mu nu)$ gives us
$
  0 = P^(lambda rho) cal(E)_(lambda rho) - 3/2 cal(E) = underbrace(g^(lambda rho) cal(E)_(lambda rho),=cal(E)) + underbrace(n^lambda n^rho cal(E)_(lambda rho),=1/2 cal(H)) - 3/2 cal(E) = 1/2 cal(H) - 1/2 cal(E).
$
Thus, $cal(E) = cal(H)$, and by @eqTraceRevEEProjected,
$
  tensor(P,+lambda,-mu) tensor(P,+rho,-nu) cal(E)_(lambda rho) = 1/2 gamma_(mu nu) cal(H).
$

To derive an evolution equation for $cal(H)$, we project $nabla^mu cal(E)_(mu nu)$ onto $n^nu$, to obtain
$
  0 &= n^nu nabla^mu cal(E)_(mu nu) = nabla^mu (cal(E)_(mu nu) n^nu) - cal(E)_(mu nu) nabla^mu n^nu\
&= nabla^mu lr((-n_mu underbrace(n^lambda cal(E)_(lambda nu) n^nu,=1/2 cal(H)) + underbrace(tensor(P,-mu,+lambda) cal(E)_(lambda nu) n^nu,=cal(M)_mu) ),size:#45%) + cal(E)_(mu nu) (K^(mu nu) + n^mu a^nu)\
&= -1/2 nabla_mu (n^mu cal(H)) + nabla^mu cal(M)_mu + underbrace(cal(E)_(mu nu) tensor(P,+mu,-lambda)tensor(P,+nu,-rho),=1/2 gamma_(lambda rho) cal(H)) K^(nu rho) + fM_mu a^mu \
&= -1/2 n^mu nabla_mu cal(H) - 1/2 underbrace((nabla_mu n^mu),=-K) cal(H) + nabla^nu (tensor(P,+mu,-nu) cal(M)_mu) + 1/2 K cal(H) + cal(M)_mu a^mu \
&= -1/2 n^mu nabla_mu cal(H) + K cal(H) + underbrace(tensor(P,+mu,-nu) nabla^nu cal(M)_mu,=mnabla^mu cal(M)_mu) + underbrace(cal(M)_mu  nabla^nu (n^mu n_nu),=cal(M)_mu a^mu) + cal(M)_mu a^mu \
&= -1/2 n^mu nabla_mu cal(H) + K cal(H) + mnabla^mu cal(M)_mu + 2 cal(M)_mu a^mu 
$
Solving for $n^mu nabla_mu cal(H)$, we find
$
  n^mu nabla_mu cal(H) = 2 K cal(H) + 2 mnabla^mu cal(M)_mu + 4 cal(M)_mu a^mu.
$
Switching to adapted coordinates and making use of 
$
  n = 1/alpha (diff_t - beta^i diff_i), quad a_mu = 1/alpha mnabla_mu alpha
$
we finally arrive at
$
  diff_t cal(H) = beta^i diff_i cal(H) + 2 alpha K cal(H) + 2 alpha mnabla^i cal(M)_i + 4 cal(M)_i mnabla^i alpha.
$

We can proceed similarly for the momentum constraint, this time projecting onto the foliation-tangent space:
$
  0 &= tensor(P,+lambda,-nu) (nabla^mu cal(E)_(mu lambda)) = nabla^mu (tensor(P,+lambda,-nu) cal(E)_(mu lambda)) - cal(E)_(mu lambda) nabla^mu tensor(P,+lambda,-nu)\
  &= nabla^mu lr((-n_mu underbrace(n^rho tensor(P,+lambda,-nu)cal(E)_(rho lambda),=cal(M)_nu) + underbrace(tensor(P,+rho,-mu)tensor(P,+lambda,-nu)cal(E)_(rho lambda),=1/2 gamma_(mu nu) cal(H))),size:#35%) - cal(E)_(mu lambda) nabla^mu (n^lambda n_nu)\
  &= - nabla^mu (n_mu cal(M)_nu) + 1/2 nabla^mu (gamma_(mu nu) cal(H)) + cal(E)_(mu lambda) (K^(mu lambda) + n^mu a^lambda) n_nu + cal(E)_(mu lambda) n^lambda (tensor(K,+mu,-nu) + n^mu a_nu)\
  &= -n^mu nabla_mu cal(M)_nu - cal(M)_nu underbrace(nabla^mu n_mu,=-K) + 1/2 mnabla_nu cal(H) + 1/2 cal(H) nabla^mu (n_mu n_nu) + 1/2 cal(H) K n_nu + cal(M)_mu a^mu  n_nu \
  & wide + tensor(K,+mu,-nu) cal(M)_mu + 1/2 cal(H) a_nu\
  &= -n^mu nabla_mu cal(M)_nu + K cal(M)_nu + 1/2 mnabla_nu cal(H) + 1/2 cal(H)a_nu + cal(M)_mu a^mu n_nu + tensor(K,+mu,-nu) cal(M)_mu + 1/2 cal(H) a_nu
$
We again want to find an expression in adapted coordinates, and then solve for $diff_t cal(M)_i$. Since $n_i = 0$, we can immediately drop the term proportional to $n_nu$ to find
$
  n^mu nabla_mu cal(M)_i = K cal(M)_i + tensor(K,+j,-i) cal(M)_j + 1/2 mnabla_i cal(H) + 1/alpha cal(H) mnabla_i alpha
$
Since $cal(M)_mu$ is a covector, we cannot simply replace the normal derivative on the left-hand side with $1/alpha (diff_t - beta^i diff_i)$ as before. Instead, let us expand the covariant derivative explicitly. For any spatial covector $X_mu$ with $n^mu X_mu = 0$, 
$
  fL_n X_nu &= n^mu nabla_mu X_nu + X_mu nabla_nu n^mu
$
implying that
$
  n^mu nabla_mu X_nu &= fL_n X_nu - X_mu nabla_nu n^mu\
  &= fL_n X_nu + X_mu (tensor(K,+mu,-nu) + n_nu a^mu)\
  &= fL_n X_nu + tensor(K,+mu,-nu) X_mu + n_nu X_mu a^mu.
$
Further using @eqTensorialityNormalTangentialLieDeriv, we can rewrite the spatial components of this expression as
$
  n^nu nabla_nu X_i = 1/alpha diff_t X_i - 1/alpha fL_beta X_i + tensor(K,+j,-i) X_j.
$<eqTangentCovectorNormalDerivative>
Inserting this result back into the original equation and expanding the Lie derivative along $beta$, we get
$
  diff_t cal(M)_i =  beta^j diff_j cal(M)_i  + cal(M)_j diff_i beta^j + alpha K cal(M)_i + 1/2 alpha mnabla_i cal(H) + cal(H) mnabla_i alpha
$ 
For quicker reference, let us summarise the results for both $diff_t cal(H)$ and $diff_t cal(M)_i$ below, making compact again the Lie derivatives along $beta$:
#bottom-number[$
  diff_t cal(H) &= cal(L)_beta cal(H) + 2 alpha K cal(H) + 2 alpha mnabla^i cal(M)_i + 4 cal(M)_i mnabla^i alpha,\
  diff_t cal(M)_i &= cal(L)_beta cal(M)_i  + alpha K cal(M)_i + 1/2 alpha mnabla_i cal(H) + cal(H) mnabla_i alpha.
$<eqConstraintEvolution>]
*Remarks:*
+ The simplest check we can carry out to verify whether this result is plausible is to see whether evolution preserves initially satisfied constraints. This is indeed the case, since all terms on the right-hand sides of both $diff_t cal(H)$ and $diff_t cal(M)_i$ are proportional to either $cal(H)$, $cal(M)_i$, or spatial derivatives thereof. Hence, if $cal(H)=cal(M)_i = 0$ initially, then the right-hand sides vanish initially. Since $cal(H) = cal(M)_i$ is itself a solution of the constraint-propagation equations, uniqueness of the corresponding initial-value problem implies that initially staisfied constraints remain satisfied throughout the evolution.

+ Although the @eqConstraintEvolution[evolution equations] look rather involved, all terms can be given a physical interpretation:
  - The terms $cal(L)_beta cal(H)$ and $cal(L)_beta cal(M)_i$ describe advection by the shift, separating coordinate motion from evolution along the normal direction. The terms proportional to $alpha K$ then modify the evolution according to the expansion or contraction of the normal congruence, decreasing the magnitude of the constraints when the congruence expands and increasing it when the congruence contracts. 

  - The terms $2 alpha mnabla^i cal(M)_i$ and $1/2 alpha mnabla_i cal(H)$ couple spatial variations of the momentum and Hamiltonian constraints, respectively. The remaining terms involving $mnabla_i alpha$ provide additional lower-order coupling when the lapse varies spatially.

=== Propagation of Constraints in GR/ADM
The time evolution of $cal(H)$ and $cal(M)_i$ can be used to derive second-order wave-like equations for $cal(H)$ and $cal(M)_i$, which reveal the propagation behaviour (or lack thereof) of constraint violations across the spacetime. In this section, we carry out such a derivation in a simplified scenario. Specifically, we make the simplifying assumption of unit lapse and vanishing shift, i.e.
$
  alpha = 1, quad beta^i = 0.
$
This choice is locally fully general: around any spacelike hypersurface, one can construct Gaussian normal coordinates for which the above holds, at least within a sufficiently small neighbourhood of the hypersurface. 

In this gauge---and replacing $diff_t$ with an overdot---@eqConstraintEvolution[the constraint evolution] breaks down to
$
  dot(cal(H)) &= 2 K cal(H) + 2 mnabla^i cal(M)_i,\
  dot(cal(M))_i &= K cal(M)_i + 1/2 mnabla_i cal(H).
$ 
We begin by taking a second time derivative of $cal(H)$. This leads to
$
  ddot(cal(H)) &= 2 dot(K) cal(H) + 2 K dot(cal(H)) + 2 diff_t (gamma^(i j) (diff_i cal(M)_j - tensor(macron(Gamma),+k,-i j) cal(M)_k))\
  &= 2 dot(K) cal(H) + 2 K dot(cal(H)) + 2dot(gamma)^(i j) mnabla_i cal(M)_j + 2 mnabla^i dot(cal(M))_i - 2 gamma^(i j) tensor(dot(macron(Gamma)),+k,-i j) cal(M)_k
$
For $dot(cal(H))$ and $dot(cal(M))_i$, we can simply insert the simplified evolution equations from above. The terms involving $dot(gamma)^(i j)$ and $tensor(dot(macron(Gamma)),+k,-i j)$ require slightly more attention. Starting with $dot(gamma)^(i j)$ and keeping in mind that $alpha = 1$, $beta^i = 0$, we use the @eqADM[ADM equation] for $diff_t gamma_(i j)$ to rewrite
$
  dot(gamma)^(i j) = - gamma^(i k) gamma^(j ell) underbrace(diff_t gamma_(k ell),=-2 K_(i j)) = 2 K^(i j).
$
For the time derivative of the connection coefficients, we observe that it is the (infinitesimal) difference of connections, and hence a tensor. In normal coordinates, it follows that
$
  tensor(dot(macron(Gamma)),+k,-i j) = 1/2 gamma^(k ell) (mnabla_i diff_t gamma_(j ell) + mnabla_j diff_t gamma_(i ell) - mnabla_ell diff_t gamma_(i j)).
$
Again using $diff_t gamma_(i j) = -2 K_(i j)$, we arrive at
$
  gamma^(i j) tensor(dot(macron(Gamma)),+k,-i j) &= - gamma^(i j) gamma^(k ell) (mnabla_i K_(j ell) + mnabla_j K_(i ell) - mnabla_ell K_(i j))\
  &= mnabla^k K - 2 mnabla_i K^(i k)
$
Inserting all this back into $ddot(cal(H))$, we arrive at
$
  ddot(cal(H)) &= 2 dot(K) cal(H) + 2 K (2 K cal(H) + 2 mnabla^i cal(M)_i) + 4 K^(i j) mnabla_i cal(M)_j + 2 mnabla^i (K cal(M)_i + 1/2 mnabla_i cal(H))\
  &wide - 2 (mnabla^k K - 2 mnabla_i K^(i k)) cal(M)_k\
  &= (2 dot(K) + 4 K^2) cal(H) + 4 (gamma^(i j) K + K^(i j)) mnabla_i cal(M)_j + cancelr(2 (mnabla^i K) cal(M)_i) + 2 K mnabla^i cal(M)_i + mnabla^i mnabla_i cal(H)\
  &wide - cancelr(2 (mnabla^k K) cal(M)) + 4 (mnabla_i K^(i k)) cal(M)_k\
  &=mnabla^i mnabla_i cal(H) + (2 dot(K) + 4K^2 )cal(H) + 6 K mnabla^i cal(M)_i + 4  mnabla_i (K^(i j) cal(M)_j) 
$
We can rearrange this into a wave equation for $cal(H)$,
$
  ddot(cal(H))-mnabla^i mnabla_i cal(H) = (2 dot(K) + 4K^2 )cal(H) + 6 K mnabla^i cal(M)_i + 4  mnabla_i (K^(i j) cal(M)_j).
$
We now carry out the same procedure for $ddot(cal(M))_i$, where we make use of the fact that $mnabla_i cal(H) = diff_i cal(H)$, since it is a scalar:
$
  ddot(cal(M))_i &= dot(K) cal(M)_i + K dot(cal(M))_i + 1/2 mnabla_i dot(cal(H))\
  &= dot(K) cal(M)_i + K(K cal(M)_i + 1/2 mnabla_i cal(H)) + 1/2 mnabla_i (2 K cal(H) + 2 mnabla^k cal(M)_k)\
  &= (dot(K) + K^2) cal(M)_i + 1/2 K mnabla_i cal(H) + mnabla_i (K cal(H)) + mnabla_i mnabla^k cal(M)_k\
$
To isolate the principal symbol, we move the second-order terms to the left. This leads us to---where for a complete summary, we repeat the equation for $cal(H)$ as well---the propagation equations
$
  ddot(cal(H))-mnabla^i mnabla_i cal(H) &= (2 dot(K) + 4K^2 )cal(H) + 6 K mnabla^i cal(M)_i + 4  mnabla_i (K^(i j) cal(M)_j),\
  ddot(cal(M))_i - mnabla_i mnabla^k cal(M)_k &= (dot(K) + K^2) cal(M)_i + 3/2 K mnabla_i cal(H) + (mnabla_i K) cal(H).
$<eqConstraintPropagationGR>
The equation for $cal(H)$ is a wave equation; the Hamiltonian constraint hence propagates hyperbolically. The lower-order terms on the right-hand side introduce coupling to $cal(M)_i$, as well as an effective mass term dependent upon the extrinsic curvature. 

The equation for the momentum constraint is similar, but has a fundamental difference. Instead of the spatial derivatives forming a Laplacian (i.e., the divergence of the gradient), instead, they constitute the gradient of the divergence. This implies that the momentum constraint does not propagate like a standard wave---let us examine its behaviour more precisely. We consider a spatial 3-covector $vX = (X_i)$ on flat spacetime, satisfying the equation
$
  ddot(X)_i - diff_i diff^k X_k = 0.
$
We employ the mode ansatz
$
  vX (t,vx) = vC e^(i (omega t - vk dot vx))
$
where $vk$ is the wave vector and $vC$ the polarisation vector. Inserting this ansatz, we find
$
  - omega^2 C_i + k_i k^k C_k = 0 quad <=> quad omega^2 vC = (vk dot vC) vk.
$
There are two distinct ways to satisfy this equation to obtain a non-trivial solution $vX$:
+ Either, the polarisation is longitudinal, $vC = C vk$ for some constant $C$, and $vk^2 = omega^2$,

+ or, the polarisation is transverse, meaning that $vk dot vC = 0$ and hence necessarily, $omega^2 = 0$.

This means that the question of whether the constraints $cal(H),cal(M)_i$ in GR propagate has an answer that is more subtle than a pure yes or no. We have just found that while the Hamiltonian constraint $cal(H)$, as well as the longitudinal modes of $cal(M)_i$---the contributions for which $mnabla^i cal(M)_i = 0$---propagate hyperbolically, the transverse modes of $cal(M)_i$ stay fixed in place. Hence, in numerical simulations, such transverse momentum constraint violations do not get propagated to the boundary, but instead remain where they are produced, and potentially grow to scales that make the simulation become unphysical and crash.
=== The Auxiliary Field $Z_mu$
To improve the stability of numerical simulations of GR, we would ideally like to modify the @eqConstraintPropagationGR[propagation system] such that all constraint-violation modes have hyperbolic propagation. However, the constraint propagation equations cannot be altered independently: they follow from the choice of evolution equations. Thus to change the constraint propagation, we must instead modify the evolution system for the dynamical variables and allow the constraint propagation to change as a consequence. In principle, such modifications could be chosen to produce a more desirable constraint subsystem, but constructing them indirectly is difficult, particularly because we would like the dynamics on the constraint manifold to remain exactly those of GR. 

In essence, since the constraints are defined as
$
  cal(H) = 2 n^mu n^nu cal(E)_(mu nu), wide cal(M)_mu = tensor(P,+lambda,-mu) n^nu cal(E)_(lambda nu),
$
they are nothing but a measure of how well certain projections of the full Einstein equations $cal(E)_(mu nu) = 0$ are satisfied. Suppose we introduce new degrees of freedom to the Einstein equations, which likewise measure the deviation from $cal(E)_(mu nu) = 0$, but whose equations of motion we can control directly. Such degrees of freedom would necessarily be related to the constraints, and could therefore provide a more direct means of controlling their propagation behaviour. 

Since we want to control four constraints, a reasonable guess is to introduce four degrees of freedom---arranged into a covector we call $Z_mu$. The motivation to pick a covector is that $cal(H)$ and $cal(M)_mu$ are really just further projections of the partial projection $cal(E)_(mu nu) n^nu$, which is a covector. 

Next, we should think about how to modify the Einstein equations, $cal(E)_(mu nu) = 0$, to incorporate additional terms involving $Z_mu$. To obtain the simplest possible equations of motion, we should only add terms linear in $Z_mu$---this has the additional effect that $Z_mu = 0$ recovers the original Einstein equations. Further, whatever we add should be symmetric in $mu nu$, and should contain derivatives of $Z_mu$ so that it is a dynamical field, rather than merely algebraically fixed. These requirements are satisfied by the terms $nabla_(\(mu)Z_(nu\))$ and $g_(mu nu) nabla^lambda Z_lambda$. Although there are more terms we could construct from the available tensors, such as $R_(mu nu) nabla^lambda Z_lambda$, the incorporation of extra curvature couplings will only make the dynamics of $Z_mu$ more involved, so we do not add any such terms. Using the two terms we have identified, we modify the Einstein equations into 
$
  cal(E)_(mu nu) + nabla_mu Z_nu + nabla_nu Z_mu - g_(mu nu) nabla^lambda Z_lambda = 0.
$<eqZ4EE>
The specific choice of relative normalisation and sign between the symmetrised gradient and divergence terms will become clear shortly; let us first review the fundamental properties of this modification.

*Remarks:*
+ This modification can be thought of as redefining the Ricci tensor,
  $
    R_(mu nu) -> tilde(R)_(mu nu) + nabla_mu Z_nu + nabla_nu Z_mu,
  $
  which causes the associated Ricci scalar to turn into
  $
    R -> tilde(R) = R + 2 nabla^lambda Z_lambda.
  $
  Hence, the Einstein tensor becomes
  $
    G_(mu nu) = R_(mu nu) - 1/2 g_(mu nu) R -> tilde(G)_(mu nu) &= tilde(R)_(mu nu) - 1/2 g_(mu nu) tilde(R) \ &= G_(mu nu) + nabla_mu Z_nu + nabla_nu Z_mu - g_(mu nu) nabla^lambda Z_nu.
  $
  Since $cal(E)_(mu nu) = G_(mu nu) + Lambda g_(mu nu) - 8pi T_(mu nu)$, this reproduces the modified Einstein equations above:
  $
    cal(E)_(mu nu) = 0 -> tilde(cal(E))_(mu nu) = cal(E)_(mu nu) + nabla_mu Z_nu + nabla_nu Z_mu - g_(mu nu) nabla^lambda Z_lambda = 0.
  $

+ Next, let us examine how the constraints $cal(H)$ and $cal(M)_mu$ are related to $Z_mu$. To do so, we simply take the corresponding projections of the modified Einstein tensor, where the $cal(E)_(mu nu)$ term reproduces the constraints and the remaining terms give us their relationship to $Z_mu$. It makes sense to also introduce names for the normal and tangential projections of $Z_mu$, writing
  $
    Theta = n^mu Z_mu, quad Z^perp_mu = tensor(P,+nu,-mu) Z_nu.
  $
  The full vector is then recovered as
  $
    Z_mu = - Theta n_mu + Z_mu^perp.
  $
  To obtain relationships between $cal(H),cal(M)_mu$ and the newly introduced variables $Theta$ and $Z_mu^perp$, we take the corresponding projections of Einstein's equations, starting with the doubly normal one: 
  #bottom-number[$
    0 &= 2 n^mu n^nu (cal(E)_(mu nu) + nabla_mu Z_nu + nabla_nu Z_mu - g_(mu nu) nabla^lambda Z_lambda)\
    &= cal(H) + 4 n^mu n^nu nabla_mu Z_nu + 2 nabla^lambda Z_lambda\
    &= cal(H) + 4 n^mu nabla_mu underbrace((n^nu Z_nu),=Theta) - 4 Z_nu underbrace(n^mu nabla_mu n^nu,=a^nu perp n^mu) -2 nabla^mu (n_mu Theta) + 2 nabla^mu Z_mu^perp\
    &= cal(H) + 4 n^mu nabla_mu Theta - 4 a^mu Z_mu^perp - 2 (nabla^mu n_mu)Theta - 2n_mu nabla^mu Theta + 2 (gamma^(mu nu) - n^mu n^nu) nabla_mu Z_nu^perp\
    &= cal(H) + 2 n^mu nabla_mu Theta - 4 a^mu Z_mu^perp + 2 K Theta + 2 mnabla^mu Z_mu^perp+ 2 a^nu Z_nu^perp\
    &= cal(H) + 2 n^mu nabla_mu Theta + 2 K Theta + 2 mnabla^mu Z_mu^perp - 2 a^mu Z_mu^perp.
  $<eqnnProjModEE>]
  Here, the raised index of $mnabla$ is to be interpreted as being raised with the metric the connection is compatible with, i.e. with $gamma^(mu nu)$. The above implies that
  $
    cal(H) = - 2 n^mu nabla_mu Theta - 2 K Theta - 2 mnabla^mu Z_mu^perp + 2 a^mu Z_mu^perp,
  $
  or in adapted coordinates, 
  $
    cal(H) = -2/alpha diff_t Theta + 2/alpha beta^i diff_i Theta - 2 K Theta - 2 mnabla^i Z_i + 2/alpha gamma^(i j) Z_i diff_j alpha.
  $
  Note that since $n_mu sim delta_mu^t$, $Z_i^perp = Z_i$. We can proceed similarly for $cal(M)_mu$:
  #bottom-number[$
    0 &= n^mu tensor(P,+nu,-lambda) (cal(E)_(mu nu) + nabla_mu Z_nu + nabla_nu Z_mu - g_(mu nu) nabla^rho Z_rho)\
    &= cal(M)_lambda + n^mu tensor(P,+nu,-lambda) nabla_mu Z_nu + n^mu tensor(P,+nu,-lambda) nabla_nu Z_mu - underbrace(n^mu gamma_(mu lambda),=0) nabla^rho Z_rho\
    &= cal(M)_lambda + n^mu nabla_mu underbrace((tensor(P,+nu,-lambda) Z_nu),=Z_lambda^perp) - (n^mu nabla_mu tensor(P,+nu,-lambda)) Z_nu + underbrace(tensor(P,+nu,-lambda) nabla_nu,=mnabla_lambda) underbrace((n^mu Z_mu),Theta) - Z_mu tensor(P,+nu,-lambda) #h(-1em) underbrace(nabla_nu n^mu,=-tensor(K,-nu,+mu) - n_nu a^mu)\
    &= cal(M)_lambda + n^mu nabla_mu Z_lambda^perp - n^mu nabla_mu (n^nu n_lambda) Z_nu + mnabla_lambda Theta + tensor(K,+mu,-lambda) Z_mu^perp\
    &= cal(M)_lambda + n^mu nabla_mu Z_lambda^perp + tensor(K,+mu,-lambda) Z_mu^perp - n_lambda a^mu Z_mu^perp - a_lambda Theta + mnabla_lambda Theta
  $<eqntProjModEE>]
  This implies
  $
    cal(M)_lambda = -n^mu nabla_mu Z_lambda^perp - tensor(K,+mu,-lambda) Z_mu^perp + n_lambda a^mu Z_mu^perp - a_lambda Theta + mnabla_lambda Theta
  $
  or in adapted coordinates, making use of @eqTangentCovectorNormalDerivative,
  $
    cal(M)_i = - 1/alpha diff_t Z_i + 1/alpha fL_beta Z_i - 2 tensor(K,+j,-i) Z_j- 1/alpha (diff_i alpha) Theta + diff_i Theta.
  $
  In summary, we have the covariant identities
  $
    cal(H) &= - 2 n^mu nabla_mu Theta - 2 K Theta - 2 mnabla^mu Z_mu^perp + 2 a^mu Z_mu^perp,\
    cal(M)_lambda &= -n^mu nabla_mu Z_lambda^perp - tensor(K,+mu,-lambda) Z_mu^perp + n_lambda a^mu Z_mu^perp - a_lambda Theta + mnabla_lambda Theta,
  $
  and in adapted coordinates, 
  $
  cal(H) &= -2/alpha diff_t Theta + 2/alpha beta^i diff_i Theta - 2 K Theta - 2 mnabla^i Z_i + 2/alpha gamma^(i j) Z_i diff_j alpha,\
  cal(M)_i &= - 1/alpha diff_t Z_i + 1/alpha fL_beta Z_i - 2 tensor(K,+j,-i) Z_j- 1/alpha (diff_i alpha) Theta + diff_i Theta.
  $
  Hence, $Theta$ and $Z_mu^perp$, or equivalently $Z_mu$, can be seen as a sort of "vector potential" for the constraints; we obtain the constraints by acting on $Z_mu$ with a linear differential operator.
  
  Note that the adapted coordinate expressions also give us the evolution equations for $Theta$ and $Z_i$---we can rearrange them for their time derivatives and obtain
  $
    diff_t Theta &= beta^i diff_i Theta - alpha K Theta - alpha mnabla^i Z_i + gamma^(i j) Z_i diff_j alpha - 1/2 alpha cal(H),\
    diff_t Z_i &= fL_beta Z_i - 2 alpha tensor(K,+j,-i) Z_j - (diff_i alpha)Theta + alpha diff_i Theta - alpha cal(M)_i.
  $
  This allows us to arrive at an important conclusion: while $(Theta,Z_i)$ form a mostly closed evolution subsystem, they are actively sourced by the constraint violations $cal(H)$ and $cal(M)_i$. 
+ Although in the previous remark, we already derived evolution equations in adapted coordinates, their character is rather opaque. Although we could technically proceed as in the previous section---take a second time derivative and replace any occuring first time derivative with the above---this would be a severely tedious undertaking. It is much easier to compute second-order evolution equations for $Z_mu$ directly from the covariant modified Einstein equations 
  $
    cal(E)_(mu nu) + nabla_mu Z_nu + nabla_nu Z_mu - g_(mu nu) nabla^lambda Z_lambda = 0
  $
  instead. In particular, this derivation will show us why the $-g_(mu  nu) nabla^lambda Z_lambda$ term is necessary; without it, the principal symbol would not be hyperbolic. 

  To proceed with the derivation, we recall that $nabla^mu cal(E)_(mu nu) = 0$. Hence, if the modified Einstein equations hold, the remaining terms together must be divergence-free as well, i.e.
  $
    0 &= nabla^mu (nabla_mu Z_nu + nabla_nu Z_mu - g_(mu nu) nabla^lambda Z_lambda)\
    &= Box_g Z_nu + nabla^mu nabla_nu Z_mu - nabla_nu nabla^mu Z_mu\
    &= Box_g Z_nu - tensor(R,+lambda,-mu,+mu,-nu) Z_lambda\
    &= Box_g Z_nu + tensor(R,+lambda,-nu) Z_lambda.
  $
  As alluded to earlier, the $-g_(mu nu) nabla^lambda Z_lambda$ term is needed---without it, we get additional second order terms because we cannot invoke the Ricci identity. 

  We conclude that $Z_mu$ satisfies a wave equation containing an additional curvature contribution which allows the components of $Z_mu$ to be "rotated into each other"---importantly, $Z_mu$ propagates hyperbolically, at the speed of light. In particular, since the constraint violations $cal(H)$ and $cal(M)_i$ derive from the action of a linear operator on $Z_mu$, this means that in evolution governed by the @eqZ4EE[modified Einstein equations], _all_ constraint violations propagate as waves. Further---since the trace-reversed Einstein equations give an expression for $R_(mu nu)$ in terms of the energy-momentum tensor and first derivatives of $Z_mu$---we can see that both the presence of matter as well as gradients in $Z_mu$ cause its components to mix and/or oscillate. 

+ Since we have introduced four additional degrees of freedom, we should ask if all of them are physical _in the framework of the generalised theory_. Within GR, of course, they are all non-physical; viewing the modified Einstein equations as a new physical system, though, we should check if $Z_mu$ has gauge freedom. Clearly---due to their linearity in $Z_mu$---the modified Einstein equations remain invariant under transformations $Z_mu -> Z_mu + X_mu$, provided that
  $
    nabla_mu X_nu + nabla_nu X_mu = g_(mu nu) nabla^lambda X_lambda.
  $
  Though this might look like some special case of the conformal Killing equation, the reality is much simpler (in $d!=2$): taking the trace yields
  $
    (d-2) nabla^lambda X_lambda = 0. 
  $
  This means that the condition breaks down to the Killing equation
  $
    nabla_mu X_nu + nabla_nu X_mu = 0.
  $
  Hence, a "physical" configuration in the sense of the modified theory fixes $Z_mu$ only up to a Killing field. The amount of physical degrees of freedom that are lost due to such transformations depends on the configuration of the spacetime. In general, the background will have no Killing vectors, and hence $Z_mu$ has no gauge freedom. However, for e.g. stationary/static or even more symmetric backgrounds, dimension of the isometry algebra is nonzero and finite, so that the moduli space---quotient of the configuration space of $Z_mu$ quotiented by the finite dimensional isometry algebra---has some of its physical configurations identified.
=== Adding Damping
Up until now, we have introduced the additional degrees of freedom $Z_mu$ and modified the Einstein equations in such a way that 
+ $Z_mu$ propagates as a wave;

+ and $cal(H)$ and $cal(M)_i$ emerge from the action of a linear first-order differential operator on $Z_mu$, in particular causing the implication $Z_mu = 0 => cal(H),cal(M)_i=0$.

Although this is already a significant improvement over unmodified GR---constraints now propagate---we would ideally introduce an additional mechanism that drives $Z_mu$ towards zero throughout its evolution. This would make small enough constraint violations decay over time, in consequence driving simulations towards the constraint manifold and hence improving stability. This is the goal of this section.

Before modifying the equations further to introduce damping, we should consider the implications for general covariance. The un-damped modified Einstein equations are fully covariant and independent of any choice of foliation. Introducing damping, however, inherently requires selecting a preferred timelike direction along which constraint violations decay---naturally provided by the unit normal $n = -alpha dt^sharp$ of a _chosen_ spacelike foliation.

Consequently, manifest four-dimensional covariance is broken. From the point of view of fundamental physics, this would be a significant drawback; however, Z4 is not intended as a new physical theory, but rather as a numerical tool developed to stabilise evolutions and drive solutions back onto the constraint manifold. Because numerical relativity fundamentally relies on a 3+1 split and thus the choice of a foliation, breaking manifest covariance for the sake of constraint damping is entirely justified, as long as the dynamics on the constraint manifold reproduce those of GR.

To introduce damping along $n$ into the wave equation satisfied by $Z_mu$, we must introduce terms of the form $n^mu nabla_mu X$, where $X$ is either $Z_mu$ itself or one of its projections. Ideally, we would like to have control over the damping of both $Theta$ and $Z_mu^perp$ individually, whence we should have multiple parameters we can control. For a damping term of the form $n^mu nabla_mu X$ to emerge in the wave equation for $Z_mu$, we have to add additional terms to the modified Einstein equations. Specifically, these terms must contain both $Z_mu$ and $n_mu$ linearly, which immediately leads us to the symmetric combinations
$
  n_mu Z_nu + n_nu Z_mu quad "and" quad g_(mu nu) n^mu Z_mu.
$
These are the only symmetric tensors linear in both $n_mu$ and $Z_mu$. We are thus motivated to introduce the tensor
$
  D_(mu nu) = -kappa_1[n_mu Z_nu + n_nu Z_mu + kappa_2 g_(mu nu) n^lambda Z_lambda]
$
to the modified Einstein equations. The specific way the tunable constants $kappa_1$ and $kappa_2$ are introduced will become clearer momentarily. We use $D_(mu nu)$ to modify the Einstein equations into
$
  cal(E)_(mu nu) + nabla_mu Z_nu + nabla_nu Z_mu - g_(mu nu) nabla^lambda Z_lambda + D_(mu nu) = 0.
$
We could check if we do in fact produce a $n^mu nabla_mu Z_nu$ term by computing the divergence of the above. Clearly, though, taking the divergence $nabla^mu D_(mu nu)$ will yield a plethera of terms that we don't really want to deal with, so let's entertain another approach. Recall that in the previous section, we derived _first-order_ equations of motion for $Theta$ and $Z_mu^perp$ by projecting the modified Einstein equations onto the normal-normal and normal-tangential directions. By repeating this, we generate additional terms from the projections of $D_(mu nu)$; to arrive at damping, we need the first-order decay structure
$
  n^mu nabla_mu Theta sim -Theta + ...,
$
and an analogous structure for $Z_mu^perp$.

Let us compute the relevant projections of $D_(mu nu)$. We start with the doubly normal projection,
$
  2 n^mu n^nu D_(mu nu) &= -kappa_1 n^mu n^nu lr([n_mu Z_nu + n_nu Z_mu + kappa_2 g_(mu nu) underbrace(n^lambda Z_lambda,=Theta)],size:#65%)\
  &= -2 kappa_1 [-Theta -Theta - kappa_2 Theta]\
  &= 2 kappa_1 (kappa_2 + 2) Theta.
$
Hence, appending this to the @eqnnProjModEE[doubly-normal-projected undamped Einstein equation], we get
$
  &&0&=cal(H) + 2 n^mu nabla_mu Theta + 2 K Theta + 2 mnabla^mu Z_mu^perp - 2 a^mu Z_mu^perp + 2 kappa_1 (kappa_2 + 2) Theta\
  <=>&quad& n^mu nabla_mu Theta&=- kappa_1 (kappa_2 + 2)Theta + K Theta - mnabla^mu Z_mu^perp + a^mu Z_mu^perp -1/2 cal(H)
$
Switching back to adapted coordinates, we get
$
  diff_t Theta = beta^i diff_i Theta - alpha kappa_1 (kappa_2 + 2)Theta + alpha K Theta - alpha mnabla^i Z_i + gamma^(i j) Z_i diff_j alpha - 1/2 alpha cal(H).
$
The normal-tangential projection is even simpler:
$
  n^mu tensor(P,+nu,-lambda) D_(mu nu) &= -kappa_1 n^mu tensor(P,+nu,-lambda)lr([n_mu Z_nu + n_nu Z_mu + kappa_2 g_(mu nu) underbrace(n^lambda Z_lambda,=Theta)],size:#65%)\
  &= -kappa_1[-Z_lambda^perp + 0 + 0]\
  &= kappa_1 Z_lambda^perp.
$
Consequently, the @eqntProjModEE[normal-tangential-projected undamped Einstein equations] turn into
$
  &&0&= cal(M)_lambda + n^mu nabla_mu Z_lambda^perp + tensor(K,+mu,-lambda) Z_mu^perp - n_lambda a^mu Z_mu^perp - a_lambda Theta + mnabla_lambda Theta\
  <=>&quad& n^mu nabla_mu Z_lambda^perp &= -kappa_1 Z_lambda^perp- tensor(K,+mu,-lambda) Z_mu^perp + n_lambda a^mu Z_mu^perp + a_lambda Theta - mnabla_lambda Theta - cal(M)_lambda
$
In adapted coordinates, this rearranges to
$
  diff_t Z_i = fL_beta Z_i - alpha kappa_1 Z_i - 2 alpha tensor(K,+j,-i) Z_j - (diff_i alpha)Theta + alpha diff_i Theta - alpha cal(M)_i
$
In summary, the new equations of motion for $Theta$ and $Z_i$ read 
$
  diff_t Theta &= beta^i diff_i Theta - alpha kappa_1 (kappa_2 + 2)Theta + alpha K Theta - alpha mnabla^i Z_i + gamma^(i j) Z_i diff_j alpha - 1/2 alpha cal(H)\
  diff_t Z_i &= fL_beta Z_i - alpha kappa_1 Z_i - 2 alpha tensor(K,+j,-i) Z_j - (diff_i alpha)Theta + alpha diff_i Theta - alpha cal(M)_i
$
Both the covariant and adapted coordinate forms of these equations show that both $Theta$ and $Z_i$ are exponentially damped, with damping factors $alpha kappa_1(kappa_2 + 2)$ and $alpha kappa_1$, respectively. 

=== Outlook: Z4c
Briefly outline what Z4c does, mention that also evolution equations are modified (doubly tangential projections), make note of trace reversal.
== #text(fill:red)[Wave Extraction & Diagnostics]
=== #text(fill:red)[Weyl Scalar $Psi_4$]
=== #text(fill:red)[Constraint Monitoring]
