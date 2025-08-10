#import "template.typ": *
#import "macros.typ": *
#import "@preview/xarrow:0.3.1": xarrow
#import "@preview/fletcher:0.5.8" as fletcher: diagram, node, edge

= Introduction

To formulate the ADM formalism, one not only requires Hamiltonian field theory and the Dirac-Bergmann algorithm, but also a means of defining a foliation of spacetime---thereby introducing a preferred time direction for evolution. This involves a precise formulation of what a foliation is, the construction and interpretation of its associated normal vector field, and, more broadly, a systematic treatment of differential geometry on submanifolds. Submanifolds carry a variety of important geometric structures: an induced metric, projected covariant derivatives, intrinsic as well as extrinsic curvature, and decompositions of the ambient curvature in terms of intrinsic and extrinsic contributions. These will all play a central role in what follows.

= Basics
== Definition of a Submanifold
In the following, $cal(M)$ will denote a ((pseudo)-Riemannian) $m$-manifold, and $x^mu$ local coordinates on it. We denote the components of the metric tensor on $cal(M)$ by $g_(mu nu) (x)$, with $mu,nu,...$ indicating coordinate indices ranging from $1,...,m$. The metric tensor defines a line element on $cal(M)$,
$
  ds^2 = g = g_(mu nu) dx^mu otimes dx^nu.
$
A _submanifold_ of $cal(M)$ is a subset $cal(S) subset cal(M)$ that is itself a manifold, say of dimension $s = dim cal(S) <= m$. Here it is understood that the topological structure on $cal(S)$ is the induced subtopology from $cal(M)$. Additional structures, such as the induced metric, connections and curvature, we will introduce gradually in the following sections.

We call the quantity
$
  codim_cal(M) cal(S) = dim cal(M) - dim cal(S) 
$
the _codimension_ of $cal(S)$. This is a simple but useful notion that conveys "how many dimensions less" the submanifold has in comparison to the ambient manifold $cal(M) supset cal(S)$, which sometimes is more important than the dimension itself. 

For example, a _hypersurface_ is, by definition, a submanifold $cal(S) subset cal(M)$ of codimension 1, and will be particularly important in the development of foliations and the ADM formalism. A simple way of generating (smooth) hypersurfaces is by considering a function $f:cal(M)->RR$ and constructing a level surface as
$
  cal(S) = f^(-1)({0}) = {p in cal(M) | f(p) = 0} subset cal(M).
$
For $cal(S)$ defined in this way to be a proper, smooth submanifold, $0$ has to be a regular value of $f$. This means that 
$
  d f|_cal(S) != 0
$
or in words, the differential (or equivalently, the gradient) of $f$ must not vanish on $cal(S)$.

More generally, given a collection of smooth functions $f_i :cal(M)->RR$ for $i=1,...,k$, the common level set
$
  cal(S) = sect.big_i f_i^(-1)({0})
$
is a submanifold of codimension $k$, provided the $f_i$ are functionally independent, $cal(S)$ is non-empty and $0$ is a regular value of the map $f= (f_1,...,f_k)$.
== Induced Metric <sectionInducedMetric>
So far, we have established that a submanifold $cal(S) subset cal(M)$ inherits its topological structure from the subspace topology induced by $cal(M)$. This makes $cal(S)$ a topological manifold. However, since we are studying submanifolds of a (pseudo-)Riemannian manifold, it should come as no surprise that we also wish for $cal(S)$ to carry a (pseudo-)Riemannian structure.

As a manifold in its own right, $cal(S)$ admits local coordinates $y^i$, where $i = 1,...,s = dim cal(S)$. In principle, endowing $cal(S)$ with a (pseudo-)Riemannian metric is straightforward: one simply picks a symmetric tensor field $gamma in T cal(S) otimes T cal(S)$ with the desired signature and defines the corresponding line element
$
  dms^2 = gamma = gamma_(i j) dy^i otimes dy^j.
$
While this makes $cal(S)$ a (pseudo-)Riemannian manifold, it does so independently of the geometry on the ambient manifold $cal(M)supset cal(S)$---it disregards the structure already present on $cal(M)$. 

Yet $cal(S)$ is fully contained in $cal(M)$; any point $p in cal(S)$ is also a point of $cal(M)$. That much is obvious---but it matters, because the metric on $cal(M)$ carries geometric information we may wish to preserve. What, after all, does a metric provide? In particular, it can be used to assign lengths to curves. If we consider a path lying in $cal(S)$, then this path is also a path in $cal(M)$. Naturally, we would want its length to be the same whether we compute it using the metric $ds^2 = g$ of the ambient manifold or the metric $dms^2 = gamma$ of the submanifold.

This requirement reflects a choice, but a canonical one: the geometry of $cal(S)$ should arise from that of $cal(M)$ by restriction, not by arbitrary redefinition. That is, the structures on $cal(M)$ should be _induced_ on $cal(S)$, not constructed from scratch. 

Requiring that all lengths of curves in $cal(S)$ agree under both metrics leads to the precise condition that
$
  dms^2 = ds^2|_cal(S).
$<inducedMetricDefn>
This is no longer just a heuristic, but an equation---it specifies how the metric on the submanifold must be related to the metric of the ambient space. 

Though precise, the definition in @inducedMetricDefn[relation] is not especially practical for computations. Every time we wish to measure the length of a path on $cal(S)$, we must return to the structure on $cal(M)$ and consult what _it_ "thinks" the length is. We would prefer to treat $cal(S)$ as a manifold in its own right, without this constant indirection via $cal(M)$. To achieve that, we need explicit knowledge of the components $gamma_(i j)$ of the _induced metric_ on $cal(S)$---the symmetric tensor satisfying
$
  gamma_(i j) dy^i otimes dy^j = dms^2 attach(=,t:!) ds^2|_cal(S) = (g_(mu nu) dx^mu otimes dx^nu)|_cal(S).
$
The remainder of this section is devoted to deriving a general formula for $gamma_(i j)$ in terms of the ambient metric $g_(mu nu)$.

To do so, observe that for any point $p in cal(S)$, we have two coordinate systems at our disposal: the intrinsic coordinates $y^i (p)$ on $cal(S)$, and the ambient coordinates $x^mu (p)$ from $cal(M)$. More formally, since $cal(S) subset cal(M)$, we obtain an _injective coordinate map_
$
  y^i |-> x^mu (y^i),
$
expressing the ambient coordinates as functions of the submanifold coordinates. This has a consequence for the differential of $x^mu$ restricted to $cal(S)$, namely
$
  dx^mu|_cal(S) = d(x^mu (y^i))|_cal(S) = (diff x^mu)/(diff y^i) dy^i =: E^mu_i dy^i,
$
where we defined the so-called _pushforward matrix_ or _vielbein_
$
  E_i^mu  = (diff x^mu)/(diff y^i). 
$<pushforwardEDefn>
Using this relation, we can rewrite the restriction of the line element as
$
   gamma_(i j) dy^i otimes dy^j &= dms^2 = ds^2|_cal(S)\
    &= (g_(mu nu) dx^mu otimes dx^nu)|_cal(S)\
    &= g_(mu nu)(x^mu (y^i)) E^mu_i dy^i otimes E^nu_j dy^j\
  &= (E^mu_i E^nu_j g_(mu nu)) dy^i otimes dy^j
$
From this, we identify the components of the induced metric as
$
  gamma_(i j) = E^mu_i E^nu_j g_(mu nu)|_cal(S) = (diff x^mu)/(diff y^i) (diff x^nu)/(diff y^j) g_(mu nu)|_cal(S).
$
This gives a general formula for the induced metric on a submanifold, expressed in terms of the ambient metric and a choice of compatible#footnote[Obviously, the coordinate charts must overlap appropriately and be smooth enough for the construction to be valid.] coordinate systems on both $cal(M)$ and $cal(S)$.

This expression for $gamma_(i j)$ guarantees that the length of any path measured intrinsically on $cal(S)$, via the line element $dms^2$, matches the length computed externally in the ambient manifold via $ds^2$. It is, of course, a coordinate-dependent representation---whereas the earlier definition in @inducedMetricDefn[equation] was fully coordinate-independent.
== Example Submanifold: $S^2 subset RR^3$
Before diving into more abstract constructions, let us ground our intuition with a concrete example of a manifold and a natural submanifold. To this end, we consider the flat Euclidean space $cal(M) = RR^3$, with Cartesian coordinates $x^mu = (x,y,z) in RR^3$ and endowed with the Euclidean metric
$
  ds^2 = dx^2 + dy^2 + dz^2 quad <=> quad g_(mu nu) = delta_(mu nu).
$<euclideanLineElement>
As our submanifold $cal(S)$, we consider the $2$-sphere $S^2$, defined by
$
  S^2 = {(x,y,z) in RR^3 mid(|) x^2 + y^2 + z^2 = 1}.
$
It is clear that this is a submanifold, as it is the preimage of ${0}$ under the map
$
  f(x,y,z) = x^2 + y^2 + z^2 - 1,
$
which has $0$ as a regular value, i.e., the gradient vanishes nowhere where $f(x,y,z) = 0$. 

The 2-sphere is most conveniently parameterised by two angular coordinates $y^i = (theta,phi) in (0,pi) times (0,2pi)$, with the embedding into the ambient space $RR^3$ given by
$
  x(theta,phi) &= sin theta cos phi,\
  y(theta,phi) &= sin theta sin phi,\
  z(theta,phi) &= cos theta.
$<embeddingS2inR3>
We now compute the induced metric on $S^2$ using two methods: first by means of the pushforward matrix, and then in a more direct way involving differentials of the above relations. 

The components of $E^mu_i = (diff x^mu)/(diff y^i)$ are given by
$
  E^x_theta &= (diff x)/(diff theta) = cos theta cos phi, &quad&& E^x_phi &= (diff x)/(diff phi) = -sin theta sin phi,\
  E^y_theta &= (diff y)/(diff theta) = cos theta sin phi, &&& E^y_phi &= (diff y)/(diff phi) = sin theta cos phi,\
  E^z_theta &= (diff z)/(diff theta) = -sin theta,&&& E^z_phi &= (diff z)/(diff phi) = 0.
$
This allows us to compute the components of the induced metric as
$
  gamma_(theta theta) &= E_theta^mu E_theta^mu delta_(mu nu) = E^x_theta 
  E^x_theta + E^y_theta E^y_theta + E^z_theta E^z_theta\
  &= cos^2 theta cos^2 phi + cos^2 theta sin^2 phi + sin^2 theta\
  &= 1,
$ 
$
  gamma_(phi phi) &= E_phi^mu E_phi^mu delta_(mu nu) = E^x_phi 
  E^x_phi + E^y_phi E^y_phi + E^z_phi E^z_phi\
  &= sin^2 theta sin^2 phi + sin^2 theta cos^2 phi\
  &= sin^2 theta,
$
$
  #h(0.9em)gamma_(theta phi) = gamma_(phi theta) &= E^x_theta E^x_phi + E^y_theta E^y_phi + E^z_theta E^z_phi\
  &= - sin theta cos theta sin phi cos phi + sin theta cos theta sin phi cos phi\
  &= 0.
$
Using this, the induced line element becomes
$
  dms^2 = gamma_(i j) dy^i otimes dy^j = dtheta otimes dtheta + sin^2 theta dphi otimes dphi
$
---an anticipated result, the standard metric on $S^2$.

We now briefly go over how one can compute the induced metric slightly more directly by considering the differentials of the @embeddingS2inR3[embeddings]. These differentials read
$
  dx|_(S^2) &= cos theta cos phi dtheta - sin theta sin phi dphi,\
  dy|_(S^2) &= cos theta sin phi dtheta + sin theta cos phi dphi,\
  dz|_(S^2) &= -sin theta d theta.
$
We can then simply insert these relations into the @euclideanLineElement[Euclidean line element] to find (writing $dx^mu otimes dx^nu$ as $dx^mu dx^nu$ for brevity)
$
  dms^2 &= ds^2|_(S^2) = dx^2|_(S^2) + dy^2|_(S^2) + dz^2|_(S^2)\
        &= (cos theta cos phi dtheta - sin theta sin phi dphi)^2\
        &quad+ (cos theta sin phi dtheta + sin theta cos phi dphi)^2\
        &quad+ (-sin theta d theta)^2\
        &= cos^2 theta cos^2 phi dtheta^2 - cancelr(2 cos theta sin theta cos phi sin phi dtheta dphi) + sin^2 theta sin^2 phi dphi^2\
        &quad+ cos^2 theta sin^2 phi dtheta^2 + cancelr(2 cos theta sin theta cos phi sin phi dtheta dphi) + sin^2 theta cos^2 phi dphi^2\
        &quad+ sin^2 theta dtheta^2\
        &= (cos^2 theta (cos^2 phi + sin^2 phi) + sin^2 theta) dtheta^2 + sin^2 theta (sin^2 phi  + cos^2 phi) dphi^2\
        &= dtheta^2 + sin^2 theta dphi^2.
$
This alternative method avoids explicitly constructing the pushforward matrix, relying instead on computing the differentials of the embedding directly and inserting them into the ambient metric---though of course, fully equivalent to the first approach. _Which_ of the two approaches is better comes down to being mostly situational, and of course, personal preference.

= Tangent Spaces, Pushforwards, Pullbacks, and Projections
== Preliminaries: Maps Between Manifolds and Their Derivatives
Before diving into the formalism of pushforwards and pullbacks, which is the goal of this section, let us first understand the geometric motivation behind them. 

Suppose we have two smooth manifolds $cal(M)$ and $cal(N)$, with dimensions $m$ and $n$, and $phi:cal(M) -> cal(N)$ a smooth map between them. This map assigns to every point $p in cal(M)$ a point $phi(p) in cal(N)$. But more than just matching up points, $phi$ also relates the local _geometry_ near those points.

To illustrate this, imagine a flexible sheet of paper $(cal(M))$ that we smoothly press onto another, possibly curved surface $cal(N)$. Every point on the paper ends up on some point of the surface, but the way the paper bends or stretches as it moves gives us more than just a positional correspondence---it tells us how _directions_ and _infinitesimal displacements_ on $cal(M)$ are transformed under the map. For example, a tiny arrow drawn on the paper (a tangent vector) will point in some new direction when the paper is curved onto the surface.

This is the essence of what the pushforward captures: it tells us how tangent vectors on $cal(M)$---representing directions of motion away from a point---are mapped to tangent vectors on $cal(N)$. Formally, the _pushforward of $phi$_, denoted by $phi_*$, is a map between tangent spaces,
$
  phi_* : T_p cal(M) -> T_(phi(p))cal(N).
$
We will now make this precise. 
\ \ 
*Definition:* (Pushforward) We begin by recalling what we want to construct: given a smooth map $phi:cal(M)-> cal(N)$, we aim to associate to any vector $X in T_p cal(M)$ a vector $phi_* X in T_(phi(p)) cal(N)$, which we call the _pushforward_ of $X$ under $phi$. This should represent the idea of how the direction $X$ "looks" after the map $phi$ has distorted the space onto $cal(N)$. 

But how do we compare directions at two different points of two different manifolds? We us the fact that vectors act on functions as differential operators---they are directional derivatives. A vector $X in T_p cal(M)$ acting on a function on $cal(M)$ tells us how it changes as we move in the direction of $X$. Similarly, $phi_* X in T_phi(p) cal(N)$ should tell us how a function on $cal(N)$ changes in the "corresponding direction" under $phi$. 

Here's the key trick: we cannot directly evaluate how $X$ acts on a function on $cal(N)$, because $X$ lives on $cal(M)$. But we _can_ make any function on $cal(N)$ into a function on $cal(M)$ by composing it with $phi$. That is, given any smooth function $f : cal(N) -> RR$, the composition $f compose phi : cal(M) -> RR$ is a new function defined on $cal(M)$ to which we _can_ apply $X$.

The function $f compose phi$ is, in essence, just $f$, but with the codomain rearranged according to $phi$. The function values don't change---only _where_ in the domain they are sampled from does. This captures the distortion $phi$ induces: how orientation, stretching, or compression affect the function's appearance when viewed on $cal(M)$. So taking $X[f compose phi]$ tells us how $f$ changes along the image of the direction $X$ under the map $phi$. 

This leads us naturally to the definition
$
  (phi_* X)[f] := X[f compose phi], quad forall f in C^infty (cal(N)).
$<pushforwardDefn>
That is, $phi_* X$ is the unique vector at $phi(p)$ action on any function $f$ on $cal(N)$ agrees with the action of $X$ on the pulled-back version of $f$. This expression captures exactly what the pushforward is doing: it tells us how to "translate" a directional derivative on $cal(M)$ to one on $cal(N)$, while preserving how it affects functions. 

This is all very abstract, so let us make it more concrete with an analogy. Let us think of $cal(M)$ as the flat Mercator map, a rectangle showing latitude and longitude lines. Then, we think of $cal(N)$ as the actual globe, a $2$-sphere. The map $phi:cal(M)->cal(N)$ we envision as the wrapping map, which tells us how to take the rectangle and bend it into a globe. Suppose $f$ assigns temperature values across the globe---a real-valued function on $cal(N)$. Then $f compose phi$ is the pulled-back function: it contains the same information, the same temperature values, but now shown on the Mercator map instead of the globe. This function $f compose phi$ lets us look at the globe's data $f$ as if it lived on the Mercator map. We then ask: how does the temperature change when moving along some vector $X$ on that map? This is what $X[f compose phi]$ computes: how the pulled-back temperature map changes as we follow $X$ on the rectangle. 

Now, on the globe, there should be a corresponding direction, $phi_* X$, along which the same change occurs---after all, we have the same function values, just arranged differently. This is the direction we have to move towards on the sphere to see the same change in temperature. This should not just happen for this one temperature map---for any function $f$ on the globe, the pushforward $phi_* X$ must give us the corresponding direction whose action on $f$ replicates the directional derivative seen on the Mercator map. And that is what the @pushforwardDefn[definition] encodes.

For explicit calculations, it is typically more convenient to move away from the abstract, coordinate-independent @pushforwardDefn[definition], and instead adopt a coordinate-based approach. To this end, let $p in cal(M)$ a point, with local coordinates $y^alpha$ around $p$, and let $x^mu$ be local coordinates on $cal(N)$ around the image point $phi(p)$. Then $phi$ induces a local map between coordinate systems,
$
  x^mu (y^alpha) :=x^mu (phi(p(y^alpha))):RR^m -> RR^n,
$
which is differentiable in the standard sense. That is, composing the chart inverse on $cal(M)$, the map $phi$, and the chart on $cal(N)$, we obtain a smooth map between subsets of Euclidean space.

To derive the coordinate expression for the pushforward, consider the coordinate bases
$
  {diff_alpha = diff/(diff y^alpha)} "of" cal(M) quad "and" quad {diff_mu = diff/(diff x^mu)} "of" cal(N).
$
A vector $X in T_p cal(M)$ can be written as $X = X^alpha diff_alpha$, and its pushforward as $phi_* X = (phi_* X)^mu diff_mu$. The @pushforwardDefn[definition] then becomes to
$
  (phi_* X)^mu diff/(diff x^mu) f(x) &= X^alpha diff/(diff y^alpha) f(x(y))\
  &= X^alpha (diff x^mu)/(diff y^alpha) diff/(diff x^mu) f(x(y)),
$ 
where we applied the standard multidimensional chain rule on the right-hand side. We may hence identify the components of the pushforward $phi_* X$ as
$
  (phi_* X)^mu = (diff x^mu)/(diff y^alpha) X^alpha
$
This shows that the pushforward acts linearly on the components of the vector, via the Jacobian matrix
$
  tensor((phi_*),+mu,-alpha) = (diff x^mu)/(diff y^alpha) = diff/(diff y^alpha) x^mu (phi(p(y^alpha))),
$
i.e.
$
  (phi_* X)^mu = tensor((phi_*),+mu,-alpha) X^alpha.
$
This expression makes the pushforward entirely concrete: it is simply a linear transformation of the vector components under the coordinate ap induced by $phi$. In terms of implementation and computation, it behaves exactly like a Jacobian---and it reproduces the familiar coordinate transformation rule for vector components when $phi:cal(M)->cal(M)$ is a diffeomorphism.
\ \
*Definition* (Canonical Pairing) For a vector $X in T_p cal(M)$ and a 1-form $omega in T_p^* cal(M)$, we define their _canonical pairing_ (sometimes called an inner product)
$
  inprod(dot,dot)_cal(M) : T_p^* cal(M) times T_p cal(M) -> RR
$
by
$
  inprod(omega, X)_cal(M) :=omega(X).
$
Here, $omega(X)$ denotes the natural pairing between a covector and a vector---i.e., the evaluation of the linear map $omega:T_p cal(M) ->RR$ on the argument $X in T_p cal(M)$. In coordinate components, this reads $omega(X) = omega_mu X^mu$. Notice that this pairing is defined independently of any metric structure.
We will use it below to define the pullback $phi^*$, a dual operation to the pushforward, acting on cotangent spaces.
\ \
*Definition* (Pullback) We define the _pullback_ associated to $phi$ as the unique linear map
$
  phi^* : T_phi(p)^* cal(N) -> T_p^* cal(M)
$
satisfying
$
  inprod(phi^* omega, X)_cal(M) = inprod(omega, phi_* X)_cal(N), quad forall X in T_p cal(M), quad omega in T_phi(p)^* cal(N).
$<pullbackDefn>
That is, the pairing between a covector and a vector is preserved under the pushforward-pullback action of $phi$. This is a natural requirement: having defined how directions deform under $phi$, we want to ensure that objects measuring those directions (i.e. 1-forms) adapt in a compatible way. Equivalently, this can be written as
$
  (phi^* omega)(X) = omega(phi_* X), quad forall X in T_p cal(M), quad omega in T_phi(p)^* cal(N).
$

Another perspective is the following: a 1-form projects a component from a vector. But if the vector is altered---say, by a pushforward---then extracting the "same" component now requires a different projection. The pullback gives us precisely this adjusted 1-form: one that reproduces the same scalar when applied to the original vector as the original 1-form applied to the pushforward.

As with the pushforward, let us now extract a coordinate-based expression for the pullback. This case turns out to be even simpler. Using the same coordinate setup---local coordinates $y^alpha$ around $p in cal(M)$, and $x^mu$ around $phi(p) in cal(N)$---we express
$
  X = X^alpha diff_alpha, quad omega = omega_mu dx^mu.
$
We begin by expanding the left-hand side of the @pullbackDefn[definition],
$
  inprod(phi^* omega, X)_cal(M) = (phi^* omega)#h(0em)_alpha X^alpha.
$
On the right-hand side, we compute the action of $omega$ on the pushforward of $X$,
$
  inprod(omega,phi_* X)_cal(N) = omega_mu (phi_* X)^mu = omega_mu tensor((phi_*),+mu,-alpha)X^alpha.
$
Equating both sides and using the arbitrariness of $X^alpha$, we identify the components of the pullback:
$
  (phi^* omega)#h(0em)_alpha = tensor((phi_*),+mu,-alpha) omega_mu = (diff x^mu)/(diff y^alpha) omega_mu. 
$
In words: the pullback transforms a 1-form using the same Jacobian matrix as the pushforward acts with on vectors, but the roles of the indices reversed---covariant versus contravariant transformation. Equivalently, the pullback acts via the transpose of the Jacobian of $phi$, in keeping with the duality between tangent and cotangent spaces.
\ \
*Remark* _Pullbacks of Vectors?_ We have seen the pushforward as a map
  $
    phi_* : T_p cal(M) -> T_(phi(p)) cal(N).
  $
  A natural question is whether this relationship can be reversed: can we pull back vectors from $T_phi(p) cal(N)$ to $T_p cal(M)$? 

  In general, the answer is no. Algebraically, this follows from the fact that $phi_*$ is not necessarily invertible. As a linear map, invertibility requires its component representation
  $
    tensor((phi_*),+mu,-alpha)  =(diff x^mu)/(diff y^alpha)
  $
  to be an invertible matrix. This is only possible if $dim cal(M) = dim cal(N)$, since otherwise the Jacobian is rectangular and hence cannot be inverted. And even when the dimensions match, the Jacobian must be non-singular (i.e. its determinant must be nonzero). Thus, $phi_*$ is only invertible when $phi$ is a diffeomorphism---a smooth bijective map with smooth inverse.

  This is the algebraic reasoning. Fortunately, there's also geometric intuition to support it. The pushforward tracks how a direction away from $p in cal(M)$ translates to a direction away from $phi(p) in cal(N)$, consistent with the local deformation $phi$ induces. In that process, we may preserve the number of independent directions, or lose some---but we can never gain new ones. The rank of a linear map $phi_*:RR^m -> RR^n$ is at most $m$, so if $n>m$, its image lies in a proper subspace of $RR^n$. In other words, a pushforward cannot invent new degrees of freedom. 

  In the best-case scenario $dim im phi_* = m$, we retain full information and the map $phi_*$ is invertible. But if information is lost---i.e. if $phi_*$ collapses multiple directions into the same image---then we can no longer recover those original directions when attempting to go back. We do not know how to "blow it back up again"---that is, how to reconstruct the full-dimesional direction space of $T_p cal(M)$. A helpful instance to picture is a pushforward collapsing $T_p cal(M)$ onto a lower-dimensional plane inside $T_phi(p) cal(N)$; there is no way to lift this back up uniquely.

  Therefore, in general, it is not possible to pull back vectors from $T_phi(p) cal(N)$ to $T_p cal(M)$. 
\ \
*Generalisation to Tensors* 
Now that we know how we should associate vectors and 1-forms/covectors on two manifolds $cal(M)$ and $cal(N)$ between which we have an association of points given by $phi: cal(M)->cal(N)$, we can generalise to tensors. Due to the difference in "functorial direction" between the pushforward and pullback, it is not possible to define a combined map for arbitrary, mixed tensors. For fully covariant or fully contravariant tensors, however, it is possible---or, if $phi$ is a diffeomorphism. This latter case we will not consider, though---this is because it essentially reduces to a change of coordinates. 

Let us begin with the pullback of covariant tensors $T in T^(0,s)_phi(p) cal(N)$. Such a tensor can be written as 
$
  T = T_(mu_1...mu_s) dx^(mu_1)otimes...otimes dx^(mu_s),
$
which allows us to naturally pair it with $s$ vectors $X_((i)) in T_phi(p) cal(N)$ to form a scalar
$
  T(X_((1)),...,X_((s))) = T_(mu_1...mu_s) X^(mu_1)_((1))... X^(mu_s)_((s)).
$
To ensure compatibility between $T$ and its pullback $phi^* T in T_p^(0,s) cal(M)$, we impose that for all vectors $X_((1)),...X_((s)) in T_p cal(M)$, we have
$
  (phi^* T)(X_((1)),...,X_((s))) = T(phi_* X_((1)),...,phi_* X_((s)))
$
In words, we require the pullback of $T$ to "have the same opinion" about a set of vectors as the original tensor has about their pushforwards---which it should, if they are to encode the "same" tensorial structure on the two spaces related by $phi$. It shouldn't matter whether we ask $phi^*T$ on the original vectors for a scalar value, or first squish the space around and then ask the same question to $T$, with the correspondents $phi_* X_((i))$ as arguments. 

Of course, this definition can also be given in terms of coordinate components. Expanding both sides of the equation, we get
$
  (phi^* T)_(alpha_1...alpha_s) X_((1))^(alpha_1) ... X_((s))^(alpha_s) = T_(mu_1...mu_s) tensor((phi_*),+mu_1,-alpha_1) X^(alpha_1)_((1)) ... tensor((phi_*),+mu_s,-alpha_s) X^(alpha_s)_((s)).
$
This identifies the components of the pullback as
$
  (phi^* T)_(alpha_1...alpha_s) = tensor((phi_*),+mu_1,-alpha_1)  ... tensor((phi_*),+mu_s,-alpha_s)T_(mu_1...mu_s).
$
Hence, the pullback simply acts multilinearly on all indices by contraction with the pushforward matrix $tensor((phi_*),+mu,-alpha)$.

We now move to the pushforward for contravariant objects, i.e. of tensors $T in T^(r,0)_p cal(M)$. These can be represented as
$
  T = T^(alpha_1...alpha_r) thin diff_alpha_1 otimes ...otimes diff_(alpha_r)
$
in terms of the coordinate basis $diff_alpha_1 = diff/(diff y^(alpha_1))$. Such a tensor is a multilinear map, taking $r$ 1-forms $omega^((i)) in T_p^* cal(M)$ and mapping them to the scalar
$
  T(omega^((1)),...,omega^((r))) &= T^(alpha_1 ... alpha_r) thin omega^((1))(diff_alpha_1)... omega^((r))(diff_alpha_r)\
  &= T^(alpha_1 ... alpha_r)thin omega^((1))_(alpha_1) ... omega_(alpha_r)^((r)).
$
We define the pushforward $phi_* T in T_phi(p)^(r,0) cal(N)$ as the unique tensor satisfying
$
  T(phi^* omega^((1)), ..., phi^* omega^((r))) = (phi_* T) (omega^((1)), ..., omega^((r))) quad forall omega^((i)) in T_phi(p)^* cal(N).
$
In words: The pushforward $phi_* T$ "says" the same about the $omega^((i))$ as the original tensor $T$ would have to say about their pullbacks. In terms of coordinates, we verify that this really is just the multilinear extension of what happens to vectors: Inserting coordinate expressions for both sides we get
$
  T^(alpha_1...alpha_r) tensor((phi_*),+mu_1,-alpha_1)omega_(mu_1)^((1)) ... tensor((phi_*),+mu_r,-alpha_r)omega_(mu_r)^((r)) = (phi_* T)^(mu_1...mu_r) omega_(mu_1)^((1)) ... omega_(mu_r)^((r)).
$
We may thus identify the components of $phi_* T$ as
$
  (phi_* T)^(mu_1...mu_r) = tensor((phi_*),+mu_1,-alpha_1) ... tensor((phi_*),+mu_r,-alpha_r) T^(alpha_1...alpha_r),
$
which is simply just the index-wise transformation with the pushforward matrix $tensor((phi_*),+mu,-alpha)$. 

== Special Case: Submanifolds and the Pushforward as Inclusion

Now that we have established some general theory of smooth maps between manifolds---along with their associated pushforwards and pullbacks---we are ready to specialise to a particularly important class: the inclusion maps of submanifolds. Given a submanifold $cal(S) subset cal(M)$, the inclusion map $iota:cal(S) arrow cal(M)$ allows us to view $cal(S)$ as embedded in the ambient space $cal(M)$

This setup will naturally recover familiar constructions from Riemannian geometry. Most notably, we will see that the induced metric on the submanifold arises as the pullback of the ambient metric under $iota$.  Framed this way, the induced geometry of submanifolds is revealed to be not an isolated trick, but part of the broader formalism of pullbacks---shedding light on how submanifold geometry fits seamlessly into the general machinery of differential geometry. 

*Definition* (Inclusion Map for Submanifolds) Let $cal(M)$ be an $m$-dimensional manifold and $cal(S) subset cal(M)$ an $s$-dimensional submanifold. We call the map
$
  iota: cal(S) -> cal(M), quad p|->p
$
the _inclusion map_ or simply _inclusion_ of $cal(S)$ into $cal(M)$. 

At first glance, this may not seem particularly interesting---it merely sends each point to itself. But despite its simplicity, $iota$ is a smooth map between manifolds, and as such, all the machinery developed previously (pushforward, pullback, tensor behaviour) can and does apply. This makes the inclusion map a powerful conceptual tool, especially when studying how geometric or tensorial structures on the ambient manifold $cal(M)$ restrict or induce structure on the submanifold $cal(S)$. So, let us now examine $iota$ through this lens.
\ \
*The Pushforward $iota_*$* Firstly, let us examine the pushforward $iota_*$. Before diving into concrete definitions, it helps to build some intuition. In general, the pushforward $phi_*$ of a smooth map $phi:cal(M)->cal(N)$ maps tangent vectors from $T_p cal(M)$ to $T_phi(p) cal(N)$, capturing how directions around $p$ deform under $phi$. 

In our case, $phi = iota$ is the inclusion of a submanifold $cal(S) subset cal(M)$. Since $iota(p) = p$, the point itself remains unchanged, and so do the local relationships between points. A tangent vector $X in T_p cal(S)$ should therefore remain to be the same geometric object under the pushforward, just now interpreted as living inside the larger ambient space $T_p cal(M)$. We are thus led to expect that $iota_*$ simply embeds $T_p cal(S)$ into $T_p cal(M)$. 

This is what we will now make precise. For that, we consider the pushforward $iota_* X$ of some vector $X in T_p cal(S)$. By definition, for any function $f in C^infty (cal(M))$, its action is given by
$
  (iota_* X)[f] = X[f compose iota].
$
Since $iota(p) = p$ for all $p in cal(S)$, the composition $f compose iota$ is simply the restriction of $f$ to the submanifold, i.e. $f|_cal(S)$. That is,
$
  (iota_* X)[f] = X[f|_cal(S)]. 
$
This shows that $iota_* X in T_p cal(M)$ acts on ambient functions $f in C^infty (cal(M))$ exactly as $X in T_p cal(S)$ acts on their restrictions to $cal(S)$. In other words, the pushforward simply embeds the vector $X$ into the ambient tangent space in a way that preserves its actions on functions.

We therefore conclude that the pushforward acts as the canonical inclusion of tangent spaces,
$
  iota_*: T_p cal(S) -> T_p cal(M), quad X|->X, quad "for all" p in cal(S).
$
That is, we naturally identify each tangent space $T_p cal(S)$ with a subspace of the ambient tangent space $T_p cal(M)$, and the map $iota_*$ serves as a pointwise injective linear embedding.

Let us now consider how the pushforward looks in coordinate components. To this end, let $p in cal(S) subset cal(M)$, and suppose $y^i$, $i = 1,...,s$ and $x^mu$, $1,...,m$, are local coordinates around $p$ on $cal(S)$ and $cal(M)$, respectively. We consider the pushforward of a vector $X in T_p cal(S)$, written in the coordinate basis $diff_i = diff/(diff y^i)$ as 
$
  X = X^i diff_i
$ 
From the general theory of pushforwards established in the previous section, its image under the inclusion $iota$ is given by
$
  iota_* X = (iota_* X)^mu diff_mu = tensor((iota_*),+mu,-i)X^i diff_mu
$
where $diff_mu = diff/(diff x^mu)$ is the coordinate basis of $T_p cal(M)$, and the pushforward matrix $tensor((iota_*),+mu,-i)$ is given by the Jacobian matrix
$
  tensor((iota_*),+mu,-i) = (diff x^mu)/(diff y^i).
$
This is the same object we previously denoted by $tensor(E,+mu,-i)$ in @pushforwardEDefn[equation], and we now return to that notation. Thus, the pushforward becomes
$
  iota_* X = (tensor(E,+mu,-i) X^i) diff_mu quad <=> (iota_* X)^mu = tensor(E,+mu,-i) X^i.
$
*The Pullback $iota^*$* We now turn to the pullback $iota^* : T_p^* cal(M) -> T_p^* cal(S)$ associated with the immersion $iota: cal(S) -> cal(M)$. By definition, for any $X in T_p cal(S)$ and $omega in T_p^* cal(M)$, it satisfies
$
  inprod(iota^* omega, X)_cal(S) = inprod(omega, iota_* X)_cal(M).
$
In terms of components, this reads
$ 
  (iota^*omega)#h(0em)_i X^i = omega_mu (iota_* X)^mu = omega_mu tensor(E,+mu,-i) X^i,
$
from which we identify the components of the pullback as
$
  (iota^* omega)#h(0em)_i = tensor(E,+mu,-i) omega_mu.
$
*The Induced Metric* As we have seen in the previous section, we can also define the pushforward and pullback for purely contra- or covariant tensors, by applying the transformation with $tensor(E,+mu,-i)$ to each index separately. This arose from the requirement of compatibility under the map between manifolds of the multilinear map a tensor defines. A tensor of particular interest to differential geometry is the metric. Since it is a $(0,2)$-tensor, it is purely covariant, and we can apply the pullback to the metric on $cal(M)$ to get a tensor on $cal(S)$. Intuition would tell us that this is the induced metric, but let us go through this calmly. The metric on $cal(M)$ is the symmetric tensor
$
  ds^2 = g = g_(mu nu) dx^mu otimes dx^nu.
$
It is a symmetric bilinear map of two vectors $A,B in T_p cal(M)$ to $RR$, with
$
  g(A,B) = g_(mu nu) dx^mu (A) dx^nu (B) = g_(mu nu) A^mu B^nu.
$
Its pullback $iota^* g$ is defined by the relationship
$
  (iota^* g)(X,Y) = g(iota_* X, iota_* Y), quad X,Y in T_p cal(S).
$<inducedMetricDefnAbstract>
Let us briefly interpret this, then we move to the coordinate expression which will match the induced metric we derived in @sectionInducedMetric. @inducedMetricDefnAbstract tells us that the the pullback $iota^* g$ is a symmetric bilinear map that simply uses the metric on $cal(M)$ to measure the pushforwards $iota_* X$ and $iota_* Y$. Recall, however, that we had previously found $iota_* X$ and $iota_* Y$ to be nothing more than the natural embeddings of $X$ and $Y$ in $T_p cal(M)$. So, in essence, the equation states that $(iota^* g)$ simply returns the same value as $g$ would, if $X$ and $Y$ were to be seen as vectors in $T_p cal(M)$, which they can as $T_p cal(S)$ is embedded in it as a subspace. This narrative of "asking the ambient metric what it thinks and reproducing that" is precisely the reasoning we used in @sectionInducedMetric to motivate its definition---but now, we have rediscovered it in a much more general context, in that of pushforwards and pullbacks associated to smooth maps between manifolds.

Going through the component expressions, we find
$
  (iota^* g)(X,Y) = g(iota_* X, iota_* Y) = g_(mu nu) tensor(E,+mu,-i) tensor(E,+nu,-j) X^i Y^j\
  = (g_(mu nu) tensor(E,+mu,-i) tensor(E,+nu,-j)) X^i X^j = gamma_(i j) X^i Y^j = gamma(X,Y).
$
This confirms explicitly that $iota^* g = gamma$; in other words, the induced metric on $cal(S)$ is simply the pullback of the metric from the ambient manifold $cal(M)$. We have come full circle: the geometric idea that guided our definition of the induced metric has now emerged naturally from algebraic considerations grounded in a broader theoretical framework. That coherence gives us confidence to admit the idea into our formal foundations---when algebra and intuition converge, we are likely on the right path.
== The Left-Inverse: Projections onto Tangent Spaces
We have previously remarked that for a general pushforward
$
  phi_* : T_p cal(M) -> T_phi(p) cal(N)
$
associated with a smooth map $phi: cal(M)-> cal(N)$, it is in general not possible to define an inverse
$
  (phi_*)#h(0em)^(-1):T_(phi(p))cal(N) -> T_p cal(M).
$
However, if $phi_*$ is injective, then it _is_ possible to define a _left-inverse_. A particularly relevant case where this holds is that of immersion maps, such as
$
  phi = iota : cal(S) -> cal(M),
$
where we have previously seen that the associated pushforward 
$
  iota_* : T_p cal(S) -> T_p cal(M)
$
is an injective linear embedding.

Let us now examine this in detail. We may define a (non-unique) left-inverse#footnote[We use the notation $(dot)^(-1)$ here to indicate just a left- and not a proper inverse.] 
$
  (iota_*)#h(0em)^(-1) : T_p cal(M) -> T_p cal(S)
$
as a linear map satisfying
$
  (iota_*)#h(0em)^(-1) compose iota_* = id_(T_p cal(S)).
$<leftInverseCondition>
Though this constrains $(iota_*)#h(0em)^(-1)$ fully on $im(iota_*)$, its action on the complement of $im(iota_*)$ remains arbitrary---hence we have no uniqueness.

 In components, using the same coordiante systems as before, let us expand both sides of this identity for a vector $X in T_p cal(S)$. On the right-hand side, we get
$
  id_(T_p cal(S))X = X = X^i diff_i = delta^i_j X^j diff_i
$
For the left-hand side, we use the component expression for the pushforward to find
$
  ((iota_*)#h(0em)^(-1) compose iota_*)(X) = (iota_*)#h(0em)^(-1) (iota_* X) = (iota_*)#h(0em)^(-1) (tensor(E,+mu,-i)X^i diff_mu) = tensor(((iota_*)#h(0em)^(-1))#h(0em),+j,-mu) tensor(E,+mu,-i) X^i diff_j.
$
where we have defined
$
  tensor(E,+i,-mu) = tensor(((iota_*)#h(0em)^(-1))#h(0em),+i,-mu)
$
as the component representation of the left-inverse.
Equating both sides allows us to express @leftInverseCondition[equation] as
$
  tensor(E,+i,-mu) tensor(E,+mu,-j) X^j diff_i = delta^i_j X^j diff_i.
$
By arbitrariness of $X$, we conclude the matrix identity
$
  tensor(E,+i,-mu) tensor(E,+mu,-j) = delta^i_j.
$<conditionLeftInversePushforward>
This is exactly the condition for $tensor(E,+i,-mu)$ to be a left-inverse of $tensor(E,+mu,-i)$ in the standard sense of matrix algebra. Notice that this implies that both the ranks of $tensor(E,+mu,-i)$ and $tensor(E,+i,-mu)$ must be $s = dim cal(S) = rank delta^i_j$. 

The left-inverse $(iota_*)#h(0em)^(-1)$ might not appear very useful, since all it does is allow us to form an identity on $T_p cal(S)$. We can, however, attempt to use it as a right-inverse, to define the map 
$
  P :=iota_* compose (iota_*)#h(0em)^(-1) : T_p cal(M) -> im (iota_*) subset T_p cal(M).
$
Clearly, the component representation of this map acting on $X in T_p cal(M)$ is
$
  P(X) = (iota_* compose (iota_*)#h(0em)^(-1)) (X^mu diff_mu) = (tensor(E,+mu,-i ) tensor(E,+i,-nu) X^nu) diff_mu = (tensor(P,+mu,-nu) X^nu) diff_mu,
$
with the components of $P$ denoted by $tensor(P,+mu,-nu) := tensor(E,+mu,-i) tensor(E,+i,-nu)$.
Since $P$ is a linear map, and
$
  P^2 = iota_* compose underbrace((iota_*)#h(0em)^(-1) compose iota_*,=id_(T_p cal(S))) compose (iota_*)#h(0em)^(-1) = iota_* compose (iota_*)#h(0em)^(-1) = P,
$
we infer that $P$ is a projection map. Moreover, since we have $rank P = rank iota_* = dim cal(S) = s$, we infer that 
$
  P:T_p cal(M) -> im(iota_*)
$
is a projector of $T_p cal(M)$ onto the embedding $im(iota_*)$ of $T_p cal(S)$ in $T_p cal(M)$. 

Recall that the left-inverse $(iota_*)#h(0em)^(-1)$ is non-unique---it depends on the choice of complement to $im(iota_*)$ during its construction. Without proof, we claim that it is always possible to choose it such that $P$ satisfies the condition
$
  g(P(X),Y) = g(X,P(Y)),quad X,Y in T_p cal(M),
$<orthoCondition>
turning $P$ into an orthogonal projection. 

In words, what $P$ tells us is the following: For any vector $X in T_p cal(M)$, the associated $P(X) in im(iota_*)$ represents the part of $X$ that aligns with the tangent space $T_p cal(S)$ of the submanifold. This construction of $P$ allows us to split any vector $X in T_p cal(M)$ into two parts: a component tangent to the submanifold $cal(S)$, and a component orthogonal to it. That is,
$
  X = P(X) + (X-P(X))
$
where $P(X) in im(iota_*)$ and $X-P(X) in im(iota_*)^perp$.

What we should take away from this section is the following. The left-inverse $(iota_*)#h(0em)^(-1)$ is not particularly interesting on its own, as it is inherently non-unique---its definition depends on an arbitrary choice of complement to $im(iota_*)$. However, it does enable the construction a projection $P:T_p cal(M) -> im(iota_*)$, which can be made into an orthogonal projection by @orthoCondition[condition], equivalent to the symmetry condition#footnote[$g(P(X),Y) = g_(mu nu) tensor(P,+mu,-lambda)X^lambda Y^nu = P_(mu nu) X^mu Y^nu = P_(nu mu) X^mu Y^nu =...$]
$
  P_(mu nu) = P_(nu mu)
$
where $P_(mu nu) = g_(mu lambda) tensor(P,+lambda,-nu)$. This condition ensures that $P$ is self-adjoint with respect to the metric $g$, and hence defines an orthogonal projection.

When $P$ is symmetric in this sense, we can interpret it as specifying an orthogonal splitting of the tangent space,
$
  T_p cal(M) = underbrace(im(iota_*),"tangent to" cal(S)) plus.circle underbrace(ker(P), "normal to" cal(S)),
$<orthonormalDecompositionTpS>
where the symbol $plus.circle$ indicates a direct sum of mutually orthogonal subspaces. The image of $iota_*$ thus spans the tanget space to the submanifold $cal(S)$ (or rather, its embedding in $T_p cal(M)$), while the kernel of $P$ corresponds to vectors in $T_p cal(M)$ orthogonal to $cal(S)$---i.e., the _normal directions_.

== Ambient Metric Decomposition and the Pullback of $P_(mu nu)$ 
In the previous section, we saw that for $P$ to define an _orthogonal_ projection, its components must be symmetric with respect to the metric $g$. Since the metric tensor itself is symmetric, this hints at the possibility of a decomposition of $g$ in terms of $P$ and a complementary projection $Q$. The goal of this section is to examine the properties and consequences of such a decomposition.

The symmetry condition for $P$ was previously written in components as
$
  P_(mu nu) = P_(nu mu), quad "with" quad P_(mu nu) = g_(mu lambda) tensor(P,+lambda,-nu),
$
but this can be expressed more naturally in coordinate-free language, which we now introduce. 

To do so, we reinterpret $P$, originally a $(1,1)$-tensor (a linear map on $T_p cal(M)$) as a $(0,2)$-tensor by lowering one index using the metric. That is, we define the bilinear form $tilde(P)$
$
  tilde(P)(X,Y) = g(P(X),Y),
$
for all $X,Y in T_p cal(M)$. In componnts, this reads
$
  tilde(P)(X,Y) = (g_(mu lambda) tensor(P,+lambda,-nu)) X^mu Y^nu 
$
so the components of $tilde(P)$ are precisely what we earlier denoted $P_(mu nu)$---the object that must be symmetric. Thus, the requirement of symmetry of $P$ turns into 
$
  tilde(P)(X,Y) = tilde(P)(Y,X), quad X,Y in T_p cal(M).
$
We will retain the tilde notation for $tilde(P)$ throughout this section for clarity. However, it should be understood that identifications between $(1,1)$- and $(0,2)$-tensors via the metric are always possible, and we may implicitly make such conversions in later sections, slightly abusing notation for brevity.

Having worked out a coordinate-independent symmetry condition for $P$, let us now approach the decomposition of the metric in terms of $tilde(P)$. We introduce it as
$
  g(X,Y) = tilde(P)(X,Y) + tilde(Q)(X,Y), quad X,Y in T_p cal(M)
$
where $tilde(Q)$, trivially given by
$
  tilde(Q)(X,Y) = g(X,Y) - tilde(P)(X,Y), quad X,Y in T_p cal(M),
$
is necessarily symmetric as well---due to symmetry of both $tilde(P)$ and $g$. The components of $tilde(Q)$ (and its associated $(1,1)$-tensor $Q$) are given by
$
  tilde(Q)_(mu nu) = g_(mu nu) - P_(mu nu) quad <=> quad tensor(Q,+mu,-nu) = g^(mu lambda) tilde(Q)_(lambda nu) = delta^mu_nu - tensor(P,+mu,-nu).
$
Notice that hence, $Q(X) = X- P(X)$ is nothing but the projection of $X$ onto the normal space $ker(P)$ in the @orthonormalDecompositionTpS[decomposition], and $tilde(Q)$ its associated bilinear form. The projections $P$ and $Q$ are orthogonally complete, as we have (this is easily verifiable by plugging in definitions)
$
  P + Q = id_(T_p cal(M)), wide P compose Q = Q compose P = 0.
$
What we learned thus far is the following: Writing the metric in terms of $tilde(P)$ and collecting the remaining parts into $tilde(Q)$ naturally decomposes it into a projection onto the embedding of $T_p cal(S)$ into $T_p cal(M)$ and its orthogonal complement---the _normal space_ $N_p cal(S) := ker(P)$.

Having decomposed the metric, we now consider its pullback onto $cal(S)$. Since we have established $iota^* g = gamma$, i.e., that the ambient metric pulls back to the induced metric, it is natural to ask: which part of $gamma$ arises from $tilde(P)$, and which from $tilde(Q)$? We approach this question from two perspectives: firstly, via a coordinate-independent formulation; then, secondly, by a component-based calculation.

By linearity of $iota^*$, and since $iota^* g = gamma$, it is sufficient to compute $iota^* tilde(P)$---this immediately gives us the pullback of $tilde(Q)$ via $iota^* tilde(Q) = gamma - iota^* tilde(P)$. We compute this now; let $X,Y in T_p cal(S)$---then, by definition,
$
  (iota^* tilde(P))(X,Y) &= tilde(P)(iota_*X,iota_*Y) = g(P(iota_* X), iota_* Y) = g lr((lr((iota_* compose underbrace((iota_*)#h(0em)^(-1) compose iota_*,=id_(T_p cal(S)))),size:#35%) X, iota_* Y ),size:#45%)\
  &= g(iota_* X, iota_* Y) = (iota^* g)(X,Y) = gamma(X,Y);\ \
  ==> wide iota^* tilde(Q) & = gamma-gamma = 0.
$
In summary, we have found
$
  gamma = iota^* g = iota^* tilde(P)
$<indMetricPullbackMetricPullbackProjection>
---hence, the pulbback $iota^* g$ depends only on the tangential part $tilde(P)$---the normal contribution $tilde(Q)$ vanishes under pullback. Notice that $rank P = rank gamma = dim cal(S)$; in light of the above, this implies that $gamma$ and $tilde(P)$ represent the same bilinear map on $T_p cal(S)$ and its embedding $im(iota_*) subset T_p cal(M)$. 

Let us now repeat the same calculation in components. We find
$
  (iota^* tilde(P))#h(0em)_(i j) &= tensor(E,+mu,-i) tensor(E,+nu,-j) P_(mu nu) = tensor(E,+mu,-i) tensor(E,+nu,-j) g_(mu lambda) tensor(P,+lambda,-nu) =  tensor(E,+mu,-i) tensor(E,+nu,-j) g_(mu lambda) tensor(E,+lambda,-k) tensor(E,+k,-nu)\
  &= tensor(E,+mu,-i) delta^k_j tensor(E,+lambda,-k) g_(mu lambda) = tensor(E,+mu,-i) tensor(E,+nu,-j) g_(mu nu) = gamma_(i j).
$
This reproduces the @indMetricPullbackMetricPullbackProjection[identity] in terms of components as
$
  gamma_(i j) = tensor(E,+mu,-i) tensor(E,+nu,-j) g_(mu nu) =  tensor(E,+mu,-i) tensor(E,+nu,-j) P_(mu nu).
$
This result offers a clear structural interpretation of the induced metric $gamma$: it is precisely the tangential part of the ambient metric $g$, isolated by projection through $P$ and realised by the pullback. The normal component $Q$ plays no role in the geometry intrinsic to the submanifold, as expected---$gamma$ contains only the information relevant to distances and angles _within_ $cal(S)$. What initially appeared as a simple construction now revelas itself as a direct manifestation of the geometry of orthogonal decomposition.
= Bundles
This section aims to introduce the notion of various types of _bundles_ one can define on a smooth manifold. Though this is not strictly necessary to study submanifolds, it seems like it would be a useful digression to prepare for the differential geometry lecture in Part III, so I will go over it briefly here.
== Vector Bundles: Intuition and Definitions
So far, when we discussed vectorial (or tensorial) objects, our expressions have been entirely pointwise. We have considered, for example, maps from $T_p cal(M)$ to $T_phi(p) cal(N)$, where $cal(N)$ is another manifold and $phi: cal(M)-> cal(N)$ is a smooth map between them. These are relations between tangent spaces at individual points.

However, since there exists a tangent space at every point $p in cal(M)$, it is natural to seek a way to _assemble_ or _bundle together_ all these tangent spaces into a single structure. Intuitively, we take the manifold $cal(M)$ and, at each point, attach the tangent space at that point. This yields a new, higher-dimensional manifold-like object that encodes all the tangent spaces and their relation to points in $cal(M)$. This construction is known as the _tangent bundle_.

Of course, this is a loose and purely intuitive description. In this section, we aim to make rigorous the idea of "attaching a vector space to each point of a manifold" by introducing the notion of a _vector bundle_. In the next section, we will see that the tangent (and cotangent) bundles are special instances of this general concept.
\ \
*Definition* (Vector Bundle) Let $cal(M)$ be a smooth manifold of dimension $m$. A _smooth real vector bundle_ of rank $n$ over $cal(M)$ is a triple $(cal(E), pi, cal(M))$, where:

+ $cal(E)$ is a smooth manifold of dimension $n+m$, called the _total space_.

+ $pi: cal(E) -> cal(M)$ is a smooth surjective map, called the _bundle projection_.

+ For each $p in cal(M)$, the _fibre_ $cal(E)_p := pi^(-1)({p})$ is equipped with the structure of a real vector space of dimension $n$.

+ We have _local triviality_, i.e. for each $p in cal(M)$ there exists an open neighbourhood $U subset cal(M)$ of $p$ and a diffeomorphism
  $
    Phi: pi^(-1) (U) -> U times RR^n
  $
  such that
    - $Phi$ is a fibre-preserving vector space isomorphism on fibres, meaning that for each $q in U$, 
    $
      Phi|_(cal(E)_q) : cal(E)_q -> {q} times RR^n tilde.eq RR^n
    $
    is a vector space isomorphism.

    - the diagram
    $
      #diagram({
        node((-1, -2), [$cal(E) supset pi^(-1)(U)$])
        node((1, -2), [$U times RR^n$])
        node((-1, 0), [$U$])
        edge((-1, -2), (-1, 0), [$pi$], label-side: right, "->")
        edge((-1, -2), (1, -2), [$Phi$], label-side: left, "->")
        edge((1, -2), (-1, 0), [$P_1$], label-side: left, "->")
      })
    $<fibrationDiagram>
    commutes, where $P_1$ is the projection onto the first component of the Cartesian product ($(a,b)|->a$)
    
Let us now go through this definition calmly, and explain the meaning and intuition behind each of the constructions separately.

Firstly, we should give a summary of what a vector bundle is supposed to be. Intuitively, a vector bundle over $cal(M)$ is a smooth family of vector spaces ${cal(E)_p}#h(0em)_(p in cal(M))$, smoothly parameterised by $cal(M)$, such that near each point, the collection of fibres looks like a product $U times RR^n$. This means we can locally identify each fibre with $RR^n$ in a way that varies smoothly with the base point, and respects the vector space operations.

Now, let us go over each part of the definition in detail.

+ _Total Space_: The total space $cal(E)$ can be viewed as the union of all the vector spaces at each point. Formally, it is the disjoint union of all fibres,
  $
    cal(E) = union.sq.big_(p in cal(M)) cal(E)_p.
  $
  This is an $(m+n)$-dimensional manifold, since the base manifold has dimension $cal(M)$, and attaching an $n$-dimensional vector space at each point increases the dimension by $n$.
  
  Very loosely, one might imagine this as dragging a window (the vector space) across a screen (the manifold) in Windows XP, with the bug that leaves behind a smeared trail of it. As you drag it, the trail being formed represents this union of all these smeared copies (counting overlapping regions appropriately), which is what is known as the total space.

+ _Bundle Projection_: A point in $cal(E)$ encodes two pieces of information: a point $p in cal(M)$, and a vector attached to it, describing a _direction from $p$_. The bundle projection's job is simple: it takes such a pair and tells us to which point on the base manifold it belongs.

  It must be surjective, so that each point $p in cal(M)$ has a corresponding vector space attached to it. If it were not surjective, there would exist some $p$ for which no vectors are associated, contradicting the idea of "attaching a vector space at each point". 

+ _Fibres_: The fibre at $p in cal(M)$ is the collection of all vectors attached or associated to $p$, i.e. its preimage under the bundle projection $pi$. More precisely,
  $
    cal(E)_p = pi^(-1)({p}).
  $
  Since we want to attach not just _any_ kind of fibre, but specifically an $n$-dimensional vector space, we impose the additional requirement that each fibre $cal(E)_p$ carries the structure of a real vector space of dimension $n$.

+ _Local Triviality_: This is likely the most convoluted part of the definition, but can also be broken down intuitively. What local triviality demands is that _locally_, in some neighbourhood $U subset cal(M)$ of $p in cal(M)$, the total space $cal(E)$ "looks like" the space $U times RR^n$. This is the simplest way of "attaching vector spaces to each point"---the Cartesian product does exactly that. In more formal terms, "looks like" is replaced by the notion of the diffeomorphism $Phi$. Since $RR^n$ is the concrete representation of the attached vector space, we would also like the fibres $cal(E)_q$, $q in U$ to map to $RR^n$ under $Phi$ in a way that respects the algebraic structure---hence the condition on $Phi|_cal(E)_q$. The requirement that the @fibrationDiagram[diagram] commutes then further ensures that the fibres get attached to the correct points on $cal(M)$. 


Now that we have defined vector bundles, let us introduce a notion that makes use of it. Specifically, we consider so-called _smooth sections_ of vector bundles. A smooth section is, intuitively speaking, the selection of one vector in the fibre $cal(E)_p$ at each $p in cal(M)$, in a way that creates a smooth surface in the total space $cal(E)$. Such a surface can be viewed as a vector field, since it maps each point on the manifold to one vector in its fibre.
\ \ 
*Definition* (Smooth Section) Let $cal(M)$ be a smooth manifold and $(cal(E),pi,cal(M))$ a smooth real vector bundle over $cal(M)$. A _smooth section_ of $cal(E)$ is a smooth map 
$
  sigma : cal(M) -> cal(E)
$
such that
$
  pi compose sigma = id_cal(M)
$
In other words, for each $p in cal(M)$, the map $sigma$ selects a vector $sigma(p) in cal(E)_p$ lying in the fibre over $p$, and this assignment varies smoothly with $p$. In particular, for $pi compose sigma$ to be the identity on $cal(M)$, $sigma$ must be injective. 

We denote the set of all smooth sections of $cal(E)$ by
$
  Gamma(cal(E)) &= {sigma:cal(M) -> cal(E) | pi compose sigma = id_cal(M), sigma "smooth"},\
  &= {sigma:cal(M)->cal(E) | sigma "smooth section on" cal(E)}.
$
This is a real vector space under pointwise addition and scalar multiplication. In the case where $cal(E) = T cal(M)$ is the tangent bundle (which we introduce in the next section), $Gamma(T cal(M))$ is the space of smooth vector fields on $cal(M)$.

We are now prepared for the next section, in which we will define the tangent and cotangent bundles $T cal(M)$ and $T^* cal(M)$, and use smooth sections to give an alternative perspective on vector fields and differential 1-forms.
== The Tangent and Cotangent Bundles
In this section, we introduce the two most important vector bundles in differential geometry: the tangent bundle $T cal(M)$ and the cotangent bundle $T^* cal(M)$. In essence, these are the special cases where one chooses the fibres $cal(E)_p$ of a vector bundle to be the (co)-tangent spaces $T_p cal(M)$ and $T_p^* cal(M)$, respectively---but let us now introduce this rigorously.
\ \
*Definition* (Tangent Bundle)
Let $cal(M)$ be a smooth manifold of dimension $m$ and denote by $T_p cal(M)$ its tangent space at any point $p in cal(M)$. Define the $2m$-dimensional total space $T cal(M)$ by the disjoint union
$
  T cal(M) = union.sq.big_(p in cal(M)) T_p cal(M).
$
We define the bundle projection
$
  pi:T cal(M) -> cal(M)
$
as the map 
$
  pi(p,v) = p,
$
choosing to write elements of $T cal(M)$ as pairs $(p,v)$ with $v in T_p cal(M)$. 
\ \
*Remarks:* A few direct consequences follow from this construction:
- _Fibres_: The fibre over a point $p in cal(M)$ is given by
  $
    (T cal(M))_p = pi^(-1) ({p}) = T_p cal(M),
  $
  meaning that we are indeed attaching the tangent space $T_p cal(M)$ to each point $p$. 

- _Local triviality and charts on $T cal(M)$_: Given a coordinate chart $(U,x^mu)$ on $cal(M)$, we can define a chart on the preimage $pi^(-1)(U) subset T cal(M)$ by identifying
  $
    pi^(-1)(U) tilde.eq U times RR^m, quad (p,V)|-> (x^mu(p), V^mu),
  $
  where $V^mu$ are the components of the tangent vector $V = V^mu diff_mu$ in the coordinate basis induced by $x^mu$. This provides $T cal(M)$ with a smooth manifold structure of dimension $2m$.

- _Zero section_: We may define the zero section
  $
    sigma_0 : cal(M) -> T cal(M), quad p |-> (p,0),
  $
  or, in coordinates, 
  $
    x^mu |-> (x^mu,0),
  $
  where $0$ denotes the zero vector in $T_p cal(M)$. This map is a smooth embedding, and its image forms a submanifold of $T cal(M)$ diffeomorphic to $cal(M)$ itself. 

  Intuitively, this makes perfect sense: if we attach a vector space to each point of the manifold in such a way that their origins coincide with the base points, then the collection of those origins sweeps out a copy of $cal(M)$ embedded inside $T cal(M)$. 

- _Vector fields as section_: Any smooth vector field on $cal(M)$, i.e. a smooth assignment $X:p|-> X_p in T_p cal(M)$, is precisely a smooth section of the tangent bundle,
  $
    X : cal(M) -> T cal(M), quad "with" quad pi compose X = id_cal(M).
  $
  That is, $X in Gamma (T cal(M))$. This is the prototypical example of a smooth section, and shows how the familiar notion of a vector field fits directly into the general formalism of vector bundles.

*Definition* (Cotangent Bundle) The definition of the cotangent bundle is analogous to that of the tangent bundle, with the distinction being that we use
$
  T^* cal(M) = union.sq.big_(p in cal(M)) T_p^* cal(M)
$
as our total space, and the bundle projection 
$
  pi: T^* cal(M) -> cal(M)
$
now given by the map
$
  pi(p,omega) = p,
$
writing elements of $T^* cal(M)$ as pairs $(p, omega)$ with $omega in T_p^* cal(M)$. The same remarks as for the tangent bundle hold, with the occasional change of terminology from vectors to $1$-forms.
== Bundle Maps and Vector Bundle Morphisms
Now that we have introduced the notion of vector bundles as well as concrete (and important) examples thereof, the tangent bundle $T cal(M)$ and the cotangent bundle $T^* cal(M)$, we can begin considering maps between them. A map between two tangent bundles becomes particularly interesting if it respects the algebraic structure of the fibres---i.e., if it maps vectors in one fibre linearly to vectors in another. Such maps we will refer to as _bundle morphisms_. We now first introduce the more general concept of a bundle map, and then impose the additional algebraic structure to define the notion of bundle morphisms.
\ \
*Definition* (Bundle Map) Let $cal(M)$, $cal(N)$ be smooth manifolds and $(cal(E),pi_cal(M), cal(M))$ as well as $(cal(F), pi_cal(N), cal(N))$ vector bundles on them. Further, let 
$
  phi: cal(M) -> cal(N)
$
be a smooth map from one manifold to the other. We call a map
$
  Phi : cal(E) -> cal(F)
$
a _bundle map_ if the following diagram commutes:
$
  #diagram({
    node((0, 0), [$cal(E)$])
    node((1.5, 0), [$cal(F)$])
    node((0, .75), [$cal(M)$])
    node((1.5, .75), [$cal(N)$])
    
    edge((0, 0), (1.5, 0), [$Phi$], label-side: left, "->") 
    edge((0, .75), (1.5, .75), [$phi$], label-side: right, "->") 
    edge((0,0), (0,.75), [$pi_cal(M)$], "->")
    edge((1.5,0), (1.5,.75), [$pi_cal(N)$], label-side: left, "->")
  })
$
This is equivalent to the condition that
$
  pi_cal(N) compose Phi = phi compose pi_cal(M).
$<bundleMapCondition>
Without further interpretation, the definition above may seem opaque---so let us walk through it and develop some intuition. Consider an element $(p, X) in cal(E)$, where $p in cal(M)$ and $X in cal(E)_p$, the fibre over $p$. We examine the two sides of the @bundleMapCondition[condition] act on such a point.

Starting with the right-hand side, we have
$
  (phi compose pi_cal(M))(p,X) = phi(pi_cal(M)(p,X)) = phi(p).
$
his means ve first project $(p,X)$ onto its base point $p in cal(M)$ via $pi_cal(M)$, then apply $phi$ to obtain a point in $cal(N)$. 

Now for the left-hand side:
$
  (pi_cal(N) compose Phi)(p,X) = pi_cal(N) (Phi(p,X)) = pi_cal(N) (q, Y) = q.
$
Here, we write $Phi(p,X)$ as some pair $(q,Y) in cal(F)$, and $pi_cal(N)$ returns the base point $q in cal(N)$. 

Thus, the commutativity condition simplifies to
$
  q = phi(p).
$
In other words, a bundle map $Phi$ must map $(p,X) in cal(E)$ to $(phi(p),Y) in cal(F)$; the base point of the image is determined entirely by the underlying map $phi$. The vector part $Y$ can be chosen freely (within the fibre over $phi(p)$), but the association between fibres is rigidly tied to that between base points.

A third, equivalent way to phrase this is to say that $Phi$ maps each fibre $cal(E)_p$ into the fibre $cal(F)_(phi(p))$. That is,
$
  Phi|_cal(E)_p : cal(E)_p -> cal(F)_phi(p).
$
The content of the definition is precisely this: a bundle map over $phi$ is one that maps vectors attached to a point $p in cal(M)$ into vectors attached to $phi(p) in cal(N)$, without violating the structure of the fibration.

Besides the structure of the fibration, a vector bundle has additional algebraic structure that one could demand to be preserved; each fibre is a vector space, and we could demand a bundle map to respect it by imposing linearity. This notion, called _bundle morphisms_, is what we now define rigorously.
\ \
*Definition* (Vector Bundle Morphism) Let $cal(M)$, $cal(N)$ be smooth manifolds, and let $(cal(E), pi_cal(M), cal(M))$ and $(cal(F), pi_cal(N), cal(N))$ be smooth real vector bundles over them. Suppose we are given a smooth map $phi:cal(M)->cal(N)$ and a bundle map $Phi:cal(E)->cal(F)$ covering $phi$, i.e.
$
  pi_cal(N) compose Phi = phi compose pi_cal(M).
$
We say that $Phi$ is a _vector bundle morphism_ (over $phi$) if, for each $p in cal(M)$, the induced map on the fibres
$
  Phi|_(cal(E)_p) : cal(E)_p -> cal(F)_phi(p)
$
is a linear map of vector spaces. 

Intuitively, a vector bundle morphism over $phi$ can be seen as a fibrewise linear transformation that "respects the base": it transforms each vector in a fibre over $p in cal(M)$ to a vector in the fibre over $phi(p) in cal(N)$, via a linear map.

While there are many abstract vector bundle morphisms on could define and study, there is a particularly natural one associated with a smooth map $phi: cal(M)->cal(N)$ that we have already encountered: the pushforward. Though we initially introduced the pushforward $phi_* : T_p cal(M) -> T_phi(p) cal(N)$ as a pointwise linear map between tangent spaces, it readily extends to a global vector bundle morphism.

To see this, consider the map
$
  Phi : T cal(M) -> T cal(N), quad (p,X) |-> (phi(p), phi_* X), quad p in cal(M), X in T_p cal(M).
$
This construction maps each element of the tangent bundle $T cal(M)$ to the tangent bundle $T cal(N)$ by pushing forward the vector $X$ and sending its base point $p$ to $phi(p)$. It is easy to check that this satisfies the bundle map condition
$
  pi_cal(N) compose Phi = phi compose pi_cal(M),
$
and that the fibrewise maps $Phi|_(T_p cal(M)) = phi_* : T_p cal(M) -> T_phi(p) cal(N)$ are linear. Hence, $Phi$ defines a (vector) bundle morphism from $T cal(M)$ to $T cal(N)$ over $phi$.
== The Normal Bundle and Orthogonal Decomposition

In this section, we explore how the tangent bundle $T cal(S)$ of a submanifold $cal(S) subset cal(M)$ can be understood as a subbundle of the restriction of the tangent bundle $T cal(M)$ to $cal(S)$. Further, we define the normal bundle $N cal(S)$ and explain how the pointwise projections $P$ and $Q$, as introduced earlier, extend naturally to smooth vector bundle morphisms between $T cal(M)|_cal(S)$, $T cal(S)$ and $N cal(S)$. While this may seem like an unnecessary abstraction at first glance, it will turn out to offer geometric clarity and prepare us for later constructions involving intrinsic curvature. 
\ \ 
*Definition* (Restriction of vector bundles to submanifolds) Let $cal(M)$ be a smooth manifold of dimension $m$, and let $cal(S) subset cal(M)$ be a submanifold of dimension $s < m$. Given a vector bundle $(cal(E), pi, cal(M))$ over $cal(M)$, its _restriction to_ $cal(S)$ is defined as the triple
$
  (cal(E)|_cal(S), pi|_cal(S), cal(S)), quad "where" cal(E)|_cal(S) := pi^(-1)(cal(S)) subset cal(E).
$
This construction simply discards all fibres of $cal(E)$ lying over points $p in cal(M) without cal(S)$, retaining only the portion of the bundle sitting above $cal(S)$. 

As a key example, consider the tangent bundle $T cal(M)$. Its restriction to $cal(S)$, denoted by $T cal(M)|_cal(S)$, consists of the collection of tangent spaces $T_p cal(M)$ for $p in cal(S)$. This restricted bundle is not the same as the tangent bundle of $cal(S)$---the fibres of $T cal(M)|_cal(S)$ are $m$-dimensional, while the fibres of $T cal(S)$ are only $s$-dimensional. 

However, the inclusion map $iota : cal(S) -> cal(M)$ induces a smooth injective bundle morphism
$
  iota_* : T cal(S) -> T cal(M)|_cal(S)
$
which embeds each fibre $T_p cal(S) subset T_p cal(M)$ for $p in cal(S)$. As a result, $T cal(S)$ is realised as a smooth vector subbundle of $T cal(M)|_cal(S)$. This perspective will be useful not only for defining and understanding the normal bundle but also for lifting previously pointwise constructions---such as the projections $P$ and $Q$---into global morphisms between bundles. 

Recall that for any submanifold $cal(S) subset cal(M)$ of a (pseudo-)Riemannian manifold $cal(M)$, we obtain a canonical orthogonal decomposition of the tangent space at each point $p in cal(S)$,
$
  T_p cal(M) = T_p cal(S) plus.circle N_p  cal(S),
$
where $T_p cal(S) subset T_p cal(M)$ is identified via the pushforward $iota_*$ of the inclusion $iota:cal(S)->cal(M)$, and the normal space $N_p cal(S)$ is defined as the orthogonal complement
$
  N_p cal(S) = (T_p cal(S))^perp subset T_p cal(M).
$
Note that $dim N_p cal(S) = codim_cal(M) cal(S)$, so that $dim T_p cal(S) + dim N_p cal(S) = dim T_p cal(M)$. 

Just as bundling the tangent spaces $T_p cal(S)$ yields the tangent bundle $T cal(S)$, we may bundle the normal spaces to obtain a smooth vector bundle over $cal(S)$. We now define this more precisely.
\ \ 
*Definition* (Normal bundle) Let $cal(S)subset cal(M)$ be an embedded submanifold of a (pseudo-)Riemannian manifold $cal(M)$. The _normal bundle_ $N cal(S)$ is the smooth vector bundle over $cal(S)$ defined by 
$
  N cal(S) := union.sq.big_(p in cal(S)) N_p cal(S), wide pi: N cal(S) -> cal(S), quad (p,X)|-> p,
$
choosing the normal spaces at each $p in cal(S)$ as the fibres.

Since the decomposition 
$
  T_p cal(M) = T_p cal(S) plus.circle N_p cal(S)
$<orthoDecompTpM>
is orthogonal and varies smoothly with $p$, the total bundle $T cal(M)|_cal(S)$ likewise decomposes as a direct sum of vector bundles,
$
  T cal(M)|_cal(S) = T cal(S) plus.circle N cal(S).
$<tangentBundleOrthoSplitting>
Here, the direct sum $plus.circle$ is understood as the fibrewise orthogonal sum within each $T_p cal(M)$. 

Let us make a final addition to the bundle perspective we have been building in this section. Specifically, let us consider how the orthogonal projections
$
  &P : T_p cal(M) -> im(P) = T_p cal(S) quad "and"\ 
  &Q: T_p cal(M) -> im(Q) = ker(P) = N_p cal(S)
$
are lifted from their pointwise definitions to smooth vector bundle morphisms. Since these projections were defined for each $p in cal(S)$ using the @orthoDecompTpM[orthogonal decomposition] and since this splitting varies smoothly over $cal(S)$, we obtain globally defined, smooth bundle morphisms 
$
  P : T cal(M)|_cal(S) -> T cal(S), quad Q : T cal(M)|_cal(S) -> N cal(S),
$
which act fibrewise as the orthogonal projections onto $T_p cal(S)$ and $N_p cal(S)$, respectively. Both $P$ and $Q$ are idempotent ($P^2 = P$, $Q^2 = Q$) and satisfy
$
  ker P = im Q quad "and" quad ker Q = im P,
$
pointwise. This allows us to view the @tangentBundleOrthoSplitting[splitting] not merely as a statement about individual tangent spaces, but as a decomposition of vector bundles, mediated by smooth projections. 
= Foliations
== Motivation and Definition of Foliations
In the previous sections, we introduced and examined submanifolds of (pseudo-)Riemannian manifolds in detail, including the tensorial structures they inherit, such as the induced metric. Rather than focusing on isolated submanifolds, we now turn our attention to decompositions of an entier manifold $cal(M)$ (or open subsets thereof) into a smooth family of non-intersecting submanifolds. These decompositions, referred to as _foliations_, naturally arise in a variety of contexts. For example, in the ADM formalism, foliating spacetime into spacelike hypersurfaces facilitates isolating a direction of dynamical evolution---typically a timelike one. In other settings, such as the study of flows or congruences, one may be interested in decompositions into curves or integral lines.

Either way, the concept of splitting a manifold into lower-dimensional submanifolds is fundamental, and we now develop it in the form of such foliations. In this section, we first provide the definition of foliations of arbitrary codimension, as well as coordinates adapted to them. We then specialise to the codimension-1 case of _hypersurface foliations_, which are particularly relevant in the ADM formalism and for which certain equations and identities take on a simpler form.
\ \
*Definition* (Foliation) Let $cal(M)$ be a smooth manifold. A _foliation of codimension $k$_ is a $k$-parameter family $lr({Sigma_(t_0)},size:#80%)#h(0em)_(t_0 in RR^k)$ of smooth, embedded submanifolds $Sigma_(t_0) subset cal(M)$ of codimension $k$ such that 
$
  Sigma_(t_0) sect Sigma_(t'_0) = nothing quad "for" quad t_0 != t'_0, quad "and" quad cal(M) = union.big_(t_0 in RR^k) Sigma_(t_0).
$
The submanifolds $Sigma_(t_0)$ are referred to as the _leaves_ of the foliation.

This definition, while clean, hides a more powerful and flexible characterisation. Since the leaves are disjoint and cover all of $cal(M)$, each point $p in cal(M)$ lies in a unique leaf $Sigma_(t_0^A)$. We can therefore associate to each point its corresponding label $t_0 = (t_0^A)$, $A = 1,...,k$, giving rise to a map
$
  t^A : cal(M) -> RR^k, quad t(p) = t_0 "such that" p in Sigma_(t_0).
$
This defines $k$ smooth scalar fields on $cal(M)$, and the leaves of the foliation may then be expressed as the family of level sets
$
  Sigma_(t_0) = t^(-1)(t_0) =  sect.big_(A=1)^k (t^A)^(-1)(t_0).
$
In order for each $Sigma_(t_0)$ to be a smooth submanifold, we require $t_0$ to be a regular value of all component maps $t^A$ of $t = (t^A)$. For this to hold for all $t_0$, we demand each $t^A$ be a _submersion_, i.e.
$
  dt^A != 0 quad "everywhere", quad forall A = 1,...,k.
$
We thus arrive at an equivalent perspective: a codimension-$k$ foliation of $cal(M)$ may be defined by a set of $k$ (functionally independent, i.e. ${dt^A}$ is linearly independent) scalar fields $t = (t^A) : cal(M) -> RR^k$, $A = 1,...,k$, whose differentials are nowhere vanishing. The intersection of their level sets then define the leaves of the foliation.

Note that because $dt^A != 0$, the map $t = (t^A)$ can be extended to a coordinate chart on $cal(M)$, where the coordinates $(t^A,y^i)$ describe both the foliation parameter and the coordinates $(y^i)$ on the leaves. The number of transverse coordiantes $y^i$ is given by $dim Sigma_t = dim cal(M) - k$. Such coordinates are called _weakly adapted_ to the foliation. In particular, fixing $t^A = t^A_0$ to some constant value $t^A_0$ yields a coordinate chart $(y^i)$ on the leaf $Sigma_(t_0)$. For this reason, the $y^i$ are referred to as _transverse coordinates_.
\ \
*Definition* (Hypersurface Foliation) We call a foliation $Sigma = {Sigma_t}$ of codimension $1$ a _hypersurface foliation_. 

A hypersurface foliation is defined by a single scalar field
$
  t : cal(M) -> RR, quad dt !=0 "everywhere",
$
which can be extended to a coordinate chart as $(t,y^i)$, $i= 1,...,dim cal(M)-1$. 

== The Normal 1-Form and Normal Vector Field <sectionFoliationTangentBundleDecomposition>

The choice of a hypersurface foliation of a manifold $cal(M)$ naturally gives rise to a host of associated mathematical objects, each playing its own role in its geometry. Among the most central of these is the normal vector field, together with its corresponding one-form. In what follows, we shall first build some intuition for why such an object may be constructed, before proceeding to define it rigorously. 

As we have seen previously, a submanifold $cal(S) subset cal(M)$ induces, at each point $p in cal(S)$, a decomposition
$
  T_p cal(M) = T_p cal(S) plus.circle N_p cal(S),
$
splitting the tangent space of $cal(M)$ into components tangent and normal to $cal(S)$. However, this decomposition is inherently local to $cal(S)$; at points $p in cal(M)$ outside of $cal(S)$, the notions $T_p cal(S)$ and $N_p cal(S)$ have no meaning. Hence one cannot define a vector field on $cal(M)$ with the property that it is normal to $cal(S)$ everywhere---this statement simply does not make sense away from $cal(S)$. 

A foliation $Sigma = {Sigma_t}#h(0em)_(t in RR)$ changes this picture dramatically. Given such a foliation, every point $p in cal(M)$ belongs to a unique leaf $Sigma_t(p)$, for a unique $t(p) in RR$. Consequently, we may write
$
  T_p cal(M) = T_p Sigma_t(p) plus.circle N_p Sigma_t(p), 
$<foliationTpMDecomp>
thereby defining a decomposition into tangent and normal directions at _every_ point of $cal(M)$. In other words, a foliation allows us to speak globally about directions that are tangent or normal to the slices $Sigma_t$, since every point of $cal(M)$ lies on precisely one such submanifold.

At this point it is natural to briefly revisit the notion of vector bundles. Using the @foliationTpMDecomp[decomposition], we may define two vector bundles over $cal(M)$ that together decompose the tangent bundle $T cal(M)$. Namely,
$
  T Sigma = union.big.sq_(p in cal(M)) T_p Sigma_t(p), wide N Sigma = union.big.sq_(p in cal(M)) N_p Sigma_t(p),
$
equipped with their canonical projections onto $cal(M)$, yield the tangent bundle and normal bundle of the foliation, respectively. It is then immediate that the tangent bundle of the manifold decomposes as
$
  T cal(M) = T Sigma plus.circle N Sigma.
$<foliationTangentBundleDecomp>
We may now also formalise the notion of a _normal vector field_. Recall that a vector field on $cal(M)$ is a smooth section of the tangent bundle $T cal(M)$. Given the @foliationTangentBundleDecomp[decomposition] we call a vector field _tangent_ to the foliation $Sigma$ if it is a smooth section of $T Sigma$, and normal to $Sigma$ if it is a smooth section of $N Sigma$. In more elementray terms, a vector field $X$ is tangent (respectively normal), if at each point $p in cal(M)$, the vector $X(p)$ lies in the subspace $T_p Sigma_t$ (respectively, $N_p Sigma_t$) associated to the unique slice $Sigma_t$ containing $p$. 

In our specific case of a codimension-one foliation, the normal spaces $N_p Sigma_t(p)$ are one-dimensional at every point. This observation is powerful. it implies that any vector field normal to the foliation can be written as a scalar multiple of a single, globally defined, nowhere-vanishing basis vector field. Concretely, we may express any such field $X$ as
$
  X = lambda thin n^sharp,
$
where $n^sharp$ denotes a chosen normal vector field, and $lambda$ is a smooth scalar function on $cal(M)$. The use of the musical isomorphism in this notation is deliberate: rather than constructing $n^sharp$ directly, it is often more natural to begin with a normal one-form $n$, and then obtain the vector field $n^sharp$ via the defining identity
$
  g(n^sharp, X) = n(X) quad forall X quad <=> quad (n^sharp)^mu = n^mu = g^(mu nu) n_nu.
$<sharpDefinition>
If the normal subspaces $N_p Sigma_t(p)$ are timelike (or spacelike) everywhere, then we may impose a canonical normalisation to fix $n^sharp$ uniquely (up to sign),
$
  g(n^sharp,n^sharp) = cases(-1\,quad&"timelike"\,,+1\,quad&"spacelike".)
$
This condition pins down a distinguished unit normal vector field. However, we have not yet addressed how to construct such a form or field in practice.

Let us now turn to this task. For a vector field $n^sharp$ to be normal to the foliation $Sigma$, it must satisfy
$
  g(n^sharp, X) = 0, quad forall X in T_p Sigma_t(p)
$
at every point $p in cal(M)$; that is, $n^sharp$ must be orthogonal to all vectors tangent to the slice $Sigma_t$ through $p$. By the defining property of the musical isomorphism, this is equivalent to requiring
$
  n(X) = 0 quad forall X in T_p Sigma_t(p)
$
i.e. the 1-form $n$ must annihilate all tangent vectors to the foliation.

This insight significantly simplifies the problem in adapted coordinates. Let $x^mu = (t,y^i)$ be a coordinate chart adapted to the foliation, so that each slice $Sigma_t$ is locally given by $t = const$, and the $y^i$ serve as coordinates on each leaf. Then the tangent space to $Sigma_t$ at $p$ is spanned by the coordinate basis vectors
$
  T_p Sigma_t(p) = span{diff_i = diff/(diff y^i) mid(|) i = 1,...,dim cal(M)-1}.
$
We thus seek a 1-form $n$ satisfying
$
  n(diff_i) = 0 quad i = 1,...,dim cal(M) -1.
$<eq1.4.16>
The standard coordinate duality relation
$
  dx^mu (diff_nu) = delta^mu_nu
$
immediately suggests a solution: the 1-form
$
  n = alpha dt
$ 
clearly satisfies @eq1.4.16[condition], as
$
  n(diff_i) = alpha dt(diff_i) = delta^t_i = 0.
$
Here, $alpha$ is a smooth scalar field on $cal(M)$, to be determined by a normalisation condition. If we require $n^sharp$ to be unit-normalised, this yields the relation
$
  pm 1 = g(n^sharp,n^sharp) = g^(mu nu) n_mu n_nu = alpha^2 g^(t t) quad <=> quad g^(t t) = (pm 1)/alpha^2.
$
Moreover, the components of $n^sharp$ can be written in terms of the metric as well,
$
  (n^sharp)#h(0em)^mu = n^mu = g^(mu nu)n_nu = alpha g^(mu nu) delta_nu^t= alpha g^(mu t).
$
In hindsight, it is only natural that the function $t:cal(M)->RR$ which defines the foliation $Sigma$ plays a central role in the construction of the normal form. After all, moving along the gradient of $t$ corresponds to moving between leaves, so it makes sense that for an arbitrary set of coordinates $macron(x)^alpha$,
$
  dt = diff_alpha t thin d macron(x)^alpha
$
somehow encodes the normal direction.

This construction of $n$ as being proportional to $dt$ was rather simple in coordinates. Keeping with the spirit of these notes, though, we should also show the equivalence
$
  dt_p (X) = 0 quad <=> quad X in T_p Sigma_t(p)
$<equivalenceCoordIndep>
in a coordinate-independent way. This is not that hard, and gives some additional geometric intuition as it makes use of the definition of tangent spaces in terms of derivatives along curves. 

_Proof_: (of @equivalenceCoordIndep[eq.]) We first prove the implication "$arrow.l.double$". Let $X in T_p Sigma_t(p)$. Then there exists a smooth curve $gamma:(-epsilon,epsilon)->Sigma_t(p) subset cal(M)$ with $gamma(0) = p$ and $dot(gamma)(0) = X$, such that for any $f in C^infty (cal(M))$,
$
  X[f] = d/(ds) f(gamma(s))|_(s=0).
$
In particular,
$
  dt(X) = X[t] = d/ds t(gamma(s)) = d/ds t(p) = 0,
$
since $gamma(s) subset Sigma_t = t^(-1)(t(p))$ implies $t(gamma(s)) = t(p)$ identically (and hence is constant).

For the converse, note that $dt != 0$ everywhere by construction, and $codim T_p Sigma_t (p) = 1$. Hence, for any vector $X in T_p cal(M) without T_p Sigma_t(p)$, we must have $dt(X) != 0$, otherwise $dt_p$ would vanish on all of $T_p cal(M)$, contradicting $dt != 0$. This implies
$
  T_p Sigma_t(p) = ker(dt_p),
$<pointwiseIdentificationkerdt>
which establishes the equivalence. #h(1fr)#Box

We can even take the pointwise @pointwiseIdentificationkerdt[identification] one step further. To this end, let us define the notion of a kernel of a differential 1-form $omega in Gamma(T^* cal(M))$. At any point $p$, we have
$
  ker(omega_p) = {X in T_p cal(M) | omega_p (X) = 0},
$
the standard definition of the kernel from linear algebra. We can extend this to the entire differential form $omega$ by attaching to each point in $cal(M)$ the corresponding subspace $ker(omega_p) subset T_p cal(M)$, yielding the subbundle
$
  ker(omega) := union.big.sq_(p in cal(M)) ker(omega_p) subset T cal(M),
$
together with the canonical projection map onto $cal(M)$. It is now straightforward to see that when choosing $omega = dt$, due to the @pointwiseIdentificationkerdt[identification] we find
$
  ker(dt) = T Sigma subset T cal(M).
$ 
The tangent bundle of a foliation $Sigma$ generated by the level sets of a scalar function $t : cal(M) -> RR$ is hence simply the kernel of the 1-form $dt$ associated to it by the exterior derivative. 

== Hypersurface-Orthogonal Distributions

In the previous section, we defined what it means for a vector field to be normal to a hypersurface foliation $Sigma$ generated by a scalar function $t in C^infty (cal(M))$ with $dt != 0$ everywhere. In particular, we constructed an explicit example of such a vector field 
$n^sharp$ by setting
$
  n = alpha dt quad "such that"quad g(n^sharp, n^sharp) = pm 1.
$
This addressed the question _"Given a hypersurface foliation, can we find a vector field that is normal to it everywhere?"_.

In this section, we turn that question around: _"Given a vector field, does there exist a hypersurface foliation to which it is everywhere normal?"_. This leads us to the notion of _hypersurface-orthogonal_ vector fields---those that are locally normal to a family of hypersurfaces. Without proof, we will also give the so-called _Frobenius condition_ which can be used to check this property directly. 
\ \
*Definition* (Hypersurface-Orthogonal Vector Fields) Let $X in Gamma(T cal(M))$ be a smooth vector field on a smooth manifold $cal(M)$. We say that $X$ is _hypersurface-orthogonal_ if, for every point $p in cal(M)$, there exists a neighbourhood $U subset cal(M)$ of $p$ and a local foliation $Sigma = {Sigma_t}#h(0em)_(t in RR)$ of $U$ with the property that
$
  g(X_p, Y) = 0 quad "for all" Y in T_p Sigma_t(p).
$
In words: $X$ is hypersurface-orthogonal if, at every point, it is orthogonal to the leaves of some local foliation.

Though geometrically intuitive, this definition of hypersurface-orthogonality is rather difficult to verify in practice. We now derive an equivalent _analytic_ criterion that allows us to check whether a given vector field $X in Gamma(T cal(M))$ is hypersurface-orthogonal.

The derivation proceeds in a few steps:
- First, recall that any local foliation $Sigma = {Sigma_t}$ defined on an open set $U subset cal(M)$ can be represented by a smooth function $phi : U -> RR$ with $dphi != 0$, whose level sets define the leaves as $Sigma_phi(p) = phi^(-1)(phi(p))$. Conversely, any such function defines a local foliation by level sets. In the previous section, we saw that the tangent spaces to the leaves are given by the pointwise kernels of the differential,
  $
    ker(d phi_p) = T_p Sigma_phi(p)
  $<kerCharacterisationTpSigmaphip>
- A vector field $X$ is orthogonal to the foliation if and only if for all $Y in T_p cal(M)$
  $
    0 = g(X_p, Y) = X_p^flat (Y) quad <=> quad Y in T_p Sigma_phi(p),
  $
  i.e. $X_p^flat$ annihilates all vectors tangent to the leaf through $p$. But by @kerCharacterisationTpSigmaphip[equation], those tangent vectors are precisely the kernel of $dphi_p$. Hence, the two kernels must agree,
  $
    ker(X_p^flat) = ker(dphi_p).
  $
- Since both $X_p^flat$ and $dphi_p$ are non-vanishing 1-forms, they must be pointwise proportional,
  $
    X_p^flat prop dphi_p.
  $
  That is, there exists a smooth function $lambda$ such that
  $
    X^flat = lambda dphi
  $
  on $U$.

We have thus shown that $X$ is hypersurface-orthogonal if and only if, for every point $p in cal(M)$, there exist smooth functions $phi,lambda in C^infty (U)$ on a neighbourhood $U subset cal(M)$ of $p$ such that
$
  X^flat = lambda dphi.
$<analyticConditionHypersurfaceOrthogonality>
This provides a concrete and (more) checkable analytic characterisation of hypersurface-orthogonality.

@analyticConditionHypersurfaceOrthogonality is a condition on the 1-form $X^flat$ associated to the vector field $X$. Conversely, one can also say that if a 1-form $omega$ can locally (say, on $U subset cal(M)$) be expressed as
$
  omega = lambda dphi, quad lambda, phi in C^infty (U),
$
then the associated vector field $X = omega^sharp$ is hypersurface-orthogonal. An immediate example is when $omega = dx^(mu_0)$ is the differential of a fixed coordinate $x^(mu_0)$, which trivially has the above form and hence generates the hypersurface-orthogonal vector field
$
  X = omega^sharp = (dx^(mu_0))^sharp = g^(mu_0 nu) diff_nu,
$
which is orthogonal to the foliation generated by fixing the coordinate $x^(mu_0)$, i.e.
$
  Sigma_(t) = {p in U subset cal(M) | x^(mu_0)(p) = t}.
$
However, there are many vector fields or 1-forms for which it is nontrivial to determine whether a representation as in @analyticConditionHypersurfaceOrthogonality[condition] exists. In such cases, the so-called _Frobenius condition_ comes into play, which we formulate (but not prove) below.
\ \
*Theorem* (Frobenius Condition) Let $cal(M)$ be a smooth manifold, $U subset cal(M)$ an open neighbourhood, and $omega$ a differential 1-form on $U$. Then, $omega$ admits a local expression of the form
$
  omega = lambda dphi
$
for some $lambda, phi in C^infty (U)$ if and only if the condition
$
  omega wedge domega = 0
$
holds.

This is now a purely computational way to verify whether a vector field $X$ is hypersurface-orthogonal---in components, this condition reads
$
  0 attach(=,t:!) X^flat wedge dX^flat = X_mu dx^mu wedge (diff_nu X_lambda) dx^nu wedge dx^lambda = X_mu diff_nu X_lambda dx^mu wedge dx^nu wedge dx^lambda,
$
which yields the standard identity
$
  X_(\[mu) nabla_nu X_(lambda\]) attach(=,t:!) 0.
$
Notice that we were able to replace partial by covariant derivatives due to the antisymmetrisation, which removes the (here torsion-free) connection components from the expression.

As announced, we will not prove the Frobenius condition here, at least not the difficult implication that $omega wedge domega = 0$ guarantees a decomposition into $omega = lambda dphi$. To make the claim somewhat more plausible, though, we can opt to briefly examine the converse direction, which is immediate from a direct computation---let us assume $omega = lambda dphi$, and plug in:
$
  omega wedge domega = lambda dphi wedge d(lambda dphi) = lambda dphi wedge dlambda wedge dphi = 0,
$
where the last equality holds by anti-symmetry of the wedge product. 


== Decomposition of the Metric under Foliation

In many applications such as ADM, it is essential to decompose a metric in terms of its contributions tangential and normal to the leaves of a hypersurface foliation. In this section, we approach this decomposition from two angles, both from a coordinate-independent and -dependent perspective. 

=== Coordinate-Independent Perspective
We begin by recalling from @sectionFoliationTangentBundleDecomposition that given a hypersurface foliation $Sigma = {Sigma_t}$ of a (pseudo-)Riemannian manifold $cal(M)$, the tangent bundle can be decomposed as
$
  T cal(M) = T Sigma plus.circle N Sigma.
$<orthoDecompFoliationTangentNormalBundles>
This decomposition is implemented by the corresponding bundle maps 
$
  P : T cal(M) -> T Sigma, quad Q : T cal(M) -> N Sigma
$
---pointwise, these are simply the orthogonal projections
$
  P_p: T_p cal(M) -> T_p Sigma_t(p),quad Q_p : T_p cal(M) -> N_p Sigma_t(p), quad p in cal(M).
$
For any vector field $X in Gamma(T cal(M))$, we have
$
  X = P(X) + Q(X)
$
where $P(X) in Gamma(T Sigma)$ and $Q(X) in Gamma(N Sigma)$. 

Alternatively, in the pointwise view, for any $X_p in T_p cal(M)$, it holds that
$
  X_p = P_p (X_p) + Q_p (X_p),
$
such that $P_p (X_p) in T_p Sigma_t(p)$ and $Q_p (X_p) in N_p Sigma_t(p)$, and hence
$
  g(P_p (X_p), Q_p (X_p)) = 0.
$
Further recall that the associated $(0,2)$-tensors of $P$ and $Q$ decompose the metric as
$
  g = P + Q.
$
Having set the stage for orthogonal decompositions of the metric, we now turn to the special case where $Sigma$ is a hypersurface foliation generated by $t in C^infty (cal(M))$ with $dt != 0$ everywhere. Our goal is to derive the orthogonal projection $Q$ onto $N Sigma$, as this will enable us to decompose the metric as $P+Q$ using $P = g - Q$.

To make this decomposition explicit, we now compute the form of $Q$ in terms of the normal 1-form $n$ we introduced in @sectionFoliationTangentBundleDecomposition. More precisely, it is given by
$
  n = alpha dt,
$
with $alpha in C^infty (cal(M))$ fixed by the condition
$
  g(n^sharp, n^sharp) = epsilon = pm 1, 
$
where the sign $epsilon$ reflects whether $n^sharp$ is spacelike or timelike (we do not consider the null case here). Further recall that $n$ annihilates vectors tangent to $Sigma$, i.e.
$
  g(n^sharp, X) = n(X) = 0 quad forall X in Gamma(T Sigma).
$
Let us consider how $n$ acts on an arbitrary vector $X in Gamma(N Sigma)$. For our case of a hypersurface foliation, $N Sigma$ is a vector bundle with one-dimensional fibres. Since $n^sharp in Gamma(N Sigma)$ is nowhere vanishing, we can express $X$ as a scalar multiple of it, i.e.
$
  X = lambda n^sharp, quad "for some" lambda in C^infty (cal(M)).
$
This allows us to derive the action of the normal 1-form $n$ on any $X in Gamma(N Sigma)$ as
$
  n(X) = lambda n(n^sharp) = lambda g(n^sharp, n^sharp) = epsilon lambda.
$
Multiplying this equation by $epsilon n^sharp$ yields
$
  epsilon n^sharp dot n(X) = epsilon^2 lambda n^sharp = X.
$
Let us briefly digest this last result. It tells us that the map
$
  epsilon n^sharp otimes n : T cal(M) -> N Sigma
$
acts as the identity when restricted to $N Sigma$. Moreover, since $n$ annihilates any vector in $Gamma(T Sigma)$ and due to the orthogonal @orthoDecompFoliationTangentNormalBundles[decomposition], we conclude that $epsilon n^sharp otimes n$ is the orthogonal projection of $T cal(M)$ onto $N Sigma$. Consequently, we have 
$
  Q = epsilon n^sharp otimes n. 
$<generalResultNormalProjector>
The metric hence decomposes as
$
  g &= P + Q\ 
  &= (g - epsilon n otimes n) + epsilon n otimes n\
  &= (g - epsilon alpha^2 dt otimes dt) + epsilon alpha^2 dt otimes dt
$
where we identify 
$
  P = g - epsilon alpha^2 dt otimes dt
$<projectorPinADM>
as the orthogonal projector onto $T Sigma$. Recall that $P$ is the tangential component of the metric, pulling back to the induced metric on $Sigma_t$ as
$
  gamma = iota^* g = iota^* P
$
due to $iota^* Q = 0$. 

=== ADM-Type Metric Decomposition in Coordinates
We again assume the (pseudo-)Riemannian manifold $cal(M)$ to be equipped with a hypersurface foliation $Sigma = {Sigma_t}$ generated by a scalar function $t in C^infty (cal(M))$ with $dt != 0$ everywhere. We use $t$ as a coordinate, extending it to a full, local coordinate system $x^mu = (t, y^i)$, $i=1,...,dim cal(M)-1$ by transverse coordinates $y^i$. Let us again consider the normal 1-form
$
  n = alpha dt,
$
subject to the condition
$
  g(n^sharp, n^sharp) = epsilon.
$<normalisationConditionAlpha>
Our goal will be to write the metric $g$ in terms of the components of $n^sharp$---i.e. of the function $alpha$ and a vector $beta^i$ we will introduce momentarily---as well as the induced metric on the leaves, which has the components
$
  gamma_(i j) := gamma(diff_i,diff_j) = g(diff_i, diff_j).
$

The vector $n^sharp$ associated to the normal 1-form can be written in terms of its coordinate components, yielding
$
  n^sharp = alpha delta^t_mu g^(mu nu) diff_nu = alpha g^(t mu) diff_mu = alpha g^(t t) diff_t + alpha g^(t i) diff_i.
$
Making use of the @normalisationConditionAlpha[normalisation condition], we find
$
  epsilon = g(n^sharp, n^sharp) = n(n^sharp) = alpha^2 g^(t t) underbrace(dt(diff_t),=1) + alpha^2 g^(t i) underbrace(dt(diff_i),=0) = alpha^2 g^(t t),
$
or equivalently,
$
  g^(t t) = epsilon/alpha^2.
$
This is our first direct relationship between $alpha$ and a component of the (inverse) metric. Inserting this back into the expansion for $n^sharp$ above, we get
$
  n^sharp = epsilon/alpha diff_t + alpha g^(t i) diff_i.
$
We now introduce the _shift vector_ $beta^i$ as
$
  beta^i = - epsilon alpha^2 g^(t i),
$
where the factor is chosen such one can extract $epsilon\/alpha$ as a common factor in $n^sharp$, i.e.
$
  n^sharp = epsilon/alpha ( diff_t - beta^i diff_i).
$
Equivalently, we may express $diff_t$ in terms of $n^sharp$ and $beta = beta^i diff_i$,
$
  diff_t = epsilon alpha n^sharp + beta.
$
This allows us to derive both $g_(t t)$ and $g_(t i)$ via
$
  g_(t t) &= g(diff_t, diff_t) = epsilon^2 alpha^2 underbrace(g(n^sharp,n^sharp),=epsilon) + 2 epsilon alpha underbrace(g(n^sharp, beta),=0) + underbrace(g(beta,beta),=gamma_(i j) beta^i beta^j)\
  &= epsilon alpha^2 + gamma_(i j) beta^i beta^j,\ \
  g_(t i) &= g(diff_t, diff_i) = epsilon alpha underbrace(g(n^sharp, diff_i),=0) + g(beta, diff_i) = gamma_(i j) beta^j
$
Here, we made use of the fact that $n^sharp$ is normal to the foliation, whereas the $diff_i$ and hence also $beta$ are tangent to it---or, algebraically,
$
  g(n^sharp, diff_i) = n(diff_i) =  alpha underbrace(dt(diff_i),=0) = 0.
$ 
In summary, the components of the metric are given by
$
  g_(t t) = epsilon alpha^2 + gamma_(i j) beta^i beta^j, wide g_(t i) = gamma_(i j) beta^j,wide g_(i j) = gamma_(i j),
$
or in block matrix form,
$
  g_(mu nu) = mat(epsilon alpha^2 + gamma_(i j) beta^i beta^j , gamma_(i j) beta^j; gamma_(i j) beta^j, gamma_(i j)).
$
The metric tensor is hence given by
$
  g = g_(mu nu) dx^mu otimes dx^nu = epsilon alpha^2 dt otimes dt + gamma_(i j) (dy^i + beta^i dt) (dy^j +  beta^j dt).
$<ADMsplitMetric>

*Remarks*
- _The Lapse_: In the literature, the function $alpha$ is known as the _lapse function_. To understand its meaning (assuming $t$ is a timelike coordinate), consider an observer moving orthogonally to the spatial slices $Sigma_t$, i.e. following the integral curves of the normal vector field $n^sharp$. Since $n^sharp$ is normalised, the tangent vector to such a path $gamma(tau)$, parameterised by proper time $tau$, satisfies
  $
    dot(gamma thin) = n^sharp.
  $
  Thus, the tangent vector is normal to the foliation, and when computing the elapsed proper time $dtau$ between nearby leaves, only the normal part of the metric contributes. Hence, for such a path, we have
  $
    dtau^2 = -g = -Q = n otimes n = alpha^2 dt^2.
  $
  This further implies
  $
    dtau = pm alpha dt.
  $
  In other words, when moving from one hypersurface $Sigma_t$ to the next along the normal flow, the lapse function $alpha$ gives the rate at which proper time $tau$ advances with respect to coordinate time $t$. The larger the lapse, the "more" proper time elapses between adjacent slices for such an observer. Naturally, if one moves between slices while also displacing in tangential directions, additional contributions appear---the lapse governs the conversion only for motion orthogonal to the foliation.

  

- _The Shift_: Given the definition of the normal 1-form,
  $
    n = alpha dt,
  $
  one might naively expect the associated vector to be proportional to $diff_t$, perhaps something like
  $
    n^sharp = epsilon/alpha diff_t. 
  $
  However, as we have seen before, the correct expression is
  $
    n^sharp = epsilon/alpha (diff_t - beta),
  $
  which includes an additional component from $beta$ along the tangential directions $diff_i$. Since $n^sharp$ is normal to the foliation, this tells us that $diff_t$ itself is _not_ normal in general. The shift vector $beta = beta^i diff_i$ encodes precisely this discrepancy---it provides the spatial correction needed to turn $diff_t$ into a vector normal to the leaves.

  Put differently, the shift quantifies the failure of $diff_t$ to be orthogonal to the foliation. It tells us how much spatial displacement occurs---within the leaves themselves---when evolving forward in coordinate time. Alternatively, we can re-express this as
  $
    diff_t = epsilon alpha n^sharp + beta,
  $
  making explicit that $diff_t$ consists of both a normal and a tangential component. From this viewpoint, the shift is simply the tangential part of the time flow vector $diff_t$.

  Motion along the flow of $diff_t$ involves keeping the coordinates $y^i$ fixed. Therefore, the choice of the coordinates $y^i$ is intimately related to the shift vector $beta$---in particular, as we will see in an upcoming example, the coordinates $y^i$ can (typically) be chosen such that the shift is zero. 

- _Induced Metric vs. Tangential Projector_: We have seen before that under the orthogonal decomposition
  $
    g = P + Q,
  $
  with $P$ the tangential and $Q$ the normal projector (or rather, its associated bilinear form), that the induced metric $gamma$ is related via
  $
    gamma = iota^* g = iota^* P + underbrace(iota^* Q,=0) = iota^* P.
  $
  In most contexts, it is fine to say that $gamma = P$, since they act the same on $T Sigma$. However, it should be kept in mind that $P$ is a bilinear form on $T cal(M)$ whereas $gamma$ is only defined on $T Sigma$. In other words, $P$ "knows" how to handle vectors with components in $N Sigma$---namely by mapping those components to zero---while $gamma$, by definition, does not.

  This manifests itself when considering the @projectorPinADM[expression] for $P$ after inserting the ADM decomposition for $g$:
  $
    P = gamma_(i j) (dy^i + beta^i dt) otimes (dy^j + beta^j dt).
  $
  This is in comparison to the induced metric on the leaves, which simply reads
  $
    gamma = gamma_(i j) dy^i otimes dy^j.
  $
  The $dt$ terms equip $P$ with the (correct) means to handle $diff_t$, which generically has both tangential ($prop beta$) and normal components. The tangential components are what $P$ needs to isolate, which is why the $dt$ bits are necessary.

  However, note here that $diff_t$ lies outside of $T Sigma$. Hence, it is a foreign object to the induced metric $gamma$, and it "does not have to know" how to deal with it---it is safe to discard the $dt$ pieces. It need not concern itself with vectors pointing "off the leaf"---that is, directions leading from one slice to another---even if they do have tangential components. The induced metric is content with only measuring vectors fully tangential to the leaf.

  If the shift is zero, though, then $P = gamma$ not just in effect but in substance. One can see this either algebraically, or intuitively as follows: when the shift vanishes, we have $diff_t prop n^sharp$, meaning $diff_t$ is fully normal. In that case, $P$ never has to process $diff_t$ at all---there is no tangential contribution to extract. Hence, no $dt$ terms need to appear, and $P$ reduces directly to $gamma$.

- _Inverse Metric_: Besides the ADM split of the metric, @ADMsplitMetric[], one frequently needs to use the inverse metric,
  $
    g^(-1) = g^(mu nu) diff_mu otimes diff_nu,
  $
  as well. Its components are given by
  $
    g^(t t) = epsilon/alpha^2, quad g^(t i) = -epsilon/alpha^2 beta^i , quad g^(i j) = gamma^(i j) + epsilon/alpha^2 beta^i beta^j,
  $
  where $gamma^(i j)$ is inverse to $gamma_(i j)$ in the sense that
  $
    gamma^(i k) gamma_(k j) = delta^i_j.
  $
  As a block matrix, this reads
  $
    g^(mu nu) = mat(epsilon\/alpha^2, -epsilon beta^i \/ alpha^2; -epsilon beta^i\/alpha^2, gamma^(i j) + epsilon beta^i beta^j \/alpha^2).
  $
  The inverse metric tensor can hence be written as
  $
    g^(-1) &= epsilon/alpha^2 (diff_t - beta^i diff_i) otimes (diff_t - beta^j diff_j) + gamma^(i j) diff_i otimes diff_j\
    &= epsilon n^sharp otimes n^sharp + gamma^(i j) diff_i otimes diff_j .
  $
  Notice that also here, we have a separation into a normal and a tangential part.

== Example: Foliation of $RR^3 without {0}$ into Spheres <exampleR3ConcentricSpheres>

At this point, we should consider an example to solidify our grasp of hypersurface foliations and the ADM decomposition of the metric. To this end, we examine a simple (but still nontrivial) foliation of the manifold $cal(M) = RR^3 without {0}$, equipped with the flat Euclidean metric 
$
  g = dx otimes dx + dy otimes dy + dz otimes dz
$
where $x,y,z$ are the Cartesian coordinates on $cal(M)$. We foliate it into origin-centered spheres of varying radii. Here, we will derive the objects we defined generally in the previous section---the normal 1-form, its associated normal vector, the lapse, the shift, and the ADM decomposition of the metric, as well as the projectors $P$ and $Q$ onto the tangent and normal bundles of the foliation.

First, let us define the foliation by introducing the function
$
  r : cal(M) -> RR, quad r(p) = sqrt(x^2 + y^2 + z^2), quad p = (x,y,z)
$
and define the leaves of the foliation $Sigma = {Sigma_r}$ to be its level sets,
$
  Sigma_r_0 := {p in cal(M) | r(p) = r_0} = r_0 S^2.
$
Notice that this does indeed define a foliation, since
$
  dr = x/r dx + y/r dy + z/r dz != 0
$
everywhere on $cal(M)$---each of the $Sigma_r subset cal(M)$ is a submanifold in its own right. 

We begin by analysing the normal $1$-form and the associated vector, which also yields an explicit expression for the lapse funtion $alpha$. In general, the normal 1-form is given by
$
  n = alpha dr.
$
Since the (inverse) metric in Cartesian coordinates is trivial, we have
$
  n^sharp = alpha dr^sharp = alpha(x/r diff_x + y/r diff_y + z/r diff_z).
$
The normalisation condition on $n$ hence yields
$
  1 &attach(=,t:!) g(n^sharp ,n^sharp) = n(n^sharp) =alpha^2 dr (x/r diff_x + y/r diff_y + z/r diff_z)\
  &= alpha^2 (x^2 /r^2 dx(diff_x) + y^2/r^2 dy(diff_y) + z^2/r^2 dz(diff_z)) \
  &= alpha^2/r^2 underbrace((x^2 + y^2 + z^2),=r^2) = alpha^2,
$
which fixes the lapse as $alpha^2 = pm 1$---we choose the positive sign, i.e. $alpha = +1$. Moreover, we now have the concrete expressions
$
  n &= dr = x/r dx + y/r dy + z/r dz,\
  n^sharp &= x/r diff_x + y/r diff_y + z/r diff_z.
$

Let us now consider the shift vector $beta$. For this, we need to introduce additional coordinates $q^i$, $i=1,2$ to amend $r$ into a full coordinate system. Naturally, a great candidate for the $q^i$ are the angles $theta,phi$ used with spherical coordinates---for now though, let us remain general. What we have now is a relationship between two sets of coordinates,
$
  cases(r &= r(x,y,z),q^i &= q^i (x,y,z),reverse:#true) quad<-->quad cases(x &= x(r,q^i), y&=y(r,q^i), z&=z(r,q^i))
$
Recall from the previous section that the shift $beta$ is given by
$
  beta = diff_r - underbrace( epsilon alpha,=1) n^sharp = diff_r - n^sharp,
$
quantifying how much tangential motion (since $beta = beta^i diff_i in Gamma(T Sigma)$) displacement along the flow of $diff_r$ entails. In terms of the Cartesian coordinate basis ${diff_x, diff_y, diff_z}$, we may express $diff_r$ as
$
  diff_r = (diff x)/(diff r) diff_x + (diff y)/(diff r) diff_y + (diff z)/(diff r) diff_z.
$
Note that this cannot be simplified further without assumptions on the coordinates $q^i$. Either way, we can now explicitly write down the shift as
$
  beta = diff_r - n^sharp = ((diff x)/(diff r) - x/r ) diff_x + ((diff y)/(diff r) - y/r ) diff_y +((diff z)/(diff r) - z/r ) diff_z.
$
At this point, we cannot really proceed without some sort of assumption on the $q^i$. So, let us do this: we require that $beta = 0$. In other words, we want the $q^i$ to be such that we have zero shift, that $diff_r$ is normal to $T Sigma$. 

This requires that each of the coefficients of $beta$ above must be zero. Notice that
$
  (diff x)/(diff r) - x/r = 0 quad <=>quad x(r,q^i) = r f_x (q^i)
$
for some function $f_x (q_i)$, and analogously for $y$ and $z$. Thus, requiring $beta = 0$ implies that the coordinates $x,y,z$ are simply parametrisations of the $2$-sphere $S^2$ scaled linearly by the radius $r$. Put differently, if we require $beta = 0$, we must parameterise each of our leaves $Sigma_r = r S^2$ with the same coordinates $q^i$, simply stretched by the sphere's radius. 

All that is left to do now is to choose one's favourite parametrisation of $S^2$. Here, we opt to use the standard coordinates $theta,phi$ which parameterise the 2-sphere of radius $r$ as
$
  x &= r sin theta cos phi,\
  y &= r sin theta sin phi,\
  z &= r cos theta.
$
We now have our lapse and shift, $alpha =1$ and $beta = 0$---all that remains to be derived are the transverse metric components
$
  gamma_(i j) = g(diff_i, diff_j), quad i,j in {theta,phi}.
$
The relevant vectors are
$
  diff_theta &= r cos theta (cos phi thin diff_x + sin phi thin diff_y) - r sin theta thin diff_z,\
  diff_phi &= r sin theta (-sin phi thin diff_x + cos phi thin diff_y).
$
A further mechanical computation reveals that
$
  gamma_(theta theta) = r^2, quad gamma_(theta phi) = 0, quad gamma_(phi phi) = r^2 sin^2 theta.
$
We are now ready to write down the full ADM decomposition of the metric. Inserting into the @ADMsplitMetric[general result], we find
$
  g &= epsilon alpha^2 dr otimes dr + gamma_(i j)(dq^i + beta^i dr)(dq^j + beta^j dr)\
  &= dr otimes dr + r^2 dtheta otimes dtheta + r^2 sin^2 theta dphi otimes dphi
$
At this point, one might reasonably point out that we've essentially just recovered the familiar expression for the flat metric in spherical coordinates through a rather elaborate detour. So far, this may seem suspiciously close to reinventing spherical coordinates. That, however, would be a disheartening conclusion after having gone through all this effort---so let us not argue that. Instead, let us appreciate having seen the ADM decomposition machinery at work in a setting where the outcome is familiar, but derived by taking the scenic route. After all, we've just learned that spherical coordinates arise naturally when seeking a zero-shift ADM decomposition of the flat Euclidean metric on concentric spheres. Which, as every self-respecting postgraduate knows, is the entire point of advanced studies: to rediscover well-known results in a more baroque fashion---ideally in a way that makes the undergrads look impressed.

Now that we have derived a particular ADM decomposition of the flat Euclidean metric on $cal(M) = RR^3 without {0}$, we can also briefly consider the bilinear forms $P$, $Q$ associated to the orthogonal decomposition
$
  T cal(M) = T Sigma oplus N Sigma.
$
The normal part of the metric, $Q$, is simply given by
$
  Q = epsilon n otimes n = dr otimes dr.
$
Correspondingly, the tangential projection $P$---or equivalently (in this case, since $beta = 0$), the induced metric $gamma$ on the leaves---is given by
$
  gamma = P = g - Q = r^2 (dtheta otimes dtheta + sin^2 theta dphi otimes dphi).
$
This is also what one would expect: an appropriately scaled metric on $S^2$. 


= Covariant Derivatives on Foliations and Submanifolds 
The goal of this section is to define how a connection on a manifold $cal(M)$ induces a connection on the leaves of a foliation $Sigma$ of $cal(M)$. To this end, we review the definition of a connection $nabla$ as well as how the conditions of vanishing torsion and metric compatibility uniquely single out the Levi-Civita connection. We then proceed to define the induced connection $mnabla$ on a foliation as the tangential projection of $nabla$, and show that if $nabla$ is of Levi-Civita type, then so is $mnabla$.
== Review: Affine Connections
Before discussing how a connection on a (pseudo-)Riemannian manifold $cal(M)$ induces a connection on the submanifolds $Sigma_(t^A) subset cal(M)$ that comprise a foliation $Sigma$, we review the definition of affine connections and recall how the conditions of metric compatibility and torsion-freeness uniquely determine the Levi-Civita connection. 

Let
$
  T^((r,s)) cal(M) = (T cal(M))^(otimes r) otimes (T^* cal(M))#h(0em)^(otimes s)
$
denote the _bundle of $(r,s)$-tensors on $cal(M)$_, where the tensor product between vector bundles is to be understood pointwise---i.e., as a tensor product of the fibres at each point. Tensor fields are smooth sections of this bundle; in other words, an $(r,s)$-tensor field $T$ is a section of $T^((r,s))cal(M)$, that is, $T in Gamma(T^((r,s)) cal(M))$. With this notational preface out of the way, we are now ready to give the definition of affine connections.
\ \
*Definition* (Affine Connection) Let $cal(M)$ be a smooth manifold, and consider a linear map
$
  nabla : Gamma(T^((r,s)) cal(M)) -> Gamma(T^((r,s+1))cal(M)) = Gamma(T^* cal(M) otimes T^((r,s))cal(M)),
$
which maps $(r,s)$-tensors to $(r,s+1)$-tensors. 

With the interpretation of a tensor as a multilinear map on a set of vectors $Y_1,...,Y_s in Gamma(T cal(M))$ and covectors $omega_1,...,omega_r in Gamma(T^* cal(M))$, we introduce the shorthand notation
$
  nabla_X T(omega_1,...,omega_r,Y_1,...,Y_r)  = (nabla T) (X,omega_1,...,omega_r,Y_1,...,Y_r),
$
which simply indicates that $nabla_X$ populates the additional vector argument introduced when passing from $T$ to $nabla T$ with the vector $X$.

We call $nabla$ an _affine connection_ or _covariant derivative_ if it satisfies:

+ Reduction to the exterior derivative on functions $phi in C^infty (cal(M))$,
  $
    nabla phi = dphi,
  $

+ Tensorial Leibniz rule,
  $
    nabla (T otimes S) = (nabla T) otimes S + T otimes (nabla S),
  $

+ and compatibility with contraction#footnote[The contraction operator $tr^i_j$ acts on a tensor by contracting its $i$-th contravariant with the $j$-th covariant index when written in components.] of a contra- and a covariant index,
  $
    nabla (tr^i_j T) = tr^i_j (nabla T).
  $

This is a rather abstract (though likely familiar, if you got this far in these notes) definition, which benefits from building some intuition. Let us make a few remarks:
- *Connection as Derivative* The first axiom is the most concrete---requiring that the connection reduce to the exterior derivative when acting on functions immediately introduces a differential structure to $nabla$. In this context, it makes sense to consider the associated map $nabla_X$ for a vector field $X in Gamma(T cal(M))$, and observe what it reduces to on a smooth function $phi in C^infty (cal(M))$. By definition,
  $
    nabla_X phi = dphi (X) = X[phi],
  $
  so the covariant derivative along $X$ coincides with the familiar directional derivative. This anchors the definition in standard calculus and motivates the interpretation of $nabla$ as a generalisation of differentiation to tensor fields. Since both $nabla$ and $nabla_X$ map tensors to tensors of the same type (modulo the extra slot), they provide a coordinate-independent extension of the partial derivative which respects the underlying bundle structure.

- *Role of the Connection on Vectors* In components associated to some coordinate system $x^mu$, we may locally expand a vector field $X in T cal(M)$ as
  $
    X = X^mu diff_mu = X^mu otimes diff_mu.
  $
  While the tensor product notation is formally unnecessary here---since we are simply multiplying smooth functions with the coordinate basis vector fields---it clarifies the use of the Leibniz rule in what follows.

  It is best to think of $nabla$ as acting on $X$ by acting separately on the component functions $X^mu in C^infty (cal(M))$ and the basis vector fields $diff_mu$, 
  $
    nabla X = nabla (X^mu otimes diff_mu) = dX^mu otimes diff_mu + X^mu otimes nabla diff_mu.
  $
  The components $X^mu$ are smooth functions, so by the first axiom of the connection, $nabla X^mu = dX^mu$. However, the axioms impose no direct condition on $nabla diff_mu$. Its role is to restore tensoriality: the first term alone is not a tensor, which we now demonstrate.

  Consider the term $dX^mu otimes diff_mu$ in components, i.e.
  $
    dX^mu otimes diff_mu = (diff_nu X^mu) dx^nu otimes diff_mu.
  $
  Under a coordinate transformation $x^mu -> y^alpha$ with Jacobians
  $
    J^mu_alpha = (diff x^mu)/(diff y^alpha), quad J^alpha_mu = (diff y^alpha)/(diff x^mu),
  $
  the transformed components of a vector field read $X^alpha = J^alpha_mu X^mu$. Differentiating yields
  $
    diff_alpha X^beta = J^nu_alpha diff_nu (J^beta_mu X^mu) = J_alpha^nu J^beta_mu diff_nu X^mu + X^mu (J^nu_alpha diff_nu J^beta_mu).
  $
  The first term transforms as expected, but the second term spoils tensoriality---it introduces an inhomogeneous piece. Therefore,
  $
    dX^mu otimes diff_mu "is not a tensor."
  $
  This is precisely the failure remedied by the $X^mu otimes nabla diff_mu$ term. Writing out the coordinate-transformed derivative in the new basis leads to
  $
    (diff_alpha X^beta) dy^alpha otimes diff_beta = [J^nu_alpha J_mu^beta diff_nu X^mu + X^mu (J^nu_alpha diff_nu J_mu^beta)] dx^nu otimes diff_mu.
  $ 
  We see that the connection term must absorb the inhomogeneous contribution. That is,
  $
    X^mu otimes nabla diff_mu = X^alpha otimes nabla diff_alpha + X^mu (J_alpha^nu diff_nu J_mu^beta) dx^nu otimes diff_mu.
  $<connectionTransformationRulePrecursor>
  The final term determines how $nabla diff_mu$ must transform in order to make $nabla X$ a genuine tensor. We will return to this structure momentarily when introducing the connection coefficients explicitly.

- *Connection Coefficients* As we have seen in the previous remark, the object $nabla diff_mu$ is of particular interest---we should write it in terms of components. To this end, let us introduce additional shorthand notation for covariant derivatives along the coordinate directions $diff_mu$, as
  $
    nabla_mu := nabla_(diff_mu).
  $
  From the derivations in the previous remark, it is clear that the object $nabla diff_nu$, in terms of coordinates, must be of the form
  $
    nabla diff_nu = (nabla_mu diff_nu) otimes dx^mu = tensor(Gamma,+lambda,-mu nu) dx^mu otimes diff_lambda, quad "for coeffiecients" quad tensor(Gamma,+lambda,-mu nu) in C^infty (cal(M)).
  $
  This implicitly defines the _connection coefficients_ $tensor(Gamma,+lambda,-mu nu)$ through
  $
    tensor(Gamma,+lambda, -mu nu) diff_lambda = nabla_mu diff_nu.
  $
  Fully explicitly, they are given by
  $
    tensor(Gamma,+lambda,-mu nu) = dx^lambda (nabla_mu diff_nu).
  $
  Using these coefficients, we can write the action of the connection on a vector $X = X^mu diff_mu in Gamma(T cal(M))$ concretely as
  $
    nabla X &= nabla (X^nu diff_nu) = dX^nu otimes diff_nu + X^nu nabla diff_nu\
    &= (diff_mu X^lambda) dx^mu otimes diff_lambda + (X^nu tensor(Gamma,+lambda,-mu nu)) dx^mu otimes diff_lambda\
    &= (diff_mu X^lambda + tensor(Gamma,+lambda,-mu nu) X^nu) dx^mu otimes diff_lambda.
  $
  Along the coordinate directions, this turns into
  $
    nabla_mu X = (diff_mu X^lambda + tensor(Gamma,+lambda,-mu nu) X^nu) diff_lambda.
  $
  This is nothing but the partial derivative of $X$, amended by a term linear in $X$ that ensures tensoriality of $nabla X$. 

- *Interpretation of the Connection* In light of our previous remark that the connection along some vector $X in Gamma(T cal(M))$ is a covariant generalisation of the directional derivative along $X$, we can give an interpretation to the equality
  $
    nabla_mu diff_nu = tensor(Gamma,+lambda,-mu nu) diff_lambda.
  $
  On the left-hand side, we compute a generalised directional derivative along the coordinate direction $diff_mu$. In other words, we ask ourselves: how does the basis vector $diff_nu$ change as we move along the direction $diff_mu$? This question is answered by the right-hand side, which tells us, for the given combination of $mu$ and $nu$, what the rate of change is, in terms of a linear combination of the coordinate basis vectors. So, the component $tensor(Gamma,+lambda,-mu nu)$ encodes information about "by how much of $diff_lambda$ does $diff_nu$ change when moving along the flow of $diff_mu$?". A connection hence imposes a relationship between vectors in infinitesimally neighbouring tangent spaces.

  A priori, there exists no canonical way to relate vectors at different points of the manifold---that is, to compare vectors in distinct tangent spaces and transport them around the manifold. While the transformation law of the connection coefficients is constrained by the requirement that $nabla T$ be tensorial for all tensors $T$, the connection itself remains arbitrary otherwise. However, by imposing additional (natural) conditions on $nabla$, one can uniquely determine a distinguished connection: the _Levi-Civita connection_. We will pursue this in the next section.

- *Coordinate Transformation Behaviour* From @connectionTransformationRulePrecursor[equation] it can be derived that under a coordinate transformation $x^mu -> y^alpha$, the connection coefficients
  $
    dx^lambda (nabla_mu diff_nu) = tensor(Gamma,+lambda,-mu nu)
  $
  must transform as
  $
    tensor(Gamma,+lambda,-mu nu) -> tensor(Gamma,+alpha,-beta gamma) = J^alpha_lambda J^mu_beta J^nu_gamma tensor(Gamma,+lambda,-mu nu) + J^alpha_lambda diff_beta J^lambda_gamma.
  $<connectionCoeffTransformRule>
  The second term is what makes the connection coefficients transform non-tensorially. Unlike the first term, which is the usual triple contraction with (inverse) Jacobians for a $(1,2)$-tensor, the inhomogeneous second term violates tensorial transformation behaviour.

- *Action on Arbitrary Tensors* Since we know how $nabla$ acts on the coordinate basis vector fields $diff_mu$, we can readily extend it to any contravariant tensor $T in Gamma(T^((r,0)) cal(M))$ by using the Leibniz rule,
  $
    nabla T &= nabla (T^(mu_1...mu_r) diff_mu_1 otimes ... otimes diff_mu_r)\
    &= dT^(mu_1 ... mu_r) otimes diff_mu_1 otimes ... otimes diff_mu_r + T^(mu_1...mu_r) (nabla diff_mu_1) otimes diff_mu_2 otimes ... otimes diff_mu_r + ...\
    &= (diff_nu T^(mu_1 ... mu_r) + tensor(Gamma,+mu_1,-nu lambda) T^(lambda mu_2 ... mu_r) + ... + tensor(Gamma,+mu_r,-nu lambda) T^(mu_1...mu_(r-1) nu)) dx^nu otimes diff_mu_1 otimes ... otimes diff_mu_r.
  $
  In essence, every contravariant index produces a term where the connection coefficients are contracted with the tensor components, as by the Leibniz rule, $nabla$ acts on each of the basis vectors once. 

  For tensors that are not purely contravariant, i.e. general tensors $T in Gamma(T^((r,s)) cal(M))$, however, the general basis expansion reads
  $
    tensor(T,+mu_1...mu_r,-nu_1...nu_s) diff_mu_1 otimes ... otimes diff_mu_r otimes dx^(nu_1) otimes ... otimes dx^(nu_s).
  $
  To be able to extend the action of $nabla$ for general tensors---as above, by means of the Leibniz rule---we must know what 
  $
    nabla dx^mu
  $
  evaluates to. This is where the (so far unused!) third axiom comes into play, the compatibility with contractions: together with the Leibniz rule, it is possible to derive an explicit expression for the above from what we already know. Consider an arbitrary vector $X  = X^mu diff_mu in Gamma(T cal(M))$ as well as an arbitrary 1-form $omega = omega_mu dx^mu in Gamma(T^* cal(M))$. Clearly, since $omega(X) in C^infty (cal(M))$, we have
  $
    tr^1_2 nabla (omega otimes X) = nabla (omega(X)) = d (omega(X)) = d (omega_mu X^mu).
  $ <eq1.5.26>
  Let us take a closer look at the left-hand side. We can derive
  $
    &tr^1_2 nabla (omega otimes X) = tr^1_2 nabla (omega_mu X^nu dx^mu otimes diff_nu)\
    &= tr^1_2 [d(omega_mu X^nu) otimes dx^mu otimes diff_nu + omega_mu X^nu ((nabla dx^mu) otimes diff_nu + dx^mu otimes nabla diff_nu)]\
    &= d(omega_mu X^mu) + omega_mu X^nu tr^1_2((nabla dx^mu) otimes diff_nu + dx^mu otimes nabla diff_nu).
  $
  Inserting this back into @eq1.5.26[equation] above and canceling the exterior derivative terms, by arbitrariness of $omega_mu$ and $X^mu$ we get
  $
    tr_2^1 [(nabla dx^mu) otimes diff_nu] = -tr^1_2 [dx^mu otimes nabla diff_nu].
  $
  Evaluating this along a coordinate direction $diff_lambda$, we find
  $
    tr^1_1 [nabla_lambda dx^mu otimes diff_nu] &= - tr^1_1 [dx^mu otimes nabla_lambda diff_nu] = -tr^1_1 [dx^mu otimes tensor(Gamma,+rho,-lambda nu )diff_rho]\
    &= - tensor(Gamma,+rho,-lambda nu)dx^mu (diff_rho) = - tensor(Gamma,+rho,-lambda nu) delta^mu_rho = - tensor(Gamma,+mu,-lambda nu).
  $
  The left hand side simplifies as
  $
    tr_1^1 [nabla_lambda dx^mu otimes diff_nu] = (nabla_lambda dx^mu)(diff_nu),
  $
  which implies that overall,
  $
    nabla_lambda dx^mu = -tensor(Gamma,+mu,-lambda nu) dx^nu.
  $
  Thus, in summary, $nabla$ acts on the coordinate basis vectors and 1-forms as
  $
    nabla diff_mu &= tensor(Gamma,+lambda,-nu mu) dx^nu otimes diff_lambda,\ 
    nabla dx^mu &= -tensor(Gamma,+mu,-nu lambda) dx^nu otimes dx^lambda,
  $<generalResultConnectionActionInComponents>
  or, along a coordinate direction $diff_mu$,
  $
    nabla_mu diff_nu &= tensor(Gamma,+lambda,-mu nu) diff_lambda,\
    nabla_mu dx^nu &= -tensor(Gamma,+nu,-mu lambda) dx^lambda.
  $
  This, together with the axiom that $nabla phi = dphi$ (or $nabla_mu phi = diff_mu phi$), allows one to evaluate the action of $nabla$ (or $nabla_mu$) on arbitrary tensors. 
  
  To reiterate: the expressions for $nabla_mu diff_nu$ and $nabla_mu dx^nu$ above inform us what the rate of change of the basis vectors $diff_nu$ and 1-forms $dx^nu$ are as one moves along the coordinate direction $diff_mu$. 

== Review: Levi-Civita Connection 
Now that we have introduced the general notion of a connection---a way of encoding the change of basis vectors and 1-forms as one moves through a manifold---we turn to the most prominent example: the _Levi-Civita connection_, which is uniquely determined by two conditions. These are:

+ Vanishing torsion,
+ Metric compatibility.
We will now introduce both conditions rigorously and examine their consequences for the connection coefficients $tensor(Gamma,+lambda,-mu nu)$, which ultimately leads to the _Christoffel Symbols_, which are the coefficients of the Levi-Civita connection. 
\ \
*Vanishing Torsion* A connection provides us with a means of comparing vectors at nearby points and describing how they change as we move infinitesimally along a given direction. Given two vector fields $X,Y in Gamma(T cal(M))$, the covariant derivatives $nabla_X Y$ and $nabla_Y X$ describe how $Y$ changes along the flow of $X$, and how $X$ changes along the flow of $Y$, respectively.

On a flat space, one can think of the vectors $X,Y,X+alpha nabla_Y X$, and $Y + alpha nabla_X Y$ as forming a parallelogram for infinitesimal $alpha$, assuming $[X,Y] = 0$. However, on a general manifold with arbitrary connection, this parallelogram may fail to close. There are two distinct reasons for this failure:
+ The flow paths of $X$ and $Y$ do not commute, i.e., following $X$ then $Y$ leads to a different point than following $Y$ then $X$. This is encoded by the Lie bracket $[X,Y]$. 
+ The change in the transported vector fields differs, i.e., $nabla_X Y != nabla_Y X$. This is an intrinsic feature of the connection.

The _torsion_ is designed to isolate the second phenomenon: it captures the failure of symmetry in the connection itself, independent of the non-commutativity of the vector fields.
\ \
*Definition* (Torsion) Let $cal(M)$ be a smooth manifold, and $nabla$ a connection on $cal(M)$. The torsion is the vector-valued map defined by
$
  T : Gamma(T cal(M)) otimes Gamma(T cal(M)), quad T(X,Y) = nabla_X Y - nabla_Y X - [X,Y].
$
The first two terms, $nabla_X Y - nabla_Y X$, measure the asymmetry of the connection, but also include contributions from the possible non-closure of paths due to non-commuting vector fields. The Lie bracket $[X,Y]$ encodes this latter effect---by subtracting it, we isolate the connection's contribution to the failure of the parallelogram to close, i.e., the torsion.

If $T= 0$, the connection is said to be _torsion-free_. To make this condition explicit, and to check its consequences on the connection coefficients, we now express the torsion in terms of components. Take $X = X^mu diff_mu$, $Y = Y^nu diff_nu$; then
$
  T(X,Y) &= X^mu nabla_mu (Y^nu diff_nu)  - (X^mu diff_mu Y^nu)diff_nu - (X <->Y)\
  &= X^mu lr((cancelr((diff_mu Y^nu) diff_nu) + Y^nu underbrace(nabla_mu diff_nu, = tensor(Gamma,+lambda,-mu nu) diff_lambda)),size:#30%) - cancelr((X^mu diff_mu Y^nu)) diff_nu - (X<->Y)\
  &= X^mu Y^nu tensor(Gamma,+lambda,-mu nu) diff_lambda - (X<->Y) \
  &= X^mu Y^nu (tensor(Gamma,+lambda,-mu nu) - tensor(Gamma,+lambda,-nu mu))
$
Note that subtracting the Lie bracket was essential to remove contributions from non-commutativity of the vector fields $X$ and $Y$. 

The torsion thus has the following properties:

+ It is antisymmetric: $T(X,Y) = -T(Y,X)$;

+ It is $C^infty (cal(M))$-linear in both arguments;
+ It can be written in terms of components as
  $
    T(X,Y) = X^mu Y^nu tensor(T,+lambda,-mu nu) diff_lambda quad "where" quad tensor(T,+lambda,-mu nu) = tensor(T,+lambda,-[mu nu]) = tensor(Gamma,+lambda,-mu nu) - tensor(Gamma,+lambda,-nu mu),
  $
  and is therefore proportional to the anti-symmetric part of the connection coefficients in the lower indices;
+ It is tensorial, since the transformation behaviour under $x^mu -> y^alpha$ is given by
  $
    tensor(T,+lambda,-mu nu) prop tensor(Gamma,+lambda,-[mu nu]) -> tensor(Gamma,+alpha,-[beta gamma]) = J^alpha_lambda J^mu_beta J^nu_gamma tensor(Gamma,+lambda,-[mu nu]) + J^alpha_lambda underbrace(diff_(\[beta) J^lambda_(gamma\]),=0),
  $
  where the inhomogeneous part cancels under anti-symmetrisation because
  $
    diff_beta J^lambda_gamma = (diff^2 x^lambda)/(diff y^beta diff y^gamma)
  $
  is symmetric in $beta$ and $gamma$. 

This fourth property ensures that $T = 0$ is a meaningful, coordinate-independent condition. If torsion were not tensorial, the expression $T(X,Y)=0$ could be true in one chart and false in another. 

We are now ready to impose the condition of vanishing torsion. By the third property above, we have
$
  T = 0 quad <=> quad tensor(T,+lambda,-mu nu) = tensor(Gamma,+lambda,-mu nu) - tensor(Gamma,+lambda,-nu mu) = 0 quad <=> quad tensor(Gamma,+lambda, -mu nu) =  tensor(Gamma,+lambda, -(mu nu)),
$
i.e., a connection is torsion-free if and only if the connection coefficients are symmetric in their lower indices.
\ \
*Metric Compatibility* Luckily, this second condition is simpler to motivate and introduce than that of vanishing torsion. Recall that any connection $nabla$ is, by definition, compatible with contractions between one co- and one contravariant index,
$
  nabla tr^i_j T = tr^i_j nabla T.
$
But as we know, contractions are not restricted to mixed index types---the metric (or its inverse) also allows contractions of two covariant or two contravariant indices. 

Coordinate-independently, these operations take the form
$
  T |-> tr^(i,j)_(1,2) (g otimes T) quad "or" quad T |-> tr^(i,j)_(1,2) (g^(-1) otimes T).
$
depending on whether one is contracting co- or contravariant indices. Perhaps more concretely, in terms of components, this becomes
$
  T^(mu_1 ... mu_r)_(nu_1 ... nu_s) |-> g_(mu_i mu_j) T^(mu_1 ... mu_r)_(nu_1 ... nu_s) quad "or" quad T^(mu_1 ... mu_r)_(nu_1 ... nu_s) |-> g^(nu_i nu_j) T^(mu_1 ... mu_r)_(nu_1 ... nu_s).
$
It is thus not far-fetched to ask that the connection be compatible with such contractions as well---after all, contraction through raising and lowering indices is a natural and frequent operation.

To formalise this, we demand that contractions involving the metric commute with covariant differentiation. Explicitly, we require#footnote[For clarity we suppress any shifts in index positions introduced by tensor products.]
$
  nabla (tr^(i,j)_(1,2) g otimes T) = tr^(i,j)_(1,2) g otimes (nabla T).
$
Using the axioms for $nabla$, we expand the left-hand side into
$
  nabla(tr^(i,j)_(1,2)  g otimes T) = tr^(i,j)_(1,2) nabla (g otimes T) = tr^(i,j)_(1,2) g otimes( nabla T) + tr^(i,j)_(1,2)(nabla g) otimes T
$
The first term matches the desired right-hand side, whereas the remaining second term vanishes if and only if $nabla g = 0$.

Hence, _metric compatibility_ is the requirement that the covariant derivative of the metric vanishes,
$
  nabla g = 0.
$
Notice that this also covers the case involving the inverse metric, since $nabla g = 0 <=> nabla g^(-1) = 0$. In components, the condition $nabla g = 0$ reads
$
  0 &attach(=,t:!) nabla g = nabla (g_(mu nu) dx^mu otimes dx^nu) = (diff_lambda g_(mu nu) - g_(rho nu) tensor(Gamma,+rho,-lambda mu) - g_(mu rho) tensor(Gamma,+rho,-lambda nu)) dx^lambda otimes dx^mu otimes dx^nu\
  &= (diff_lambda g_(mu nu) - Gamma_(nu lambda mu) - Gamma_(mu lambda nu)) dx^lambda otimes x^mu otimes dx^nu,
$
where we defined $Gamma_(lambda mu nu)= g_(lambda rho) tensor(Gamma,+rho,-mu nu)$. Making use of the torsion-freeness---i.e., symmetry of the connection coefficients in the last two indices---this is equivalent to
$
  g_(mu nu,lambda) = Gamma_(mu nu lambda) + Gamma_(nu mu lambda),
$<metricGradientAndConnection>
with $g_(mu nu,lambda) = diff_lambda g_(mu nu)$. This is a direct relationship between the derivatives of the metric and the connection coefficients.

Since we are attempting to single out a specific connection by imposing constraints, the next step is to solve the equation above for $Gamma_(mu nu lambda)$, from which the explicit form of the connection coefficients can be extracted using the inverse metric. This derivation is largely an exercise in index manipulation and combining permutations of the above equation to isolate the connection. The result is the expression for the _Christoffel symbols_, 
$
  tensor(Gamma,+lambda,-mu nu) = 1/2 g^(lambda rho) [g_(rho mu,nu) + g_(rho nu, mu) - g_(mu nu, rho)].
$<ChristoffelSymbols>
The fact that this is the correct expression can also be checked by inserting @metricGradientAndConnection[equation] to replace the expressions $g_(mu nu, lambda)$ with the appropriate index permutations and verifying that the terms cancel correctly. This can be done in a few lines.
\ \
*Remark* The @ChristoffelSymbols[Christoffel symbols] above are the coefficients of the so-called _Levi-Civita connection_, which is the _unique_ connection that is both torsion-free and metric-compatible.

== Induced Connection on Foliations <sectionDefInducedConnection>

The goal of this section is to define an affine connection on the submanifolds constituting a foliation $Sigma$ of a (pseudo-)Riemannian manifold $cal(M)$. Just as with the ambient metric and the choice of induced metric on submanifolds, one is, in principle, free to assign any connection to these submanifolds. However, if $cal(M)$ is already equipped with an affine connection $nabla$, it is natural to study a connection _induced_ by this ambient structure, rather than constructing one from scratch. We begin by developing an intuitive, geometric picture of how such an induced connection should behave, before proceeding to formalise it in a precise definition.

As we have a seen, an affine connection is fully characterised by its action on vector fields, encoded via the connection coefficients $tensor(Gamma,+lambda,-mu nu)$ appearing in
$
  nabla diff_nu = tensor(Gamma,+lambda,-mu nu) dx^mu otimes diff_lambda.
$
Its action on functions $phi in C^infty (cal(M))$, 1-forms, and general tensors is then determined entirely by the axioms of the connection. For example,
$
  nabla phi = dphi, quad nabla dx^nu = - tensor(Gamma,+nu,-mu lambda) dx^mu otimes dx^lambda.
$
In other words, specifying how $nabla$ acts on vector fields suffices to determine its behaviour on all tensor fields. 

Thus, the task of defining a geometrically meaningful connection on a foliation $Sigma= {Sigma_(t^A)}$ reduces to defining how it acts on vector fields tangent to the leaves, and extending that action to arbitrary tensors via the standard axioms of an affine connection.

Therefore, given a connection $nabla$ on a manifold $cal(M)$ endowed with a foliation $Sigma = {Sigma_(t^A)}$, to define an induced connection $mnabla$ on each leaf $Sigma_(t^A)$, we must provide a prescription for
$
  mnabla_X Y quad "for any" quad X,Y in Gamma(T Sigma),
$
which relates back to $nabla$ in a geometrically meaningful way. 

Let us now build some geometric intuition. Recall that the expression $nabla_X Y$ describes how a vector field $Y in Gamma(T Sigma)$ changes as it is transported along the flow of $X$. Further, the pushforward $iota_*$ of the inclusion map $iota$ embeds $T Sigma$ into the ambient tangent bundle $T cal(M)$, allowing us to view $X$ and $Y$ as vector fields living in $T cal(M)$ via their pushforwards $iota_* X$ and $iota_* Y$. 

Since we want the connection $mnabla$ on $Sigma$ to "mimic" the ambient connection $nabla$ on $cal(M)$ as closely as possible, a natural first attempt to define $mnabla$ would be
$
  mnabla_X Y = nabla_(iota_* X)(iota_*Y).
$
In words, this says that $mnabla$ transports $Y$ along $X$ in precisely the same way that the ambient connection $nabla$ transports the embedded version of $Y$ along that of $X$ within $T cal(M)$.

This is a good starting point---essentially all of the geometric structure of $nabla$ is being transferred to $mnabla$. However, this definition has a fundamental flaw. It is subtle but crucial: an induced connection on the foliation $Sigma$ must restrict to a connection on each individual leaf $Sigma_(t^A)$, which is a map
$
  mnabla : Gamma(T^((r,s))Sigma_(t^A)) -> Gamma(T^((r,s+1)) Sigma_(t^A)), quad T |-> mnabla T,
$
and, when acting on (and along) vector fields, is given more concretely by a map
$
  mnabla : Gamma(T Sigma_(t^A)) times Gamma(T Sigma_(t^A)) -> Gamma(T Sigma_(t^A)), quad (X,Y)|-> mnabla_X Y.
$
In short, the connection must send vector fields tangent to the leaves to other vector fields tangent to the leaves. But $nabla_(iota_* X) (iota_* Y)$, although well-defined in $Gamma(T cal(M))$, need not necessarily lie in the subbundle $T Sigma$. The ambient connection $nabla$ is under no obligation to preserve tangency to the leaves---it can easily produce components orthogonal to them when transporting vectors through the submanifolds.

We are therefore forced to modify our first attempt so as to eliminate any normal components that may arise. For this purpose, recall the left-inverse $(iota_*)#h(0em)^(-1): T cal(M) -> T Sigma$ introduced earlier. While there exist infintely many such left-inverses, we singled out a unique one by requiring that the projection $P = iota_* compose (iota_*)#h(0em)^(-1)$ be orthogonal with respect to the ambient metric. This construction gives us precisely the tool we need: $(iota_*)#h(0em)^(-1)$ acts as the identity on $T Sigma$, while annihilating vectors in the normal bundle $N Sigma$; that is,
$
  ker((iota_*)#h(0em)^(-1)) = N Sigma.
$
Thus, applying $(iota_*)#h(0em)^(-1)$ to our first attempt yields a vector field in $Gamma(T Sigma)$, with all normal components stripped away. This motivates the corrected definition
$
  mnabla_X Y = (iota_*)#h(0em)^(-1) nabla_(iota_* X) (iota_* Y).
$
This expression is admittedly cumbersome to read, but its geometric interpretation is clear. To compute $mnabla_X Y$ for $X,Y in Gamma(T Sigma)$, we push both vector fields forward (i.e., embed) into the ambient tangent bundle, transport one along the other with the ambient connection, and then project the result back onto $T Sigma$ using the orthogonal "projection" $(iota_*)#h(0em)^(-1)$. This removes any component normal to the foliation that may have been introduced by $nabla$, while remaining true to the tangent contributions.

Let us now formalise this geometric construction as a rigorous definition.
\ \
*Definition* (Induced Connection on a Foliation) Let $cal(M)$ be a smooth manifold, $Sigma = {Sigma_(t^A)}$ a foliation of $cal(M)$, and 
$
  nabla : Gamma(T^((r,s))cal(M))->Gamma(T^((r,s+1))cal(M))
$
a connection on $cal(M)$. We call the linear map
$
  mnabla : Gamma(T^((r,s))Sigma) -> Gamma(T^((r,s+1)) Sigma)
$
the _induced connection on $Sigma$_ (associated to $nabla$) if

+ it reduces to the exterior derivative on functions $phi in C^infty (cal(M))$,
  $
    mnabla phi = dphi,
  $

+ acts on vector fields $X,Y in Gamma(T Sigma)$ as 
  $
    mnabla_X Y = (iota_*)^(-1) (nabla_(iota_* X) (iota_* Y))
  $
  where $(iota_*)^(-1):T cal(M) -> T Sigma$ is the left-inverse of $iota_*$ such that
  $
    ker((iota_*)^(-1)) = N Sigma,
  $
  or equivalently, $P = iota_* compose (iota_*)^(-1)$ is the orthogonal projection $P : T cal(M) -> T Sigma$;

+ satisfies a tensorial Leibniz rule,
  $
    mnabla (T otimes S) = (mnabla T) otimes S + T otimes (mnabla S),
  $
  on arbitrary tensors $T,S$;

+ and is compatible with contractions,
  $
    mnabla (tr^i_j T) = tr^i_j (mnabla T),
  $
  for any tensor $T$.

*Remarks:*
- The axioms 1., 3. and 4. guarantee that $mnabla$ is itself an affine connection. 

- Since the connection coefficients $tensor(macron(Gamma),+k,-i j)$ in coordinates $(t^A,y^i)$ adapted to $Sigma$ are defined by the action of the connection on basis vectors, axiom 2. allows for their explicit computation in terms of the coefficients $tensor(Gamma,+lambda,-mu nu)$ of the ambient connection $nabla$, the pushforward matrix $E^mu_i = (diff x^mu)/(diff y^i)$ and its orthogonal left-inverse $E^i_mu$. Concretely, we may derive
  $
    tensor(macron(Gamma),+k,-i j)diff_k &= mnabla_i diff_j = (iota_*)^(-1) (nabla_(iota_* diff_i) (iota_* diff_j))\ &= (iota_*)^(-1) (nabla_(E^mu_i diff_mu) (E^nu_j diff_nu))\
    &= E^mu_i (iota_*)^(-1) (nabla_mu (E^nu_j diff_nu))\ &= E^mu_i (iota_*)^(-1) lr(((diff_mu E^nu_j)diff_nu + E^nu_j underbrace(nabla_mu diff_nu,=tensor(Gamma,+lambda,-mu nu) diff_lambda)),size:#30%)\
    &= E^mu_i ((diff_mu E^lambda_j) + E^nu_j tensor(Gamma,+lambda,-mu nu)) underbrace((iota_*)^(-1)(diff_lambda),=E_lambda^k diff_k)\
    &= (E^k_lambda E^mu_i E^nu_j tensor(Gamma,+lambda,-mu nu) + E^k_lambda diff_i E^lambda_j) diff_k
  $
  from which we conclude
  $
    tensor(macron(Gamma),+k,-i j) = E^k_lambda E^mu_i E^nu_j tensor(Gamma,+lambda,-mu nu) + E^k_lambda diff_i E^lambda_j.
  $<projectionConnectionCoeffs>
  This matches the structure of the @connectionCoeffTransformRule[coordinate transformation behaviour] of the connection. 
- Since $mnabla$ is an affine connection connection, the above relationship that
  $
    mnabla_i diff_j = tensor(macron(Gamma),+k,-i j) diff_k
  $
  on vectors extends to 1-forms by
  $
    mnabla_i dy^j = - tensor(macron(Gamma), +j,-i k) dy^k
  $
  by the @generalResultConnectionActionInComponents[general result]. Further, the linearity and Leibniz rule that hold for $mnabla$ allow for the extension onto any tensor $T$.
== Proof: Induced Connection from Levi-Civita Connections
The previous section introduced the notion of a connection on the submanifolds of a foliation, induced from a connection on the ambient manifold. This construction defines a specific connection on each submanifold, but it is not immediately clear which properties of the ambient connection are preserved, or whether the induced connection might acquire new properties absent from the ambient. Such questions typically require case-by-case analysis.

There is, however, an important special case, which we explore in this section. Recall that the Levi-Civita is uniquely characterised by two conditions: vanishing torsion and compatibility with the metric. We will show that if the ambient connection satisfies these conditions, then the induced connection does as well---in this sense, they are inherited. It follows that the connection induced by a Levi-Civita connection is itself Levi-Civita. In particular, this yields an alternative to the @projectionConnectionCoeffs[projection of connection coefficients], as the coefficients of the induced connection can now be computed directly from linear combinations of partial derivatives of the induced metric. 

In the following, let $cal(M)$ denote a (pseudo-)Riemannian manifold with metric tensor $g in Gamma(T^((0,2))cal(M))$, $Sigma = {Sigma_(t^A)}$ a foliation of $cal(M)$, and $gamma = iota^* g$ the induced metric on $Sigma_(t^A)$. Further, let $nabla$ be a connection on $cal(M)$ and $mnabla$ the connection on $Sigma$ induced by $nabla$. 
\ \
*Vanishing Torsion* Suppose $nabla$ has vanishing torsion. Concretely, this means that
$
  T(X,Y) = nabla_X Y - nabla_Y X - [X,Y] = 0 quad forall X,Y in Gamma(T cal(M)),
$
or equivalently,
$
  nabla_X Y - nabla_Y X = [X,Y]
$

The torsion of the induced connection, evaluated on vector fields $X,Y in Gamma(T Sigma)$ is given by
$
  macron(T)(X,Y) &= mnabla_X Y - mnabla_Y X - [X,Y] \ 
  &= (iota_*)#h(0em)^(-1)[nabla_(iota_* X) (iota_*Y) - nabla_(iota_* Y) (iota_*X)] - [X,Y]\
  &= underbrace((iota_*)#h(0em)^(-1)[iota_* X, iota_* Y],=[X,Y] quad (*)) - [X,Y] = 0.
$
Thus, as a consequence of the vanishing torsion of $nabla$, the induced connection $mnabla$ is torsion-free as well. The identity $(*)$ originates from
$
  [iota_* X, iota_* Y] = iota_* [X,Y],
$
which can be shown coordinate-independently using the definition of the commutator,
$
  [X,Y][f]= X[Y[f]]-Y[X[f]], quad f in C^infty (cal(M)),
$
and the fact that $(iota_* X)[f] = X[f compose iota]$.

Alternatively, one can show this on the level of components as well. The ambient connection being torsion free implies that its connection coefficients $tensor(Gamma,+lambda,-mu nu)$ are symmetric in the lower indices. The induced connection coefficients, given by
$
tensor(macron(Gamma),+k,-i j) = E^k_lambda E^mu_i E^nu_j tensor(Gamma,+lambda,-mu nu) + E^k_lambda diff_i E^lambda_j
$
(cf. @connectionCoeffTransformRule[eq.]), are then symmetric in the lower indices $i,j$ as well---the first term is symmetric due to the symmetry of $tensor(Gamma,+lambda,-mu nu)$, and the second because of
$
  diff_i E^lambda_j = (diff^2 x^lambda)/(diff y^i diff y^j),
$
which is symmetric as well. 

We have now given a coordinate-free and a component-based proof of the fact that an induced connection inherits torsion-freeness from the ambient connection. It remains to show that metric compatibility is preserved as well, allowing us to conclude that an ambient Levi-Civita connection induces a Levi-Civita connection on the submanifolds of a foliation. 
\ \
*Metric Compatibility* To show that metric compatibility of $nabla$ is inherited to $mnabla$, let us first derive a useful identity---a product rule for the inner product induced by the metric $g$. Observe that for $X,Y,Z in Gamma(T cal(M))$, we have
$
  &nabla_X g(Y,Z) = nabla_X tr_(1,2)^(3,4) (g otimes Y otimes Z)\
   &= tr^(3,4)_(1,2) [(nabla_X g) otimes Y otimes Z + g otimes (nabla_X Y) otimes Z + g otimes Y otimes (nabla_X Z)]\
   &= (nabla_X g) (Y,Z) + g(nabla_X Y, Z) + g(Y, nabla_X Z),
$<proofDerivationOfMetric>
which is equivalent to
$
  nabla_X g(Y,Z) - (nabla_X g)(Y,Z) = g(nabla_X Y, Z) + g(Y, nabla_X Z).
$<metricConnectionIdentity>
In particular, a connection $nabla$ is compatible with a metric $g$ if and only if the first term vanishes, i.e. if
$
  nabla_X g(Y,Z) = g(nabla_X Y, Z) + g(Y, nabla_X Z) quad forall X,Y,Z in Gamma(T cal(M)).
$<compatibleMetricConnectionIdentity>
We will now assume that this holds for $nabla$ and $g$, and show that as a consequence, it is also true for the induced connection $mnabla$ and metric $gamma$. This proceeds as follows, for $X,Y,Z in Gamma(T Sigma)$:
$
  mnabla_X gamma(Y,Z) - (mnabla_X gamma)(Y,Z) &attach(=,t:1.) gamma(mnabla_X Y, Z) + (Y <-> Z)\
  &attach(=,t:2.) (iota^* g) ((iota_*)#h(0em)^(-1) nabla_(iota_* X) (iota_* Y), Z) + (Y <-> Z)\
  &attach(=,t:3.) g lr((underbrace(iota_* compose(iota_*)#h(0em)^(-1),=P) nabla_(iota_* X) (iota_* Y), iota_* Z),size:#35%)+ (Y <-> Z)\
  &attach(=,t:4.) g lr((nabla_(iota_* X) (iota_* Y), underbrace(P iota_* Z,=iota_* Z)),size:#35%) + (Y <-> Z)\
  &attach(=,t:5.) nabla_(iota_* X) g(iota_* Y, iota_* Z)\
  &attach(=,t:6.) mnabla_X gamma(Y,Z)

$<inducedMetricCompatibilityDerivation>
which shows that
$
  (mnabla_X gamma) (Y,Z) = 0,
$
i.e. that $mnabla gamma = 0$ and hence the induced connection is compatible with the induced metric. Since the algebra is quite dense here, let us outline the individual steps in the @inducedMetricCompatibilityDerivation[derivation]: 

+ Make use of @metricConnectionIdentity[identity];
+ Insert definitions, $gamma = iota^* g$ and $mnabla_X Y = (iota_*)#h(0em)^(-1) nabla_(iota_*X)(iota_*Y)$;
+ Employ $(iota^* g)(X,Y) = g(iota_* X, iota_* Y)$
+ Apply orthogonality of $P$, i.e. $g(P X,Y) = g(X,P Y)$;
+ Use @compatibleMetricConnectionIdentity[identity]---here, the assumption that $nabla$ is metric-compatible enters;
+ Employ $(iota^* g)(X,Y) = g(iota_* X, iota_* Y)$ as well as $nabla_(iota_* X) phi = (iota_* X)[phi] = X[phi] = mnabla_X phi$ on functions $phi in C^infty (cal(M))$.

Though somewhat algebra-heavy and light on geometric intuition, the above derivations yield a powerful result: The induced connection $mnabla$ on the foliation $Sigma$---induced by the Levi-Civita connection $nabla$ on $cal(M)$---is of Levi-Civita type as well. In particular, it follows that the induced connection components are given by the Christoffel symbol expression for $gamma$, i.e.,
$
  tensor(macron(Gamma),+k,-i j) = 1/2 gamma^(k ell)(gamma_(ell i, j) + gamma_(ell j, i) - gamma_(i j, ell)).
$<christoffelSymbolsInduced>
This is a direct construction of the connection coefficients that makes no reference to the ambient connection.
== Example: Induced Connection on the Foliation of $RR^3 without {0}$ into Spheres <sectionInducedConnectionExample>
To see the machinery of induced connections as well as the inheritance of the Levi-Civita propery in action, in this section, we reconsider the example of $cal(M) = RR^3 without {0}$, foliated into concentric origin-centered spheres from @exampleR3ConcentricSpheres. Let us briefly reestablish the setting.

On $RR^3 without {0}$, in Cartesian coordinates $x^mu = (x,y,z)$, the Euclidean metric reads
$
  g = g_(mu nu) dx^mu otimes dx^nu = dx otimes dx + dy otimes dy + dz otimes dz.
$
We foliate $RR^3 without {0}$ by introducing the function
$
  r : cal(M) -> RR, quad r(p) = sqrt(x^2 + y^2 + z^2), quad p = (x,y,z)
$
and define the leaves of the foliation $Sigma = {Sigma_r}$ to be its level sets,
$
  Sigma_r_0 := {p in cal(M) | r(p) = r_0} = r_0 S^2.
$
In spherical coordinates $(r,y^i) = (r,theta,phi)$ on $RR^3 without {0}$, the induced metric on a leaf $Sigma_r$ was derived to be
$
  gamma = gamma_(i j) dq^i otimes dq^j = r^2 dtheta otimes dtheta + r^2 sin^2 theta dphi otimes dphi,
$
with nonzero components $gamma_(theta theta) = r^2$, $gamma_(phi phi) = r^2 sin^2 theta$. 

We equip the ambient manifold $RR^3 without {0}$ with the Levi-Civita connection, which in Cartesian coordinates has vanishing coefficients, i.e.
$
  tensor(Gamma,+lambda,-mu nu) = 1/2 g^(lambda rho) (g_(rho mu, nu) + g_(rho nu,mu) - g_(mu nu, rho)) = 0,
$
since the metric components $g_(x x) = g_(y y) = g_(z z) = 1$ are constant.

We now have two formulae to compute the components $tensor(macron(Gamma),+k,-i j)$ of the connection $mnabla$ induced by $nabla$ at our disposal; we can either use the @projectionConnectionCoeffs[projection formula], or calculate them directly using the @christoffelSymbolsInduced[Christoffel expression]. We now proceed to evaluate both.

The evaluation of the @projectionConnectionCoeffs[projection formula] requires us to compute the components of the pushforward matrix
$
  E^mu_i = (diff x^mu)/(diff q^i)
$
as well as its left-inverse $E^i_mu$ subject to the condition that $tensor(P,+mu,-nu) = E^mu_i E^i_nu$ is orthogonal. The components of the pushforward matrix read
$
  E^x_theta &= r cos theta cos phi, &quad&& E^x_phi &= -r sin theta sin phi,\
  E^y_theta &= r cos theta sin phi, &&& E^y_phi &= r sin theta cos phi,\
  E^z_theta &= -r sin theta,&&& E^z_phi &= 0.
$
The orthogonal left-inverse can be computed using
$
  E^i_mu = gamma^(i j) g_(mu nu) E_j^nu,
$
as then
$
  E^i_mu E^mu_k = gamma^(i j) underbrace(g_(mu nu) E^nu_j E^mu_k,=gamma_(j k)) = gamma^(i j) gamma_(j k) = delta^i_k,
$
ensures the left-inverse property and 
$
  P_(mu nu) = g_(mu lambda) tensor(P,+lambda,-nu) = g_(mu lambda) E^lambda_i E^i_nu = g_(mu lambda) E^lambda_i gamma^(i j) g_(nu rho) E^rho_j
$
symmetry of $P_(mu nu)$ which is equivalent to orthogonality. 

The components of the orthogonal left-inverse read
$
  E^theta_x  &= 1/r cos theta cos phi, &quad&& E_x^phi &= -1/(r sin theta) sin phi,\
  E^theta_y &= 1/r cos theta sin phi, &quad&& E^y_phi &= 1/(r sin theta) cos phi,\
  E^theta_z &= -1/r sin theta, &&& E^z_phi &= 0.
$
Since the ambient connection coefficients vanish in Cartesian coordinates, the induced components reduce to the inhomogeneous term in @projectionConnectionCoeffs[], that is,
$
  tensor(macron(Gamma),+k,-i j) = E^k_lambda diff_i E^lambda_j.
$
While entirely mechanical, the computation of all components is somewhat laborious. We'll carry out one as an illustration and state the remaining results without derivation. As our example, we choose to calculate
$
  tensor(macron(Gamma),+theta,-phi phi) &= E^theta_lambda diff_phi E^lambda_phi =  E^theta_x diff_phi E^x_phi + E^theta_y diff_phi E^y_phi \
  &= (1/r cos theta cos phi) diff_phi (-r sin theta sin phi) + (1/r cos theta sin phi) diff_phi (r sin theta cos phi)\
  &= - cos theta sin theta cos^2 phi - cos theta sin theta sin^2 phi= - cos theta sin theta \
  &= - 1/2 sin(2 theta)
$
The remaining non-vanishing connection coefficients are given by
$
  tensor(macron(Gamma),+phi,-theta phi) = tensor(macron(Gamma),+phi,-phi theta) = cot theta,
$
where $cot theta = (cos theta)/(sin theta)$. 

We now proceed to recompute these coefficients using the standard Christoffel formula. To do so, it is usually advantageous to first compute the components
$
  tensor(macron(Gamma),-k i j) = 1/2(g_(k i, j) + g_(k j,i) - g_(i j, k))
$
and to then raise the first index with the inverse metric. Before wildly starting to evaluate all possible index combinations, we should first examine the nature of the components $gamma_(i j)$ of the induced metric. There is only one component, $gamma_(phi phi) = r^2 sin^2 theta$, which is non-constant#footnote[Recall that on a leaf, $r$ is constant.] and hence has a chance of contributing to $macron(Gamma)_(k i j)$. It only depends on $theta$, and hence the only non-zero component of the partial gradient of $gamma_(i j)$ is
$
  diff_theta gamma_(phi phi) = 2 r^2 sin theta cos theta.
$
Consequently, a coefficient $macron(Gamma)_(k i j)$ can only be non-zero if one of the indices is $theta$ and the other two are $phi$. Due to symmetry in the last two indices, this leaves us with two options,
$
  macron(Gamma)_(theta phi phi) &= 1/2 lr((underbrace(gamma_(theta phi, phi),=0) + underbrace(gamma_(theta phi, phi),=0) - gamma_(phi phi, theta)),size:#35%) = - r^2 sin theta cos theta,\
  macron(Gamma)_(phi theta phi) = macron(Gamma)_(phi phi theta) &= 1/2 lr((gamma_(phi phi, theta) + underbrace( gamma_(phi theta,phi),=0) - underbrace(gamma_(phi theta, phi),=0)),size:#35%) =  r^2 sin theta cos theta.
$
Given the components of the inverse metric,
$
  gamma^(theta theta) = 1/(r^2), quad gamma^(phi phi) = 1/(r^2 sin^2 theta),
$
we can compute the connection coefficients as
$
  tensor(macron(Gamma),+theta,-phi phi) &= g^(theta theta) macron(Gamma)_(theta phi phi) = -sin theta cos theta = -1/2 sin(2 theta),\
  tensor(macron(Gamma),+phi,-theta phi) = tensor(macron(Gamma),+phi, -phi theta) &= g^(phi phi)macron(Gamma)_(phi phi theta) = (cos theta)/(sin theta) = cot theta.
$
This reproduces exactly the same result as the projection formula, as we would expect from our derivations made for general manifolds---it is nevertheless satisfying to see that the abstract machinery does indeed work when applied to concrete examples.

= Curvature
In this section, we examine the different kinds of curvature that arise in the study of submanifolds and foliations. We begin by reviewing the definition of the Riemann curvature tensor and its contractions, which will allow us to define both the _ambient curvature_ of a manifold $cal(M)$ and the _instrinsic curvature_ of the leaves of a foliation $Sigma$ on it. The ambient curvature is defined via the ambient connection $nabla$ on $cal(M)$, whereas the intrinsic curvature is derived from the induced connection $nabla$, which---as established in the previous section---is the tangential projection of $nabla$ onto $Sigma$. 

This immediately suggests a relationship between the ambient and intrinsic curvatures. One might naïvely expect, by analogy  with the connections, that the intrinsic curvature is simply the projection of the ambient one. A simple counterexample will demonstrate that this cannot be entire picture. This leads us to the notion of _extrinsic curvature_, which---roughly speaking---captures the normal component of the ambient connection $nabla$ that is discarded when passing to $mnabla$. 

This will prepare us for the derivation of the Gauss equation in the next section, which relates the intrinsic, extrinsic and projected ambient curvatures in a precise and elegant way.
== Intrinsic Curvature
=== Curvature of Manifolds: the Riemann Tensor
On a space like $RR^2$ equipped with the Euclidean metric,
$
  g = dx otimes dx + dy otimes dy,
$
the Christoffel symbols of the Levi-Civita connection vanish in Cartesian coordinates. This means that the covariant derivatives of coordinate basis vectors vanish,
$
  nabla_mu diff_nu = 0.
$<basisParallel>
Informally speaking, the basis vectors are constant with respect to the connection---vectors with this property are called _parallel_. 

By contrast, consider the $2$-sphere $S^2$ with its standard round metric,
$
  g = dtheta otimes dtheta + sin^2 theta dphi otimes dphi.
$
Here, the Christoffel symbols do not vanish, and more fundamentally, there is no coordinate system in which the basis vectors are parallel, i.e., in which @basisParallel[equation] holds. This reflects a key geometric difference: the tangent spaces on $S^2$ "tilt" as one moves across the surface, and vectors transported between them must adjust accordingly. 

This geometric tilt can be detected via _parallel transport around a closed loop_. In $RR^2$, transporting a vector around any closed path will return it unchanged to its starting point. On $S^2$, however, the result typically differs: the transported vector may fail to return aligned with the original. This discrepancy encodes the _curvature_ of the manifold. Let us now formalise this by introducing the Riemann curvature tensor.
\ \
*Definition* (Riemann Curvature Tensor) Let $cal(M)$ be a smooth manifold equipped with an affine connection $nabla$, and let $X,Y,Z in Gamma(T cal(M))$ be vector fields. The _Riemann curvature tensor_ is the map 
$
  R : Gamma(T cal(M)) times Gamma(T cal(M)) times Gamma(T cal(M)) -> Gamma(T cal(M)),
$
defined by
$
  R(X,Y)Z = nabla_X nabla_Y Z - nabla_Y nabla_X Z - nabla_[X,Y]Z.
$
The first two terms compare the changes $Z$ undergoes when transported along $Y$ and then $X$, versus along $X$ and then $Y$, respecitvely---they measure the failure of the transported versions of $Z$ to align when moving around a parallelogram. More algebraically speaking, their difference measures the non-commutativity of covariant derivatives. However, since the flows of $X$ and $Y$ do not necessarily form a closed parallelogram, the discrepancy induced by the transport along the gap, i.e. the commutator $[X,Y]$, must also be accounted for. This is implemented by the final term, $nabla_[X,Y]Z$.

We define the components of the Riemann tensor through
$
  R(diff_mu, diff_nu) diff_lambda = tensor(R,+rho, -lambda mu nu) diff_rho,
$<RXYonBasisVectors>
such that on arbitrary vector fields $X,Y,Z in Gamma (T cal(M))$ we have
$
  R(X,Y)Z = tensor(R,+rho,-lambda mu nu) X^mu Y^nu Z^lambda.
$
In particular, along coordinate directions $diff_mu$, we have
$
  R(diff_mu, diff_nu) Z = [nabla_mu, nabla_nu] Z = tensor(R,+rho,-lambda mu nu) Z^lambda diff_rho
$<RiemannTensorComponentDef>
---the $nabla_[diff_mu, diff_nu]$ term vanishes since $[diff_mu,diff_nu] = 0$.

Notice that for fixed $X,Y$, the expression $R(X,Y)$ can naturally be viewed as a linear operator on $Gamma(T cal(M))$. Beyond that, its action can be extended to arbitrary tensor fields, by defining
$
  R:Gamma(T cal(M)) times Gamma(T cal(M)) times Gamma(T^((r,s))cal(M)) -> Gamma(T^((r,s)) cal(M)),
$
$
  R(X,Y)T = nabla_X nabla_Y T - nabla_Y nabla_X T - nabla_[X,Y] T.
$
Let us now examine the consequences of this definition in closer detail. 
\ \
*Remarks*
- *Annihilation of Functions* In particular, on functions $phi in C^infty (cal(M))$, we have 
  $
    R(X,Y)phi &= nabla_X nabla_Y phi - nabla_Y nabla_X phi - nabla_[X,Y] phi\
              &= underbrace(X[Y[phi]] - Y[X[phi]],= [X,Y][phi]) - [X,Y][phi]\
              &= 0.
  $
  That is, $R(X,Y)$ annihilates functions. 

- *Derivation Property* The linear operator $R(X,Y)$ on $Gamma(T^((r,s))cal(M))$ is a _derivation_: for arbitrary tensor fields $T$ and $S$, it satisfies the Leibniz rule
  $
    R(X,Y)(T otimes S) = (R(X,Y) T) otimes S + T otimes (R(X,Y)S).
  $<RXYLeibnizRule>
  To see this, note that the curvature operator $R(X,Y)$ can be written in the compact form
  $
    R(X,Y) = [nabla_X, nabla_Y] - nabla_[X,Y].
  $
  The second term $nabla_[X,Y]$ is plainly a derivation, as it is a covariant derivative. It therefore suffices to verify the Leibniz property for the commutator $[nabla_X, nabla_Y]$. Consider $nabla_X nabla_Y (T otimes S)$. Using the product rule for $nabla$, we compute
  $
    &nabla_X nabla_Y (T otimes S) = nabla_X ((nabla_Y T) otimes S + T otimes nabla_Y S)\
    &= (nabla_X nabla_Y T) otimes S + underbrace((nabla_Y T) otimes (nabla_X S) + (nabla_X T) otimes (nabla_Y S),"symmetric in" X "and" Y) + T otimes (nabla_X nabla_Y S),
  $
  Note that the middle two terms together are symmetric in $X$ and $Y$, and hence cancel upon antisymmetrisation. Therefore, 
  $
    [nabla_X, nabla_Y] (T otimes S) = ([nabla_X, nabla_Y]T) otimes S + T otimes ([nabla_X, nabla_Y] S).
  $
  Combining this with the Leibniz property of $nabla_[X,Y]$, we conclude that $R(X,Y)$ is indeed a derivation. 

- *Action on 1-Forms* The derivation property of $R(X,Y)$ allows the explicit computation of its action on 1-forms, and thus the extension to arbitrary tensors. Given a 1-form $omega in Gamma(T^* cal(M))$ and a vector field $Z in Gamma(T cal(M))$, we have
  $
    0 &= R(X,Y) omega(Z) = R(X,Y) (tr^1_1 omega otimes Z) = tr_1^1 R(X,Y) (omega otimes Z)\
    &= tr_1^1 [(R(X,Y) omega) otimes Z + omega otimes (R(X,Y)Z)]\
    &= (R(X,Y)omega)(Z) + omega(R(X,Y)Z).
  $
  This demonstrates that the $1$-form $R(X,Y)omega$ is given by
  $
    (R(X,Y)omega)(Z) = - omega(R(X,Y)Z),
  $
  which characterises the action of curvature on basis $1$-forms as
  $
    R(diff_mu,diff_nu) dx^rho = -tensor(R,+rho,-lambda mu nu) dx^lambda.
  $
  This is the dual to the action on basis vectors, @RXYonBasisVectors[identity].
- *Action on Arbitrary Tensors*
  The action of $R(X,Y)$ on functions, basis vectors $diff_mu$ and basis 1-forms $dx^mu$, as well as the @RXYLeibnizRule[Leibniz rule] allow for the evaluation of its action on arbitrary tensors. For example, on a (0,2)-tensor $T in Gamma(T^((0,2))cal(M))$, we have
  $
    R(diff_mu,diff_nu) T &= R(diff_mu, diff_nu) (T_(rho sigma) dx^rho otimes dx^sigma)\ &= T_(rho sigma) (R(diff_mu, diff_nu)dx^rho) otimes dx^sigma + dx^rho otimes (R(diff_mu, diff_nu))dx^sigma)\
    &= T_(rho sigma) (-tensor(R,+rho,-lambda mu nu) dx^lambda otimes dx^sigma - tensor(R,+sigma,-lambda mu nu) dx^rho otimes dx^sigma)\
    &= (- T_(lambda sigma) tensor(R,+lambda,-rho mu nu) - T_(rho lambda) tensor(R,+lambda,-sigma mu nu)) dx^rho otimes dx^sigma.
  $
  In general, the action of $R$ on tensor components involves contracting each index of the tensor with either the first or second slot of $tensor(R,+rho,-lambda mu nu)$, depending on the variance: contravariant indices contract against the first index $rho$, acquiring a minus sign, while covariant indices contract against the second index $lambda$ with a positive sign.

- *Relationship to Connection Coefficients* Since the Riemann tensor is defined through the connection, there exists a relationship between its components and the connection coefficients appearing in 
  $
    nabla_mu diff_nu = tensor(Gamma,+lambda,-mu nu) diff_lambda.
  $
  We may derive this relationship from @RiemannTensorComponentDef[equation]. Before anti-symmetrisation, the second covariant derivatives that appear read
  $
    nabla_mu nabla_nu diff_lambda &= nabla_mu (tensor(Gamma,+rho,-nu lambda) diff_rho) = (diff_mu tensor(Gamma,+rho,-nu lambda)) diff_rho + tensor(Gamma,+rho,-nu lambda) nabla_mu diff_rho\
    &= (diff_mu tensor(Gamma,+rho,-nu lambda) + tensor(Gamma,+rho,-mu sigma) tensor(Gamma,+sigma,-nu lambda) ) diff_rho.
  $
  Anti-symmetrising this expression in $mu$ and $nu$ leads us to
  $
    R(diff_mu, diff_nu) diff_lambda = [nabla_mu,nabla_nu] diff_lambda = (diff_mu tensor(Gamma,+rho,-nu lambda) - diff_nu tensor(Gamma,+rho,-mu lambda) + tensor(Gamma,+rho,-mu sigma) tensor(Gamma,+sigma,-nu lambda) - tensor(Gamma,+rho,-nu sigma) tensor(Gamma,+sigma,-mu lambda)) diff_rho,
  $
  from which we identify the components of the Riemann tensor as
  $
    tensor(R,+rho,-lambda mu nu) = diff_mu tensor(Gamma,+rho,-nu lambda) - diff_nu tensor(Gamma,+rho,-mu lambda) + tensor(Gamma,+rho,-mu sigma) tensor(Gamma,+sigma,-nu lambda) - tensor(Gamma,+rho,-nu sigma) tensor(Gamma,+sigma,-mu lambda).
  $<riemannTensorComponents>
=== Symmetries and Contractions of the Riemann Tensor
The Riemann curvature tensor introduced in the previous section exhibits a variety of intrinsic algebraic symmetries and contraction properties that depend on the underlying connection's characteristics. In this section, we summarise and derive these symmetries, distinguishing those that hold purely by definition from those that arise under additional assumptions such as metric compatibility and the absence of torsion. This will enable us to identify a canonical contraction of the Riemann tensor---the Ricci tensor---which is unique up to sign for the Levi-Civita connection, along with its trace, the _scalar curvature_.

To discuss algebraic symmetries, i.e. relations between permutations of the tensor slots, we first introduce the fully covariant $(0,4)$-tensor associated to the Riemann curvature. This is defined by lowering the upper index of the $(1,3)$-curvature operator using the metric by defining
$
  R(W,Z,X,Y) := g(R(X,Y)Z,W).
$
In components, this corresponds to
$
  R_(rho lambda mu nu) W^rho Z^lambda X^mu Y^nu = g_(rho sigma) (tensor(R,+rho,-lambda mu nu) X^mu Y^nu Z^lambda) W^sigma,
$
or equivalently,
$
  R_(rho lambda mu nu) = g_(rho sigma) tensor(R,+rho,-lambda mu nu).
$
Thus the fully covariant tensor arises simply by lowering the vector index of the curvature operator via the metric---in this sense, one could also write
$
  R(dot, Z,X,Y) = (R(X,Y)Z)^flat.
$
We are now ready to discuss the symmetries of $R(W,Z,X,Y)$: 

+ *Anti-Symmetry in Second Pair* By the anti-symmetry of the curvature operator, $R(X,Y)=-R(Y,X)$, which follows directly from the definition and the skew-symmetry of the commutator, 
  $
    R(X,Y) = [nabla_X,nabla_Y] - nabla_[X,Y] = -R(Y,X),
  $
  we obtain anti-symmetry in the last two slots of the fully covariant tensor, i.e.
  $
    R(W,Z,X,Y) = - R(W,Z,Y,X).
  $
  This holds for any connection. In components, we have
  $
    R_(rho sigma mu nu) = -R_(rho sigma nu mu) quad <=> quad tensor(R,+rho,-sigma mu nu) = -tensor(R,+rho,-sigma nu mu).
  $

+ *Anti-Symmetry in the First Pair* If the connection is metric-compatible (i.e. $nabla g =0$), the curvature operator acts on the metric as $R(X,Y)g = 0$.
  By the derivation propery of $R(X,Y)$, for any vector fields $W,Z in Gamma(T cal(M))$, we have#footnote[This makes use of the property that $D g(X,Y) = (D g)(X,Y) + g(D X,Y) + g(X,D Y)$ for any derivation $D$. This is a generalisation of what we have shown before in @proofDerivationOfMetric[eq.] for the particular case of $D = nabla$.]
  $
    0 &= R(X,Y) g(Z,W)\ &= underbrace((R(X,Y)g),=0)(Z,W) + underbrace(g(R(X,Y)Z,W),=R(W,Z,X,Y)) + underbrace(g(R(X,Y)W,Z),=R(Z,W,X,Y))\
    &=R(W,Z,X,Y) + R(Z,W,X,Y)
    $
  where the first equality with zero holds because $R(X,Y)$ is acting on a function. The above can be rearranged for 
  $
    R(W,Z,X,Y) = - R(Z,W,X,Y),
  $
  which proves anti-symmetry in the first pair of arguments of the covariant tensor. In components, this reads
  $
    R_(rho sigma mu nu) = -R_(sigma rho mu nu).
  $
+ *Bianchi Identity* We first derive an identity relating the sum of the cyclic permutations of $R(X,Y)Z$ with $X,Y,Z in Gamma(T cal(M))$ and the torsion tensor. This identity will then produce a symmetry in the case of zero torsion. Note that for the torsion, we have 
  $
    nabla_X Y - nabla_Y X = T(X,Y) + [X,Y]
  $
  using which we can begin writing out
  $
    &R(X,Y)Z + R(Y,Z)X + R(Z,X)Y\
    &= mhighlight(nabla_X nabla_Y Z) - nabla_Y nabla_X Z - nabla_[X,Y]Z\
    &quad+ nabla_Y nabla_Z X - nabla_Z nabla_Y X mhighlight(- nabla_[Y,Z] X)\
    &quad+ nabla_Z nabla_X Y mhighlight(- nabla_X nabla_Z Y) - nabla_[Z,X] Y\
    &= mhighlight(nabla_X (nabla_Y Z - nabla_Z Y) - nabla_[Y,Z] X) + limits(#scale(150%, rotate(180deg,[$arrow.cw$])))_(\ X\,Y\,Z)\
    &= nabla_X T(Y,Z) + underbrace(nabla_X [Y,Z] - nabla_[Y,Z] X,=T(X,[Y,Z]) + [X,[Y,Z]])\
    &= nabla_X T(Y,Z) + T(X,[Y,Z]) + [X,[Y,Z]] + limits(#scale(150%, rotate(180deg,[$arrow.cw$])))_(\ X\,Y\,Z)\
  $
  When expanding the sum over cyclic permutations of $X,Y,Z$ indicated by $#rotate(180deg,[$arrow.cw$])$, the last term drops out by the Jacobi identity for the commutator,
  $
    [X,[Y,Z]] + [Y,[Z,X]] + [Z,[X,Y]] = 0.
  $
  We hence arrive at the identity
  $
    R(X,Y)Z + R(Y,Z)X + R(Z,X)Y = nabla_X T(Y,Z) + T(X,[Y,Z]) + limits(#scale(150%, rotate(180deg,[$arrow.cw$])))_(\ X\,Y\,Z).
  $
  In the case of vanishing torsion, $T=0$, this turns into the symmetry
  $
    R(X,Y)Z + R(Y,Z)X + R(Z,X)Y = 0,
  $
  and equivalently,
  $
    R(W,Z,X,Y) + R(W,X,Y,Z) + R(W,Y,Z,X) = 0.
  $
  or in terms of components,
  $
    tensor(R,+rho,-[sigma mu nu]) = tensor(R, -rho, -[sigma mu nu]) = 0.
  $
  Here, the cyclic permutation of the last three indices/slots is proportional to their anti-symmetrisation since we have anti-symmetry in the last two indices.

+ *Symmetry in First and Second Pair* For a connection that is both metric-compatible and has vanishing torsion---i.e., a Levi-Civita connection---the above symmetries imply a further symmetry,
  $
    R(W,Z,X,Y) = R(X,Y,W,Z).
  $
  That is, the expression is symmetric under the exchange of the first and second pair of slots. This is shown by repeatedly applying the symmetries 1., 2. and 3. to the left-hand side:
  #bottom-number($
    R(W,Z,X,Y) &attach(=,t:2.) -R(Z,W,X,Y)\
    &attach(=,t:3.) R(Z,X,Y,W) + R(Z,Y,W,X)\
    &attach(=,t:2.) -R(X,Z,Y,W) - R(Y,Z,W,X)\
    &attach(=,t:3.) R(X,Y,W,Z) + R(X,W,Z,Y) + R(Y,W,X,Z) + underbrace(R(Y,X,Z,W),attach(=,t:"1. & 2.")R(X,Y,W,Z))\
    &attach(=,t:2.) 2R(X,Y,W,Z) - R(W,X,Z,Y) - R(W,Y,X,Z)\
    &attach(=,t:3.) 2R(X,Y,W,Z) + R(W,Z,Y,X)\
    &attach(=,t:1.) 2R(X,Y,W,Z) - R(W,Z,X,Y).
  $)
  The claim now follows from adding $R(W,Z,X,Y)$ to both sides. In terms of the tensor components, the symmetry reads
  $
    R_(rho sigma mu nu) = R_(mu nu rho sigma). 
  $
  This is certainly not the most elegant or efficient way to establish this symmetry. However, since the result follows directly from the symmetries imposed by metric compatibility and vanishing torsion, and offers little in terms of geometric insight, a somewhat brute-force algebraic proof suffices for our purposes.

The Riemann curvature, being a rank 4 tensor, allows for contractions between its slots. In principle, for 4 slots, there are twelve combinations for contractions---however, symmetries reduce these significantly. For a metric-compatible connection, we have anti-symmetry in the first and second pairs. This means that the contractions over these pairs vanish,
$
  tr^1_2 R(dot,dot,X,Y) = tr^1_2 R(W,Z,dot,dot) = 0.
$
Hence, for a contraction not to vanish, it must contract over one slot in the first and one in the second pair. Without loss of generality, this contraction can be performed over the first an third slot, as all others are related by signs (due to anti-symmetry in the first and second pair). This contraction defines a new rank 2 tensor called the _Ricci tensor_, denoted by#footnote[Here, the contraction over two covariant slots is to be interpreted as with respect to the inverse metric.]
$
  Ric(X,Y) = tr^1_3 R(dot,X,dot,Y).
$
In components, this reads
$
  Ric(X,Y)= R_(mu nu) X^mu Y^nu,
$
where the Ricci tensor components $R_(mu nu)$ emerge as a contraction of the Riemann tensor,
$
  R_(mu nu) = tensor(R,+lambda,-mu lambda nu).
$
If further, the connection is torsion-free, we have symmetry of the covariant Riemann tensor under the exchange of the first and second pairs, making the Ricci tensor symmetric, i.e.,
$
  Ric(X,Y) = Ric(Y,X).
$
At this point, we can perform a second contraction to obtain the so-called Ricci scalar
$
  cal(R) = tr^1_2 Ric(dot,dot) = g^(mu nu) R_(mu nu).
$
For a Levi-Civita connection, it is the unique scalar contraction of the Riemann tensor#footnote[The unique scalar linear in the Riemann tensor. Of course, there are other contractions possible at higher orders.].
=== Curvature of the Induced Connection
Previously, we have defined the induced connection $mnabla$ on the submanifolds of a foliation $Sigma$ of a smooth manifold $cal(M)$ equipped with the ambient connection $nabla$. The ambient connection gives rise to the Riemann curvature tensor
$
  R : Gamma(T cal(M)) times Gamma(T cal(M)) times Gamma(T cal(M)) -> Gamma(T cal(M)),\ 
$
with
$
  R(X,Y)Z = [nabla_X,nabla_Y] Z - nabla_[X,Y] Z.
$
In this context, it is referred to as the _ambient curvature_. 

Since each of the foliation's submanifolds are equipped with the induced connection $mnabla$---the projection of the ambient connection $nabla$ onto the tangent bundle $T Sigma$---we can define a Riemann curvature tensor with respect to it on each of the leaves of the foliation, leading to
$
  macron(R) : Gamma(T Sigma) times Gamma(T Sigma) times Gamma(T Sigma) -> Gamma(T Sigma)
$
with
$
  macron(R)(X,Y)Z = [mnabla_X,mnabla_Y] Z - mnabla_[X,Y]Z.
$
We refer to $macron(R)$ as the _intrinsic (Riemann) curvature_ of the foliation $Sigma$. On any individual leaf $Sigma_(t^A)$, $macron(R)$ is simply the Riemann tensor of the induced Levi-Civita connection---that is, the curvature one would assign having access only to the intrinsic geometry, i.e., without any knowledge of the ambient geometry.

Since both metric compatibility and absence of torsion are properties $mnabla$ inherits from $nabla$, the fully covariant tensors $R(W,Z,X,Y)$ as well as $macron(R)(W,Z,X,Y)$ exhibit the same algebraic symmetries (cf. previous section). In the case that $nabla$ is the Levi-Civita connection associated to the metric $g$ on $cal(M)$, then $mnabla$ is the Levi-Civita connection associated to the induced metric $gamma = iota^* g$ on the leaves of the foliation. This makes the components of $macron(R)$ computable entirely from the components of the induced metric, due to the relationships
$
  tensor(macron(R),+k,-ell i j) &= diff_i tensor(macron(Gamma),+k,-j ell) - diff_j tensor(macron(Gamma),+k,-i ell) + tensor(macron(Gamma),+k,-i m) tensor(macron(Gamma),+m,-j ell) - tensor(macron(Gamma),+k,-j m) tensor(macron(Gamma),+m,-i ell),\
    tensor(macron(Gamma),+k,-i j) &= 1/2 gamma^(k ell)(gamma_(ell i, j) + gamma_(ell j, i) - gamma_(i j, ell))
$
we have established in preceding sections. We may define an intrinsic Ricci curvature tensor as
$
  mRic(X,Y) := tr^1_3 macron(R)(dot, X, dot, Y),
$
and an intrinsic scalar curvature,
$
  macron(cal(R)) := tr^1_2 mRic(dot,dot) = gamma^(i j) macron(R)_(i j),
$
where $macron(R)_(i j)$ are the components of the induced Ricci tensor, given by
$
  macron(R)_(i j) = tensor(macron(R),+k,-i k j).
$

=== Intrinsic vs Projected Ambient Curvature <intrinsicVsProjectedAmbient>
A natural question following the above introduction of the intrinsic curvature $macron(R)$ on a foliation $Sigma$ is how it relates to the ambient curvature $R$. Since $mnabla$ is defined as the projection of $nabla$ onto $T Sigma$, one might expect a similar relationship between the intrinsic and ambient curvature---something like
$
  macron(R) = P R|_(T Sigma),
$<naiveGuessAmbientIntrinsicCurvature>
where $P$ denotes the orthogonal projector onto $T Sigma$.

This would appear to follow the same pattern we observed for the connection: the action of $R(X,Y)Z$ is projected onto the tangent bundle $T Sigma$, and since $mnabla$ acts on $Sigma$, the result must also lie in $T Sigma$. That is, one might guess
$
  macron(R)(X,Y)Z = P(R(X,Y)Z), quad X,Y,Z in Gamma(T Sigma) subset Gamma(T cal(M)).
$
However, this is _not_ the case. There are (at least) two ways to see this---first, through a concrete counterexample; and second, via an algebraic derivation that, while less geometrically intuitive, revelas the deeper structure behind the failure of this naive guess. Let us consider both perspectives in turn.

- *Curved Submanifolds of Flat Manifolds* Let us reconsider a recurring example from these notes: The foliation by origin-centered spheres of the ambient manifold $RR^3 without {0}$, equipped with the Euclidean metric
  $
    g &= dx otimes dx + dy otimes dy + dz otimes dz\
      &= dr otimes dr + r^2 (dtheta otimes dtheta + sin^2 theta dphi otimes dphi),
  $
  ---where $x^mu = (x,y,z)$ are Cartesian and $y^alpha = (r,y^i) = (r,theta,phi)$ spherical coordinates. The connection under consideration is the Levi-Civita connection associated with this metric. Our goal is to compute and compare the ambient and induced Riemann tensor components.

  Observe that, in Cartesian coordinates, the Levi-Civita has vanishing coefficients, i.e.
  $
    tensor(Gamma,+lambda,-mu nu) = 0,
  $<cartesianChristoffelSymbols>
  as the metric components $g_(mu nu) = delta_(mu nu)$ are constant. From this, it follows immediately that the Riemann curvature tensor also vanishes in these coordinates,
  $
    tensor(R,+rho,-lambda mu nu) = 0,
  $
  due to the @riemannTensorComponents[formula]. In contrast to @cartesianChristoffelSymbols[equation], the above equation is tensorial---that is, it holds in any coordinate system---and we conclude that the ambient manifold is flat:
  $
    R(X,Y)Z = 0,quad forall X,Y,Z in Gamma(T cal(M)),
  $
  If our earlier naive @naiveGuessAmbientIntrinsicCurvature[guess] were correct, this would imply that
  $
    macron(R)(X,Y)Z = P(R(X,Y)Z) = 0,
  $
  and thus that the intrinsic curvature vanishes as well.

  This already casts doubt on our naive guess---after all, the word "sphere" does not readily evoke "flatness". Let us indulge this suspicion and compute the components of the intrinsic curvature explicitly to see where the discrepancy lies. We already derived the induced connection coefficients for this foliation in @sectionInducedConnectionExample, arriving at the following non-zero Christoffel symbols:
  $
      tensor(macron(Gamma),+theta,-phi phi) -sin theta cos theta, wide
  tensor(macron(Gamma),+phi,-theta phi) = tensor(macron(Gamma),+phi, -phi theta) = cot theta.
  $
  At first glance, computing the components of a rank-4 tensor like the Riemann curvature might seem daunting. However, in two dimensions---and with the Levi-Civita connection---the symmetries of the Riemann tensor drastically reduce the number of independent components. Specifically, the antisymmetry in both the first and second pair of indices,
  $
    R_(rho sigma mu nu) = -R_(rho sigma nu mu) = -R_(sigma rho mu nu),
  $
  implies that the indices $theta$ and $phi$ must each appear exactly once in both pairs. All valid permutations are then related by symmetry. Moreover, since the metric is diagonal, we only need to compute a single nontrivial component, say $tensor(R,+theta,-phi theta phi)$. We proceed by applying the standard formula:
  $
    tensor(macron(R),+theta,-phi theta phi) &= diff_theta tensor(macron(Gamma),+theta,-phi phi) - diff_phi underbrace(tensor(macron(Gamma),+theta,-theta phi),=0) + underbrace(tensor(macron(Gamma),+theta,-theta i) tensor(macron(Gamma),+i,-phi phi),=0) - tensor(macron(Gamma),+theta,-phi i) tensor(macron(Gamma),+i,-theta phi)\
    &= - diff_theta (sin theta cos theta )- tensor(macron(Gamma),+theta,-phi phi) tensor(macron(Gamma),+phi,-theta phi)\
    &= -cos^2 theta + sin^2 theta + underbrace(sin theta cos theta cot theta,=cos^2 theta)\
    &= sin^2 theta.
  $
  This is very clearly _not_ zero. We have thus found a counterexample to @naiveGuessAmbientIntrinsicCurvature[our naive guess]: the intrinsic curvature does not, in general, arise from a simple projection of the ambient curvature.

  This is a good point to take a step back and generalise the insight, in order to build further intuition for why our guess cannot be correct. What we have done is the following: we took a flat manifold, $RR^3 without {0}$, and foliated it into surfaces that are scaled copies of the 2-sphere. Intuitively, spheres possess curvature---this is evident from the fact that their normal vector field varies as one moves along their surface. Our @naiveGuessAmbientIntrinsicCurvature[guess], however, attempted to capture something quite different: it projected the ambient curvature tensor (which vanishes in this case) onto the tangent bundle of the foliation (where it still vanishes). The projection $P R|_(T Sigma)$ captures only the part of the _ambient_ curvature that is tangential to the foliation; it entirely neglects how the surface itself bends within the ambient space. In other words, this projection measures the curvature of the background in which the leaves of the foliation live, but not how those leaves curve within it. The normal field vector plays no role in this projection. Hence, while $P R|_(T Sigma)$ may contribute to the intrinsic curvature $macron(R)$, it clearly does not suffice to determine it completely: the way in which the surface curves relative to the background also generates intrinsic curvature.

- *Algebraic Argument* In @sectionDefInducedConnection, we introduced the induced connection on a foliation $Sigma$ by defining its action on vector fields as
  $
    mnabla_X Y = (iota_*)^(-1) nabla_(iota_* X) (iota_* Y),quad X,Y in Gamma(T Sigma),
  $
  where $iota_*$ is the pushforward of the inclusion map $iota:Sigma -> cal(M)$, and $(iota_*)^(-1)$ is the unique left-inverse with the property that
  $
    P = iota_* compose (iota_*)^(-1)
  $<eq7168>
  is the orthogonal projection from $T cal(M)$ onto $T Sigma$. This is a very precise definition, as it distinguishes $T Sigma$ and $im(iota_*) subset T cal(M)$ as separate objects.
  
  In practice, however---particularly when working with foliations---it is often more convenient to treat $T Sigma$ as a proper subbundle of $T cal(M)$, which we may do via the linear embedding map $iota_*$. From this perspective, the induced connection takes a simpler and more direct form:
  $
    mnabla_X Y = P nabla_X Y, quad X,Y in T Sigma.
  $<altDefnInducedConnection>
  Since $T Sigma subset T cal(M)$, both $X$ and $Y$ are valid inputs for the ambient connection $nabla$. The appearance of the full projector $P$ on the right-hand side is then a result of pushing forward the image of $(iota_*)^(-1)$ into $T cal(M)$, which gives rise to the @eq7168[combination]. We will adopt this more algebraic, embedded viewpoint for the remainder of the discussion, as it makes many derivations more transparent: @altDefnInducedConnection[equation] makes clear that the induced connection is simply the projection of the ambient connection onto the tangent bundle of the foliation.

  With this notational preface in place, we are now ready to give an algebraic argument for why the intrinsic curvature $macron(R)$ cannot, in general, be written as the orthogonal projection of the ambient curvature $R$ onto $T Sigma$.
  
  Inserting into the definition of $macron(R)$, we find that for $X,Y,Z in Gamma(T Sigma)$,
  $
    macron(R)(X,Y)Z &= [mnabla_X, mnabla_Y] Z - mnabla_[X,Y] Z\
    &= [P nabla_X, P nabla_Y] Z - P nabla_[X,Y] Z.
  $
  The projection of the ambient curvature, on the other hand, reads
  $
    P R(X,Y)Z = P [nabla_X, nabla_Y] Z - P nabla_[X,Y] Z,
  $
  so the difference lies in the first term. Let us now expand the commutator appearing in $macron(R)$ as
  $
    [P nabla_X, P nabla_Y] Z &= P nabla_X (P nabla_Y Z) - (X<->Y)\
    &= underbrace(P^2,=P) nabla_X nabla_Y Z + P (nabla_X P) nabla_Y Z - (X<->Y)\
    &= P [nabla_X, nabla_Y] Z + P lr(((nabla_X P)nabla_Y Z - (nabla_Y P)nabla_X Z), size:#130%)
  $
  We can see that there is an additional term involving covariant derivatives of the projector $P$---we hence obtain the identity 
  $
    macron(R)(X,Y)Z = P R(X,Y)Z + P lr(((nabla_X P)nabla_Y Z - (nabla_Y P)nabla_X Z), size:#130%).
  $<precursorGaussEqn>
  This demonstrates that the intrinsic curvature contains more than just the projection of the ambient curvature. The second term involves the covariant derivative of the projector $P$ and encodes how $T Sigma$ varies under parallel transport within the ambient manifold. It is, in fact, the algebraic seed of the _Gauss equation_, which we will derive in the following sections.

  However, to properly interpret this additional term, we must first introduce the concept of _extrinsic curvature_. As it stands, the expression
  $
    P lr(((nabla_X P) nabla_Y Z - (nabla_Y P) nabla_X Z),size:#130%)
  $<extraTermForIntrinsicCurvature>
  clearly introduces a dependence on how the submanifold is situated within the ambient space---it involves the derivative of the projector onto the tangent bundle, and hence reflects how that bundle varies from point to point. This suggests that the _shape_ of the submanifold is encoded here---the part that we found to be missing in the concrete example above. However, the geometric content of this term is not transparent in its current algebraic form. The introduction of extrinsic curvature will allow us to isolate and interpret this shape-dependence more clearly, ultimately leading to the _Gauss equation_, which links intrinsic, extrinsic and ambient curvature in a precise way.
== Extrinsic Curvature
=== Algebraic Motivation for Extrinsic Curvature
In this section, we introduce the concept of _extrinsic curvature_. Its connection to the @extraTermForIntrinsicCurvature[additional term] appearing in the expression for $macron(R)$ will not be immediate---we will establish that link later. For now, we treat it as a concept in its own right, motivated by the structure behind the ambient and induced connections.

Recall from the previous section that the induced connection can be written as
$
  mnabla_X Y = P nabla_X Y,
$
when we regard $T Sigma subset T cal(M)$ as an actual subbundle via the embedding $iota_*$. This form makes it clear that $mnabla$ is simply the tangential projection of the ambient connection. Crucially, this means that certain components of $nabla_X Y$ are _discarded_ in the process.

Which components are lost? The short answer is: the _normal_ ones---but this deserves a more precise formulation.

To that end, recall that the metric $g$ on $cal(M)$ allows a decomposition into tangential and normal projectors,
$
  g = P + Q,
$
when viewed as a $(0,2)$-tensor. Interpreted as $(1,1)$-tensors, this becomes the orthogonal identity decomposition
$
  id = P + Q,
$
where
$
  id = delta^mu_nu diff_mu otimes dx^nu, quad P = tensor(P,+mu,-nu) diff_mu otimes dx^nu,quad Q = tensor(Q,+mu,-nu) diff_mu otimes dx^nu.
$
This follows from raising an index of the metric, 
$
  g^(mu lambda) g_(lambda nu) = delta^mu_nu.
$
Now consider the action of the ambient connection on two tangent vector fields $X,Y in Gamma(T Sigma)$. Since $id = P + Q$, we may decompose
$
  nabla_X Y = P nabla_X Y + Q nabla_X Y = mnabla_X Y + Q nabla_X Y.
$
This naturally leads to the definition of the _extrinsic curvature_ (also called the _second fundamental form_#footnote[The first fundamental form is the induced metric, but I will not use this terminology in these notes.]) as
$
  K(X,Y) := Q nabla_X Y.
$<defnExtrinsicCurvature>
It captures the normal part of the ambient covariant derivative of one tangent vector along another. We can thus write the clean decomposition
$
  nabla_X Y = mnabla_X Y + K(X,Y),
$<relationAmbientInducedConnectionAndExtrCurv>
which splits the ambient connection into its tangential and normal components relative to the foliation $Sigma$.

This relationship can also be turned around:
$
  mnabla_X Y = nabla_X Y - K(X,Y).
$
While this may appear trivial at first glance, it offers a useful geometric perspective: $K(X,Y)$ is the _correction_ necessary to make $nabla_X Y$ tangent to $Sigma$. Such a correction becomes necessary when the tangent plane "tilts"---for instance, on a sphere, moving from one point to another causes the tangent plane to rotate. The full derivative $nabla_X Y$ will generally not remain in the new tangent plane, and thus must have parts of it "chopped off" to become tangent. That "chopped-off" part is $K(X,Y)$. It hence encodes the change in orientation of the tangent spaces---or equivalently, the change in the normal spaces---as one moves along the submanifold. This makes it clear why $K$ is rightly referred to as the _extrinsic curvature_: it measures how the submanifold is shaped _within_ the ambient manifold. This is exactly what we found to be missing when guessing that $macron(R) = P R|_(T Sigma)$: a measure of how the leaves themselves curve on top of the background. How $K(X,Y)$ relates to $macron(R)$ precisely remains to be established---but it certainly feels like a step in the right direction.

This idea is also reflected in the alternate expression
$
  K(X,Y) = Q nabla_X Y = nabla_X (Q Y) - (nabla_X Q) Y.
$
Since $Q Y= 0$ for $Y in Gamma(T Sigma)$, this reduces to
$
  K(X,Y) = -(nabla_X Q) Y.
$
This version highlights another aspect of the same idea: extrinsic curvature measures how the normal projection operator changes along $X$, when acting on $Y$. In other words, it tracks how the normal bundle twists and turns along the leaves of the foliation---consistent with our earlier interpretation in terms of tilting tangent spaces. 


=== Geometric Interpretation via the Normal Vector <extCurvatureNormalVector>
The definition of extrinsic curvature we presented in the previous section is valid for any foliation $Sigma$ of a smooth manifold $cal(M)$, regardless of the codimension of the leaves. In this section, we specialise to a hypersurface foliation, where the expression simplifies considerably.

In the case where $Sigma = Sigma_t$ is a hypersurface foliation generated by a scalar function $t in C^infty (cal(M))$, the normal bundle has one-dimensional fibres. From @sectionFoliationTangentBundleDecomposition, we recall that in this scenario, we are provided with a normal vector field $n^sharp$, associated with the normal 1-form $n = alpha dt$, which satisfies the relations
$
  g(n^sharp,n^sharp) = epsilon = pm 1, quad g(n^sharp,X) = 0,quad forall X in T Sigma.
$
Thus any normal vector field $N in Gamma(N Sigma)$ can be written as a scalar multiple of $n^sharp$, 
$
  N = lambda n^sharp, quad lambda in C^infty (cal(M)).
$
This means that for any hypersurface foliation, we may express $K(X,Y)$---because it is normal to $Sigma$---as
$
  K(X,Y) = k(X,Y)n^sharp,
$
where
$
  k(X,Y) := epsilon g(K(X,Y), n^sharp)
$
is the proportionality factor $lambda$ from the previous equation. The map $k(X,Y)$ is a scalar-valued bilinear map,
$
  k : Gamma(T Sigma) times Gamma(T Sigma) -> RR,
$
and in components, it can be written as
$
  k = k_(mu nu) dx^mu otimes dx^nu.
$
Since $k$ only takes vectors tangent to the foliation as inputs, in adapted coordinates $(t,y^i)$, we can also write it in terms of the cobasis associated with the transverse coordinates $y^i$ as
$
  k = k_(i j) dy^i otimes dy^j = k_(i j) E_mu^i E_nu^j dx^mu otimes dx^nu,
$
implying that $k$ is its own pullback; we have $k = iota^* k$. Its action on vector fields $X = X^mu diff_mu = X^i diff_i$ and $Y = Y^mu diff_mu = Y^i diff_i$, both in $Gamma(T Sigma)$, is given by
$
  k(X,Y) = k_(mu nu) X^mu X^nu = k_(i j) X^i X^j,
$
making use of the fact that $T Sigma subset T cal(M)$ is a proper subspace via the embedding $iota_*$.

We now examine the proportionality factor $k(X,Y)$ in more detail. Explicitly---by inserting the definition of extrinsic curvature $K(X,Y)$---it expands to
$
  k(X,Y) = epsilon g(K(X,Y),n^sharp) = epsilon g (Q nabla_X Y, n^sharp) = epsilon g(nabla_X Y, Q n^sharp) = epsilon g(nabla_X Y, n^sharp),
$
since $Q$ is orthogonal and acts on $n^sharp$ as the identity. For a metric-compatible connection (i.e., when $nabla g = 0$, an assumption we make from hereon out), we may further rewrite this as
$
  k(X,Y) = epsilon g(nabla_X Y, n^sharp) = epsilon g(Y, -nabla_X n^sharp).
$
Upon defining the _shape operator_#footnote[A priori, it is unclear that the domain is $Gamma(T Sigma)$ and not $Gamma(T cal(M))$. However, due to normalisation of $n^sharp$, we have $0 = nabla_X epsilon = nabla_X g(n^sharp,n^sharp) = -2g(S(X),n^sharp)$ which shows that $S(X)$ has no normal component.] $S:Gamma(T Sigma)-> Gamma(T Sigma)$
$
  S(X) = -nabla_X n^sharp,
$
this becomes 
$
  k(X,Y) = epsilon g(S(X), Y).
$
Explicitly, the components of $k$ are given by
$
  k_(mu nu) = -epsilon gamma_(nu lambda) nabla_mu n^lambda = -epsilon nabla_mu n_nu = -epsilon nabla_mu (alpha delta^t_nu),
$
or equivalently,
$
  k_(i j) = - epsilon E^mu_i E^nu_j nabla_mu n_nu,
$
with respect to the transverse coordinates $y^i$.

We have now derived a rich set of equations. Let us conclude by reflecting on the geometric meaning behind these expressions. The central identity
$
  k(X,Y) = epsilon g(S(X),Y), quad "where" quad S(X) = -nabla_X n^sharp,
$
and the relationship
$
  K(X,Y) = k(X,Y) n^sharp,
$
show that the extrinsic curvature is closely linked to the change in the normal vector field as one moves along the foliation. The shape operator $S(X)$ encodes how the normal vector evolves when traversing a leaf, revealing the curvature of the hypersurface within the ambient manifold. The scalar $k(X,Y)$ then measures the projection of this change along a tangent vector $Y$, providing a clear geometric interpretation of the curvature in the case of hypersurface foliations.

In summary, the extrinsic curvature provides a measure of how the leaves are embedded and deformed within the ambient manifold, with $k(X,Y)$ quantifying the degree of this deformation. In the specific case of hypersurfaces, this bending is elegantly described by the change of the normal vector field $n^sharp$, as the normal bundle has one-dimensional fibres, simplifying the geometric interpretation of the curvature.
=== Symmetry of the Extrinsic Curvature
We have previously observed that, under certain conditions such as metric compatibility or vanishing torsion, the Riemann curvature tensor acquires additional symmetries. A similar result holds for the extrinsic curvature $K(X,Y)$, which becomes a symmetric bilinear form when the connection is torsion-free.

Recall that the torsion tensor is defined as
$
  T(X,Y) = nabla_X Y - nabla_Y X - [X,Y], quad X,Y in Gamma(T cal(M)),
$
such that the requirement of vanishing torsion is equivalent to
$
  nabla_X Y = nabla_Y X + [X,Y].
$
Moreover, since the Lie bracket is a closed operation on the tangent spaces of (sub)manifolds, for $X,Y in T Sigma$ we also have $[X,Y]in T Sigma$.

With this in mind, it is rather straightforward to show that the extrinsic curvature tensor is symmetric. We simply derive
$
  K(X,Y) = Q nabla_X Y = underbrace(Q nabla_Y X,=K(Y,X)) + underbrace(Q [X,Y],=0) = K(Y,X).
$
Here, we made use of the fact that $Q$ is the normal projection onto $N Sigma$, whereas $[X,Y] in T Sigma$, so that $Q$ annihilates it. We have thus shown that the intrinsic curvature is symmetric, i.e.,
$
  K(X,Y) = K(Y,X).
$
For the case of a hypersurface foliation, we also obtain symmetry for the scalar factor $k(X,Y)$ in the expression $K(X,Y) = k(X,Y)n^sharp$. Specifically, we find
$
  k(X,Y) = k(Y,X),
$
showing that the scalar function $k(X,Y)$, which encodes the magnitude of the extrinsic curvature, is also symmetric.

=== Kosmann-Type Formula
In the case that the ambient manifold $cal(M)$ is endowed with a metric-compatible and torsion-free connection $nabla$, and $Sigma$ is a hypersurface foliation of it, one can derive a further identity relating the magnitude factor $k(X,Y)$ of the extrinsic curvature to the Lie derivative $cal(L)_(n^sharp)$ of the metric along the normal flow generated by $n^sharp$. 

To derive this, we begin by expanding
$
  (cal(L)_(n^sharp) g)(X,Y) = n^sharp [g(X,Y)] - g([n^sharp,X],Y) - g(X,[n^sharp,Y]),quad X,Y in Gamma(T Sigma),
$
where we made use of the fact that $cal(L)$ fulfils a Leibniz rule, that it acts as a directional derivative on scalars, and $cal(L)_X Y =[X,Y]$ on vectors. We now treat the individual terms separately.

For the first term, we may replace the directional with a covariant derivative, which after applying the product rule and metric compatibility yields
$
  n^sharp [g(X,Y)] = nabla_(n^sharp) g(X,Y) = g(nabla_(n^sharp)X,Y) + g(X, nabla_(n^sharp) Y).
$
The second (and third term, analogously), may be rewritten by making use of the identity
$
  [X,Y] = nabla_X Y - nabla_Y X
$
that holds for torsion-free connections. Inserting it, we find
$
  g([n^sharp,X],Y) = g(nabla_(n^sharp) X, Y) - g(nabla_X n^sharp, Y).
$
Putting everything back together, we arrive at
$
  (cal(L)_(n^sharp) g)(X,Y) &= cancelr(g(nabla_(n^sharp) X, Y)) - cancelr(g(nabla_(n^sharp) X, Y)) + underbrace(g(nabla_X n^sharp, Y),=-epsilon k(X,Y)) + (X<->Y)\
  &= -epsilon k(X,Y) - epsilon k(Y,X)\
  &= -2epsilon k(X,Y),
$
implying the identity
$
  k(X,Y) = -epsilon/2 (cal(L)_(n^sharp) g)(X,Y).
$

This gives a direct relationship between the extrinsic curvature and the metric; more precisely, it tells us that its magnitude corresponds to the change of the metric along the normal flow of the foliation.

=== Example: Extrinsic Curvature of the Foliation of $RR^3 without {0}$ into Spheres 

= Gauss-Codazzi-Mainardi Equations
Now that we established the concepts of ambient, intrinsic and extrinsic curvature, we can finally derive relationships between them. 
== Gauss Equation
In this section, we derive the Gauss equation, which we already touched upon very lightly at the end of @intrinsicVsProjectedAmbient. We first treat the general case, and then specialise to the case of a hypersurface foliation, where terms simplify somewhat. In the following, $cal(M)$ will denote a (pseudo-)Riemannian manifold equipped with a metric $g$, an affine connection $nabla$, and a foliation $Sigma$. The leaves of the foliation are endowed with the induced metric $gamma = iota^* g$ as well as the induced connection $mnabla$. Moreover, if not stated otherwise, $X,Y,Z,W$ will denote foliation-tangent vector fields in $Gamma(T Sigma)$.
=== General Case
To derive the Gauss equation, let us first collect some definitions and identities we have derived in the previous sections. The starting point of the derivation will be the @precursorGaussEqn[equation], 
$
  macron(R)(X,Y)Z = P R(X,Y)Z + P lr(((nabla_X P)nabla_Y Z - (nabla_Y P)nabla_X Z), size:#130%),
$<precursorGaussEqn2>
Further, we will need the definition and alternate characterisation of the extrinsic curvature,
$
  K(X,Y) = Q nabla_X Y = -(nabla_X Q) Y.
$
Lastly, the decomposition of the identity into the orthogonal projectors $P$ and $Q$,
$
  id = P + Q,
$
will be of use as well. 

The derivation is best performed using the fully covariant Riemann tensor
$
  macron(R)(W,Z,X,Y) = g(macron(R)(X,Y)Z,W).
$
We begin by inserting the @precursorGaussEqn2[precursor] into this, finding
$
  &#h(-3em) macron(R)(W,Z,X,Y) = g(macron(R)(X,Y)Z,W)\
  &= g lr((P lr((R(X,Y)Z + (nabla_X P)nabla_Y Z - (nabla_Y P) nabla_X Z), size: #135%), W),size:#135%)\
  &= g lr((R(X,Y)Z + (nabla_X P) nabla_Y Z - (nabla_Y P) nabla_X Z, underbrace(P W,=W)),size:#55%)\
  &= R(W,Z,X,Y) + lr((g lr(((nabla_X P) nabla_Y Z, W),size:#135%) - (X <-> Y)),size: #135%)
$<eq815>
The first term is clear; it is nothing but the ambient Riemann tensor applied to the vector tuple $(W,Z,X,Y)$. The tangential projection $P$ drops out due to the fact that it is orthogonal and $P W = W$, as $W in Gamma(T Sigma)$. The second term requires a more detailed examination. It is an anti-symmetrisation in $X,Y$, so it is sufficient to consider only
$
  g lr(((nabla_X P) nabla_Y Z, W), size:#135%)
$
and perform the anti-symmetrisation afterwards. The treatment of this expression requires a series of operations, the steps behind which we list for clarity below:
$
  g lr(((nabla_X P) nabla_Y Z, W), size:#135%) &attach(=,t:1.) -g lr(((nabla_X Q) nabla_Y Z, W), size:#135%)\
  &attach(=,t:2.) - g lr(( nabla_X underbrace((Q nabla_Y Z),=K(Y,Z)), W), size:#37%) + underbrace(g(Q nabla_X nabla_Y Z, W),=g lr((nabla_X nabla_Y, underbrace(Q Z,=0)),size:#35%) = 0)\
  &attach(=,t:3.) -nabla_X underbrace(g lr(( underbrace(K(Y,Z),in Gamma(N Sigma)), underbrace(W,in Gamma(T Sigma))),size:#40%),=0) + g(K(Y,Z),nabla_X W) \
  &attach(=,t:4.) g(K(Y,Z),K(X,W)).
$
The steps are as follows:

+ Replace $P$ using $P = id - thin Q$ and note that
  $
    nabla_X P = underbrace(nabla_X id,=0) - nabla_X Q = -nabla_X Q
  $
+ Integrate first argument by parts,
  $
    (nabla_X Q)nabla_Y Z = nabla_X (Q nabla_Y Z)- Q nabla_X nabla_Y Z
  $
+ Integrate by parts using
  $
    g(nabla_X A,B) = nabla_X g(A,B) - g(A, nabla_X B)
  $
  This step assumes metric compatibility of the connection.

+ First term vanishes since $K(Y,Z) perp W$. In the second term, the first argument lies in $Gamma(N Sigma)$ such that only the normal part $Q nabla_X W = K(X,W)$ of the second argument contributes. 

Inserting this partial result back into @eq815[equation] above and carrying out the anti-symmetrisation yields the _Gauss equation_ 
$
  macron(R)(W,Z,X,Y) = R(W,Z,X,Y) + g(K(X,W),K(Y,Z)) - g(K(Y,W),K(X,Z)).#h(2em)
$<gaussEqn>
This is the result we anticipated in @intrinsicVsProjectedAmbient, now fully worked out in terms of the extrinsic curvature. It relates the intrinsic curvature $R$ of the leaves of the foliation to the (pullback of) the ambient curvature $R$ and the extrinsic curvature $K$.
\ \
*Remarks:*
- *Metric Compatibility* The derivation of the Gauss equation made use of the assumption that the connection is metric-compatible. This reduces its generality, but also introduces anti-symmetry in the $W,Z$-pair for both $macron(R)(W,Z,X,Y)$ and $R(W,Z,X,Y)$. Consequently, the terms involving the extrinsic curvature must satisfy this anti-symmetry as well---to remain consistent. It is straightforward to see that this is the case; the exchange $(W<->Z)$ yields a negative sign.

- *Vanishing Torsion* Though we did not have to assume vanishing torsion in the derivation, we may add it as a further requirement. This introduces further symmetries of both the ambient and intrinsic curvature tensors, namely symmetry under the exchange of the first and second pair (since the connection is also metric-compatible), 
  $
    R(W,Z,X,Y) = R(X,Y,Z,W),
  $
  as well as the Bianchi identity
  $
    R(W,Z,X,Y) + R(W,X,Y,Z) + R(W,Y,Z,X) = 0.
  $
  The extrinsic curvature terms must satisfy these as well to remain consistent. Let us verify that this is the case for the former symmetry. Recall that for a torsion-free connection, $K$ is symmetric, i.e.
  $
    K(X,Y) = K(Y,X). 
  $
  This implies that
  #bottom-number($
    macron(R)(W,Z,X,Y) - R(W,Z,X,Y) &= g(K(X,W),K(Y,Z)) - g(K(Y,W),K(X,Z))\
    &= g(K(W,X),K(Z,Y)) - g(K(W,Y),K(Z,X))\
    &= macron(R)(X,Y,W,Z) - R(X,Y,W,Z),
  $)
  as required.
=== Hypersurface Foliation
The Gauss @gaussEqn[equation] can be specialised to hypersurface foliations by making use of the relationship 
$
  K(X,Y) = k(X,Y) n^sharp
$
that we derived in @extCurvatureNormalVector. Inserting this into the Gauss equation yields
$
  &macron(R)(W,Z,X,Y)\ 
  &= R(W,Z,X,Y) + g(k(X,W)n^sharp, k(Y,Z) n^sharp) - g(k(Y,W)n^sharp, K(X,Z)n^sharp)\
  &= R(W,Z,X,Y) + underbrace(g(n^sharp,n^sharp),=epsilon)(k(X,W) k(Y,Z) - k(Y,W)k(X,Z))\
  &= R(W,Z,X,Y) + epsilon(k(X,W) k(Y,Z) - k(Y,W)k(X,Z))
$
We have thus found the Gauss equation for hypersurface foliations,
$
  macron(R)(W,Z,X,Y) = R(W,Z,X,Y) + epsilon(k(X,W)k(Y,Z)-k(Y,W)k(X,Z)).
$<gaussEqnHypersurface>
In terms of components in the adapted coordinates $(t,y^i)$, this reads
$
  macron(R)_(k ell i j) = R_(k ell i j) + epsilon (k_(i k) k_(j ell) - k_(j k)k_(i ell)),
$
where 
$
  R_(k ell i j) = (iota^* R)_(k ell i j) = E^rho_k E^sigma_ell E^mu_i E^nu_j R_(rho sigma mu nu)
$
are the components of the pushforward of the ambient Riemann tensor. 
=== Example: Foliation of $RR^3 without {0}$ into Spheres 

== Codazzi-Mainardi Equation
=== General Case
=== Hypersurface Foliation
=== Example: Foliation of $RR^3 without {0}$ into Spheres 



= Ricci-Voss Identity
