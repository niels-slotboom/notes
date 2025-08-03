#import "template.typ": *
#import "macros.typ": *
#import "@preview/xarrow:0.3.1": xarrow

= Hamiltonian Mechanics
== Legendre Transform
=== One-Dimensional Case
Say we have a function $cal(L)(q,dot(q))$. Our goal is to encode the information it contains in terms of slopes with respect to $dot(q)$ in place of this function depending on $dot(q)$: this is done using a _Legendre transformation_. In the following, we suppress the dependence on $q$ since the Legendre transformation acts only on the $dot(q)$-dependence---the $q$-dependence is treated parametrically. The Legendre transformation does the following: Take a slope $p$, and let $cal(H)(p)$ be such that 
$
  p dot(q) - cal(H) (p) <= cal(L)(dot(q))
$
everywhere, with equality in exactly one point $dot(q)$ (This gives implicit assumptions on $cal(L)$, e.g. that it is strictly convex). This is equivalent to 
$
  p dot(q) - cal(L) (dot(q)) <= cal(H)(p). 
$
Since we have equality in exactly one point $dot(q)$ and $<$ everywhere else, the point of equality is a maximum (in the case of differentiable functions and that point not being at infinity), so  we must have
$
  cal(H)(p) = max_(dot(q)) [p dot(q) - cal(L)(dot(q))]. 
$<maximumCondition>
The maximality condition hence demands that
$
  (diff cal(L))/(diff dot(q)) (dot(q)(p)) = p,
$
where $dot(q)(p)$ is the point at which the maximum @maximumCondition[] is attained, i.e.
$
  dot(q)(p) = op("arg max", limits: #true)_(dot(q)) [p dot(q) - cal(L)(dot(q))].
$ 
For a function $cal(L)$ which is convex in $dot(q)$, this is an invertible expression, which implicitly defines $dot(q) = dot(q)(p)$ and hence
$
  cal(H)(p) = p dot(q)(p) - cal(L)(dot(q)(p)).
$<LegendreTransformExpression>
Now, one might wonder why $cal(H)(p)$, the Legendre transform of $cal(L)(dot(q))$, is of any interest. The answer to this question is: It contains the same information as $cal(L)(dot(q))$, just phrased differently: since $p = p(dot(q))$ is invertible, we may just as well write @LegendreTransformExpression[eq.] as 
$
  cal(L)(dot(q)) = p(dot(q)) dot(q) - cal(H)(p(dot(q))),
$
which recovers the original $cal(L)(dot(q))$. This implies that no information was lost in transitioning to $cal(H)(p)$, since we can recover the full $cal(L)(dot(q))$ from it.
=== Higher, Finite-Dimensional Case
In higher dimensions, say with $q = q^mu$, $dot(q) = dot(q)^mu$, the condition we start with is
$
  p_mu dot(q)^mu - cal(H)(p) <= cal(L)(dot(q)).
$
with the assumption that we have equality in exactly one point. Here, the $p_mu dot(q)^mu$ generates a hyperplane in place of just a line. The condition that we have equality in exactly one point $dot(q)$ again translates to
$
  cal(H)(p) = max_dot(q) [p dot(q) - cal(L)(dot(q))].
$
This time, the maximum is attained where
$
  (diff cal(L))/(diff dot(q)^mu) (dot(q)(p)) = p_mu.
$<maximumConditionHigherDim>
Under the condition that this relationship is invertible, which is equivalent to the condition 
$
  det((diff^2 cal(L))/(diff dot(q)^mu diff dot(q)^nu)) != 0
$
on the Jacobian of @maximumConditionHigherDim[relation] (or equivalently, the Hessian of $cal(L)$ with respect to $dot(q)$), we can express
$
  dot(q)^mu = dot(q)^mu (p).
$
Consequently, we have
$
  cal(H)(p) = p_mu dot(q)^mu (p) - cal(L)(dot(q)(p)).
$<multiDimHamiltonianAsLegendreTransform>
== Equations of Motion from Lagrangian Mechanics
In Lagrangian mechanics, we work with a Lagrangian $cal(L)(q, dot(q))$, where $q = q^mu$ are the coordinates of a particle. The Lagrangian defines an action functional via
$
S[q(tau)] = integral cal(L)(q, dot(q)) dtau,
$
where $q(tau)$ is a parametrisation of the path. The equations of motion are obtained by varying $q -> q + delta q$ and requiring $S -> S$, i.e. that $delta S = 0$. This leads to the Euler-Lagrange equations
$
(diff cal(L))/(diff q^mu) = d/(dtau) (diff cal(L))/(diff dot(q)^mu). 
$<ELeqns>
These are typically second-order differential equations, since $(diff cal(L))/(diff dot(q)^mu)$ usually depends on $dot(q)$.

We can now use these equations to derive the canonical equations of motion for the variables $q^mu$ and $p_mu$ associated with a Hamiltonian $cal(H)$. Define the canonical momentum as
$
p_mu = (diff cal(L))/(diff dot(q)^mu).
$
Inserting this into @ELeqns yields
$
dot(p)_mu = (diff cal(L))/(diff q^mu) = (diff)/(diff q^mu) [p_nu dot(q)^nu - cal(H)(q, p)] = - (diff cal(H))/(diff q^mu),
$
which is the equation of motion for $p_mu$.

To derive the equation for $q^mu$, we begin with the definition of the Hamiltonian as a Legendre transform,
$
cal(H)(q, p) = p_mu dot(q)^mu - cal(L)(q, dot(q)). 
$<multiDimHamiltonianAsLegendreTransform>
Taking the derivative with respect to $p_mu$, we obtain
$
(diff cal(H))/(diff p_mu) = (diff)/(diff p_mu) [p_nu dot(q)^nu - cal(L)(q, dot(q))] = dot(q)^mu + underbrace( p_nu (diff dot(q)^nu)/(diff p_mu) - (diff cal(L))/(diff dot(q)^nu) (diff dot(q)^nu)/(diff p_mu),=0) = dot(q)^mu,
$
where the last equality follows from the identity $p_nu = (diff cal(L))/(diff dot(q)^nu)$.

In summary, we have derived the canonical equations
$
dot(q)^mu = (diff cal(H))/(diff p_mu), quad dot(p)_mu = - (diff cal(H))/(diff q^mu),
$<canonicalEqns>
which are a set of first-order differential equations, in contrast to the second-order Euler-Lagrange @ELeqns[equations]. This makes them typically more convenient for numerical integration.
== Poisson Brackets
An important concept in the Hamiltonian formalism is that of Poisson brackets. For any two observables $A,B$, functions of $q,p$, i.e.
$
  A = A(q,p), quad B = B(q,p),
$
their _Poisson bracket_ is defined as
$
  {A,B} = (diff A)/(diff q^mu) (diff B)/(diff p_mu) - (diff B)/(diff q^mu) (diff A)/(diff p_mu).
$
Before we establish why this is useful, let us compute the Poisson brackets of the fundamental observables $q^mu$ and $p_mu$. This is rather straightforward:
$
  {q^mu, q^nu} &= (diff q^mu)/(diff q^lambda) underbrace((diff q^nu)/(diff p_lambda),=0) - (diff q^nu)/(diff q^lambda) underbrace((diff q^mu)/(diff p_lambda),=0) = 0,\
  {p_mu, p_nu} &= underbrace((diff p_mu)/(diff q^lambda),=0) (diff p_nu)/(diff p_lambda) - underbrace((diff p_mu)/(diff q^lambda),=0) (diff p_nu)/(diff p_lambda) = 0,\
  {q^mu, p_nu} &= underbrace((diff q^mu)/(diff q^lambda),=delta^mu_lambda) underbrace((diff p_nu)/(diff p_lambda),=delta^lambda_nu) - underbrace((diff p_nu)/(diff q^lambda) (diff q^mu)/(diff p_lambda),=0) = delta^mu_nu.
$
To be able to reference this more compactly later on, we repeat this result:
$
  {q^mu, q^nu} = {p_mu, p_nu} = 0, quad {q^mu, p_nu} = delta^mu_nu.
$<PoissonBracketqp>
We are now ready to see why the notion of the Poisson bracket is useful. To this end, we consider a particular observable---the Hamiltonian $cal(H)(q,p)$, and see how its Poisson bracket with another observable $A(q,p)$ behaves. To this end we recall the canonical @canonicalEqns[equations] and compute
$
  {A, cal(H)} = (diff A)/(diff q^lambda) underbrace((diff cal(H))/(diff p_lambda),=dot(q)^lambda) - underbrace((diff cal(H))/(diff q^lambda),=-dot(p)_lambda) (diff A)/(diff p_lambda) = (diff A)/(diff q^lambda) (d q^lambda)/(d tau) + (diff A)/(diff p_lambda) (d p_lambda)/(d tau) = dA/dtau.
$
This is a profound result, put more succinctly as
$
  dA/dtau = {A, cal(H)}.
$
Put into words, this means the following: the evolution of any observable $A(q(tau), p(tau))$ along the particle's path is given by its Poisson bracket with the Hamiltonian (as long as the particle obeys the equations of motion). In particular, for the observables $q^mu$ and $p_mu$ we confirm
$
  {q^mu, cal(H)} = delta^mu_lambda (diff cal(H))/(diff p_lambda) = dot(q)^mu, quad {p_mu, cal(H)} = - (diff cal(H))/(diff q^lambda) delta^lambda_mu = dot(p)_mu.
$
We have hence established the Poisson bracket ${dot, cal(H)}$ to be a valuable tool in obtaining the first-order evolution equation for any observable. One particular result that follows from the anti-symmetry of the Poisson bracket is that 
$
  dot(cal(H)) = {cal(H), cal(H)} = 0, 
$
implying that along any physical path, the Hamiltonian is conserved.
== Free Non-Relativistic Particle
We now have the necessary theory to begin looking at some examples. In this section, we consider a free particle---in the next, we will put it into a quadratic potential, turning the problem into a harmonic oscillator. A non-relativistic free particle has the Lagrangian
$
  cal(L) = 1/2 m  dot(q)^mu dot(q)_mu, 
$
assuming a flat Euclidean space (for simplicity, let us not try to do things with Minkowski-space for now. We will revisit that in Hamiltonian Field Theory). Here, $m$ denotes the mass of the particle. We want to first derive the Hamiltonian, and then the concrete form of the equations of motion for the free particle. For this, we need the conjugate momenta, which we compute as
$
  p_mu = (diff cal(L))/(diff dot(q)^mu) = m dot(q)_mu quad <=> quad dot(q)_mu = 1/m p_mu
$
This is the expected formula from classical mechanics: momentum equals mass times velocity. Using this momentum, we can compute the Hamiltonian,
$
  cal(H)(q,p) = p_mu dot(q)^mu - cal(L) (q, dot(q)) = 1/m p_mu p^mu - 1/(2 m) p_mu p^mu = 1/(2m) p_mu p^mu.
$
This is also the anticipated result---the Hamiltonian is the total energy of the system, which we know from classical mechanics to be $E = 1/2 m v^2 = p^2\/2m$. Let us now move to the equations of motion. The canonical equations can straightforwardly be derived as
$
  dot(q)^mu = (diff cal(H))/(diff p_mu) = 1/m p^mu, wide dot(p)_mu = -(diff cal(H))/(diff q^mu) = 0.
$
The latter equation has the simple interpretation that the momentum does not change. This makes sense for a free particle, since changes in momentum are caused by forces, of which there act none. The former equation is simply restating the relationship between velocity and momentum we had already found before.
== Harmonic Oscillator
To make our previous example a bit more exciting, we can add a potential term. In general, the Lagrangian of a particle in a potential $V$ has the form
$
  cal(L) = 1/2 m dot(q)^mu dot(q)_mu - V(q).
$
A harmonic oscillator has a particularly simple potential, with the full Lagrangian reading 
$
  cal(L) = 1/2 m dot(q)^mu dot(q)_mu - 1/2 m omega^2 q^mu q_mu.
$
The conjugate momenta evaluate to the same as in the free case,
$
  p_mu = (diff cal(L))/(diff dot(q)^mu) = m dot(q)_mu,
$
since no additional dependence on $dot(q)^mu$ is introduced into the Lagrangian. We can again compute the Hamiltonian, finding
$
  cal(H) = p_mu dot(q)^mu - cal(L) &= 1/m p_mu p^mu - 1/(2m) p_mu p^mu + 1/2 m omega^2 q_mu q^mu\ &= 1/(2m) p_mu p^mu + 1/2 m omega^2 q_mu q^mu
$ 
This should look familiar, from e.g. the Hamiltonian operator of a harmonic oscillator in quantum mechanics, where
$
  hat(cal(H)) = hat(p)^2/(2m) + 1/2 m omega^2 hat(q)^2
$
Returning to Hamiltonian mechanics, let us now compute the canonical equations. This is also rather straightforward, leading to
$
  dot(q)^mu = (diff cal(H))/(diff p_mu) = 1/m p^mu, wide dot(p)_mu = - (diff cal(H))/(diff q^mu) = m omega^2 q_mu.
$
These equations have by oscillatory solutions,
$
  q^mu (tau) &= A^mu cos(omega tau) + B^mu sin(omega tau),\
  p^mu (tau) &=  m omega [-A^mu sin(omega tau) + B^mu cos(omega tau)],
$
where $A^mu, B^mu in RR^n$ are arbitrary constant vectors determined by initial conditions.
= Hamiltonian Field Theory: Basics

== Review: Lagrangian Field Theory
Before constructing the Hamiltonian formalism for classical field theories, we first review the Lagrangian approach. As in mechanics, the Legendre transform will then allow us to transition to the Hamiltonian picture. This process introduces some complications, which we will address in subsequent sections.

Let us begin by recalling the essential elements of Lagrangian mechanics, and then draw analogies with Lagrangian field theory.

In Lagrangian mechanics, we consider a one-dimensional parameter $tau$ living in a connected subset of $RR$. The dynamical variables are the coordinate embeddings $x^mu (tau)$, and the Lagrangian $cal(L)$ is a function of these variables and their derivatives,
$
  cal(L) = cal(L)(x,dot(x),...),
$
typically involving only first-order derivatives. From this, we define the action functional,
$
  S[x] = integral dtau thin cal(L)(x(tau), dot(x)(tau)).
$
Infinitesimally varying the path as $x^mu -> x^mu + delta x^mu$ results in a variation of the action,
$
  delta S = integral dtau thin delta cal(L)(x,dot(x)),
$
where the variation of $cal(L)$ with respect to $x^mu (tau)$ is given by
$
  delta cal(L) = [(diff cal(L))/(diff x^mu) - d/dtau (diff cal(L))/(diff dot(x)^mu)] delta x^mu.
$
Requiring the action to remain invariant under all variations of the path leads to the Euler-Lagrange equations of motion
$
    (diff cal(L))/(diff x^mu) - d/dtau (diff cal(L))/(diff dot(x)^mu) = 0.
$

We now turn to field theory. A field is, loosely speaking, a function defined on a spacetime manifold $M$. Whether the spacetime is flat or curved will not concern us here, as our focus is the structure of the formalism rather than its geometric details. Locally, the manifold provides coordinate charts, and fields become coordinate-dependent functions.

In contrast to mechanics, where the parameter $tau$ was one-dimensional, field theory takes its parameters from the manifold itself: the space of independent variables is $N$-dimensional (for an $N$-dimensional manifold). The dynamical objects are the fields we denote as $Phi^A (x)$, where $A$ labels field components (or species). 

This setup leads to several important differences:
+ *Multiple Derivatives*: Since the fields depend on spacetime coordinates $x^mu$, their derivatives $diff_mu Phi^A (x)$ are now gradients rather than time derivatives. Depending on the theory, higher-order derivatives may also appear.

+ *Infinite Degrees of Freedom*: In field theory, each field $Phi^A (x)$ is defined over all of spacetime, which means specifying its value at every point $x in M$. This leads to an uncountably infinite number of degrees of freedom (This was also the case in classical mechanics when formulated for continuous systems, though over a lower-dimensional base space.).


The Lagrangian $cal(L)$ must now be a scalar function of the fields and their derivatives,
$
  cal(L) = cal(L)(Phi^A (x), diff_mu Phi^A(x), ...). 
$
Accordingly, the action becomes an integral over the full manifold,
$
  S[Phi^A] = integral_M d^N x thin cal(L)(Phi^A, diff_mu Phi^A,...).
$
As before, we determine the equations of motion by requiring stationarity of the action under infinitesimal variations $Phi^A -> Phi^A + delta Phi^A$. For simplicity, let us assume dependence of $cal(L)$ on only $Phi^A$ and its gradient $diff_mu Phi^A$. We thus vary 
$
  Phi^A &-> Phi^A + delta Phi^A,\
  diff_mu Phi^A &-> diff_mu Phi^A + diff_mu delta Phi^A.
$
Under such a variation, the action changes as $S -> S + delta S$, where 
$
  delta S = integral_M d^N x [(diff cal(L))/(diff Phi^A) delta Phi^A + (diff cal(L))/(diff(diff_mu Phi^A)) diff_mu delta Phi^A ].
$
Using integration by parts on the second term and assuming that $delta Phi^A (x)$ vanishes on the boundary (or that the manifold is compact without boundary), we obtain
$
  delta S = integral_M d^N x [(diff cal(L))/(diff Phi^A) - diff_mu ((diff cal(L))/(diff(diff_mu Phi^A)))] delta Phi^A.
$
Since $delta Phi^A (x)$ is arbitrary (within suitable regularity conditions), stationarity $delta S = 0$ implies the Euler-Lagrange equations for fields,
$
  (diff cal(L))/(diff Phi^A) - diff_mu ((diff cal(L))/(diff(diff_mu Phi^A))) = 0.
$
This equation is often written in terms of the _functional derivative_ of the action#footnote[More generally, this contains derivatives of the integrand $cal(L)$ with respect to all derivatives of $Phi^A$, not just $Phi^A$ and $diff_mu Phi^A$.],
$
  (delta S)/(delta Phi^A (x)) := [(diff cal(L))/(diff Phi^A) - diff_mu ((diff cal(L))/(diff(diff_mu Phi^A)))](Phi^A (x), diff_mu Phi^A (x)).
$
We will consider the functional derivative in more detail in the next section.
== Functional Derivatives
Let $Phi^A (x)$ be a field configuration over a spacetime manifold $M$, and let $S[Phi^A]$ be a functional---that is, a map from functions to real numbers, typically given as an integral over $M$ of a local function of $Phi^A$ and its derivatives,
$
  S[Phi^A] = integral_M d^N x cal(L)(Phi^A (x), diff_mu Phi^A (x),...).
$
The _functional derivative_ of $S$ with respect to $Phi^A (x)$, denoted by $(delta S) / (delta Phi^A (x))$, is defined by the condition that for any test variation $delta Phi^A (x)$ with compact support, 
$
  delta S = integral_M d^N x (delta S)/(delta Phi^A (x)) delta Phi^A (x).
$<defFunctionalDerivative>
In practice, when $cal(L)$ depends on $Phi^A$ and its first derivatives $diff_mu Phi^A$, the functional derivative takes the familiar Euler-Lagrange form,
$
  (delta S)/(delta Phi^A (x)) = [(diff cal(L))/(diff Phi^A (x)) - diff_mu ((diff cal(L))/(diff_mu Phi^A))](x),
$
evaluated at the point $x$. This is the same expression that appears in the Euler-Lagrange equation for fields. 

Conceptually, the functional derivative measures the sensitivity of the action $S[Phi^A]$ to infinitesimal local changes in the field. It plays the role of a gradient in infinite-dimensional function space, and it vanishes precisely on field configurations that extremise the action.

The @defFunctionalDerivative[definition] tells us that in order to compute the functional derivative of an action, we may compute its variation $delta S$ by the standard rules of variational calculus, collect terms and integrate by parts such that we isolate a $delta Phi^A (x)$ factor, and subsequently read it off.

Alternatively, we can derive some computation rules for the functional derivative. Firstly, note that $Phi^A (x)$ can be seen as an integral functional when expressing it as
$
  Phi^A (x) = integral_M d^N y thin delta (x-y) delta^(A B) Phi^B (y).
$
Using this, the variation becomes
$
  delta Phi^A (x) = integral_M d^N y thin delta^((N-1))(x-y) delta^(A B) delta Phi^B (y),
$
from which we can read off that
$
  (delta Phi^A (x))/(delta Phi^B (y)) = delta^((N-1))(x-y) delta^(A B)
$
This makes intuitive sense: the field value $Phi^A (x)$ depends on $Phi^B (y)$ if and only if $A = B$ and $x = y$---otherwise, field values are independent of each other. 

Using 
$
  diff_mu Phi^A (x) &= integral_M d^N y thin diff_mu delta^((N-1))(x-y) delta^(A B) Phi^B (y),
$
we can further conclude
$
  delta (diff_mu Phi^A (x)) = integral_M d^N y thin diff_mu delta^((N-1))(x-y) delta^(A B) delta Phi^B (y).
$
This allows us to read off
$
  delta/(delta Phi^B (y)) diff_mu Phi^A (x) = diff_mu delta (x-y) delta^(A B) = diff_mu delta/(delta Phi^B (y)) Phi^A (x).
$
We have hence shown that the partial and functional derivatives commute. 

Now that we know how the functional derivative acts on the fields and their derivatives, we can also consider how it acts on the Lagrangian itself. For this, we need to further impose a product and chain rule on it, which leads us to the following:
$
  delta/(delta Phi^B (y)) cal(L)(Phi^A (x), diff_mu Phi^A (x)) &= (diff cal(L))/(diff Phi^A)  (delta Phi^A (x))/(delta Phi^B (y)) + (diff cal(L))/(diff (diff_mu Phi^A)) (delta (diff_mu Phi^A) (x))/(delta Phi^B (y))\
  &= (diff cal(L))/(diff Phi^A) delta^(A B) delta^((N-1))(x-y) + (diff cal(L))/(diff (diff_mu Phi^A)) delta^(A B) diff_mu delta^((N-1))(x-y)\
  &= [(diff cal(L))/(diff Phi^A) - diff_mu (diff cal(L))/(diff(diff_mu Phi^A))](x) thin delta^(A B) delta^((N-1))(x-y).
$
This is consistent with 
$
  (delta S)/(delta Phi^B (y)) &= integral_M d^N x delta/(delta Phi^B (y)) cal(L)(Phi^A (x), diff_mu Phi^A (x))\ 
  &= integral_M d^N x [(diff cal(L))/(diff Phi^A) - diff_mu (diff cal(L))/(diff(diff_mu Phi^A))](x) thin delta^(A B) delta^((N-1))(x-y)\
  &= [(diff cal(L))/(diff Phi^B) - diff_mu (diff cal(L))/(diff(diff_mu Phi^B))] (y).
$
So, we now know how the functional derivatives acts on functionals as well as (local) functions of the fields, which can be viewed as functionals as well (in appropriate distribution spaces).

All of this has been very informal, and I'd need to take a proper course on calculus of variations to provide a better insight. But this was more of an interlude so it doesn't really matter.

== Legendre Transformation and the Hamiltonian (Density)
The Legendre transformation of a field theory Lagrangian proceeds very similarly to how it works in mechanics. We start from a Lagrangian density,
$
  cal(L) = cal(L) (Phi^A (x), diff_mu Phi^A (x)),
$
and define conjugate momenta as
$
  Pi_A (x) = (diff cal(L))/(diff (diff_0 Phi^A))(Phi^A (x), diff_mu Phi^A (x)).
$<canonicalMomentaFields>
Here, $diff_0 = diff/(diff x^0)$ denotes the derivative with respect to one chosen coordinate $x^0$ (not necessarily timelike), such that the full set of coordinates is given by $x^mu = (x^0, x^i)$. In mechanics, we never had to choose a parameter/direction for the Legendre transform, since the parameter space is one-dimensional, just a single parameter $tau$, so there's no freedom in choosing a Legendre direction. As we will see later on, the Hamiltonian that comes out of the Legendre transform allows one to evolve arbitrary observables along $x^0$. Because of this, typically, one chooses $x^0$ as a timelike coordinate (e.g. the parameter $t$ from a chosen foliation of a spacetime $Sigma_t$, or a time coordinate in cartesian Minkowski-space coordinates), which then yields time-evolution equations. One could, however---just as well---evolve along a spacelike or even null direction. Time-evolution is particularly useful in practice, as it directly encodes causality and allows for the imposition of initial conditions. Before moving on, let us briefly remark that this choice of coordinate to take the Legendre transformation with respect to explicitly breaks Lorentz- or diffeomorphism invariance (depending on context). This is a tradeoff one has to make to get access to evolution equations.

The equation @canonicalMomentaFields give a relationship between the field momenta $Pi_A (x)$ and the velocities $diff_0 Phi^A (x)$. To make the Legendre transform invertible, this relationship must be invertible too, requiring that
$
  det ((diff^2 cal(L))/(diff (diff_0 Phi^A) diff (diff_0 Phi^B))) != 0 quad forall x in M.
$
This is equivalent to the variable transformation @canonicalMomentaFields[] having nonsingular Jacobian everywhere, which in turn is equivalent to $Pi_A (x)$ and $diff_0 Phi^A (x)$ having a bijective relationship. The Legendre transform of the Lagrangian (read: the Hamiltonian (density)) takes the form
$
  cal(H) (Phi^A (x), diff_i Phi^A (x), Pi_A (x)) = Pi_A (x) diff_0 Phi^A (x) - cal(L)(Phi^A (x), diff_mu Phi^A (x)),
$
or more abridgedly, 
$
  cal(H)(Phi^A, diff_i Phi^A, Pi_A) = Pi_A diff_0 Phi^A - cal(L) (Phi^A, diff_mu Phi^A).
$<hamiltonianDensityDef>
Notably, the derivatives $diff_i Phi^A$ still appear as arguments and are not replaced by momenta. This reflects the fact that this canonical transformation is only partial, i.e. only in the $x^0$-direction, and not in all spacetime directions. If one decides to Legendre-transform with respect to _all_ directions, then one arrives at the _de Donder-Weyl_ covariant Legendre transform.

== Equations of Motion from Lagrangian Field Theory
Now that we have a Hamiltonian density, let us derive the equations of motion for it. This proceeds in a manner similar to the calculation in mechanics. We begin with the Euler-Lagrange equations and substitute the expression for the momentum as well as the Lagrangian. More concretely,
$
  0 &= (diff cal(L))/(diff Phi^A) - diff_mu (diff cal(L))/(diff (diff_mu Phi^A)) = (diff cal(L))/(diff Phi^A) - diff_0 Pi_A (x) - diff_i (diff cal(L))/(diff (diff_i Phi^A (x)))\
  &= -(diff cal(H))/(diff Phi^A) - diff_0 Pi_A (x) + diff_i (diff cal(H))/(diff(diff_i Phi^A))
$
which implies that
$
  diff_0 Pi_A (x) = -[(diff cal(H))/(diff Phi^A) - diff_i (diff cal(H))/(diff (diff_i Phi^A))].
$
To find the equations for $Phi^A$, we simply take the derivative of @hamiltonianDensityDef with respect to $Pi_A$, which yields
$
  (diff cal(H))/(diff Pi_A) &= diff_0 Phi^A + underbrace(Pi_B (diff( diff_0 Phi^B))/(diff Pi_A) - (diff cal(L))/diff(diff_0 Phi^B) (diff (diff_0 Phi^B))/(diff Pi_A),=0)\
  &= diff_0 Phi^A.
$
In summary, we have found the canonical equations for the fields and their momenta,
$
  diff_0 Phi^A = (diff cal(H))/(diff Pi_A), wide diff_0 Pi_A = -[(diff cal(H))/(diff Phi^A) - diff_i (diff cal(H))/(diff (diff_i Phi^A))].
$
This is a good time to introduce the Hamiltonian (not density), which is defined as
$
  H(x^0) = integral_(Sigma_(x^0)) d^(N-1) x thin cal(H)(Phi^A (x), diff_i Phi^A (x), Pi_A (x)).
$
Here, $Sigma_(x^0)$ denotes the family of spacelike hypersurfaces with $diff_(x^0)$ as a normal vector (a foliation, probably somewhat wrong the way I put it but we will get to that in Cambridge I suppose). In terms of this, we can write the equations above as
$
  diff_0 Phi^A (x) = (delta H)/(delta Pi_A (x)), quad diff_0 Pi_A = - (delta H)/(delta Phi^A),
$<canonicalEqnsFieldTheory>
which in this form are much more analogous to those we found in mechanics (cf. @canonicalEqns)
== Hamiltonian Mechanics as a 1D Field Theory
Before moving on to Poisson brackets in field theories, we should review the connection between Hamiltonian mechanics and Hamiltonian field theory. This is because there is an important observation to make: Hamiltonian mechanics is a special case of Hamiltonian field theory, where the parameter space is one-dimensional---with $tau$ as the only parameter.

This has a couple of consequences, one of which we already touched upon briefly when computing the Legendre transform for the Lagrangian of a field theory. If there is only one parameter, there is only one choice of what direction to perform the transform with respect to. Another important observation is that of initial conditions. The phase space is essentially the space of solutions. It is parameterised by initial conditions (under sufficient assumptions that they define unique solutions to the equations of motion), meaning that the dimension of the phase space is the dimension of the space of allowed initial conditions. Initial conditions are specified on a surface of codimension 1, i.e. of one dimension less than the parameter space. In the case of mechanics, with one parameter $tau$, this means that initial conditions are prescribed on dimension-0 surfaces, also known as _points_. This makes the space of initial conditions and thus phase space _finite-dimensional_. As soon as we have a proper field theory, with $>= 2$ parameters or coordinates, the dimension of the initial condition surface is $>= 1$. The space of, say, $L^2$ functions on such a surface is infinite-dimensional, making the phase space infinite-dimensional as well.

With this perspective in mind, we are ready to extend familiar concepts like Poisson brackets to the field-theoretic setting, where they retain their role as generators of time evolution—albeit now in infinite-dimensional phase space.

== Poisson Brackets
Recall that in mechanics, for two observables $f(q,p)$ and $g(q,p)$, the Poisson bracket is defined as
$
  {f,g} = (diff f)/(diff q^mu) (diff g)/(diff p_mu) - (diff g)/(diff q^mu) (diff f)/(diff p_mu).
$
To generalise this to field theory, we should first identify which structural properties are essential and must be preserved:
- *Antisymmetry and bilinearity*: The Poisson bracket is an antisymmetric bilinear operation on the space of observables. This algebraic structure should carry over to the field-theoretic setting.
- *Dependence on canonical variables*: The bracket encodes how observables vary with respect to the canonical variables $q^mu$ and $p_mu$. In field theory, the canonical variables are the fields $Phi^A (x)$ and their conjugate momenta $Pi_A (x)$, so the bracket must similarly reflect how observables depend on these.
- *Sensitivity to local structure*: Observables in field theory are typically functionals that may depend on the fields, their momenta, and possibly spatial derivatives such as $diff_i Phi^A$. To capture this dependence, the role of partial derivatives in mechanics should be played by functional derivatives in the generalisation.

Putting all of this together, and sprinkling in an integral for good measure (haha get it, measure), we get to the field-theory definition of the Poisson bracket for two observable functionals $F[Phi^A, diff_i Phi^A, Pi_A]$ and $G[Phi^A, diff_i Phi^A, Pi_A]$ given by
$
  F(x^0) &= integral_(Sigma_(x^0)) d^(N-1) x thin cal(F) (Phi^A (x), diff_i Phi^A (x), Pi_A (x)),\
  G(x^0) &= integral_(Sigma_(x^0)) d^(N-1) x thin cal(G) (Phi^A (x), diff_i Phi^A (x), Pi_A (x)),
$
then the Poisson bracket reads
$
  {F,G} = integral_(Sigma_(x^0)) d^(N-1) x [(delta F)/(delta Phi^A (x)) (delta G)/(delta Pi_A (x)) - (delta G)/(delta Phi^A (x)) (delta F)/(delta Pi_A (x))].
$
The Poisson bracket ${F,G}$ is again a functional of the canonical variables with a residual dependence on $x^0$. The functional derivative here is to be interpreted as
$
  (delta F)/(delta Phi^A (x)) = [(diff cal(F))/(diff Phi^A) - diff_i (diff cal(F))/(diff(diff_i Phi^A))](x)
$<eq264>
etc.
Let us now compute the fundamental Poisson brackets of the canonical variables, with the standard reinterpretation in the functional sense:
$
  Phi^A (x^0, x^i) = integral_(Sigma_(x^0)) d^(N-1)y thin Phi^A (x^0, y^i) delta^((N-1)) (x^i - y^i),\
  Pi_A (x^0, x^i) = integral_(Sigma_(x^0)) d^(N-1)y thin Pi_A (x^0, y^i) delta^((N-1)) (x^i - y^i).
$<eq265>
With these, we find
$
  {Phi^A, Phi^B} = {Pi_A, Pi_B} = 0
$
which follows from the $Pi^A$-independence of $Phi^B$ and vice versa. The interesting thing happens with the remaining Poisson bracket, where we can make use of
$
  (delta Phi^A (x^0,x^i))/(delta Phi^B (x^0,y^i)) = delta^A_B delta^((N-1))(x^i - y^i),quad (delta Pi_A (x^0,x^i))/(delta Pi_B (x^0,y^i)) = delta_A^B delta^((N-1))(x^i - y^i), 
$
and zero for $delta Phi\/delta Pi$ and $delta Pi \/ delta Phi$ (this follows from @eq264[eq.] and @eq265[eq.]). This allows us to compute
$
  {Phi^A (x^0, x^i), Pi_B (x^0, y^i)} &= integral_(Sigma_(x^0)) d^(N-1)z thin [(delta Phi^A (x^0,x^i))/(delta Phi^C (x^0,z^i)) (delta Pi_B (x^0,y^i))/(delta Pi_C (x^0,z^i)) - (delta Pi_B)/(delta Phi^C) (delta Phi^A)/(delta Pi_C)]\
  &= integral_(Sigma_(x^0)) d^(N-1)z thin delta^A_C delta^((N-1))(x^i-z^i) delta_B^C delta^((N-1))(y^i-z^i)\
  &= delta^A_B delta^((N-1))(x^i - y^i).
$
This is analogous to the Poisson bracket in mechanics. 

More generally, we can recover classical mechanics by considering the special case $N=1$, where $x^0 = tau$ and there are no spatial coordinates $x^i$. In this case, $N-1 = 0$, so the integrals over the spatial slice and delta functions disappear, and the Poisson bracket becomes purely algebraic. The functional derivative with respect to $Phi^A$ and $Pi_A$ reduce to ordinary partial derivatives, since the terms involving spatial derivatives (such as the second term in @eq264[equation]) vanish. Identifying $q^A (tau) = Phi^A (x^0)$ and $p_A (tau) = Pi_A (tau)$, we recover the usual Poisson bracket structure of classical mechanics.

So, we would also suspect the Poisson bracket of an observable with the Hamiltonian to correspond to a time evolution. Let us investigate this now. To this end, we simply compute
$
  {F,H} &= integral_(Sigma_(x^0)) d^(N-1)x lr([(delta F)/(delta Phi^A (x)) underbrace((delta H)/(delta Pi_A (x)),=diff_0 Phi^A (x)) - underbrace((delta H)/(delta Phi^A (x)),=-diff_0 Pi_A (x)) (delta F)/(delta Pi_A (x))],size:#55%) \
  &= integral_(Sigma_(x^0)) d^(N-1)x [(delta F)/(delta Phi^A (x)) diff_0 Phi^A (x) + (delta F)/(delta Pi_A (x)) diff_0 Pi_A (x)]\
  &= dF/dx^0.
$ 
The last step might seem a little bit abrupt, so let us go through it more carefully. Evaluating the total derivative with respect to $x^0$, making use of the symmetry of partial derivatives, we get
#[
#set math.equation(number-align:bottom)
$
  dF/dx^0 &= integral_(Sigma_(x^0)) d^(N-1)x d/dx^0 cal(F)(Phi^A (x), diff_i Phi^A (x), Pi_A (x), diff_i Pi_A (x))\
  &= integral_(Sigma_(x^0)) d^(N-1)x [(diff cal(F))/(diff Phi^A) diff_0 Phi^A + (diff cal(F))/(diff (diff_i Phi^A))diff_0 diff_i Phi^A + (diff cal(F))/(diff Phi^A) diff_0 Phi^A + (diff cal(F))/(diff (diff_i Pi_A))diff_0 diff_i Pi_A]\
  &= integral_(Sigma_(x^0)) d^(N-1)x [[(diff cal(F))/(diff Phi^A) - diff_i (diff cal(F))/(diff(diff_i Phi^A))] diff_0 Phi^A + [(diff cal(F))/(diff Pi_A) - diff_i (diff cal(F))/(diff(diff_i Pi_A))]diff_0 Pi_A]\
  &= integral_(Sigma_(x^0)) d^(N-1)x [(delta F)/(delta Phi^A) diff_0 Phi^A + (delta F)/(delta Pi_A) diff_0 Pi_A].
$
]
This completes the analogy to mechanics.

= Hamiltonian Field Theory: Examples and Constraints
In this section, we apply the formalism developed previously to specific examples of field theories. We begin with a real scalar field theory and then move on to Maxwell theory. For each case, we start from the Lagrangian, derive the canonical variables, and compute the corresponding equations of motion.

The scalar field example proceeds smoothly and serves to illustrate the mechanics of the formalism in a straightforward setting. In contrast, Maxwell theory reveals a subtlety: the Legendre transformation becomes singular due to the presence of gauge freedom. This breakdown motivates the introduction of constraints, which will extend the canonical formalism to accommodate gauge-invariant systems.
== Naive Example: Real Scalar Field Theory
We begin with one of the simplest field theories: a massive real scalar field in flat Minkowski spacetime. The Lagrangian density is given by
$
  cal(L) = -1/2 eta^(mu nu) diff_mu phi.alt diff_nu phi.alt - 1/2 m^2 phi.alt^2,
$ 
where we adopt the metric signature used by adults---mostly plus, of course ($eta_(mu nu) = diag{-1,+1,...,+1}$). The associated action reads
$
  S = integral d^N x thin cal(L),
$ 
and its variation under $phi.alt->phi.alt+delta phi.alt$ is
$
  delta S = integral d^N x lr([underbrace(diff_mu diff^mu, = Box) phi.alt - m^2 phi.alt], size:#35%) delta phi.alt 
$
From this, we extract the Euler-Lagrange equation of motion,
$
  (Box - m^2) phi.alt = 0.
$<KleinGordonEqn>
While this covariant equation is elegant, it is second-order in time and does not make the dynamics or initial value structure immediately transparent. Moreover, it is ill-suited for direct numerical implementation. To address this, we now transition to the Hamiltonian formalism, which, though non-covariant, naturally casts the dynamics in a first-order form better adapted to time evolution.

To derive the Hamiltonian, we need to perform a Legendre transform. To this end, we derive the canonical momenta, 
$
  pi (x) = (diff cal(L))/(diff (diff_0 phi.alt))(x) = -diff^0 phi.alt(x) = diff_0 phi.alt(x).
$
It is pretty clear already that this is an invertible relationship, but to be sure we can compute the Hessian (in this case, a $1 times 1$ matrix) as well. This yields
$
  (diff^2 cal(L))/(diff (diff_0 phi.alt)^2) = 1
$
which is clearly non-singular. Thus, the Legendre-transform of $cal(L)$ (aka the Hamiltonian density) is well-defined and given by
$
  cal(H) = pi diff_0 phi.alt - cal(L) &= pi^2 + 1/2 eta^(mu nu) diff_mu phi.alt diff_nu  phi.alt + 1/2 m^2 phi.alt^2\
  &= 1/2 pi^2 + 1/2 delta^(i j) diff_i phi.alt diff_j phi.alt + 1/2 m^2 phi.alt^2.
$
The Hamiltonian (non-density) is simply the integral over a spatial slice,
$
  H(t) = integral_(RR^3) d^3x thin cal(H)(x).
$
From it, we can derive the equations of motion for the canonical variables $phi.alt$ and $pi$. These read
$
  diff_0 phi.alt &= (delta H)/(delta pi) = (diff cal(H))/(diff pi) - diff_i (diff cal(H))/(diff (diff_i pi)) = pi,\
  diff_0 pi &= - (delta H)/(delta phi.alt) = delta^(i j) diff_i diff_j phi.alt -m^2 phi.alt.
$
More compactly, we can write
$
  diff_0 phi.alt = pi, quad diff_0 pi =  Delta phi.alt - m^2 phi.alt
$<canonicalEoMFreeMassiveScalar>
Clearly, these equations can be combined to find
$
  diff_0^2 phi.alt = diff_0 pi = Delta phi.alt - m^2 phi.alt,
$
which can be rearranged to recover the Euler-Lagrange equation
$
  (Box - m^2) phi.alt = 0.
$
This system is unconstrained and features a simple phase space: one real degree of freedom per point in space $x^i$, evolving in the direction $x^0$, described by the field $phi.alt (x^0,x^i)$ and its conjugate momentum $pi(x)$. The Hamiltonian equations of motion form a well-posed Cauchy problem, making the theory suitable for both analytical and numerical treatment---initial data specified on a spatial slice determines the evolution uniquely. Moreover, this canonical structure forms the basis for canonical quantisation, where one promotes $phi.alt$ and $pi$ to operators and imposes equal-time commutation relations,
$
  [hat(phi.alt)(x^0,x^i), hat(pi)(x^0, y^i)] = i delta^((N-1)) (x^i - y^i),
$
leading to a quantised scalar field theory in the Schrödinger or Fock picture. In contrast, our next example---Maxwell theory---features gauge freedom, which renders the Legendre transform singular and demands a more careful analysis involving constraints.
== First Attempt at Maxwell Theory
We now attempt to repeat the procedure from the previous section to Maxwell theory in Minkowski space. The fields involved in Maxwell theory are the components of a 1-form,
$
  A = A_mu dx^mu,
$
which define a field strength via the exterior derivative,
$
  F = dA = 1/2 F_(mu nu) dx^mu wedge dx^nu,
$
where 
$
  F_(mu nu) = diff_mu A_nu - diff_nu A_mu.
$
Though it is by no means necessary to use the language of differential forms here, it gives us the advantage of making gauge invariance obliviously obvious. A gauge transformation is a transformation of $A$ as
$
  A -> A + dPsi,
$
where $Psi$ is a 0-form (aka a scalar). On the level of the components, this reads
$
  A_mu -> A_mu + diff_mu Psi.
$
The reason this makes gauge invariance of $F$ obvious is because of the exactness of the exterior derivative, $d^2 = 0$. Making use of that, we immediately deduce the field strength's gauge transformation behaviour,
$
  F = dA -> d(A + dPsi) = dA = F.
$
Hence, $F$, as well as any quantity constructed from it, is gauge invariant. In particular, we can define an action as
$
  S[A] = integral F wedge star F = integral d^N x thin [-1/4 F_(mu nu) F^(mu nu)]
$
where $star F$ is the Hodge dual of $F$. We won't get any further into differential form language here, so we won't concern ourselves with the Hodge dual anymore. This was just a neat way to briefly introduce some very basic differential form theory from what I learned the past half year or so. Either way, from now on, we just continue with the Lagrangian density
$
  cal(L) = - 1/4 F_(mu nu) F^(mu nu).
$
To derive the Euler-Lagrange equations of motion, it is most convenient to derive the variation of the action under $A_mu -> A_mu + delta A_mu$. This yields
$
  delta S &= -integral d^N x thin 1/2 F^(mu nu) delta F_(mu nu) = -integral d^N x thin F^(mu nu) diff_(\[mu) delta A_(nu\])\ &= -integral d^N x thin F^(mu nu) diff_mu delta A_nu = integral d^N x [diff_mu F^(mu nu)] delta A_nu. 
$<ELeqnsDerivMaxwell>
The requirement of stationarity, $delta S = 0$, thus gives rise to the Euler-Lagrange equations of motion,
$
  diff_mu F^(mu nu) = 0.
$
In the following, it will be useful to decompose the action in terms of the time-space and space-space components of the field strength (there are no time-time components due to anti-symmetry). More concretely, we expand
$
  cal(L) = -1/4 F_(mu nu) F^(mu nu) = -1/2 F_(0 i) F^(0 i) - 1/4 F_(i k) F^(i k).
$
Notably, time-derivatives of the field $A_mu$ only appear in the first term. This allows us to compute the conjugate momenta. For the spatial components, we find
$
  pi^i = (diff cal(L))/(diff (diff_0 A_i)) = - F^(0 i) = F_(0i) = diff_0 A_i - diff_i A_0 = E^i,
$
where $E^i$ denotes the $i$-th component of the electric field. For the time component, $pi^0$, we encounter something unexpected (or expected, if you've given it some thought before). Namely, following the definition, we find that
$
  pi^0 = (diff cal(L))/(diff (diff_0 A_0)) = 0
$
identically since $F_00 = 0$ is the only component that could contain any dependence on $diff_0 A_0$, but is zero due to antisymmetry. This means the $pi^0 $ momentum is constrained and clearly leads to a non-invertible relationship between the velocities $diff_0 A_mu$ and the canonical momenta $pi^mu$. This is also clear from
$
  (diff^2 cal(L))/(diff(diff_0 A_mu) diff(diff_0 A_nu)) = diag{0,1,...,1} = cases(0\,quad&mu\,nu=0\,,1\, &mu = nu != 0.)
$ 
This is a rank $N-1$ matrix which has zero determinant. This reflects the fact that the velocity $diff_0 A_0$ does not appear in the Lagrangian, and hence its conjugate momentum vanishes identically. In conclusion, we are unable to perform a non-singular Legendre transform.

Let us now reinterpret this algebraic obstruction from a more physical point of view. From the expression $pi^i = F_(0 i) = diff_0 A_i - diff_i A_0$, 
we see that the relationship between the momenta $pi^i$ and the velocities $diff_0 A_i$ is invertible, provided $A_0$ is treated as a given background field. In contrast, the time component yields $pi^0 = 0$, since no term in the Lagrangian depends on $diff_0 A_0$. This means that $A_0$ does not appear with any time derivatives in the action and hence does not appear with any time derivatives in the equations of motion either. In other words, it has no dynamics of its own.

Nonetheless, $A_0$ can still influence the evolution fo the spatial components $A_i$ through its spatial derivatives $diff_i A_0$, which enter the equations of motion via the expression for $pi^i$. Thus, while $A_0$ cannot evolve independently, it acts as a source for the dynamics of the spatial sector.

This non-dynamical character allows us to remove $A_0$ entirely by a gauge choice. For example, if we take 
$
  Psi = integral dx^0 A_0,
$
then under the gauge transformation $A-> A-dPsi$, we obtain
$
  A_0 -> A_0 - diff_0 Psi = 0.
$
This is often referred to as the _temporal gauge_. Of course, this is not the only way to fix the gauge redundancy; other choices such as the Coulomb gauge, $diff_i A^i = 0$, or Lorenz gauge, $diff_mu A^mu = 0$, are equally valid. In the end, they simply relate the non-dynamical component $A_0$ to the dynamical $A_i$. 

From the variation of the action in @ELeqnsDerivMaxwell[], we can now distinguish two qualitatively different types of equation. A variation with $delta A_0$= 0, $delta A_i != 0$, yields 
$
  0 = diff_mu F^(mu i) = diff_0 F^(0 i) + diff_k F^(k i) = - diff_0 pi^i + diff_k F^(k i),
$
which is a genuine evolution equation, since it involves a time derivative of the momentum. In contrast, a variation with $delta A_0 != 0$, $delta A_i = 0$, gives
$
  0 = diff_mu F^(mu 0) = diff_k F^(k 0) = diff_k pi^k.
$
This equation contains no time derivatives at all---it is a constraint rather than a dynamical law. Recognising $pi^k = E^k$, we identify this constraint as Gauss' law.

In summary, we attempted a Legendre transformation of the Maxwell Lagrangian to obtain a Hamiltonian formulation. This revealed that the momentum $pi^0$ is constrained to vanish identically, indicating that $A_0$ is a non-dynamical field. This is reflected both in the structure of the Lagrangian and in the equations of motion. In addition to the primary constraint $pi^0 = 0$, we also uncovered a secondary constraint $diff_i pi^i =0$, corresponding to Gauss' law. Both constraints arise naturally from the underlying gauge redundancy, and their proper treatment requires a more systematic theory of constrained Hamiltonian systems, to which we now turn.


== Introduction to Constraints in Mechanics
In this subsection, we introduce the concept of constraints in classical systems. We begin with mechanical systems, where they are easier to motivate and analyse, and then generalise the formalism to field theories with an $N$-dimensional parameter space. This general framework will then allow us to revisit the situation we encountered earlier with Maxwell theory from a broader perspective and develop the Dirac-Bergmann algorithm.

=== Constraints in Lagrangian Mechanics
Consider a mechanical system with a configuration described by coordinates $q^mu (tau)$, governed by a Lagrangian of the form
$
  cal(L) = cal(L)(q,dot(q)).
$
Now suppose that the system is subject to an additional condition---such as being confined to a lower-dimensional submanifold of the ambient configuration space, or having fixed components of its velocity. Such conditions are formalised as constraints of the form
$
 F(q,dot(q)) = 0,
$
where $F$ is a function of the coordinates and velocities. 

To make this idea concrete, consider the case of a free particle in $RR^2$, with coordinates $q^mu = (x,y)$ and a free Lagrangian
$
  cal(L)(x,y,dot(x), dot(y)) = 1/2 dot(x)^2 + 1/2 dot(y)^2.
$
A possible constraint could be that the motion is restricted to the unit circle, i.e.,
$
 F(x,y, dot(x), dot(y)) = x^2 + y^2 - 1 = 0.
$
This condition restricts the admissible trajectories in configuration space. A function like $F$, which defines such a restriction, is called a _constraint function_, and the condition itself is called a _constraint_. 

While defining a constraint is straightforward---one simply specifies a function $F(q,dot(q))$ and demands it to vanish---the more subtle question is how such a condition affects the equations of motion. This leads to the formalism of _Lagrange multipliers_, which we now develop.

Let us first try to get an intuitive picture. You may have noticed something already: above, we explicitly took a _free_ Lagrangian---only to then tie our particle to a leash and drag it around the unit circle. The only truly "free" solution, in the sense of satisfying the unconstrained Euler-Lagrange equations
$
  ddot(x) = ddot(y) = 0,
$
would be a stationary particle. But motion along a circle is never inertial---it is always accelerating. So the free dynamics and the constrained dynamics are clearly not the same.

The takeaway is this: _the equations of motion get modified_. Something must supply the force that keeps the particle on the constraint surface $F=0$. These additional terms act exactly so as to cancel the parts of the dynamics that would otherwise pull the system away from $F=0$, while leaving all motion within the constraint surface untouched.

On a superficial level, this is straightforward. If we want new terms in the equations of motion, we add new terms to the Lagrangian. That much is clear. But what is not yet clear is which terms to add. We want just enough to enforce the constraint, and nothing more---no overcorrections, no spurious effects. More formally, we look for paths $q(tau)$ that minimise the action functional
$
  S[q] = integral dtau thin cal(L)(q,dot(q)),
$
while simultaneously being subject to the constraint. Informally, we can write this as
$
  q(tau) = op("arg min", limits:#true)_(f(q,dot(q))=0) S[q].
$
We aren't minimising the functional $S$ on the whole configuration space anymore, only on a restricted subset.

At this point, our intuition has taken us as far as it can. Time to be precise---and introduce the method of _Lagrange multipliers_. For this, consider the following preposterous idea: what if we could turn the constraint equation
$
 F(q,dot(q)) = 0
$<mechanicsConstraintEqn>
into an equation of motion itself? Then it would automatically be satisfied whenever the equations of motion hold. That sounds promising---but clearly cannot work _as is_. So let us reflect on what it actually means to be an equation of motion, or just a plain ol' equation. The answer is kind of obvious from the previous section, but still worth reiterating: equations of motion arise from varying the degrees of freedom in the action and requiring
$
  delta S = 0.
$
This gives a collection of conditions which we call the equations of motion. Now that we reminded ourselves of what that means, we can continue on our "fantasy journey" of turning the constraint @mechanicsConstraintEqn[] into one. 

This means we must somehow vary a degree of freedom in such a way that the constraint @mechanicsConstraintEqn[] pops out. Varying with respect to $q$, that _might_ be possible in some cases---but there is a much simpler trick. Just introduce a new degree of freedom---let us call it $lambda(tau)$. Then consider an action of the form
$
  S_C [q;lambda] = integral dtau thin lambda(tau)F(q (tau),dot(q)(tau)).
$
Varying this degree of freedom as $lambda -> lambda + delta lambda$ gives
$
  delta S_C = integral dtau thin F(q,dot(q)) delta lambda,
$
and requiring this to vanish for arbitrary $delta lambda(tau)$ gives back the constraint $F=0$. Voilà---We have turned it into an equation of motion for this new action $S_C$.

Naturally, $S_C$ carries none of the physics of the original system. But since $S$ does not depend on $lambda$, varying the total action 
$
  S_T = S - S_C = integral dtau [cal(L)(q,dot(q)) - lambda F(q,dot(q))] equiv integral dtau thin cal(L)_T (q, dot(q); lambda).
$
with respect to $lambda$ gives
$
  delta S_T =0 quad iff quad delta S_C = 0 quad iff quad F(q,dot(q)) = 0.
$
So the constraint is now _built into_ the variational principle of $S_T$---just as we wanted. The minus sign here is pure convention; it will become more natural later. As a side note: the symbol $cal(L)_tau (q,dot(q);lambda)$ we introduced in passing above will from now on be referred to as the _total Lagrangian_.

Clearly, modifying the action in this way affects the equations of motion for $q$ as well. After all, we are introducing a new term that depends on $q$ and $dot(q)$. A straightforward calculation shows that under $q-> q+ delta q$, the variation of the total action becomes
$
  delta S_T = integral dtau [(delta cal(L))/(delta q^mu) - lambda (delta F)/(delta q^mu)] delta q^mu.
$
Writing out the full equations of motion, we find
$
  [(diff cal(L))/(diff q^mu) - d/dtau (diff cal(L))/(diff dot(q)^mu)](tau) = lambda(tau) [(diff F)/(diff q^mu) - d/dtau (diff F)/(diff dot(q)^mu)](tau).
$
This is just the regular Euler-Lagrange expression on the left-hand side, with a correction term on the right that enforces the constraint. That correction represents the _constraint force_ needed to keep the system on the surface defined by $F=0$. Also, the conventional minus sign becomes of use here: it allows us to take the constraint force to the right-hand side without introducing additional signs.

The new degree of freedom $lambda(tau)$ we have introduced is not a physical coordinate like $q(tau)$, but rather a _Lagrange multiplier_. It enforces a constraint by adjusting itself to ensure $F(q,dot(q))=0$ throughout the motion. We should note here that it is not dynamical---the full action does not incorporate any of its derivatives, only $lambda$ itself. 

Whether one must explicitly solve for $lambda(tau)$ depends on the problem. Sometimes, $lambda$ can be eliminated by combining the equations of motion with the constraints and its derivatives, effectively reducing the system to unconstrained variables. In other cases, $lambda$ carries meaningful information and must be kept as an auxiliary field to correctly describe the dynamics.

Either way, the introduction of $lambda$ provides a powerful and systematic way to incorporate constraints directly into the variational principle, allowing us to handle complex constrained systems with the familiar tools of Lagrangian and Hamiltonian mechanics.
=== Constraints in Hamiltonian Mechanics
The first sections of these notes aim to develop Hamiltonian field theory from the ground up. So far, we have laid out the foundations and encountered our first major complication: a singular relationship between velocities and momenta in Maxwell theory. This led us into this detour on constrained dynamics. However, up to this point, our discussion of constraints has remained entirely within the framework of _Lagrangian_ mechanics.

In this section, we begin the translation of that framework into the Hamiltonian picture. When the relationship between velocities $dot(q)^mu$ and canonical momenta
$
  p_mu = (diff cal(L))/(diff dot(q)^mu)
$
is non-singular, the Legendre transform can be performed directly, and the treatment of constraints carries over in a rather straightforward way. That will be the focus here.

The singular case is more delicate. There, the Legendre transform fails to be invertible, and constraints arise intrinsically from the degeneracy of the Lagrangian's Hessian with respect to the velocities. We will return to that scenario later, in the context of field theories, when we introduce the Dirac-Bergmann algorithm.

Recall from the previous section that starting from a Lagrangian $cal(L)(q,dot(q))$ and a constraint
$
  F(q,dot(q)) = 0,
$
we may define a constrained _total Lagrangian_ by
$
  cal(L)_T (q,dot(q); lambda) = cal(L)(q,dot(q)) - lambda F(q,dot(q)),
$
where $lambda = lambda(tau)$ is a newly introduced degree of freedom enforcing the constraint. This total Lagrangian describes the dynamics of the constrained system, and may itself admit a Legendre transform, yielding the constrained Hamiltonian $cal(H)_T (p,q;lambda)$. This construction, as always, presupposes that the relation
$
  p_mu = (diff cal(L)_T)/(diff dot(q)^mu) = (diff cal(L))/(diff dot(q)^mu) - lambda (diff F)/(diff dot(q)^mu)
$<momentaConstraints>
between the momenta and the velocities is invertible. Since $lambda$ appears without derivatives, it does not contribute a momentum and is non-dynamical; we choose not to Legendre-transform with respect to its velocity.

Applying the standard Legendre transform yields
$
  cal(H)_T &= p_mu dot(q)^mu - cal(L)_T = [(diff cal(L))/(diff dot(q)^mu) dot(q)^mu -lambda (diff F)/(diff dot(q)^mu) dot(q)^mu] - [cal(L) - lambda F]\
  &= cal(H) - lambda cal(F),
$
where $cal(H)$ is the Legendre transform of the unconstrained Lagrangian, i.e. the unconstrained Hamiltonian
$
  cal(H) = (diff cal(L))/(diff dot(q)^mu) dot(q^mu) - cal(L),
$<LTforUnconstrainedLagrangian>
and $cal(F)$ is the Legendre transform of the constraint function $F$,
$
  cal(F) = (diff F)/(diff dot(q)^mu) dot(q)^mu - F
$<LTforConstraint>
As a side note (hopefully clarifying rather than confusing), we now have three momentum-like quantities in play:
- the _total_ canonical momentum from $cal(L)_T$,
  $ 
    p_mu = (diff cal(L)_T)/(diff dot(q)^mu),
  $
- the _unconstrained_ momentum, arising from $cal(L)$ alone,
  $
    p_mu^"UC" = (diff cal(L))/(diff dot(q)^mu),
  $
- and the constraint momentum, defined by
  $
    p_mu^C = (diff F)/(diff dot(q)^mu).
  $
From @momentaConstraints[equation], we find that they are related via
$
  p_mu = p_mu^"UC" - lambda p_mu^C.
$
This matters because the Legendre transforms $cal(L) -> cal(H)$ and $F-> cal(F)$ are constructed with respect to $p^"UC"$ and $p^C$, respectively. By the above relationship between the momenta, this causes the full Legendre transform $cal(L)_T -> cal(H)_T$ to be with respect to the total momentum $p_mu$, as required.
=== Example: Constrained Dynamics in $RR^2$ <constrainedMechanicsExample>
In this section, we apply our theory of Lagrange multipliers to the free action
$
  S[q] = integral dtau [1/2 dot(x)^2 + 1/2 dot(y)^2]
$
with $q^mu = (x,y) in RR^2$, subject to the constraint
$
  0 =F(x,y, dot(x), dot(y)) = x^2 + y^2 - 1.
$
We will first write down the total Lagrangian and its associated equations of motion, then compute the total Hamiltonian and derive its canonical equations. After that, we will analyse these results in detail---the constraint force will turn out to have a particularly intuitive interpretation. Finally, we will solve the equations explicitly, obtaining both the particle's trajectory $q(tau)$ and the corresponding Lagrange multiplier $lambda(tau)$. 

We begin with the total action, whose integrand is the total Lagrangian. From the general formalism developed earlier, we know that it is given by
$
  S_T [q;lambda] &= integral dtau [cal(L)(q,dot(q)) - lambda F(q,dot(q))]\
  &= integral dtau [1/2 dot(x)^2 + 1/2 dot(y)^2 - lambda (x^2 + y^2 - 1)].
$
The total Lagrangian is therefore
$
  cal(L)_T (q,dot(q);lambda) = 1/2 dot(x)^2 + 1/2 dot(y)^2 - lambda (x^2 + y^2 - 1).
$
Under a variation $x->x+delta x$, the action $S_T$ changes as
$
  delta S_T = integral dtau [dot(x)delta dot(x) - 2 lambda x delta x]= - integral dtau [ddot(x) + 2 lambda x]delta x,
$
where we have integrated by parts and discarded boundary terms. By symmetry under $x<->y$, the full Euler-Lagrange equations of motion are
$
  ddot(x) = - 2lambda x, wide ddot(y) = - 2lambda y.
$<mechanicsConstraintExampleELeqns>
Before analysing these equations, let us also compute the Hamiltonian and the associated canonical equations. The canonical momenta, derived from $cal(L)_T$, are
$
  p_x = (diff cal(L)_T)/(diff dot(x)) = dot(x), wide p_y = (diff cal(L)_T)/(diff dot(y)) = dot(y)
$
which is clearly an invertible relationship. Consequently, the total Hamiltonian is
$
  cal(H)_T (p,q;lambda) &= p_mu dot(q)^mu - cal(L)_T (q,dot(q)(p);lambda)\
  &= p_x^2 + p_y^2 - [1/2 p_x^2 + 1/2 p_y^2 - lambda (x^2 + y^2 - 1)]\
  &= 1/2 p_x^2 + 1/2 p_y^2 + lambda (x^2 + y^2 -1).
$
The associated canonical equations are
$
  dot(x) &= (diff cal(H))/(diff p_x) = p_x, &wide&& dot(y) &= (diff cal(H))/(diff p_y) = p_y,\
  dot(p)_x &= -(diff cal(H))/(diff x) = - 2lambda x, &wide&& dot(p)_y &= - (diff cal(H))/(diff y)= -2lambda y.
$<mechanicsConstraintExampleCanonicalEqns>
These are clearly consistent with the Euler-Lagrange @mechanicsConstraintExampleELeqns[equations] derived from the Lagrangian formulation. 

Let us now examine the physical content of the equations of motion, particularly the nature of the constraint force. Using vector notation, we may rewrite the Euler-Lagrange equations compactly as
$
  ddot(q)^mu = -2lambda q^mu.
$<mechanicsConstraintExampleELeqn>
This shows that the acceleration is proportional to the position vector, with proportionality factor $-2lambda$. Thus, the force is directed radially---towards the origin if $lambda>0$, or away from it if $lambda < 0$ (though only the former is consistent with circular motion on the constraint surface). 

This behaviour aligns with the geometry of the problem: the constraint surface is the unit circle in configuration space, or more precisely, a submanifold of $RR^2 times RR^2$, parameterised by $(q,dot(q))$, defined by the condition $x^2 + y^2 = 1$. In order to remain on that circle, the particle must experience a force that points normal to the surface at all times. The normal vector to the circle at any point is proportional to the position vector itself, which is exactly the direction of the constraint force appearing on the right-hand side of @mechanicsConstraintExampleELeqn[equation]. 

The magnitude of this foce is determined by the Lagrange multiplier $lambda$. Physically, $lambda$ must take a value that ensures the centripetal acceleration precisely matches the kinematics of the motion. If $lambda$ is too small, the particle spirals outward; if too large, it collapses inward. It is thus intuitively clear that only the correct value of $lambda(tau)$ maintains the particle on a stable circular path.

We now turn to solving the equations explicitly and determining the form of $lambda(tau)$. We consider three different aspects: the purely theoretical existence-of-solutions, a geometric derivation of $lambda$ which will allow the easy solution of the full equations of motion, and lastly a shortcut in which we explicitly parameterise the constraint surface to reduce the system to a simpler one where the constraints are implemented automatically.

+ Structurally, equation @mechanicsConstraintExampleELeqn[] resembles a 2D time-independent Schrödinger equation with potential $V(tau) = 2lambda(tau)$. As a differential equation, solutions always exist for arbitrary $lambda(tau)$, but these trajectories may not satisfy the constraint $F=0$. Only if $lambda(tau)$ is determined as the correct Lagrange multiplier does a given solution remain confined to the circle. So while equation @mechanicsConstraintExampleELeqn[] admits solutions for any function $lambda$, only for the physically meaningful one the system is actually constained. 
 
+ As a second point of view, we can derive solutions of @mechanicsConstraintExampleELeqn[equation] that actually satisfy the constraint by first computing $lambda(tau)$. The Lagrange multiplier turns out to be very simple for this system: a non-negative constant. Intuitively this makes sense---after all, we are restricting a free particle to a circle. This means that there should be no change in velocity _along_ the circle, meaning that the circular motion has constant angular velocity and hence the centripetal force must have a constant magnitude as well. The Lagrange multiplier $lambda$ constitutes that magnitude, and hence, it should be constant. Moreover, it cannot be negative, as otherwise the particle would be accelerated away from the origin, which cannot lead to bound circular motion.

 Though these are compelling arguments, let us make this mathematically rigorous as well. First of all, we note that we can write the constraint $F=0$ as
  $
    q dot q = |q|^2 = 1.
  $ 
  Differentiating once with respect to $tau$ gives us the orthogonality relation
  $
    dot(q) dot q = 0,
  $
  and a second differentiation leads to
  $
    ddot(q) dot q + |dot(q)|^2 = 0.
  $
  Inserting the equation of motion @mechanicsConstraintExampleELeqn[] into this further leads to 
  $
    -2 lambda underbrace(|q|^2,=1) + |dot(q)|^2 quad <==> quad lambda = 1/2|dot(q)|^2.
  $
  This is already an explicit expression for $lambda$. Without much more effort, we can even show that this is constant---let us do this as well. We simply differentiate $lambda$ to find
  $
    dot(lambda) = d/dtau 1/2|dot(q)|^2 = ddot(q) dot dot(q) = -2lambda underbrace(q dot dot(q),=0) = 0.
  $
  This tells us that $lambda = const$, drastically simplifying the remaining @mechanicsConstraintExampleELeqn[equation] into one solvable in closed form. It is also easy to see now that $lambda < 0$ leads to exponential/hyperbolic solutions, which cannot be constrained---we must hence have $lambda>=0$. Note here that, however, $lambda$ can still not be chosen freely. It depends on the initial velocity of the particle, the magnitude of which is preserved in the circular motion, as we will see momentarily. Nonetheless, for $lambda>0$ constant we can write down the solutions as
  $
    q(tau) = q_0 cos sqrt(lambda) tau + 1/sqrt(lambda) dot(q)_0 sin sqrt(lambda) tau,
  $
  where $q_0$ and $dot(q)_0$ are the values of $q(tau)$ and $dot(q)(tau)$ at $tau = 0$, respectively. By our derivations above, they must satisfy $dot(q)_0 dot q_0 = 0$. Moreover, to satisfy the constraint, it needs to hold that
  $
    1 = |q(tau)|^2  = |q_0|^2 cos^2 sqrt(lambda) tau + (|dot(q)_0|^2)/lambda sin^2 sqrt(lambda)tau,
  $
  which can only be valid if $|q_0|^2 = 1$ and $lambda = |dot(q)_0|^2$. We have thus completely solved the constrained dynamics, arriving at the solution
  $
    q(tau) = q_0 cos (|dot(q)_0|tau) + dot(q)_0/(|dot(q)_0|) sin(|dot(q)_0|tau).
  $<constrainedCircleSolutions>
  
  Lastly, for $lambda=0$, we find the constant solutions $q(tau) = q_0 in S^1 subset RR^2$. 


+ Another way to approach the problem is to move to new dynamical variables that automatically encode the constraint. This is not always as straightforward, particularly if $F$ involves velocities, but nonetheless is a powerful trick to keep in mind. In our case, where the constraint is
  $
    x^2 + y^2 = 1,
  $
  we can reduce the system to a single dynamical variable $phi(tau)$, related to $x(tau)$ and $y(tau)$ via
  $
    x(tau) = cos phi(tau), quad y(tau) = sin(tau). 
  $
  By the trigonometric identity $sin^2 + cos^2 = 1$, this reparametrisation always satisfies the constraint, and the system formulated in terms of $phi$ becomes unconstrained. In particular, we have
  $
    dot(x) = - dot(phi) sin phi, quad dot(y) = dot(phi) cos phi,
  $
  implying that the unconstrained free Lagrangian becomes
  $
    cal(L) = 1/2 dot(x)^2 + 1/2 dot(y)^2 = 1/2 dot(phi)^2 sin^2 phi + 1/2 dot(phi)^2 cos^2 phi = 1/2 dot(phi)^2.
  $
  This is the simplest non-trivial system one can encounter: A free particle in one dimension, with equation of motion
  $
    ddot(phi) = 0
  $
  and solutions 
  $
    phi(tau) = dot(phi)_0 tau + phi_0.
  $
  In terms of the variables $x(tau)$ and $y(tau)$, the solutions read
  $
    x(tau) = cos(dot(phi)_0 tau + phi_0), wide y(tau) = sin(dot(phi)_0 tau + phi_0), 
  $
  which are precisely the solutions found earlier in @constrainedCircleSolutions[equation]. 
== Constraints in Field Theories
=== Generalising Mechanical Constraints to Field Theories
In the previous section, we constructed the mathematical tool of Lagrange multipliers to deal with constraints in mechanical systems. In this section, we generalise these results to field theories, again starting from the Lagrangian picture and transitioning to the Hamiltonian formalism for non-singular Lagrangians. The case of singular Lagrangians will be treated separately in the section after the next.

Recall that in field theory, the spacetime Manifold $M$ serves as the parameter space for the dynamical variables. It is locally coordinatised by $x^mu$, and the fields (i.e. dynamical variables) are functions
$
  Phi^A (x) in C^infty (M).
$
The Lagrangian density $cal(L)$ is a local function of the fields and their derivatives,
$
  cal(L) = cal(L)(Phi^A (x), diff_mu Phi^A (x)). 
$
The corresponding unconstrained action is given by
$
  S[Phi^A] = integral_M d^N x thin cal(L)(Phi^A, diff_mu Phi^A).
$
Constraints generalise naturally in this setting: as in mechanics, they take the form of equations that must be satisfied on shell, but now depend on both the fields and their derivatives,
$
  F(Phi^A (x), diff_mu Phi^A (x)) = 0,
$<fieldTheoryConstraintCondition>
for $F$ a function of the fields and their derivatives. 

To enforce this constraint, we follow the same procedure as in mechanics: introduce a new, non-dynamical degree of freedom $lambda(x)$ acting as a Lagrange multiplier. Define the total action as
$
  S_T [Phi^A] = integral_M d^N x thin [cal(L)(Phi^A, diff_mu Phi^A) - lambda F(Phi^A, diff_mu Phi^A)].
$
This promotes the constraint @fieldTheoryConstraintCondition[] to an equation of motion---varying $lambda -> lambda + delta lambda$ arbitrarily, while keeping the fields $Phi^A$ fixed, leads to 
$
  0 attach(=,t:!) delta S_T = integral_M d^N x thin F delta lambda,
$
which requires 
$
  F = 0
$
must hold identically.

As in the mechanical case, the equations of motion for the fields $Phi^A$ are also modified. Under a variation with respect to the fields, $Phi^A -> Phi^A + delta Phi^A$, the total action varies as
$
  0 attach(=,t:!) delta S_T = integral_M d^N x [(delta cal(L))/(delta Phi^A) - lambda (delta F)/(delta Phi^A)]. 
$
We can read off the modified equations of motion to be
$
  [(diff cal(L))/(diff Phi^A) - diff_mu (diff cal(L))/(diff (diff_mu Phi^A))](x) = lambda(x) [(diff F)/(diff Phi^A) - diff_mu (diff F)/(diff (diff_mu Phi^A))](x).
$
As before, the left-hand side gives the standard field equations for the unconstrained system, while the right-hand side encodes the "normal direction" to the constraint surface#footnote[More precisely, this is variation orthogonal to the surface of admissible configurations], scaled pointwise by the multiplier $lambda(x)$ to ensure the appropriate magnitude of the constraint force.

Of course, the Lagrangian picture can also be translated to the Hamiltonian setting. We pursue this for non-singular Lagrangians here, and will return to the singular case---which introduces additional constraints---in a later section. This translation casts the constrained dynamics in terms of phase space variables and prepares the ground for a systematic treatment of constraints, where their classification and algebraic structure will become central in the Dirac-Bergmann formalism.

Suppose that the total Lagrangian has a non-singular Hessian with respect to its velocities,
$
  det((diff^2 cal(L)_T)/(diff(diff_0 Phi^A) diff(diff_0 Phi^B))) != 0.
$
This gives an inveritble relationship between the velocities and the momenta,
$
  Pi_A &= (diff cal(L)_T)/(diff (diff_0 Phi^A)) = (diff cal(L))/(diff(diff_0 Phi^A)) - lambda (diff F)/(diff(diff_0 Phi^A)).
$
The Hamiltonian density associated with $cal(L)_T$ is hence
$
  cal(H)_T &= Pi_A diff_0 Phi^A - cal(L)_T\
  &= [(diff cal(L))/(diff(diff_0 Phi^A)) diff_0 Phi^A - cal(L)] - lambda [(diff F)/(diff (diff_0 Phi^A)) diff_0 Phi^A - F]\
  &= cal(H) - lambda cal(F),
$
where $cal(H)(Phi^A, Pi_A)$ is the unconstrained Hamiltonian associated with $cal(L)$ and $cal(F) (Phi^A, Pi_A)$ is the Legendre transform of the constraint function $F$.

As a last step, we can compute the canonical equations of motion of $cal(H)_T$ as
$
  diff_0 Phi^A &= (delta cal(H)_T)/(delta Pi_A) = (delta cal(H))/(delta Pi_A) - lambda (delta cal(F))/(delta Pi_A),\
  diff_0 Pi_A &= - (delta cal(H)_T)/(delta Phi^A) = - [(delta cal(H))/(delta Phi^A) - lambda (delta cal(F))/(delta Phi^A)].
$
These equations extend the unconstrained canonical equations  @canonicalEqnsFieldTheory[] by terms involving the transformed constraint function $cal(F)$. 

=== Example: Constrained Biscalar Field Theory
==== Lagrangian Picture
Now that we have developed the general formalism for constraints in field theory, let us consider a concrete example. The natural field-theoretic analogue of the mechanical system studied in @constrainedMechanicsExample is a theory of two free real scalar fields, $phi_1$ and $phi_2$, in Minkowski spacetime, subject to a constraint. The unconstrained free Lagrangian density is
$
  cal(L) (phi_i, diff_mu phi_i) = -1/2 eta^(mu nu) diff_mu phi_1 diff_nu phi_1 - 1/2 eta^(mu nu) diff_mu phi_2 diff_nu phi_2,
$
where we omit the inclusion of mass terms for simplicity. As a constraint, we impose the field-theoretic analogue of motion restricted to a circle,
$
  0 = F(phi_1, phi_2) = phi_1^2 + phi_2^2 - 1.
$
This is not meant to be physically deep (though one could view it as a complex scalar field constrained to take values on the unit circle), but it serves as a useful toy model to explore constrained dynamics in a field theory.

To enforce the constraint, we introduce a Lagrange multiplier field $lambda(x)$ to promote $F=0$ to an equation of motion. The total Lagrangian becomes
$
  cal(L)_T (phi_i, diff_mu phi_i;lambda) = -1/2 eta^(mu nu) diff_mu phi_1 diff_nu phi_1 - 1/2 eta^(mu nu) diff_mu phi_2 diff_nu phi_2 - lambda(phi_1^2 + phi_2^2 - 1).
$
Varying the total action with respect to $phi_1$ and $phi_2$, we obtain the modified Euler-Lagrange equations,
$
  Box phi_i = 2 lambda phi_i.
$<constrainedScalarEoM>
On the lef-hand side, we recover the standard wave operator; the right-hand side encodes the constraint force. In this case, it acts as a restoring force, normal to the constraint surface $phi_1^2 + phi_2^2 = 1$ at each spacetime point.

This constraint surface is an infinite-dimensional generalisation of a circle: a submanifold of the field configuration space, formally something like
$
  {(phi_1,phi_2) in L^2 (RR^(N-1)) times L^2 (RR^(N-1)) | phi_1(x)^2 + phi_2(x)^2 = 1 "for all" x in RR^(N-1)},
$
or a more regularised version admitting weak derivatives.

There are two natural perspectives on what the equations of motion @constrainedScalarEoM[] tell us:

+ *Geometric Constraint Force*: The term $-2lambda phi_i$ is the analogue of a centripetal force: it keeps the field values $(phi_1(x),phi_2(x))$ on the unit circle at each spacetime point. At any $x$, it points in the "normal direction" to the constraint surface and adjusts its magnitude via $lambda(x)$ to maintain the constraint dynamically. This is the direct, though abstract, geometric interpretation, mirroring the mechanical case.

+ *Effective mass interpretation*: Alternatively, we can rearrange the equation as
  $
    (Box - m^2 (x))phi_i = 0, quad"with"quad m^2(x) = 2lambda(x).
  $
  This resembles a Klein-Gordon equation with a spacetime-dependent mass. While a bit abstract to visualise, this perspective highlights how the constraint enforces field behaviour by adjusting the local "mass" to keep the fields bounded to the constraint surface. Here, $lambda(x)$ plays the role of a dynamically adjusting mass profile.

==== Hamiltonian Picture
Let us now switch to the Hamiltonian picture, as it illustrates a couple of things nicely which we will see again in the Dirac-Bergmann algorithm. 

To begin, we compute the conjugate momenta
$
  pi_i = (diff cal(L)_T)/(diff (diff_0 phi_i)) = -eta^(0 mu) diff_mu phi_i = diff_0 phi_i.
$
Using these, the total Hamiltonian density can be determined as
$
  cal(H)_T (phi_i, pi_i, diff_k phi_i) &= pi_i diff_0 phi_i - cal(L)_T (phi_i, diff_mu phi_i)\
  &= sum_(i = 1,2)[1/2 (pi_i)^2 + 1/2 (nabla phi_i)^2] + lambda(phi_1^2 + phi_2^2 - 1) 
$
With the Hamiltonian $H_T$ defined as
$
  H_T = integral d^3 x thin cal(H)_T (phi_i (x), pi_i (x), nabla phi_i (x))
$
this leads to the canonical equations of motion
$
  diff_0 phi_i = (delta H_T)/(delta pi_i) = pi_i, quad diff_0 pi_i = -(delta H_T)/(delta phi_i) = Delta phi_i - 2 lambda phi_i
$<canonicalEquationsConstrainedFT>
By inserting the former equation into the latter, it is easily verified that the second-order Euler-Lagrange @constrainedScalarEoM[equation] is recovered. These equations describe the evolution of a field configuration through phase space, i.e. the space of functions $(phi_i,pi_i)$ defined over a spatial hypersurface at fixed $x^0$. This space may be modelled as the product $L^2(RR^(N-1)) times L^2 (RR^(N-1))$, provided we choose appropriate boundary conditions. Put differently, the canonical equations are evolution equations in the vector space $L^2 (RR^(N-1)) times L^2 (RR^(N-1))$.

Similarly, after transforming to the Hamiltonian picture, the constraints now also live in phase space---in the previous sections, where we developed the general theory, we saw that the constraint function $F$ also undergoes a Legendre transformation,
$
  F (phi_i, diff_mu phi_i) -> cal(F)(phi_i, pi_i, diff_k phi_i) = (diff F)/(diff (diff_0 phi_i))diff_0 phi_i - F.
$
The general form reflects that constraints, in the Hamiltonian formalism, are functions on phase space involving $phi$, $pi$ and possibly spatial derivatives thereof. In our case, the constraint is purely positional and so the Legendre transform is trivial, $cal(F) = -F$. 

With the upcoming Dirac-Bergmann algorithm in mind, let us investigate  second-class constraints. Suppose we know that the primary constraint $cal(F) = 0$ is satisfied at some initial time $x_0^0$. This is nothing outlandish---if we seek solutoins that respect the constraint, the initial data must lie on the constraint surface. To remain on that surface, the constraint must be preserved in time. That is, we must demand
$
  diff_0 cal(F) = 0
$
for all configurations such that $cal(F) = 0$. In general, this condition corresponds to the requirement that the Poisson bracket ${cal(F), cal(H)_T}$ vanishes, but here we can compute the time derivative directly as
$
  0 attach(=,t:!) diff_0 cal(F) = 2 phi_i diff_0 phi_i = 2 phi_i pi_i.
$
This leads to a secondary constraint,
$
  0 = cal(G)(phi_i, pi_i) = phi_i pi_i.
$
To ensure this constraint also holds for all times, we now demand $diff_0 cal(G) = 0$, or more concretely
$
  0 = diff_0 cal(G) = (diff_0 phi_i) pi_i + phi_i (diff_0 pi_i).
$
Substituting the canonical equations of motion (cf.  @canonicalEquationsConstrainedFT[]), we find
$
  0 = pi_i pi_i + phi_i [Delta phi_i - 2 lambda phi_i] = [pi_i pi_i + phi_i Delta phi_i] - 2 lambda phi_i phi_i.
$
Using the primary constraint $cal(F) = 0 => phi_1^2 + phi_2^2 = 1$, we can solve for the Lagrange multiplier $lambda$, finding
$
  lambda = 1/2 [pi_1^2 + phi_1 Delta phi_1 + pi_2^2 + phi_2 Delta phi_2].
$
Thus, the consistency condition for the secondary constraint does not yield a new equation to impose, but instead determines the Lagrange multiplier $lambda$ dynamically. This is characteristic of second-class constraints and will appear again in more general form when we turn to the Dirac-Bergmann algorithm.
== The Dirac-Bergmann Algorithm
As we have seen in the case of Maxwell theory, field theories with singular Lagrangians are not esoteric mathematical curiosities---they occur in physically relevant settings. In such cases, the relationship between momenta and velocities becomes non-invertible, making a straightforward Legendre transformation impossible. This, however, is not the end of the Hamiltonian story. A Hamiltonian formulation of Maxwell theory does exist---it simply requires a more intricate construction, via the so-called Dirac–Bergmann algorithm.

While one might feel this effort is unnecessary---after all, Maxwell’s equations themselves provide a consistent first-order evolution system---certain theories, like General Relativity, force the issue. In those contexts, transitioning to the Hamiltonian picture becomes essential for deriving evolution equations. This motivates us to generalise the Hamiltonian formalism to accommodate singular Lagrangians. The Dirac–Bergmann algorithm provides the necessary framework to do just that.

We will first cover the algorithm in mechanics, where the phase space is finite-dimensional and easier to visualise (though one might argue that singular Lagrangians are more esoteric there). Then we generalise to field theories, and finally treat Maxwell theory in the way it deserves---by extracting a clean Hamiltonian field theory from it.
=== Construction of the Algorithm in Mechanics
The setting of the Dirac-Bergmann algorithm is as follows. Suppose we have coordinates $q^mu (tau)$, where $tau$ parameterises a one-dimensional base manifold (e.g. time), and a Lagrangian
$
  cal(L) = cal(L)(q,dot(q)).
$
Suppose further that $cal(L)$ is singular, i.e. that the Hessian determinant of the Lagrangian with respect to velocities
$
  det ((diff^2 cal(L))/(diff dot(q)^mu diff dot(q)^nu)) = 0
$<mechanicsHessianDB>
vanishes at some configurations and for certain values of $tau$. In other words, the relationship between velocities and momenta, given by
$
  p_mu = (diff cal(L))/(diff dot(q)^mu),
$
is not invertible. 

This non-invertibility implies that certain combinations of the momenta are not independent but instead satisfy _constraints_, which we may---without loss of generality---express as
$
  phi.alt_a (q,p) = 0,
$
where the index $a$ labels the constraints, running from $1$ to the number $m$ of such relations. This number equals the number of configuration variables minus the rank of the Hessian above. Note that the constraint functions $phi.alt_a$ may depend on both the coordinates $q$ and momenta $p$, but not on the velocities $dot(q)$.

The geometric picture of this setup is as follows: With singular Lagrangians, we are no longer free to explore the entire phase space $Gamma = {(q,p)}$. Each constraint restricts us to a submanifold $tilde(Gamma)$ of codimension $m$---we lose one dimension per independent constraint.

Let us consider how such constraints can be enforced. We can mimic what happens when constraints are introduced at the Lagrangian level by adding non-dynamical degrees of freedom to the system, acting as Lagrange multipliers. Specifically, we add a term $lambda^a (tau) phi.alt_a (q,p)$ to the Hamiltonian. This introduces the necessary constraint forces into the system, though with undetermined multipliers. So, how do we determine them? 

A natural starting point is to require the initial condition to satisfy all the constraints $phi.alt_a$, i.e. to choose it on the constrained submanifold $tilde(Gamma)$. Beyond that, it is only reasonable to demand that these constraints remain valid under time evolution: we should impose $d/dtau phi.alt_a = 0$ initially. The same goes for the second derivative, the third, and so on.

This might sound alarming---infinitely many constraints?! But not to worry: at some point, usually sooner rather than later, these conditions begin to close among themselves, turning into combinations of constraints already imposed, or leading to equations that determine the $lambda^a$. In effect, one always finds a finite, closed set of constraints. The procedure for systematically identiying them---laid out here, so far, in slightly informal terms---is known as the Dirac-Bergmann algorithm. We shall now embark on its formal introduction. 

The algorithm is best broken down into steps: 
+ *Canonical Momenta:* We begin by computing the Legendre relations between the momenta and velocities,
  $
    p_mu = (diff cal(L))/(diff dot(q)^mu).
  $
  Let $N$ be the amount of coordinates, i.e. $mu=1,...,N$. If the Hessian matrix
  $
    W_(mu nu) := (diff^2 cal(L))/(diff dot(q)^mu diff dot(q)^nu)
  $
  is not invertible, i.e. if
  $
    op("rank") W equiv r < N,
  $
  then the map $dot(q)^mu |-> p_mu$ is not globally invertible. 

  Accordingly, we split the momenta into two subsets:
  - *Regular Sector:* momenta $p_alpha$ for which the relation $dot(q)^alpha <-> p_alpha$ can be inverted. This may require a reparametrisation of the coordinates $q^mu$ to separate regular from singular directions cleanly.
  
  - *Singular Sector:* momenta $p_a$ that cannot be solved for in terms of velocities. This gives rise to _primary constraints_,
    $
      phi.alt_a^((0))(q,p) = 0.
    $
  This split is often straightforward in practice. In many physical systems, the constraints manifest as simple conditions such as a momentum vanishing identically.

  The set of primary constraints $phi.alt_a^((0))$ defines the _primary constraint surface_ in phase space, specifically
  $
    Gamma^((0)) equiv {(q,p) in Gamma mid(|) phi.alt_a^((0)) (q,p) = 0}.
  $

+ *Total Hamiltonian:* On the regular sector, where the Legendre transform is valid, we define the _canonical Hamiltonian_
  $
    cal(H)_C = p_alpha dot(q)^alpha - cal(L).
  $
  We then promote the primary constraints $phi.alt_a^((0)) = 0$ to equations of motion via Lagrange multipliers $lambda^a (tau)$. This yields the _total Hamiltonian_
  $
    cal(H)_T equiv cal(H)_C + lambda^a phi.alt_a^((0)). 
  $
  The time evolution of the canonical variables is governed by
  $
    dot(q)^mu = {q^mu, cal(H)_T} = (diff cal(H)_T)/(diff p_mu), wide dot(p)_mu = {p_mu, cal(H)_T} = - (diff cal(H)_T)/(diff q^mu).
  $
  This defines off-shell evolution---i.e. before all constraints have been enforced.

+ *Constraint Consistency and Secondary Constraints:* To ensure that the primary constraints are preserved under time evolution, we impose 
  $
    d/dtau phi.alt_a^((0)) (q,p) approx 0.
  $
  Here, $approx$ denotes _weak equality_, i.e. equality only on the primary constraint surface $Gamma^((0))$.
  The consequences fall into three categories:
  - *Fixing a Lagrange multiplier:* the consistency condition determines a $lambda^a$. In that case, we simply solve for it and proceed.
  - *Closure of a constraint chain:* the time derivative of a constraint reduces to a combination of existing constraints,
  $
    d/dtau phi.alt_a^((0)) = f(phi^((0))),
  $
  with $f(0)=0$. Then no new constraint arises.
  - *New (secondary) constraint:* if neither of the above applies, the consistency condition yields a genuinely new constraint:
  $
    phi.alt_a^((1))(q,p) equiv d/dtau phi.alt_a^((0))(q,p) approx 0.
  $
  In case secondary constraints do arise, we further define the secondary constraint surface as
  $
    Gamma^((1)) &= {(q,p) in Gamma^((0)) mid(|) phi.alt_a^((1)) approx 0}\
              &= {(q,p) in Gamma mid(|) phi.alt_a^((0)) = 0, med phi.alt_a^((1)) approx 0}.
  $

+ *Higher-Level Constraints:* 
  We iterate the same procedure: demand that secondary constraints be preserved under time evolution, potentially generating _tertiary constraints_ $phi.alt_a^((2))$, and so on. Each step defines a corresponding surface,
  $
    Gamma^((i+1)) equiv {(q,p) in Gamma^((i)) mid(|) phi.alt_a^((i+1)) (q,p) approx 0}.
  $
  The process terminates when no new constraints arise.

At the end of this procedure, we obtain:
- A full set of constraints:
  $
    cal(C) = {phi.alt_a^((0)), phi.alt_a^((1)), phi.alt_a^((2)),...},
  $
- A nested sequence of constraint surfaces...
  $
    Gamma supset Gamma^((0)) supset Gamma^((1)) supset Gamma^((2)) supset ...,
  $
- ...whose intersection defines the _final constraint surface_: 
  $
    Gamma^infty equiv sect.big_i Gamma^((i)) = {(q,p) in Gamma mid(|) phi.alt_a^((i))(q,p) = 0, med forall phi.alt_a^((i)) in cal(C)}. 
  $
By construction, all Lagrange multipliers have now either been fixed or remain completely undetermined, in which case they generate gauge freedom. These unfixed multipliers do not affect evolution on $Gamma^infty$. 

And now the key result: If the initial data satisfy all constraints $phi.alt_a^((i)) = 0$, then the evolution under $cal(H)_T$ will preserve them. That is, the trajectory remains entirely within $Gamma^infty$, and all constraints are automatically maintained over time.

It is rather simple to see that this is the case: the way we constructed $Gamma^infty$ directly implies that on this surface, any derivative of the primary constraints with respect to $tau$, regardless of order (including the "zeroth derivative", the constraint itself), is zero. This can only be the case if $phi.alt_a^((0)) = 0$ for all $tau$, and hence the (primary) constraints are always satisfied if they are satisfied initially. 

There are two more questions we should at least mention here. Firstly, why do we only include constraint terms of the primary constraints in the Lagrangian, but not of higher order constraints? This can be answered as follows. The primary constraints are those that arise from the singularity of the Lagrangian. Not all constraints are independent---rather, some combinations of them must vanish. This tells us the following: Our evolution cannot happen in the whole (naive) phase space $Gamma = {(q,p)}$, but only a subset thereof. Not all choices of momentum values make sense, as not all of them can be inverted for their velocity configurations. In other words, not all choices of momenta even have a solution of the equations of motion $q(tau)$ that could produce those momenta. This means, parts of the unconstrained phase space $Gamma$ are _unphysical_---they don't really exist. So, to have something meaningful, we need to promote the primary constraints.

Secondary and higher order constraints are of a fundamentally different nature, though. Rather than being a strict requirement from the singular Legendre map, they simply arise from wanting the primary constraints to _remain_ fulfilled if they are fulfilled initially. If there are such constraints, that tells us that even if a phase space position $(q,p)$ satisfies the primary constraints, i.e. lies on $Gamma^((0))$, it _might_ still drift off that hypersurface through evolution---hence isn't a physically meaningful configuration either. A particle's path cannot be physically meaningful in one instant and become unreason in the second. We thus don't only want solutions that are physically meaningful initially, but that remain physically meaningful---or in more precise terms, we only care about phase space trajectories that satisfy the primary constraints at all times. This simply reduces the relevant portions of phase space to $Gamma^infty$. This submanifold constitutes all phase space configurations for which the evolution preserves the primary constraints for all $tau$. 

At this point, notice the fundamental difference between primary and higher order constraints: primary constraints need to be satisfied for a point $(q,p)$ to be physically meaningful at all---thus must be part of the evolution and hence the Hamiltonian. Secondary constraints, however, only further "weed out" the trajectories that still fail to "stay on course", so to speak. 

Then, the second question: We constructed a mathematically meaningful way of extracting Hamiltonian evolution from a singular Lagrangian in a rather canonical way. But, how can we trust that we did not change the physical content of the theory while doing it? Why does the new formulation describe the same physics? _This is something I would like to maybe get into at one point, not now---I'm not sure if I could properly prove this as of right now. But, I feel like this question deserves to at least be asked here, since it is rather fundamental. An intuitive answer was partially given by answering the first question, but a rigorous answer remains to be provided._

It is rather simple to see that this is the case: the way we constructed $Gamma^infty$ directly implies that on this surface, any derivative of the primary constraints with respect to $tau$, regardless of order (including the constraints themselves), is zero. This implies that $phi.alt_a^((0))=0$ for all $tau$, i.e. the primary constraints vanish identically along the evolution, and hence are always satisfied if they are satisfied initially.

There are two more questions we should at least mention here. Firstly, why do we only include constraint terms of the primary constraints in the Lagrangian, but not of higher order constraints? This can be answered as follows. The primary constraints are those that arise from the degeneracy of the Lagrangian. Not all momenta are independent---rather, some combinations of them must vanish. This tells us the following: our evolution cannot happen in the whole (naive) phase space $Gamma = {(q,p)}$, but only a subset thereof. Not all choices of momentum make values make sense, as not all of them can be inverted for their velocity configurations (notice that the converse situation always makes sense---any velocities always lead to valid momenta). In other words, not all choices of momenta even correspond to consistent solutions $q(tau)$ of the equations of motion that could produce those momenta. This means, parts of the unconstrained phase space $Gamma$ are unphysical---they do not really exist. So, to have something meaningful, we need to incorporate the primary constraints directly into the Hamiltonian through Lagrange multipliers.

Secondary and higher order constraints are of a fundamentally different nature. Rather than being a strict requirement from the Legendre map, they simply arise from wantin the primary constraints to _remain_ fulfilled if they are fulfilled initially. If there are such constraints, that tells us that even if a phase space position $(q,p)$ satisfies the primary constraints, i.e. lies on $Gamma^((0))$, it _might_ still drift off that hypersurface through evolution---hence it isn't a physically meaningful trajectory either. A particle's path cannot be physically meaningful in one instant and beocome unphysical in the next. We thus don't only want solutions that are physically meaningful initially, but that remain physically meaningful---more precisely, we only care about phase space trajectories that satisfy the primary constraints at all times. These are precisely those trajectories that satisfy all higher order constraints as well, which restricts physically admissible trajectories to lie within $Gamma^infty$. This submanifold constitutes all phase space configurations for which the evolution preserves the primary constraints for all $tau$. 

At this point, notice the fundamental difference between primary and higher order constraints: primary constraints need to be satisfied for a point $(q,p)$ to be physically meaningful at all---thus must be part of the evolution and hence included in the Hamiltonian. Secondary constraints, however, only further "weed out" the trajectories that still fail to remain within the constraint surface under time evolution, even when the primary constraints are promoted. 

Then, the second question: We constructed a mathematically meaningful way of extracting Hamiltonian evolution from a singular Lagrangian in a rather canonical way. But, how can we trust that we did not change the physical content of the theory while doing it? Why does the new formulation describe the same physics? _This is something I would like to get into at some point in the future, not now---I'm not sure if I could properly prove this as of right now. But, I feel like this question deserves to be at least asked here, since it is rather fundamental. An intuitive answer was partially given by answering the first question, but a rigorous answer remains to be provided._
=== Example: Mechanical System with a Singular Lagrangian
To elucidate the somewhat abstractly introduced Dirac-Bergmann algorithm, let us explore a simple mechanical system with a singular Lagrangian. We consider a two-dimensional system with only one dynamical degree of freedom, defined by
$
  cal(L) = 1/2 dot(x)^2  - V(x,y).
$
Technically speaking, $y(tau)$ is an externally defined control function. We will not concern ourselves with that here and still treat it as a coordinate $y(tau)$.

Let us begin by reviewing the Euler-Lagrange equations to get a sense of the theory. The equation of motion for $x(tau)$ reads
$
  ddot(x) = -diff_x V(x,y).
$
Treating $y(tau)$ the same way gives another equation, 
$
  0 = -diff_y V(x,y)
$
This tells us two things: on one hand, $y(tau)$ is completely arbitrary---on the other, the potential $V$ must be independent of it. In other words, consistency of the Euler-Lagrange equations forces $diff_y V = 0$, so our Lagrangian does not involve $y$ at all. We are thus describing a one-dimensional particle in a potential, with an additional nondynamical function $y(tau)$ that has no effect on the dynamics of $x$. 

Regardless of the questionable physical relevance of such a system, it is clear that the Lagrangian is singular. To see this, compute the Hessian determinant with respect to the velocities,
$
  det [mat(diff_(dot(x)) diff_dot(x), diff_(dot(x)) diff_dot(y);diff_(dot(y)) diff_dot(x), diff_(dot(y)) diff_dot(y)) cal(L)] = det mat(1,0;0,0) = 0.
$
This means that the Legendre transformation cannot be invertible.

We now now compute the canonical momenta and identify the primary constraints. For $x$, we have
$
  p_x = (diff cal(L))/(diff dot(x)) = dot(x).
$
For $y$, we get
$
  p_y = (diff cal(L))/(diff dot(y)) = 0.
$
This tells us that $p_y = 0$ is a primary constraint. The regular sector consists of the pair $(x,p_x)$, and the singular sector of $(y,p_y)$ with the constraint
$
  phi.alt_1^((0))(q,p) = p_y = 0.
$
This defines the primary constraint surface
$
  Gamma^((0)) = {(x,y,p_x,p_y) | p_y = 0}.
$

Next, we construct the canonical and total Hamiltonians. The former is simply given by
$
  cal(H)_C = p_x dot(x) - cal(L) = 1/2 p_x^2 + V(x).
$
To obtain the total Hamiltonian, we add a Lagrange multiplier $lambda(tau)$ enforcing the constraint,
$
  cal(H)_T = cal(H)_C + lambda^a phi.alt_a^((0)) = 1/2 p_x^2 + V(x) + lambda p_y
$
The canonical equations of motion are follow as
$
  dot(x) &= (diff cal(H)_T)/(diff p_x) = p_x, &wide&& dot(p)_x &= -(diff cal(H)_T)/(diff x) = -diff_x V(x),\ \
  dot(y) &= (diff cal(H)_T)/(diff p_y) = lambda, &&& dot(p)_y &= -(diff cal(H)_T)/(diff y) = 0. 
$
In particular, the multiplier $lambda$ governs the evolution of $y(tau)$ via $dot(y) = lambda$. 

We then move to secondary constraints. Demanding consistency of $phi.alt_1^((0))$ yields
$
  0 attach(=,t:!) d/dtau phi.alt_1^((0)) (q,p) = dot(p)_y = 0,
$
which is trivially satisfied by the equations of motion. Hence, we do not get any secondary or higher order constraints, and the multiplier $lambda$ remains arbitrary. By the relation $dot(y)(tau) = lambda(tau)$ we get from the equations of motion, this means that $y(tau)$ is arbitrary as well, reflecting our findings from the Lagrangian formulation. 

So, while this may not be the most interesting system physically, it illustrates some key features of the Dirac-Bergmann algorithm in a simple mechanical context: the appearance of primary constraints, arbitrainess of certain degrees of freedom, and the absence of further consistency conditions.

=== Example II: Mechanical System with a Gauge Potential

Things become more interesting in field theory, where derivatives with respect to parameters beyond time can give rise to richer constraint structures. To prepare for that, we will generalise the Dirac-Bergmann algorithm to field-theoretic systems in the next section. Before doing so, however, we examine a second mechanical example---slightly more intricate than the previous---to further develop our intuition and illustrate the algorithm in action once more. 

The Lagrangian of this system is given by
$
  cal(L) = 1/2 (dot(x) - y)^2.
$
It describes one dynamical variable, $x(tau)$, and one non-dynamical variable, $y(tau)$. Structurally, this resembles a kinetic term modified by a gauge potential.

Let us begin with the Lagrangian formulation. The Euler-Lagrange equations read
$
  0 &= underbrace((diff cal(L))/(diff x),=0) - d/dtau (diff cal(L))/(diff dot(x)) = ddot(x) - dot(y) wide&& ==> wide ddot(x) &= dot(y)\
  0 &= (diff cal(L))/(diff y) - d/dtau underbrace((diff cal(L))/(diff dot(y)),=0) = -(dot(x)-y) wide&&==> wide dot(x) &= y
$
Since the second equation implies the first, the equations of motion reduce to $dot(x) = y$, with $y(tau)$ arbitrary. The solution is simply
$
  x(tau) = c + integral dtau thin y(tau) 
$
where $c$ is an integration constant determined by initial data.

Now we turn to the Hamiltonian picture. The canonical momenta are
$
  p_x &= (diff cal(L))/(diff dot(x)) = dot(x) - y quad ==> quad dot(x) = p_x + y\
  p_y &= (diff cal(L))/(diff dot(y)) = 0.
$
We obtain a primary constraint, $phi.alt_1^((0)) = p_y = 0$, and the Legendre transform with respect to $dot(x)$ yields the canonical Hamiltonian
$
  cal(H)_C = p_x dot(x) - cal(L) = 1/2 p_x^2 + y p_x
$
Enforcing the constraint $phi.alt_1^((0))=1$ via a Lagrange multiplier $lambda(tau)$ leads to the total Hamiltonian
$
  cal(H)_T = 1/2 p_x^2 + y p_x + lambda p_y.
$
This Hamiltonian governs the evolution through phase space. The corresponding canonical equations are
$
  dot(x) &= (diff cal(H))/(diff p_x) = p_x + y, &wide&& dot(p)_x &= -(diff cal(H))/(diff x) = 0,\
  dot(y) &= (diff cal(H))/(diff p_y) = lambda, &&& dot(p)_y &= - (diff cal(H))/(diff y) = -p_x. 
$
These will guide our analysis of consistency conditions.

Demanding preservation of the primary constraint $phi.alt_1^((1))=0$ under time evolution yields
$
  0 attach(=,t:!) d/dtau phi.alt_1^((0)) = dot(p)_y = - p_x
$
by making use of the canonical equations derived above. This is a genuinely new constraint---it is neither a combination of known constraints nor does it fix $lambda$. Thus, we have to introduce the secondary constraint
$
  phi.alt_1^((1)) = p_x = 0.
$
At this point, one might start to feel a bit suspicious. These two conditions limit us to the secondary constraint surface
$
  Gamma^((1)) = {(x,y,p_x,p_y) mid(|) p_x = p_y = 0}.
$
While this may seem overly restrictive at first glance, it is in fact consistent with the Lagrangian solution: $p_x = 0$ implies $dot(x) = y$, in perfect agreement with the Euler-Lagrange equation. 

There are no tertiary constraints---the condition 
$
  0 attach(=,t:!) d/dtau phi.alt_1^((1)) = dot(p)_x = 0
$
is automatically satisfied according to the canonical equations. Notice that the Lagrange multiplier $lambda$ remained undetermined. 

In summary, we have now constructed the constraint surface $Gamma^infty = Gamma^((1))$, on which the evolution remains constrained. Hence, such solutions always satisfy $p_x = p_y = 0$, reducing the canonical equations to 
$
  dot(x) = y, quad dot(y) = lambda,
$
with the other two equations being fulfilled trivially. Since $lambda$ was found to be arbitrary, this reproduces the Lagrangian picture exactly. Compared to the previous example, this system features a primary _and_ a secondary constraint, giving a hint of more intricate behaviour---but still falling short of the complexity encountered in field theory.
=== Dirac-Bergmann in Field Theory
Without motivation, here we simply state the procedure behind the Dirac-Bergmann algorithm for field theories. As a starting point, suppose we have dynamical variables (fields) $Phi^A (x)$, which take values on an $N$-dimensional spacetime manifold with local coordinates $x^mu$. We further consider a Lagrangian density
$
  cal(L) = cal(L)(Phi^A, diff_mu Phi^A)
$
which defines an action given by
$
  S[Phi] = integral d^N x thin cal(L)(Phi^A (x), diff_mu Phi^A (x)). 
$
We work under the assumption that the Lagrangian is singular, i.e. that the velocity Hessian
$
  W_(A B) = (diff^2 cal(L))/(diff (diff_0 Phi^A) diff(diff_0 Phi^B))
$
does not have full rank and hence zero determinant. The procedure then goes as follows:

+ *Canonical Momenta:* As a first step, the Legendre relations between the momenta and velocities,
  $
    Pi_A = (diff cal(L))/(diff(diff_0 Phi^A))
  $
  are to be computed. Using $r equiv op("rank") W_(A B) < N$, we can split the momenta into two subsets:
  - *Regular Sector:* Momenta $Pi_cal(A)$ for which the relation $diff_0 Phi^cal(A) <-> Pi_cal(A)$ can be inverted. This may require a reparametrisation of the fields $Phi^A$ to separate regular from singular directions cleanly.

  - *Singular Sector:* Momenta $Pi_a$ that cannot be solved for in terms of the velocities $diff_0 Phi^A$. This gives rise to _primary constraints_,
    $
      phi.alt_a^((0)) (Phi^A, Pi_A) = 0
    $
    In many cases, these manifest directly as identities such as $Pi_a = 0$. 

  Just as in mechanics, the set of primary constraints $phi.alt_a^((0))$ defines the _primary constraint surface_ in phase space, given by
  $
    Gamma^((0)) equiv {(Phi^A, Pi_A) in Gamma mid(|) phi.alt_a^((0))(Phi^A, Pi_A) = 0}.
  $
  Notably, in contrary to mechanics, this surface (as well as the space it $Gamma$ it lives in), is infinite-dimensional.  

+ *Canonical and Total Hamiltonian:* With respect to the regular sector, the Legendre transformation can be performed, yielding the _canonical Hamiltonian density_ $cal(H)_C$, given by
  $
    cal(H)_C = Pi_cal(A) diff_0 Phi^cal(A) - cal(L).
  $
  To enforce the primary constraints, we introduce Lagrange multipliers $lambda^a (x)$ and use them to promote the constraints to equations of motion by defining the _total Hamiltonian density_
  $
    cal(H)_T = cal(H)_C + lambda^a phi.alt^((0))_a.
  $
  From this, we further define the total Hamiltonian as
  $
    H_T = integral d^(N-1)x thin cal(H)(Phi^A, Pi_A)
  $

  This Hamiltonian is what we use to define off-shell evolution; the equations of motion read
  $
    diff_0 Phi^A = {Phi^A, H_T} = (delta cal(H)_T)/(delta Pi_A), wide quad diff_0 Pi_A = {Pi_A, H_T} = - (delta cal(H)_T)/(delta Phi^A).
  $

+ *Constraint Consistency and Secondary Constraints:*  
  We demand constraints to remain satisfied if they are satisfied initially. More concretely, we require
  $
     diff_0 phi.alt_a^((0)) (Phi^A, Pi_A) attach(approx,t:!) 0.
  $<DBconsistencyCondition>
  where $approx$ indicates weak equality, i.e. equality only on $Gamma^((0))$. There are three possible outcomes for this:

  + *Determination of a Multiplier:* Such a condition can lead to a relationship between a multiplier $lambda^a$ and the canonical variables. In this case, we solve for it and may ignore the condition.

  + *Closure:* When the expression $diff_0 phi.alt_a^((0))$ can be written as a function of known constraints,
    $
      diff_0 phi.alt_a^((0)) = f(phi.alt_a^((0)))
    $
    ---with $f(0)=0$---the condition is automatically satisfied. We may thus ignore it.

  + *Emergence of a Secondary Constraint* If neither of the first two outcomes apply, then the condition @DBconsistencyCondition[] constitutes a new _secondary constraint_ $phi.alt_a^((1)) approx 0$, where
    $
      phi.alt_a^((1)) (Phi^A, Pi_A) = diff_0 phi.alt^((0)) (Phi^A, Pi_A).
    $
  If there are secondary constraints, they lead to the definition of the _secondary constraint surface_
  $
    Gamma^((1)) = {(Phi^A, Pi_A) in Gamma^((0)) mid(|) phi.alt_a^((1))(Phi^A, Pi_A) approx 0}.   
  $

+ *Higher-Order Constraints:* The previous step is iterated: demand that secondary constraints be preserved under time evolution, potentially generating _tertiary constraints_ $phi.alt_a^((2))$, and so on. At each step, define the corresponding surface,
  $
    Gamma^((i+1)) = {(Phi^A, Pi_A) in Gamma^((i)) mid(|) phi.alt_a^((i+1))(Phi^A, Pi_A) approx 0}
  $
  The procedure is stopped when an iteration does not introduce new constraints.

Just as in mechanics, completion of this procedure leads to:
- A full set of constraints:
  $
    cal(C) = {phi.alt_a^((0)), phi.alt_a^((1)), phi.alt_a^((2)),...}
  $
- A nested sequence of (infinite dimensional) constraint surfaces,
  $
    Gamma supset Gamma^((0)) supset Gamma^((1)) supset ... supset Gamma^infty,
  $
- where $Gamma^infty$ is the _final constraint surface_, defined as
  $
    Gamma^infty = sect.big_i Gamma^((i)) = {(Phi^A,Pi_A) in Gamma | phi.alt_a^((i)) (Phi^A, Pi_A) = 0, med forall phi.alt_a^((i)) in cal(C)}.
  $

The consequences are the same as in mechanics: By construction, all Lagrange multipliers have now either been fixed or remain completely undetermined, in which case they generate gauge freedom. These unfixed multipliers do not affect evolution on $Gamma^infty$. Further, starting from an initial configuration $(Phi^A,Pi_A) in Gamma^infty$ and evolving along the flow field generated by $cal(H)_T$, the constraints are preserved and we are in a formulation that is physically equivalent to the Lagrangian picture. 

Put more concisely,
$
  (Phi^A (0, x^i), Pi_A (0,x^i)) in Gamma^infty, wide diff_0 Phi^A = (delta cal(H)_T)/(delta Pi_A), wide diff_0 Pi_A = - (delta cal(H)_T)/(delta Phi^A)
$
implies that
$
  (diff cal(L))/(diff Phi^A) - diff_mu (diff cal(L))/(diff (diff_mu Phi^A)) = 0.
$
This construction is what will finally allow us to derive the Hamiltonian formulation of Maxwell theory in the next section.

=== Example: Maxwell Theory, Revisited Properly
We have endeavoured plenty to get to this point, but we can now do it: Treat Maxwell theory in the Hamiltonian picture, with the Dirac-Bergmann algorithm. For completeness, let us start from the beginning. In the mostly plus signature of the Minkowskian metric, the Maxwell Lagrangian density reads
$
  cal(L) = -1/4 F_(mu nu) F^(mu nu),
$
where $F_(mu nu)$ is the field strength associated with the 1-form $A = A_mu dx^mu$,
$
  F_(mu nu) = diff_mu A_nu - diff_nu A_mu.
$
To make the subsequent analysis simpler, we provide a number of alternative but equivalent expressions for this Lagrangian:
$
  cal(L) &= -1/2 F_(0 i) F^(0 i) - 1/4 F_(i k) F^(i k)\
  &= 1/2 E_i E_i - 1/4 epsilon_(i k l) B_l epsilon_(i k p) B_p\
   &= 1/2 (vE^2 - vB^2)
$
Here we made use of the expressions for the electric field components $E_i$ and magnetic field components $B_i$, which are given by
$
  E_i = diff_0 A_i - diff_i A_0  = F_(0 i), wide epsilon_(i j k) B_k = diff_i A_j - diff_j A_i  = F_(i j),
$
as well as the identity
$
  epsilon_(i j k) epsilon_(i j l) = 2 delta_(k l).
$

For inclusiveness, we restate the Euler-Lagrange equations for Maxwell theory, which read
$
  diff_mu F^(mu nu) = 0.
$
Dissecting this a bit further, we find two types of equations:
$
  0 = diff_i F^(i 0) = nabla dot vE,
$
which is simply Gauss' law, and 
$
  0 = diff_0 F^(0 k) + diff_i F^(i k) = diff_0 E_k + epsilon_(i k j) diff_i B_j = (diff_0 vE - nabla times vB)_k.
$
In the following, we will rediscover these equations in the pursuit of the Dirac-Bergmann algorithm. 

+ *Canonical Momenta:* Without further ado, let us proceed with its procedure. We begin by evaluating the Legendre map to derive the momenta. For the spatial momenta, we get
  $
    Pi_i = (diff cal(L))/(diff (diff_0 A_i)) = -(diff F_(0 i))/(diff (diff_0 A_i)) F^(0 i) =  - F^(0 i) = E_i.
  $
  In words: the momenta corresponding to the spatial velocities $diff_0 A_i$ are the electric field components $E_i$. We may rearrange this equation to find the velocity in terms of the momentum,
  $
    E_i = F_(0 i) = diff_0 A_i - diff_i A_0 quad <=> quad diff_0 A_i = E_i + diff_i A_0 
  $
  For the temporal momentum, we find
  $
    Pi_0 = (diff cal(L))/(diff(diff_0 A_0)) = 0.
  $
  To summarise, we have found three regular/invertible momenta, $E_i$, and a singular momentum $Pi_0$. The latter gives us a primary constraint,
  $
    phi.alt_1^((0)) = Pi_0 = 0.
  $

+ *Canonical and Total Hamiltonian Densities:* Performing the partial Legendre transform with respect to the regular sector yields the canonical Hamiltonian density
  $
    cal(H)_C &= E_i diff_0 A_i - cal(L) = E_i (E_i + diff_i A_0) - 1/2 (vE^2 - vB^2)\
    &= 1/2 (vE^2 + vB^2) + E_i diff_i A_0 
  $
  We are now ready to define the total Hamiltonian density by introducing a Lagrange multiplier $lambda(x)$ to enforce the constraint $phi.alt_a^((0)) = Pi_0 = 0$. It reads
  $
    cal(H)_T = cal(H)_C + lambda Pi_0  =1/2 (vE^2 + vB^2) + E_i diff_i A_0 + lambda Pi_0. 
  $
  The associated Hamiltonian is given by
  $
    H_T = integral d^(N-1)x thin cal(H)_T = integral d^(N-1)x [1/2 (vE^2 + vB^2) - A_0 (nabla dot vE) + lambda Pi_0],
  $
  where we integrated by parts, discarding vanishing boundary terms. The reason for doing so will become evident later.
+ *Secondary Constraints:* Demanding consistency of $Pi_0 = 0$ leads to
  $
    0 attach(=,t:!) diff_0 Pi_0 = - (delta H_T)/(delta A_0) = nabla dot vE,
  $
  which is simply Gauss' law. This constitutes a secondary constraint,
  $
    phi.alt_1^((1)) = nabla dot vE = 0.
  $
  Observe here that $H_T$ contains the term $A_0 (nabla dot vE)$, which is of the form of a constraint term for $phi.alt_1^((1))$---with $A_0$ acting as the Lagrange multiplier enforcing it. 

+ *Higher Order Constraints:* Imposing that $phi.alt_1^((1))$ be consistent as well, we find
  $
    0 &attach(=,t:!) diff_0 diff_i E_i = diff_i (diff_0 E_i) = -diff_i (delta H_T)/(delta A_i) \
    &= -diff_i lr([underbrace((diff cal(H)_T)/(diff A_i),=0) - diff_k (diff cal(H)_T)/(diff(diff_k A_i))],size:#56%) = diff_i diff_k F^(i k) = 0
  $
  This means that the constraint chain closes, and we do not get any tertiary constraints. 

Notice that we have not found any conditions on the two Lagrangian multipliers $A_0$ and $lambda$. They hence remain free, and we may choose them to be whatever we please to make calculations easier. 

Let us now compute the canonical equations of motion. For the spatial components these read
$
  diff_0 A_i = (delta H_T)/(delta E_i) = E_i + diff_i A_0, wide diff_0 E_i = - (delta H_T)/(delta A_i) = -lr([underbrace((diff cal(H)_T)/(diff A_i),=0) - diff_k (diff cal(H)_T)/(diff(diff_k A_i))],size:#56%) = diff_k F^(k i)
$
In index-free notation, these are simply the equations
$
  diff_0 vA = vE + nabla A_0 quad <=> quad vE = diff_0 vA - nabla A_0,
$
$
  diff_0 vE = nabla times vB,
$
which we already know and love from Maxwell theory in its original formulation. 

The canonical equations for the temporal components read
$
  diff_0 A_0 = (delta cal(H)_T)/(delta Pi_0) = lambda, wide diff_0 Pi_0 = -(delta cal(H)_T)/(delta A_0) = diff_i E_i
$
For configurations on $Gamma^infty$, the right-hand side of the latter equation vanishes and hence $Pi_0 = 0$ is preserved automatically. The former show that $lambda$ governs the (otherwise unconstrained) time evolution of $A_0$. 
=== Revisiting Poisson Brackets
Let us revisit the definition of Poisson brackets in field theory from first principles and establish a set of practical computation rules.

Consider a Hamiltonian system with canonical fields $Phi^A (x)$ and their conjugate momenta $Pi_A (x)$, defined on a spatial hypersurface of dimension $N-1$. For two differentiable local functionals $cal(A)[Phi^A, Pi_A]$ and $cal(B)[Phi^A, Pi_A]$, the Poisson bracket is defined as
$
  {cal(A), cal(B)} = integral d^(N-1)x [(delta cal(A))/(delta Phi^A (x)) (delta cal(B))/(delta Pi_A (x)) - (delta cal(B))/(delta Phi^A (x)) (delta cal(A))/(delta Pi_A (x))],
$
This generalises the finite-dimensional Poisson bracket to infinite dimensional field theory, replacing partial with functional derivatives and summations with spatial integrals.

From this definition, the equal-time fundamental Poisson brackets follow as#footnote[Here we  consider $Phi^A$ and $Pi_B$ as functions of $x^i$ only, at fixed $x^0$. We have the implicit convention that all participating functions are evaluated at the same value of $x^0$.]
$
  {Phi^A (x), Pi_B (y)} &= delta^A_B delta^(N-1)(x-y),\ {Phi^A (x),Phi^B (y)} &= {Pi_A (x), Pi_B (y)} = 0.
$<basicPoissonBrackets>
There are a few computation rules that one can further derive to make calculations with Poisson brackets easier. Let us now go through these.

- *Leibniz Rule* The Poisson bracket satisfies a Leibniz rule,
  $
    {cal(A) cal(B), cal(C)} = cal(A){cal(B),cal(C)} + {A,cal(C)}cal(B).
  $
  This is rather simple to see, as by the Leibniz rule for the functional derivative we find
  $
    {cal(A) cal(B), cal(C)} &= integral d^(N-1)x [(delta (cal(A B)))/(delta Phi^A (x)) (delta cal(C))/(delta Pi_A (x)) - (delta cal(C))/(delta Phi^A (x)) (delta (cal(A B)))/(delta Pi_A (x))]\
    &= integral d^(N-1)x [((delta cal(A))/(delta Phi^A (x)) cal(B) + cal(A) (delta cal(B))/(delta Phi^A (x)))(delta cal(C))/(delta Pi_A (x))\
    &#h(5.5em)- (delta cal(C))/(delta Phi^A (x))((delta cal(A))/(delta Pi_A (x)) cal(B) + cal(A) (delta cal(B))/(delta Pi_A (x)))]\
    &= integral d^(N-1)x [(delta cal(A))/(delta Phi^A (x)) (delta cal(C))/(delta Pi_A (x)) - (delta cal(C))/(delta Phi^A (x)) (delta cal(A))/(delta Pi_A (x))] cal(B) \
    &quad + integral d^(N-1)x thin cal(A)[(delta cal(B))/(delta Phi^A (x)) (delta cal(C))/(delta Pi_A (x)) - (delta cal(C))/(delta Phi^A (x)) (delta cal(B))/(delta Pi_A (x))]\
    &= {cal(A),cal(C)}cal(B) + cal(A){cal(B),cal(C)},
  $
  where we made use of the fact that $cal(A)$ and $cal(B)$ do not depend on the integration variable $x^i$ and can hence be pulled out of the integrals. Naturally, by the anti-symmetry of the Poisson bracket, this identity generalises to the second argument,
  $
    {cal(A),cal(B C)} = cal(B){cal(A),cal(C)} + {cal(A),cal(B)}cal(C).
  $ 
  Notably, this is the same structure as one finds for commutators---as is a necessity for canonical quantisation.

- *Power Rule* From the Leibniz rule, the power rule
  $
    {cal(A)^n, cal(B)} = n cal(A)^(n-1){cal(A), cal(B)}, quad n>= 2.
  $
  follows. This can be shown by induction, starting with the base case $n=2$:
  $
    {cal(A)^2,cal(B)} = {cal(A A),cal(B)} = cal(A){cal(A),cal(B)} + {cal(A),cal(B)}cal(A) = 2cal(A){cal(A),cal(B)}.
  $
  We can then perform the induction step as
  $
    {cal(A)^(n+1),cal(B)} &= {cal(A)cal(A)^n,cal(B)} = cal(A){cal(A)^n, cal(B)} + {cal(A),cal(B)} cal(A)^n\
    &= n cal(A)cal(A)^(n-1){cal(A),cal(B)} + {cal(A),cal(B)}cal(A)^n = (n+1)cal(A)^n {cal(A),cal(B)}.
  $

- *Interaction with Derivatives* Suppose we have a $y^i$-dependent local functional $cal(A)(y)$. We may still compute Poisson brackets such as
  $
    {cal(A)(y), cal(B)} = integral d^(N-1)x [(delta cal(A)(y))/(delta Phi^A (x)) (delta cal(B))/(delta Pi_A (x)) - (delta cal(B))/(delta Phi^A (x)) (delta cal(A)(y))/(delta Pi_A (x))]
  $
  (and in fact, we already did this when computing @basicPoissonBrackets[identities]). This bracket ${cal(A)(y),cal(B)}$ is, in particular, a function of $y^i$. Thus, we can take derivatives with respect to it---let us see how they interact with the bracket. Making use of the fact that partial and functional derivatives commute, we derive
  $
    diff_(y^i) {cal(A)(y),cal(B)} &= integral d^(N-1)x thin diff_(y^i)[(delta cal(A)(y))/(delta Phi^A (x)) (delta cal(B))/(delta Pi_A (x)) - (delta cal(B))/(delta Phi^A (x)) (delta cal(A)(y))/(delta Pi_A (x))]\
    &= integral d^(N-1)x thin [(delta (diff_(y^i)cal(A)(y)))/(delta Phi^A (x)) (delta cal(B))/(delta Pi_A (x)) - (delta cal(B))/(delta Phi^A (x)) (delta (diff_(y^i)cal(A)(y)))/(delta Pi_A (x))]\
    &= {diff_i cal(A)(y), cal(B)}. 
  $

- *Interaction with Integrals:* By the same reasoning as for derivatives---because integration and functional derivation can be exchanged---it can be shown that 
  $
    integral d^(N-1)y {cal(A)(y),cal(B)} = {integral d^(N-1)y thin cal(A)(y), cal(B)}. 
  $
  This is particularly useful for brackets with the Hamiltonian, as one has
  $
    {cal(A), H} = {cal(A), integral d^(N-1)y thin cal(H)(y)} = integral d^(N-1)y thin {cal(A),cal(H)(y)}
  $
  which turns calculations with the Hamiltonian into calculations with the Hamiltonian density which are typically simpler to handle.

For convenience, we summarise the above results in a cheat sheet:\
#box(width:1fr, stroke:1.5pt, inset:2em, outset:-1em)[
  $
    {cal(A)cal(B),cal(C)} = cal(A){cal(B),cal(C)} + {cal(A),cal(C)}cal(B),
  $
  $
    {cal(A)^n,cal(B)} = n cal(A)^(n-1) {cal(A),cal(B)},
  $
  $
    {diff_i cal(A)(y),cal(B)} = diff_(y^i) {cal(A)(y), cal(B)},
  $
  $
    {integral d^(N-1)y thin cal(A)(y),cal(B)} = integral d^(N-1)y thin{cal(A)(y), cal(B)}
  $
  ]\ \
To get more familiar with these calculation rules, let us reconsider the free scalar field theory
$
  cal(L) = -eta^(mu nu) diff_mu phi diff_nu phi -1/2 m^2 phi^2,
$
for which we have previously derived the canonical momentum $pi(x) = diff_0 phi(x)$ as well as the Hamiltonian
$
  cal(H)(x) = 1/2 pi^2 + 1/2 delta^(i j) diff_i phi diff_j phi + 1/2 m^2 phi^2. 
$
We now recompute the canonical equations of motion using Poisson brackets, purely algebraically---using the rules above. This leads to
$
  &diff_0 phi (x) = {phi(x), H} = integral d^(N-1)y {phi(x),cal(H)(y)}\
  &= integral d^(N-1)y lr([1/2{phi(x), pi(y)^2} + 1/2 delta^(i j) underbrace({phi(x),diff_i phi(y) diff_j phi(y)},=0) + 1/2 m^2 underbrace({phi(x), phi(y)^2},=0)],size:#80%)\
  &= integral d^(N-1)y thin pi (y) {phi(x),pi(y)} = integral d^(N-1)y thin pi(y) delta^((N-1))(x-y) = pi(x),
$
which is the correct equation of motion for $phi$. Let us repeat the same for $pi$:
$
  &diff_0 pi(x) = {pi(x), H} = integral d^(N-1)y {pi(x), cal(H)(y)}\
  &= integral d^(N-1)y lr([1/2underbrace({pi(x), pi(y)^2},=0) + 1/2 delta^(i j) {pi(x),diff_i phi(y) diff_j phi(y)} + 1/2 m^2 {pi(x), phi(y)^2}],size:#80%)\
  &= integral d^(N-1)y lr([delta^(i j) diff_(y^i) phi(y) underbrace({pi(x),diff_j phi(y)},=-diff_(y^j) delta^((N-1))(x-y)) + m^2 phi(y)underbrace({pi(x),phi(y)},=-delta^((N-1))(x-y))],size:#60%)\
  &= integral d^(N-1)y [-delta^(i j) diff_(y^i) phi (y) diff_(y^j) delta^((N-1))(x-y) - m^2 phi(y) delta^((N-1))(x-y) ]\
  &= integral d^(N-1)y [delta^(i j)diff_(y^i)diff_(y^j) phi(y) - m^2 phi(y)]delta^((N-1))(x-y)\
  &= Delta phi(x) - m^2 phi(x)
$
which also reproduces our result from @canonicalEoMFreeMassiveScalar[equation]. Notice that in these derivations we never had to make explicit use of the Poisson bracket's definition---it was possible to simply use the established rules of calculation as well as the fundamental brackets of the canonical variables. 

In this sense, the Hamiltonian formalism can be viewed as _generated by_ the fundamental relation
$
  {Phi^A (x), Pi_B (y)} = delta^A_B delta^((N-1))  (x-y),
$
together with the dynamical prescription
$
  diff_0 cal(A) = {cal(A), H}
$
for any observable $cal(A)$, and the algebraic calculation rules derived above. 
