#import "template.typ": *
#import "macros.typ": *
#import "@preview/xarrow:0.3.1": xarrow
#import "@preview/fletcher:0.4.5" as fletcher: diagram, node, edge
#import "@preview/cetz:0.5.1": canvas, draw
#import "@preview/cetz-plot:0.1.4": plot


= Numerical Methods
== Finite Difference Stencils
=== First Derivatives
Two-point forward difference stencil with error:
$
  f'(x) &= (f(x+epsilon)-f(x))/epsilon \ & quad- epsilon/2 f''(x) + cal(O)(epsilon²)
$
Three-point central difference stencil with error: 
$
  f'(x) &= (f(x+epsilon) - f(x-epsilon))/(2epsilon) \
  & quad - epsilon^2/6 f^((3))(x) + cal(O)(epsilon^3)
$<eqThreePointCentralDifference>
Five-point central difference stencil with error:
$
  f'(x) &= (-f(x+2epsilon) + 8f(x+epsilon) - 8f(x-epsilon) +f(x-2epsilon))/(12 epsilon)\
  &quad + epsilon^4/30 f^((5))(x) + cal(O)(epsilon^6)
$<eq2ndDerivStencil5pt>
=== Second Derivatives
==== Single-Variable
Three-point stencil with error:
$
  f''(x) &= (f(x+epsilon) - 2f(x) + f(x-epsilon))/epsilon^2\ &quad - epsilon^2/12 f^((4))(x) + cal(O)(epsilon^4)
$<eq2ndDerivStencil3pt>
Five-point stencil with error:
$
  f''(x) &= (-f(x+2epsilon) + 16 f(x+epsilon) - 30 f(x) + 16 f(x-
  epsilon) - f(x-2epsilon))/(12 epsilon^2)\ &quad- epsilon^4/90 f^((6))(x) + cal(O)(epsilon^6)
$
==== Mixed Second Derivative
With a grid spacing of $epsilon_x,epsilon_y>0$ in the $x$- and $y$-directions, respectively, we have
$
  (diff^2 f)/(diff x diff y) = (f_(++) + f_(--) - f_(+-) - f_(-+))/(4 epsilon_x epsilon_y) - 1/6(epsilon_x^2 (diff^4 f)/(diff x^3 diff y) + epsilon_y^2 (diff^4 f)/(diff x diff y^3)) + fO(epsilon^4)
$<eqMixedDerivative>
Here,
$
  f_(pm_1 pm_2) = f(x pm_1 epsilon_x, y pm_2 epsilon_y)
$
Note that @eqMixedDerivative[the expression] corresponds to taking twice the @eqThreePointCentralDifference[central difference] of $f$, once with respect to $x$ and once with respect to $y$.
=== Laplacian
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
$<eqIsotropicLaplacianStencil>
== Runge-Kutta <sectRK>
=== General Structure
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
$<eqRKApprox>
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
One way of determining optimal coefficients $vA,vb,vc$ is to expand the @eqRKApprox[approximation] around $epsilon=0$ and comparing coefficients to the Taylor expansion of $y(t_0+epsilon)$, which using $y'=F(t,y(t))$ as well as the notation
$
  F = F(t_0,y_0), quad F_t = (diff_t F)(t_0,y_0), quad F_y = (diff_y F)(t_0,y_) quad "etc."
$
can be written explicitly as
#bottom-number[$
  y(t_0+epsilon) &= y_0 + epsilon F + epsilon^2/2 (F_t + F F_y)\
  &quad + epsilon^3/6 (F_(t t) + 2 F F_(t y) + F_t F_y + F F_y^2 + F^2 F_(y y))\
  & quad + epsilon^4/24 (F_(t t t) + 3 F F_(t t y) + 3 F_t F_(t y) + 5 F F_y F_(t y) + 3 F^2 F_(t y y) + F_(t t) F_y\ 
  & #h(3.5em)+3 F F_t F_(y y) + F_t F_y^2 + F F_y^3 + 4 F^2 F_y F_(y y) + F^3 F_(y y y))\
  &quad+fO(epsilon^5)
$<eqExplicitExpansionRK>]
As is evident from the number of terms quickly growing with the order in $epsilon$, the strategy of comparing coefficients quickly becomes unmanageable for higher-order RK integrators. Further, we can see that it is highly non-linear in derivatives of $F$; this explains why we need to use the $k_i$ with $i<j$ in the evaluation of $k_j$. Without this chaining, we could not reproduce the nonlinearity in $F$ required to match the explicit @eqExplicitExpansionRK[expression].
=== RK1 (is just Euler)
As a first illustration of how the strategy of reading off coefficients works, let us indulge in the lowest-order case of $s=1$. Although _technically_, $s=0$ would be the lowest-order case, and $y(t) approx y_0$ is _technically_ an approximation of the solution to $y'=F(t,y)$, it is neither a good nor interesting one, hence the decision to start at $s=1$. We start by writing down the $k$'s, of which there is only one---namely
$
  k_1 = F(t_0,y_0) = F.
$
The RK method then prescribes the approximation
$
  y(t_0 + epsilon) approx y_0 + epsilon b_1 k_1 = y_0 + epsilon b_1 F 
$
For this to match @eqExplicitExpansionRK[expansion] up to order $fO(epsilon)$, we need to set $b_1 = 1$. Thus, RK1 is nothing but
$
  y(t_0+epsilon) = y_0 + epsilon F(t_0,y_0) + fO(epsilon^2),
$
which we immediately recognise as the explicit Euler method. 
=== RK2 (gets interesting)
Moving up to the $s=2$ case, we can employ the same strategy to obtain a (hopefully) more precise integration scheme. Here, the $k$'s, along with their expansions in $epsilon$, read
$
  k_1 &= F(t_0,y_0) = F,\
  k_2 &= F(t_0+epsilon c_2, y_0 + epsilon a_(2 1) k_1) &&= F + epsilon c_2 F_t + epsilon a_(2 1) k_1 F_y + fO(epsilon^2)\
  &&&= F+ epsilon (c_2 F_t + a_(2 1) F F_y) + fO(epsilon^2).  
$
Again assembling the RK-approximation, we find
$
  y(t_0 + epsilon) = y_0 + epsilon (b_1 + b_2) F + epsilon^2 (b_2 c_2 F_t + b_2 a_(2 1) F F_y) + fO(epsilon^3) 
$
Clearly, $b_1 + b_2 = 1$ and $b_2 != 0$ must hold for this to have a chance of matching @eqExplicitExpansionRK[expansion] up to $fO(epsilon^2)$. Further, we need $b_2 c_2 = b_2 a_(2 1) = 1/2$, which upon introduction of a parameter $alpha in [1/2,infty)$ allows us to parameterise
$
  b_1 = 1-alpha, quad b_2 = alpha, quad c_2 = 1/(2 alpha), quad a_(2 1) = 1/(2 alpha).
$
Using this parametrisation, the RK2 integration procedure reads
$
  k_1 = F(t_0,y_0), quad k_2 = F(t_0 + epsilon/(2 alpha), y_0 + epsilon/(2 alpha) k_1),
$
and
$ 
  y(t_0 + epsilon) = y_0 + epsilon((1-alpha) k_1 + alpha k_2) + fO(epsilon^3).
$
Two common choices for the free parameter $alpha$ are 
- $alpha = 1/2$ (Heun's method): This leads to the approximation scheme
  $
    y(t_0+epsilon) = y_0 + epsilon/2 (k_1 + k_2) + fO(epsilon^3)
  $
  with
  $
    k_1 = F(t_0,y_0), quad k_2 = F(t_0 + epsilon, y_0 + epsilon k_1).
  $

- $alpha=1$ (Explicit Midpoint Method): In this case, we get
  $
    y(t_0 + epsilon) = y_0 + epsilon k_2 + fO(epsilon^3),
  $
  where
  $  
    k_1 = F(t_0,y_0), quad k_2 = F(t_0 + epsilon/2,y_0 + epsilon/2 k_1).
  $
Although both these cases---and for that matter, any $alpha in [1/2,infty)$ has the same _order_ error, $fO(epsilon^3)$, different values for $alpha$ lead to different forms of the error, which might---depending on the ODE in question---lead to a better constant factor.

=== RK3 (gets ugly)
For $s=3$, the $k$'s and their expansions around $epsilon = 0$ read
$
  k_1 &= F(t_0,y_0) = F,\
  k_2 &= F(t_0 + epsilon c_2, y_0 + epsilon a_(2 1) k_1)\
    &= F + epsilon (c_2 F_t + a_(2 1) F F_y) + epsilon^2/2 (c_2^2 F_(t t) + 2 c_2 a_(2 1) F F_(t y) + a_(2 1)^2 F^2 F_(y y)) + fO(epsilon^3),\
  k_3 &= F(t_0 + epsilon c_3, y_0 + epsilon (a_(3 1) k_1 + a_(3 2) k _2))\
  &= F + epsilon (c_3 F_t + (a_(3 1) + a_(3 2)) F F_y) + epsilon^2/2 (2 a_(3 2) c_2 F_t F_y + 2 a_32 a_21 F F_y^2\
  &quad + c_3^2 F_(t t) + 2 c_3 (a_31 + a_32) F F_(t y) + (a_31 + a_32)^2 F^2 F_(y y)) + fO(epsilon^3)
$
Next, we evaluate the RK3 integration scheme
$
  y(t_0 + epsilon) approx y_0 + epsilon (b_1 k_1 + b_2 k_2 + b_3 k_3)
$
order by order in $epsilon$, comparing to @eqExplicitExpansionRK[the explicit expansion] on the left-hand side:
$
  "at" fO(epsilon): &quad& F &attach(=,t:!) b_1 F + b_2 F + b_3 F
$
which requires
$
  b_1 + b_2 + b_3 = 1 quad <=> quad b_1 = 1-b_2-b_3,
$
for some $b_2,b_3 != 0$. Moving on, we get
$
  "at" fO(epsilon^2): quad epsilon/2 (F_t + F F_y) attach(=,t:!) epsilon^2(b_2 (c_2 F_t + a_21 F F_y) + b_3 (c_3 F_t + (a_31 + a_32) F F_y))
$
From this, we read off the requirements
$
  c_2 = a_21, quad c_3 = a_31 + a_32, quad b_2 c_2 + b_3 c_3 = 1/2.
$
Note that this is not the full solution space; it does however simplify things moving forward. Concretely, it makes the right-hand side simpler; at $fO(epsilon^3)$, we find
$
  &epsilon^3/6 (F_(t t) + 2 F F_(t y) + F_t F_y + F F_y^2 + F^2 F_(y y))\  &attach(=,t:!) epsilon^3/2 [b_2 (c_2^2 F_(t t) + 2 c_2 a_(2 1) F F_(t y) + a_(2 1)^2 F^2 F_(y y)) + b_3 (2 a_(3 2) c_2 F_t F_y + 2 a_32 a_21 F F_y^2\
  &quad + c_3^2 F_(t t) + 2 c_3 (a_31 + a_32) F F_(t y) + (a_31 + a_32)^2 F^2 F_(y y))]\
  &= epsilon^3/2 [b_2 c_2^2 (F_(t t) + 2 F F_(t y) + F^2 F_(y y)) + b_3 (2 a_32 c_2 (F_t F_y + F F_y^2) + c_3^2(F_(t t) + 2 F F_(t y) + F^2 F_(y y)))]
$
This implies
$
  b_2 c_2^2 + b_3 c_3^2 = 1/3 quad "and" quad 2 b_3 a_32 c_2 = 1/3. 
$
At this point, it is most convenient to introduce $c_2 = alpha$ and $c_3 = beta$ as freely choosable constants. In doing so, we get the system of equations
$
  b_1 = 1-b_2-b_3, quad alpha b_2 + beta b_3 = 1/2, quad alpha^2 b_2 + beta^2 b_3 = 1/3, quad a_32 = 1/(6 alpha b_3).
$
This is solved by
$
  b_1 &= -(3 alpha (1 - 2 beta) + 3 beta-2)/(6 alpha beta), &quad&& b_2 &= (2-3beta)/(6 alpha(alpha-beta)),\  b_3 &= -(2-3 alpha)/(6 beta (alpha-beta)),&&quad& a_32 &= -(beta (alpha - beta))/(alpha (2 - 3 alpha)).  
$
The full RK3 integration procedure thus reads
#bottom-number[$
  y(t_0 + epsilon) approx y_0 + epsilon(-(3 alpha (1 - 2 beta) + 3 beta-2)/(6 alpha beta) k_1 + (2-3beta)/(6 alpha(alpha-beta)) k_2 -(2-3 alpha)/(6 beta (alpha-beta)) k_3) + fO(epsilon^4) \ \ \
$]
where
$
  k_1 &= F(t_0,y_0),\
  k_2 &= F(t_0 + epsilon alpha, y_0 + epsilon alpha k_1)\
  k_3 &= F(t_0 + epsilon beta, y_0 + epsilon[(beta + beta(alpha-beta)/alpha(2-3alpha))k_1 - beta(alpha-beta)/alpha(2-3alpha) k_2]).
$
A common choice for the quadrature points is $alpha = 1/2$, $beta = 1$. This leads to
$
  b_1 &= 1/6, &&quad& b_2 &= 2/3, &&quad& b_3 &= 1/6,\
  c_1 &= 0, &&& c_2 &= 1/2, &&& c_3&=1\
  a_21 &= 1/2, &&&  a_31 &= -1 &&& a_32 &= 2.
$
When inserted, this yields the scheme
$
  y(t_0 + epsilon) = y_0 + epsilon(1/6 k_1 + 2/3 k_2 + 1/6 k_3) + fO(epsilon^4)
$
with
$
  k_1 &= F(t_0,y_0),\
  k_2 &= F(t_0 + epsilon/2,y_0 + epsilon/2 k_1),\
  k_3 &= F(t_0 +epsilon, y_0 - epsilon k_1 + 2 epsilon k_2)
$
== Courant-Friedrichs-Lewy (CFL) Conditions for PDEs
=== Example: Parabolic Heat Equation
In this section, we derive the CFL conditions for various discretisations of the parabolic heat equation,
$
  diff_t phi.alt = alpha Delta phi.alt
$<eqHeatEqn>
where $alpha>0$ is the diffusion coefficient and $phi.alt = phi.alt(t,vx)$, $t in RR$, $vx in RR^d$ is the dynamical variable.
==== 7-Point Laplacian Explicit Euler CFL Condition
For this first derivation, we consider the simplest discretisation of @eqHeatEqn, taking the forward difference for the left-hand side, and the Laplacian stencil built from the @eq2ndDerivStencil3pt[second derivative stencil] with equal grid spacing $Delta x$ for all axes. In $d$ dimensions, this amounts to
$
  (phi.alt(t + Delta t,vx) - phi.alt(t,vx))/(Delta t) = alpha sum_(i=1)^d (phi.alt(t,vx + ve_i Delta x ) - 2 phi.alt(t,vx) + phi.alt(t,vx-ve_i Delta x))/(Delta x^2),\
$
or when rearranged for $phi.alt(t+Delta t,vx)$, 
$
  phi.alt(t+Delta t,vx) = phi.alt (t,vx) + C sum_(i=1)^d [phi.alt(t,vx + ve_i Delta x) - 2 phi.alt(t,vx) + phi.alt(t,vx-ve_i Delta x)].
$<eqDiscHeatEqn>
where we defined $C:= (alpha Delta t)/(Delta x^2)$. The goal is to ensure $C$ is such that the discrete linear operator acting on the right-hand side has no exponentially growing modes, i.e. that all its eigenvalues are non-negative. Since any discrete function $psi(vx) = phi.alt(t,vx)$ can be decomposed into functions of the form $e^(i vk dot vx)$, we make the ansatz for a mode
$
  phi.alt(t,vx) = G^(t\/Delta x) e^(i vk dot vx)
$
where $G$ is the growth factor per timestep, fixed by the @eqDiscHeatEqn[discretised equation]. For the evolution to be stable, we must have $|G| <= 1$---let us examine how these modes evolve. Inserting into the left-hand side, we get
$
  phi.alt(t,vx) = G^((t\/Delta x) + 1) e^(i vk dot vx) = G phi.alt(t,vx);
$
hence the name _growth factor_. For the right-hand side, we expand
$
  &phi.alt (t,vx) + C sum_(i=1)^d [phi.alt(t,vx + ve_i Delta x) - 2 phi.alt(t,vx) + phi.alt(t,vx-ve_i Delta x)]\
  &= phi.alt(t,vx) + C G^(t\/Delta x) sum_(i = 1)^d [e^(i vk dot vx)e^(i vk dot ve_i Delta x) - 2 e^(i vk dot vx) + e^(i vk dot vx)e^(-i vk dot ve_i Delta x)]\
  &= phi.alt(t,vx) lr(( 1 + C sum_(i=1) underbrace([e^(i k_i Delta x) - 2 + e^(-i k_i Delta x)],= -4 sin^2(k_i Delta x))),size: #65%)\
  &= phi.alt(t,vx)(1-4C sum_(i=1)^d sin^2 (k_i Delta x)).
$
Equating both sides and dividing by $phi.alt(t,vx)!=0$, we establish
$
  G = 1- 4C sum_(i=1)^d sin^2 (k_i Delta x)
$
The condition that $|G| <= 1$ is satisfied if $G<=1 $ and $G>=-1$. The former is always satisfied since $sin^2 >= 0$, but the latter introduces constraints on $C$---and by that, on the relationship between $Delta t$, $Delta x$ and $alpha$. Concretely, the right-hand side attains a global minimum where $k_i = pi/(2 Delta x)$, $i=1,...,d$. With this value for $vk$, we obtain
$
  -1 attach(<=,t:!) G = 1-4 C d quad <=> quad 1/(2 d) >= C = (alpha Delta t)/(Delta x^2).
$
We can turn this into a condition on the timestep $Delta t$, which leads us to the CFL condition for this discrete stepping operator,
$
  Delta t <= (Delta x^2)/(2 alpha d).
$
==== 7-Point Laplacian Implicit Euler CFL Condition
In the implicit case, one replaces the forward difference in time with the backwards difference, or equivalently, evaluates the right-hand side of the equation at $t + Delta t$ instead of $t$. That is, the heat equation is discretised as
#bottom-number[$
  (phi.alt(t + Delta t,vx) - phi.alt(t,vx))/(Delta t) = alpha sum_(i=1)^d (phi.alt(t + Delta t,vx + ve_i Delta x ) - 2 phi.alt(t+ Delta t,vx) + phi.alt(t+ Delta t,vx-ve_i Delta x))/(Delta x^2).\ \ \
$]
Inserting the ansatz $phi.alt(t,vx) = G^(t\/Delta t) e^(i vk dot vx)$ and identifying terms proceeds largely identically to the explicit case, with the difference that the Laplacian stencil terms now carry an additional factor of $G$ as they are evaluated at $t+Delta t$ instead of $t$. This has as consequence that now,
$
  G = 1- 4 G C sum_(i = 1)^d sin^2 (k_i Delta x\/2),
$
which rearranges to
$
  G = 1/(1+ 4 C sum_(i=1)^d sin^2 (k_i Delta x\/2)).
$
The implications of this result are profound. As long as $C >= 0$, we have $|G| <= 1$---this means that the evolution is unconditionally stable; it allows for an arbitrarily large timestep $Delta t$.

This, however, does not mean that one should opt to just do one's entire simulation in a single step. Although the discrete evolution is _stable_ for any $Delta t > 0$, the error of a step with the given choices of stencils is $fO(Delta t + Delta x^2)$. Increasing $Delta t$ hence also increases the error---just as increasing $Delta x$ in the explicit case (to allow for larger $Delta t$) can destroy the evolution due to a lack of resolution. Because of this, the initially perceived superiority of implicit over explicit time-stepping becomes more subtle; if error requirements pose harsher constraints on $Delta t$ than the explicit CFL condition does, then explicit stepping is preferable over implicit as it is cheaper computationally.

==== 27-point Isotropic Laplacian Explicit Euler CFL Condition
We now move on to derive the CFL condition for the explicit Euler case where on the right-hand side, the @eqIsotropicLaplacianStencil[stencil] is employed. Concretely, we thus consider the stepping scheme
#bottom-number[$
  phi.alt(t+Delta t,vx) = phi.alt(t,vx) + C (-164 dot "(center)" + 30 dot "(faces)" - 2 dot "(edges)" + 1 dot "(corners)"), \ \ \
$<eqA.3.13>]
where $C = (alpha Delta t)/(26 Delta x^2)$. Similar steps to the previous two derivations, together with the fact that the stencil applied to a constant yields 0, leads to
#bottom-number[$
  G(vk) &= 1 -4 C (30 sum_(i=1)^3 sin^2 (k_i Delta x \/2) - 2 sum_(i<j) [sin^2((k_i+k_j) Delta x \/2) + sin^2 ((k_i-k_j) Delta x \/2)]\
  &wide wide quad + sum_(sigma_y,sigma_z in {pm 1}) sin^2 ((k_x + sigma_y k_y + sigma_z k_z) Delta x\/2) )
$]
Denoting the large parentheses by $S(vk)$, we obtain the stability condition
$
  -1 <= 1 - 4C S(vk) quad <=>quad S(vk) <= 1/(2C).
$<eqA.3.15>
We are thus left to find the maximum of $S(vk)$. Introducing the auxiliary variables
$
  xi_i = (k_i Delta x)/2, quad u_i = sin^2(xi_i),
$
we can rewrite the sum over the faces as
$
  sum_(i=1)^3 sin^2(k_i Delta x\/2) = sum_i u_i.
$
For the sum over the edges, we make use of
$
  sin^2(xi_i+xi_j)+sin^2(xi_i-xi_j) &= 2 (sin^2 xi_i cos^2 xi_j + cos^2 xi_i sin^2 xi_j) \
  &= 2(u_i (1-u_j) + (1-u_i)u_j)\
  &= 2(u_i + u_j - 2 u_i u_j)
$
so that
$
  sum_(i<j) [sin^2((k_i+k_j) Delta x \/2) + sin^2 ((k_i-k_j) Delta x \/2)] &= 2 sum_(i < j) [u_i + u_j - 2u_i u_j]\
  &= 4sum_i u_i - 4 sum_(i < j) u_i u_j.
$
For the corner terms, we consider
$
  &sin(xi_x + sigma_y xi_y + sigma_z xi_z) = sin(xi_x + sigma_y xi_y) cos(xi_z) + sigma_z cos(xi_x + sigma_y xi_y) sin(xi_z)\
  &= sin(xi_x) cos(xi_y) cos(xi_z)+ sigma_y cos(xi_x) sin(xi_y) cos(xi_z)\
  &quad + sigma_z cos(xi_x)cos(xi_y) sin(xi_z) - sigma_y sigma_z sin(xi_x) sin(xi_y) sin(xi_z)
$
After squaring, this turns into a horrible mess, but since we are summing over all values of $sigma_y,sigma_z in {pm 1}$, any term linear in one of the $sigma$ will cancel in the sum. That is, only the squares of the individual summands appear in the final result, so that
$
  &#h(-2em)sum_(sigma_y,sigma_z in {pm 1}) sin^2 ((k_x + sigma_y k_y + sigma_z k_z) Delta x\/2)\ &= 4 (sin^2 xi_x cos^2 xi_y cos^2 xi_z + cos^2 xi_x sin^2 xi_y cos^2 xi_z\
  & wide + cos^2 xi_x cos^2 xi_y sin^2 xi_z - sin^2 xi_x sin^2 xi_y sin^2 xi_z)\
  &= 4(u_x (1-u_y)(1-u_z) + (1-u_x)u_y (1-u_z)\ 
  &wide + (1-u_x)(1-u_y)u_z + u_x u_y u_z)\
  &= 4(sum_i u_i - 2 sum_(i < j) u_i u_j + 4 u_x u_y u_z)
$
Hence, $S(vk)$, now as a function of $vu = (u_x,u_y,u_z)$ reads
$
  S(vu) &= 30 sum_i u_i - 2 (4 sum_i u_i - cancelr(4 sum_(i < j) u_i u_j)) + 4(sum_i u_i - cancelr(2 sum_(i<j) u_i u_j) + 4 u_x u_y u_z)\
  &= 26 (u_x + u_y + u_z) + 16 u_x u_y u_z.
$
Since $u_i in [0,1]$, $S$ clearly takes its maximum where $u_x = u_y = u_z = 1$, yielding
$
  max_([0,1]^3) S(vu) = 94.
$
Thus, the @eqA.3.15[condition] is turned into
$
  94 <= 1/(2C) quad <=> quad (alpha Delta t)/(26 Delta x^2) = C <= 1/188
$
Rearranging for $Delta t$ yields the final condition
$
  Delta t <= 13/94 (Delta x^2)/(alpha).
$
Note that 
$
  13 / 94 approx 0.13892...,
$
making the 27-point stencil CFL condition slightly more restrictive than for the 7-point case, where the factor in front of $Delta x^2 \/ alpha$ is
$
  1/(2d) = 1/6 approx 0.166 overline(6). 
$
This trade-off is typically worth it, given that the isotropy of the error is much higher. 
==== 27-point Isotropic Laplacian Implicit Euler CFL Condition
In the implicit case of the preceding section, where the right-hand side of @eqA.3.13 is evaluated at $t+Delta t$ instead, we obtain
$
  G = 1- 4 G C S(vk).
$
Solving for $G$, we find
$
  G = 1/(1+4 C S(vk)).
$
As we have shown before,
$
  S(vk(vu)) = 26 (u_x + u_y + u_z) + 16 u_x u_y u_z in [0,94],
$
so that $|G|<= 1$ is always satisfied and hence, the implicit Euler procedure remains unconditionally stable even with the isotropic stencil for the Laplacian.
=== Example: Hyperbolic Wave Equation
The derivation of CFL conditions for the wave equation,
$
  Box phi.alt = 0 quad <=> quad diff_t^2 phi.alt = c^2 Delta phi.alt,
$
where $c$ is the propagation velocity. The simplest discretisation for this equation (and the only we will consider here) is
#bottom-number[$
  (phi.alt(t+Delta t,vx) - 2 phi.alt(t,vx) + phi.alt(t-Delta t,vx))/(Delta t^2) = c^2 sum_i (phi.alt(t,vx+ve_i Delta x) -2phi.alt(t,vx) + phi.alt(t,vx-ve_i Delta x))/(Delta x^2)\ \ 
$]
This form, as is, is unsuited for a time-stepping procedure, but we can derive stability conditions from it. Again employing the ansatz $phi.alt(t,vx) = G^(t\/Delta t) e^(i vk dot vx)$ for an amplification number $G$, we can rewrite the above as
$
  (G - 2 + G^(-1)) phi.alt(t,vx) &=  C^2  sum_i (e^(i k_i Delta x) - 2 + e^(-i k_i Delta x))phi.alt(t,vx)\
  &= -4C^2 sum_i sin^2 (k_i Delta x\/2) phi.alt(t,vx)\
  &=: -C^2 S(vk) phi.alt(t,vx).
$
Rearranging and multiplying by $G$ leads to
$
  G^2 - 2(1-1/2 C^2 S(vk)) G + 1 = 0
$<eqA.3.34>
The allowed amplification factors are the two roots $G_pm$ of this polynomial in $G$. Since the constant term is $1$, they must satisfy $G_+ G_- = 1$, whence we must have $|G_pm| = 1$ to satisfy the no-growth condition $|G|<=1$. 

For a general polynomial of this form,
$
  G^2 - 2beta G + 1 = 0,quad beta in RR,
$
the roots read
$
  G_pm = beta pm sqrt(beta^2 - 1).
$
If $beta^2 > 1$, we have strictly real and distinct roots, of which at least one has a magnitude greater than 1. Thus, we must have $beta^2 <= 1$, or equivalently, $-1<=beta<=1$. Inserting the $beta$ from @eqA.3.34, this turns into
$
  -1 <= 1-1/2 C^2 S(vk) <= 1.
$
The right-hand bound holds trivially; the left-hand bound requires
$
   C^2 S(vk) <= 4 quad => quad C^2 S_max <= 4
$
Since $S_max = 4d$, this yields the final CFL condition,
$
  (c^2 Delta t^2)/(Delta x^2) = C^2 <= 1/d quad <=> quad Delta t <= (Delta x)/(c sqrt(d)).
$
This has a very nice physical interpretation: the distance travelled by a wave within one timestep, $c Delta t$, must not exceed a value proportional to the grid spacing, $Delta x$. Since the proportionality factor of $1\/sqrt(d)$ is less than 1, this means that the numerical domain of dependence is contained in the physical domain of dependence.
=== Stability Regions of RK Schemes
In this section, we derive and study the stability regions of differen RK schemes applied to a PDE of the form
$
  diff_t phi.alt = F(phi.alt).
$
In a sense, this will be a generalisation of the analysis we did in the preceding sections, where we derived the CFL conditions for the heat and wave equations. The main advantage to the approach in this section is that we will derive stability conditions for the _time integrator_ separately (at least in the linear case), where the only input from the right-hand side are the eigenvalues (or state-dependent growth factors when linearising) of the right-hand side. This separation is useful because frequently in numerical solvers for PDEs, one wants to replace the time integrator or RHS implementation while leaving the other untouched, whence it is more convenient to know the stability conditions of both components individually to then be able to bring them together into a final stability condition.
==== Linear Equations
We first consider the case where $F$ is linear,
$
  diff_t phi.alt = F phi.alt.
$
Assuming $F$ admits a complete set of spatial eigenmodes $psi_n$ with corresponding eigenvalues $lambda_n$---as is typical for spatial discretisations of linear operators---we can expand any solution as
$
  phi.alt(t,vx) = sum_n c_n (t) psi_n (vx).
$
Substituting this expression into the linear PDE yields decoupled ODEs for each mode coefficient:
$
  dot(c)_n (t) = lambda_n c_n (t)
$
Because the modes evolve independently, the stability of the full system reduces to studying the linear test equation
$
  diff_t phi.alt = lambda phi.alt
$<eqRKStabRegTestEqn>
where $lambda in CC$ represents an eigenvalue of the spatial operator $F$. A Runge-Kutta scheme applied to a PDE will be stable only if it produces a non-amplifying update for every eigenvalue $lambda$ in the spectrum of $F$.

Due to the simplicity of this test equation, we can evaluate the full update performed by an RK step explicitly. Let us begin by looking at RK1/Forward Euler. Written in full, RK1 reads
$
  k_1 &= F(phi.alt(t_0)),\
phi.alt(t_0 + Delta t) &= phi.alt (t_0) + Delta t dot k_1.
$
Inserting $k_1$ into the definition, and making use of $F phi.alt = lambda phi.alt$, we get
$
  phi.alt(t_0 + Delta t) = (1 + lambda Delta t) phi.alt(t_0).
$
Introducing the complex variable $z = lambda Delta t$, we can write this as 
$
  phi.alt(t_0 + Delta t) = p_1 (z) phi.alt(t_0), quad p_1 (z) = 1 + z.
$
  The no-growth condition now simply demands that $|p_1(z)|<=1$, which defines a region in the complex plane for $z$ where the RK1 integration scheme is stable. Because of the relationship $z= lambda Delta t$, stability hence depends on the spectrum of $F$ through $lambda$, and on the size of the time step through the factor $Delta t$; a smaller $Delta t$ generally improves stability. We will look at the _stability region_ ${z in CC : |p_1(z)| <= 1}$ later on, when we have derived the polynomials of higher-order schemes to compare it to.

Moving on to RK2, taking e.g. the Heun scheme, the update expands to
$
  k_1 &= F(phi.alt(t_0)) = lambda phi.alt(t_0),\
  k_2 &= F(phi.alt(t_0) + Delta t k_1) = lambda phi.alt(t_0) + Delta t lambda^2 phi.alt(t_0),\
  phi.alt(t_0 + Delta t) &= phi.alt(t_0) + (Delta t)/2 (k_1 + k_2)\
  &= phi.alt(t_0) + (Delta t)/2 (lambda + lambda + Delta t lambda^2) phi.alt(t_0)\
  &= (1+ lambda Delta t + 1/2 (lambda Delta t)^2) phi.alt(t_0).
$
This implies the RK2 stability condition
$
  |p_2 (z)| <= 1, quad p_2(z) := 1 + z + 1/2 z^2.
$
Though already now, one might start to see a pattern, let us briefly also derive $p_3(z)$, the growth factor polynomial for RK3. Expanding the RK3 update yields
$
  k_1 &= F(phi.alt(t_0)) = lambda phi.alt(t_0),\
  k_2 &= F(phi.alt(t_0) + (Delta t)/2 k_1) = lambda phi.alt(t_0) + (Delta t)/2 lambda^2 phi.alt(t_0),\
  &= (lambda + (Delta t)/2 lambda^2) phi.alt(t_0),\
  k_3 &= F(phi.alt(t_0) - Delta t k_1 + 2 Delta t k_2) \
  &= lambda phi.alt(t_0) - Delta t lambda^2 phi.alt(t_0) +2 Delta t lambda^2 phi.alt(t_0) + (Delta t)^2 lambda^3 phi.alt(t_0)\
  &= (lambda + Delta t lambda^2 + (Delta t)^2 lambda^3) phi.alt(t_0),\
  phi.alt(t_0 + Delta t) &= phi.alt(t_0) + (Delta t)/6 (k_1 + 4 k_2 + k_3)\
  &= (1+(Delta t)/6 (lambda + 4(lambda + (Delta t)/2 lambda^2)) + (lambda + Delta t lambda^2 + (Delta t)^2 lambda^3)) phi.alt(t_0)\
  &= (1 +  lambda Delta t + 1/2 (lambda Delta t )^2 + 1/6 (lambda Delta t)^3) phi.alt(t_0).
$
Hence, at third order, we get the stability condition
$
  |p_3(z)|<=1, quad p_3 (z) = 1 + z + 1/2 z^2 + 1/6 z^3.
$
At this point, the pattern is unmistakable: for an $n$-th order RK method with $n$ stages, the stability polynomial $p_n (z)$ is preciely the $n$-th degree Taylor polynomial of the exponential function,
$
  p_n (z) = sum_(k=0)^n z^k/k!.
$
This is not a coincidence. The exact solution to the linear test equation over one time step is 
$
  phi.alt(t_0 + Delta t) = e^(lambda Delta t) phi.alt(t_0) = e^z phi.alt(t_0).
$
By definition, an $n$-th order numerical scheme must match the Taylor expansion of the exact evolution operator $e^z$ up to $fO(z^n)$. Furthermore, an $m$-stage explicit RK scheme performs $m$ nested evaluations of $F$, so its amplification factor $p(z)$ is necessarily a polynomial in $z$ of degree $<=m$.

When $m=n$ (which Butcher's barriers permit for $n<=4$), the polynomial degree bound matches the order-matching condition, which fixes $p_n (z)$ uniquely to the partial sum of the exponential series. For $n>=5$, where maintaining order $n$ requires $m>n$ stages, the polynomial $p(z)$ contains higher-degree terms $fO(z^(n+1))$ that depend on the specific choice of Butcher tableau, though the first $n$ terms always remain $1 + z + ... + z^n\/n!$.

$
#canvas({
  plot.plot(
    size: (10, 10),
    x-domain: (-4, 4),
    y-domain: (-4, 4),
    x-label: $"Re" z$,
    y-label: $"Im" z$,
    x-grid: true,
    y-grid: true,
    y-equal: "x",
    x-tick-step: 1,
    y-tick-step: 1,
    axis-style: "school-book",
    {
      plot.add-contour(
        z: 1,
        ((x, y) => calc.pow(x+1,2) + calc.pow(y,2)),
        x-samples: 50,
        y-samples: 50,
        x-domain: (-4, 4),
        y-domain: (-4, 4),
        fill: false,
        label: "RK1",
        style: (stroke:(thickness: 1.6pt, paint: orange))
      )
      plot.add-contour(
        z: 1,
        ((x, y) => calc.pow((x*x - y*y)/2 + x + 1,2) + calc.pow((x+1)*y,2)),
        x-samples: 50,
        y-samples: 50,
        x-domain: (-3, 3),
        y-domain: (-3, 3),
        fill: false,
        label: "RK2",
        style: (stroke:(thickness: 1.6pt, paint: green))
      )
      plot.add-contour(
        z: 1,
        ((x, y) => calc.pow(1 + x + 1/2 * (x*x - y*y) + 1/6  * (x*x*x - 3 * x * y * y),2) + calc.pow((1+x)* y + 1/6 *(3 * x * x* y - y*y*y),2)),
        x-samples: 50,
        y-samples: 50,
        x-domain: (-3, 3),
        y-domain: (-3, 3),
        fill: false,
        label: "RK3",
        style: (stroke:(thickness: 1.6pt, paint: blue))
      )
      plot.add-contour(
        z: 1,
        ((x, y) => calc.pow((1+x)*y + 1/6 * (3* x*x* y - y*y*y) + 1/24 *(4 * x*x*x* y - 4* x* y*y*y),2) + calc.pow(1 + x + 1/2 * (x*x - y*y) + 1/6 * (x*x*x - 3 * x * y * y) + 1/24 * (x*x*x*x - 6* x*x* y*y + y*y*y*y),2)),
        x-samples: 50,
        y-samples: 50,
        x-domain: (-3, 3),
        y-domain: (-3, 3),
        fill: false,
        label: "RK4",
        style: (stroke:(thickness: 1.6pt, paint: purple))
      )
    }
  )
})
$
In the plot above, the stability regions are the interiors of the lines above. We can see that the higher the order of the method, the larger the region generally becomes. Many important linear PDEs have negative real or purely imaginary eigenvalues, so it makes sense to know the largest possible value of $z$ for which the scheme is still stable. For reference, we provide a table of these values below:
$
#table(
  columns: (auto, auto, auto),
  inset: 10pt,
  align: center,
  stroke: none,
  table.header(
    [*Method*], [*Min. $z$ along $RR_(<=0)$, $z_min^"Re"$*], [*Max. $|z|$ along $i RR$, $z_max^"Im"$*],
    "RK1/Euler", $-2$, $0$,
    "RK2", $-2$,$0$,
    "RK3", $approx -2.5127$, $sqrt(3) approx 1.732$,
    "RK4", $approx -2.7853$, $2sqrt(2)approx 2.8284$
  ),
)
$
There is an important observation to be made about the shape of the stability regions of the higher-order schemes RK3, RK4 and beyond. Recalling that the growth factor of the exact polynomial is $e^z$, it is clear that we can only have stability if $"Re"thin z <=0$. However, both RK3 and RK4 incorporate regions with $"Re"thin z$ in their stability regions, producing stable evolution for modes which physically speaking would be unstable. Though this is often irrelevant as the RHS operators of interest typically do not have unstable eigenmodes, it is worth keeping in mind.

To see how these affect integration schemes for different operators, let us reconsider the heat and wave equations. The main input into the RK stability considerations are the eigenvalues of the right-hand side operator when put into first-order form. The heat equation,
$
  diff_t phi.alt = Delta phi.alt,
$
is already in first-order form; the wave equation needs to be reexpressed slightly. By introducing the additional variable $Pi = diff_t phi.alt$, we can write it as
$
  diff_t vec(phi.alt,Pi) = vec(Pi,Delta phi.alt) = mat(0,1;Delta,0) vec(phi.alt,Pi),
$
so that the RHS operator is the matrix operator $mat(0,1;Delta,0)$. 

We first consider RK stability by only discretising time, i.e. considering $Delta$ a continuum operator. Starting with the heat equation, we recall that on $RR^3$, we can decompose any function as a (continuous) linear combination of the mode functions $phi.alt(vx) = e^(i vk dot vx)$. These are eigenfunctions of the Laplacian, since
$
  Delta phi.alt(vx) = - vk^2 e^(i vk dot vx) = - vk^2 phi.alt(vx).
$
Hence, the eigenvalues of the heat equation's RHS are of the form $lambda_vk^Delta  = -vk^2 in RR_(<=0)$, so that the timestep must satisfy
$
  z_min^"Re" <= lambda_vk^Delta Delta t quad <=> quad Delta t <= z_min^"Re"/lambda_vk^Delta = -z_min^"Re"/vk^2 > 0.
$
Given an initial condition containing mode contributions of arbitrarily large $|vk|$, this requires $Delta t = 0$. If however, we limit our initial conditions to have contributions only up to some fixed maximum $|vk|$, a finite non-zero value of $Delta t$ can provide stability across all of them. In essence, due to $z = lambda Delta t$, we can use $Delta t$ to rescale the relevant portion of the spectrum to lie inside the stability region.

We now move on to the spatial-continuous wave equation, whose eigenvalues are derived as follows:
$
  0=(mat(0,1;Delta,0)-lambda I)vec(phi.alt,Pi) = mat(-lambda,1;Delta,-lambda)vec(phi.alt,Pi) = vec(Pi - lambda phi.alt, Delta phi.alt - lambda Pi).
$
Combining these two linear equations, we obtain
$
  Delta phi.alt = lambda^2 phi.alt, quad Pi = lambda phi.alt,
$
implying that the eigenvalues for the wave equation RHS are the square roots of those of the Laplacian, meaning that
$
  lambda_vk^"wave" = pm sqrt(lambda_vk^Delta) pm i|vk|, quad vk in RR^3,
$
with corresponding eigenfunction
$
  phi.alt(vx) = e^(i vk dot vx), quad Pi(vx) = pm i|vk|e^(i vk dot vx).
$
In particular, this means that the eigenvalues of Since RK1 and RK2 do not contain any nonzero part of the imaginary axis, these schemes are unconditionally unstable. Both RK3 and RK4 contain parts of the imaginary axis, and limit $Delta t$ by
$
  Delta t <= (|z_max^"Im"|)/(|vk|) = cases(sqrt(3)/(|vk|)quad&"for RK3"\,,(2sqrt(2))/(|vk|)quad&"for RK4".)
$
Again, an initial condition containing contributions at arbitrarily large $|vk|$ would require $Delta t = 0$. For frequency-limited initial data, however, finite non-zero values for $Delta t$ exist for which RK$>=$3 is stable.

Recalling the derivations of the CFL conditions, we note that there, the concept of a "maximum allowed wavenumber" $|vk|$ to make it possible for $Delta t$ to be non-zero. On one hand, this has to do with the fact that there, we were working with a discretised grid, which automatically sets an upper limit for the wave number---any mode above the Nyquist frequency aliases onto a mode below it. However, there is another effect at play as well: discretising an operator changes its spectrum.

Hence, we are led to compute the spectrum of discretised versions of the Laplacian, so that we can compare to the continuum case.
+  The most basic discretisation of the Laplacian is
  $
    Delta^"disc"_1 = 1/(Delta x^2) sum_(i=1)^d (E_i^(Delta x) - 2 I + E_i^(-Delta x)),
  $
  where $E_i^a f(vx) = f(vx + a ve_i)$ and $ve_i$ is the unit vector in the $i$-th direction. Clearly, $e^(i vk dot vx)$ are eigenfunctions of $E_i^a$ and $I$, and hence also of $Delta_1^"disc"$. Since any function on $RR^d$ can be written as a (continuous) linear combination of these eigenmodes, they form a complete basis. Concretely, the eigenvalues appear as
  $
    lambda_vk^(Delta,1) e^(i vk dot vx) &= Delta_1^"disc" e^(i vk dot vx) = 1/(Delta x^2) sum_(i = 1)^d (e^(i k_i Delta x) - 2 + e^(-i k_i Delta x))e^(i vk dot vx)\
    &= -4/(Delta x^2)(sum_(i=1)^d sin^2 ((k_i Delta x)/2)) e^(i vk dot vx),
  $
  whence
  $
    lambda_vk^(Delta,1) = -4/(Delta x^2) sum_(i = 1)^d sin^2 ((k_i Delta x)/2)
  $
  For small $k_i Delta x$, we can approximate $sin(k_i Delta x\/2) approx k_i Delta x\/2$, whence
  $
    lambda_vk^(Delta,1) approx -4/(Delta x^2)sum_(i=1) ((k_i Delta x)/2)^2 = - vk^2,
  $
  reproducing the true spectrum. However, the closer the period of the spatial oscillation gets to the grid spacing, the more the spectral points of the continuous and discretised operators diverge. Moreover, unlike $-vk^2$, which grows more and more negative for larger wave vectors, the discretised spectrum is bounded and periodic in the components of $vk$. This is simply due to aliasing of higher onto lower frequency modes; the relevant portion of spatial frequency space is within the Brillouin zone $[-(pi)/(Delta x), (pi)/(Delta x)]^d$. 
  
  The largest magnitude that $lambda_vk^(Delta,1)$ can attain is $(4 d)/(Delta x^2)$. Putting this together with e.g. RK1 reproduces the heat equation's CFL condition 
  $
    Delta t <= min_vk lr(|(z_min^("Re"))/(lambda_vk^(Delta,1))|) = 2/((4d)/ (Delta x^2)) = (Delta x^2)/(2d)
  $<eqGeneralRKStabilityCondition>
  that we have derived before. However, this more general approach of deriving stability conditions---split into determining stability regions of the time integrator and the calculating the eigenvalues of the right-hand side operator---is a lot more versatile and insightful; we could easily replace the time integrator by RK2 or higher, and all that would change is an already computed value of $z_min^"Re"$. 

+ Alternative to changing the time integrator, we can also pick a different discretisation of the right-hand side operator. Let us take a look at a more precise stencil, the Laplacian built from the higher-order stencil @eq2ndDerivStencil5pt:
  $
    Delta_2^"disc" = 1/(12 Delta x^2) sum_(i=1)^d (-E_i^(2 Delta x) + 16 E_i^(Delta x) - 30 I + 16 E_i^(-Delta x) - E_i^(-2 Delta x)),
  $
  Inserting the eigenmodes $e^(i vk dot vx)$, we get
  $
    lambda_vk^(Delta,2) e^(i vk dot vx) = Delta_2^"disc"e^(i vk dot vx) &= 1/(12 Delta x^2)sum_(i=1)^d (-e^(i 2k_i Delta x) + 16 e^(i k_i Delta x) - 30 + 16 e^(-i k_i Delta x) - e^(-i 2 k_i Delta x)) e^(i vk dot vx)\
    &= -1/(3 Delta x^2) sum_(i = 1)^d (16 sin^2 ((k_i Delta x)/2) - sin^2(k_i Delta x)),
  $
  so that
  $
    lambda_vk^(Delta,2) = -1/(3 Delta x^2)sum_(i=1)^d (16 sin^2 ((k_i Delta x)/2)-sin^2(k_i Delta x))
  $
  Again, assuming $k_i Delta x$ to be small, we can approximate
  $
    lambda_vk^(Delta,2) approx -1/(3 Delta x^2) sum_(i=1)^d (16 ((k_i Delta x)/2)^2 - (k_i Delta x)^2) = -1/(3 Delta x^2) sum_(i = 1)^d (k_i Delta x)^2 = - vk^2. wide
  $
  This confirms that we again reproduce the continuum spectrum for small enough $|vk|$, though now with higher accuracy (as we will see in a plot below).

$
#canvas({
  import draw: *

  plot.plot(
    size: (12, 8),
    x-label: $k dot Delta x$,
    y-label: $-lambda_k^Delta dot Delta x^2$,
    x-min: -4,
    x-max: 5,
    y-min: 0,
    y-max: 6.05,
    x-tick-step: calc.pi/2,
    y-tick-step: 1,
    x-grid: true,
    y-grid: true,
    axis-style: "school-book",
    legend: "inner-south-east",
    legend-style: (
      stroke: none,
      fill: none,
      item-spacing: 1em, // Controls vertical gap between entries
    ),
    {
      // Exact eigenvalues
      plot.add(
        k => k*k,
        domain: (-4, 4),
        label: "continuous",
        style: (stroke: (paint: blue, thickness: 1.5pt))
      )

      // ngrow=1 discretised eigenvalues
      plot.add(
        k => 4*calc.pow(calc.sin(k/2),2),
        domain: (-4, 4),
        label: "discretised, width = 3",
        style: (stroke: (paint: red, thickness: 1.2pt, dash: "dashed"))
      )

      // Backward Euler
      plot.add(
        k => 1/3 * (16*calc.pow(calc.sin(k/2),2) - calc.pow(calc.sin(k),2)),
        domain: (-4, 4),
        label: "discretised, width = 5",
        style: (stroke: (paint: green.darken(20%), thickness: 1.2pt, dash: "dotted"))
      )
    }
  )
})
$
In this figure, the one-dimensional comparison between the continuous eigenvalues and the two discretised spectra is shown. We can see that around $k=0$, both discretised spectra match the graph of the continuum spectrum well, but diverge moving away from $k=0$. The divergence is much faster with the lower-order stencil; the 5-wide stencil resolves higher-frequency modes better. This, however, comes at the cost of having a larger maximum value for $lambda_vk^(Delta,2)$ than for $lambda_vk^(Delta,1)$, which reduces the maximum admissible timestep due to @eqGeneralRKStabilityCondition.

==== #text(fill:red)[Non-Linear Equations]
== Boundary Conditions and Grid Stability
=== Sommerfeld Radiation Boundaries
In this section, we derive _Sommerfeld radiation boundary conditions_ for the wave equation on flat Minkowski space, 
$
  diff_t^2 phi.alt = c^2 Delta phi.alt.
$
Although these boundary conditions cannot be applied directly to the more complex BSSN system, it serves as a useful entry point for more advanced radiation boundaries that "absorb" all outgoing radiation, and produce no incoming radiation. 

To distinguish what is "outgoing" and "incoming", we assume our source/region of interest to be located near the origin, and the boundaries to be far enough away from the source that we can reasonably approximate the wavefront of any radiation to be a sphere centered around the origin. Whenever origin-centered spheres come up, a reasonable choice is spherical coordinates $(t,r,theta,phi)$. Since we assume the wavefronts to be spherical for large enough $r$, the field $phi.alt$ becomes independent of $theta$ and $phi$, so
$
  phi.alt = phi.alt(t,r).
$
This turns the wave equation into
$
  diff_t^2 phi.alt = c^2/r^2 diff_r (r^2 diff_r phi.alt).
$
We rewrite this using $phi.alt = u\/r$, yielding
$
  1/r diff_t^2 u = c^2/r^2 diff_r (r diff_r u - u) = c^2/r diff_r^2 u,
$
or equivalently,
$
  diff_t^2 u = c^2 diff_r^2 u,
$
meaning that $u(t,r)$ satisfies a one-dimensional wave equation. Its general solution is a superposition of in- and outgoing waves,
$
  u(t,r) = f_+(r - c t) + f_-(r + c t),
$
for two functions $f_pm:RR->RR$. Factoring the one-dimensional wave operator as 
$
  Box = diff_t^2 - c^2 diff_r^2 = (diff_t - c diff_r)(diff_t + c diff_r) =: X_- X_+,  
$
we can see that $X_+ = diff_t + c diff_r$ annihilates the outgoing wave $f_+$, while $X_- = diff_t - c diff_r$ annihilates the incoming wave $f_-$.

Let us consider how $X_+$ acts on the incoming wave $f_-$. We get
$
  X_+ f_- (r+ c t) = c f'(r + c t) + c f' (r+c t) =2 c f'(r+c t),
$
whence
$
  X_+ u = underbrace(X_+ f_+,=0) + X_+ f_- = 2 c f'_-(r+c t).
$
This means that if we require $X_+ u = 0$, the only "incoming" component we have is an irrelevant constant, $f_- = const$. As we have just seen, the wave operator factors into $Box = X_- X_+$, meaning that if $X_+ u = 0$, then also $X_- X_+ u = 0$. Thus, requiring $X_+ u = 0$ achieves our two goals; $u$ is both a solution to the wave equation _and_ it has no incoming components. 

We can hence, at our boundary away from the origin, require $u$ to satisfy as boundary condition the equation
$
  X_+ u = 0 quad <=> quad diff_t u + c diff_r u = 0. 
$
By construction, this is compatible with the wave equation, and further ensures that no incoming waves enter the domain through the boundary---exactly what we want. We are left to translate this into a boundary condition for our original field $phi.alt = u\/r$. This is done simply by inserting $u = r phi.alt$ into the above, yielding
$
  X_+ (r phi.alt) = 0 quad <=>& quad& diff_t (r phi.alt) + c diff_r (r phi.alt) &= 0\
  <=> && r diff_t phi.alt + c r diff_r phi.alt + c phi.alt &= 0.
$
After dividing both sides by $r$, we finally arrive at the _Sommerfeld radiation boundary condition_
$
  diff_t phi.alt + c diff_r phi.alt + c/r phi.alt = 0. 
$<eqSommerfeldRadBC>
Having now derived this boundary condition, let us consider some more practical aspects of how to implement it in a numerical simulation. In most numerical simulations, one breaks down the second-order in time wave equation into two coupled first-order in time equations for $phi.alt$ and its momentum $pi = diff_t phi.alt$, which explicitly read
$
  diff_t phi.alt &= pi,\
  diff_t pi &= c^2 Delta phi.alt.
$
Both $phi.alt$ and $pi$ are dynamical variables of the problem. We note that only $phi.alt$ appears with spatial derivatives, so that the boundary conditions are only used for filling its exterior boundary ghost cells. By using the definition of $pi$, we may turn the @eqSommerfeldRadBC[Sommerfeld radiation boundary condition] into
$
  pi + c diff_r phi.alt + c/r phi.alt = 0 quad <=> quad diff_r phi.alt = -pi/c - phi.alt/r.
$<eqB.6.13>
This specifies a first derviative of $phi.alt$ which---in typical domains---can be used to isolate the normal derivative of $phi.alt$ needed to fill exterior boundary ghost cells.

If one is running the simulation with a spherical boundary and works in spherial coordinates, then $diff_r$ is the normal direction to the boundary, and there is nothing left to work out. However, most simulations are carried out on a Cartesian grid, where the boundary is a box. This has two consequences: computing $diff_r$ is not as straightforward, and the normal derivative of a boundary is one of the coordinate derviatives $diff_x$, $diff_y$ or $diff_z$ which we need to solve for. 

Luckily, these issues are remedied rather easily. A quick calculation reveals that
$
  diff_r = x/r diff_x + y/r diff_y + z/r diff_z 
$
where of course, $r(x,y,z) = sqrt(x^2 + y^2 + z^2)$. Without loss of generality, let us consider a face of the box-shaped domain's exterior boundary whose normal is $pm diff_x$---that is, a portion of the boundary parallel to the $y z$-plane. To assert its boundary conditions, we need to solve for the normal derivative $diff_x phi.alt$. Inserting the expanded expression for $diff_r$ into @eqB.6.13 allows us to do this as
$
  &&x/r diff_x phi.alt + y/r diff_y phi.alt + z/r diff_z phi.alt = -pi/c -phi.alt/r \ \ 
  <==>  &wide& diff_x phi.alt = -y/x diff_y phi.alt - z/x diff_z phi.alt - r/(c x) pi - 1/x phi.alt.
$
In this form, the boundary condition can be implemented directly. However, it is worth noting that when replacing the boundary-tangent derivatives $diff_y$ and $diff_z$ with finite difference stencils, some numerical error is introduced, which causes the boundary conditions to become satisfied less precisely the larger the angle between $diff_r$ and $diff_x$.

As an additional remark, we should note that the above form assumes $x!= 0$ at the boundary (as well as $y!= 0$ and $z!= 0$ at the respective other boundaries). In practice, this is no restriction, since we have built up this entire discussion on the assumption that the boundaries are far enough from the origin that outgoing waves can be approximated as being emanated from a point source at the origin.

Lastly, we should consider what the boundary conditions mean for $pi$. Starting again from @eqB.6.13, taking a time derivative and applying the PDE $diff_t pi = c^2 Delta phi.alt$ leads us to
$
  diff_r pi = diff_t diff_r phi.alt = - (diff_t pi)/c - (diff_t phi.alt)/r = -c Delta phi.alt - pi/r.
$
Thus, the full set of Sommerfeld boundary conditions reads
$
  diff_r phi.alt &= -pi/c -phi.alt/r,\
  diff_r pi &= -c Delta phi.alt - pi/r.
$<eqCompleteSommerfeldRadial>
Written in terms of cartesian coordinates, assuming that $diff_x$ is the face normal, we have
$
  diff_x phi.alt &= -y/x diff_y phi.alt - z/x diff_z phi.alt - r/(c x)pi - 1/x phi.alt,\
  diff_x pi &= -y/x diff_y pi - z/x diff_z pi - r/x Delta phi.alt - 1/x pi.
$<eqCompleteSommerfeldCartesian>
In numerical implementations on Cartesian grids, it is sometimes useful to simplify this by making the assumption that the waves are incident normal to the face, rather than along $diff_r$. This simplification amounts to setting $r=x$, and neglecting all transverse derivatives, leading to
$
  diff_x phi.alt &= -1/c pi - 1/x phi.alt,\
  diff_x pi &= -c Delta phi.alt - 1/x pi.
$<eqSimplifiedSommerfeld>
At $y=z=0$, this still matches the original Sommerfeld boundary conditions exactly, and remains a valid approximation near these face centers. the more oblique the angle however, the more of the radial wave gets reflected. In some cases, however, this might be an acceptable trade-off compared to the simplification of the boundary conditions one gets. 

To see this, let us first consider how we would discretise the simplified @eqSimplifiedSommerfeld[boundary condition formulation] to $fO(Delta x^2)$ for a cell-centered grid with a halo region of width 1, where the ghost cell index is $-1$, the first interior cell is at $0$, and the boundary is located at $-1/2$. The left-hand sides are straightforwardly discretised to $fO(Delta x^2)$, since the finite difference between $-1$ and $0$ is naturally located at $-1/2$ to order $fO(Delta x^2)$:
$
  diff_x phi.alt|_(-1/2) &= (phi.alt_(-1)-phi.alt_0)/(Delta x) + fO(Delta x^2),\
  diff_x pi|_(-1/2) &= (pi_(-1)-pi_0)/(Delta x) + fO(Delta x^2).
$
On the right-hand sides, we have an issue: we only know $phi.alt|_0$ and $pi|_0$, which approximate $phi.alt|_(-1/2)$ and $pi|_(-1/2)$ to order $fO(Delta x)$, which is too low. Higher accuracy is achieved by using interpolation,
$
  phi.alt_(-1/2) &= (phi.alt_(-1) + phi.alt_0)/(2) + fO(Delta x^2),\
  pi_(-1/2) &= (pi_(-1) + pi_0)/(2) + fO(Delta x^2).
$
This turns the discretised boundary conditions into an implicit (but linear) system in $(phi.alt_(-1), pi_(-1))$, which we will have to solve later.

One term remains that we have not yet adressed---the Laplacian. Since it consists of second derivatives, it is much harder to colocate on the boundary at $-1/2$---luckily, we do not have to. This is because, when multiplying through by the $Delta x$ factor from the first derivative of the left-hand side, it obtains an additional factor of $Delta x$. This increases the error order by one, so it is sufficient to have
$
  Delta phi.alt|_(-1/2) = Delta phi.alt|_0 + fO(Delta x).
$
The $Delta phi.alt|_0$ term, we separate into a normal and a transverse part,
$
  Delta phi.alt|_0 = (phi.alt_(-1) - 2 phi.alt_0 + phi.alt_1)/(Delta x^2) + Delta_perp phi.alt|_0 + fO(Delta x^2),
$
with $Delta_perp phi.alt = diff_y^2 phi.alt + diff_z^2 phi.alt$ the transverse Laplacian that we can evaluate entirely from interior cells, to order $fO(Delta x^2)$. This makes it so that
$
  Delta phi.alt|_(-1/2) = (phi.alt_(-1) - 2 phi.alt_0 + phi.alt_1)/(Delta x^2) + Delta_perp phi.alt|_0 + fO(Delta x).
$
Let us now insert everything back into the continuous equation, keeping track of error orders:
$
  (phi.alt_(-1)-phi.alt_0)/(Delta x) + fO(Delta x^2) &= -(pi_(-1)+pi_0)/(2c) - (phi.alt_(-1) + phi.alt_0)/(2x) + fO(Delta x^2)\
  (pi_(-1)-pi_0)/(Delta x) + fO(Delta x^2) &= -c (phi.alt_(-1) - 2 phi.alt_0 + phi.alt_1)/(Delta x^2) - c Delta_perp phi.alt|_0 - (pi_(-1)+pi_0)/(2x)+ fO(Delta x).
$
Multiplying through by $Delta x$ and collecting error terms, we get
$
  phi.alt_(-1)-phi.alt_0 &= -(Delta x)/(2c)(pi_(-1)+pi_0) - (Delta x)/(2x)(phi.alt_(-1) + phi.alt_0) + fO(Delta x^3),\
  pi_(-1)-pi_0 &= -c/(Delta x) (phi.alt_(-1) - 2 phi.alt_0 + phi.alt_1) - c Delta x (Delta_perp phi.alt|_0) - (Delta x)/(2x)(pi_(-1)+pi_0)+ fO(Delta x^2).
$
As stated previously, this is a linear system in the unknowns $(phi.alt_(-1),pi_(-1))$. To solve it, we first have to separate evaluations at $-1$ from those at $0$ and $1$, yielding
#top-number[$
  (1+(Delta x)/(2x))phi.alt_(-1) + (Delta x)/(2c) pi_(-1) &= (1-(Delta x)/(2x)) phi.alt_0 - (Delta x)/(2c) pi_0 + fO(Delta x^2),\
  c/(Delta x) phi.alt_(-1) + (1+(Delta x)/(2x))pi_(-1) &= (1- (Delta x)/(2x))pi_0 - c Delta x (Delta_perp phi.alt|_0) + c/(Delta x) (2 phi.alt_0 - phi.alt_1) + fO(Delta x^2).
$<eqSimplifiedSommerfeldBCBoundaryUpdate>]
Summarising the right-hand sides into the vector $vb = (b_phi.alt,b_pi)$ with
$
  b_phi.alt &= (1-(Delta x)/(2x)) phi.alt_0 - (Delta x)/(2c) pi_0,\
  b_pi &= (1- (Delta x)/(2x))pi_0 - c Delta x (Delta_perp phi.alt|_0) + (2 c)/(Delta x) (phi.alt_0 - 1/2 phi.alt_1),
$
as well as introducing the matrix
$
  vA = mat(alpha, beta; beta^(-1), alpha), quad alpha = 1+(Delta x)/(2x), quad beta = (Delta x)/(2c)
$
we can write the @eqSimplifiedSommerfeldBCBoundaryUpdate[system] as
$
  vA vx = vb + fO(Delta x^2) quad <=>quad vx = vA^(-1) vb + fO(Delta x^2),
$
with $vx = (phi.alt_(-1),pi_(-1))$ and
$
  vA^(-1) = 1/(alpha^2 - 1) mat(alpha,-beta;-beta^(-1),alpha).
$
We can now also go back to the @eqCompleteSommerfeldCartesian[full Sommerfeld boundary conditions]. Reintroducing the transverse derivatives and unsetting $r=x$, the matrix $vA$ turns into
$
  vA = mat(alpha, r/x beta; r/x beta^(-1), alpha),
$
and 
$
  vb &-> vb - Delta x (y/x diff_y vb|_0 + z/x diff_z vb|_0),\
$
=== Kreiss-Oliger Dissipation
When evolving non-linear PDEs, the non-linear terms can introduce higher-frequency components through the mechanism we discussed in @remarkSpectralMethods. Although there, we were considering spectral methods---where the mechanism presents itself most clearly---it is irrespective of the field representation used; in particular, it is also present when working with discretised field values. In that specific case, the high-frequency modes produced by non-linear terms may exceed the spatial grid cutoff set by the Nyquist limit. Such modes alias back noto lower-frequency modes, introducing unphysical "energy" that can cause a simulation to become unstable and diverge. 

Clearly, one should---as a first step---choose the resolution of the grid fine enough so that all physical modes one expects to be present can be resolved, i.e. do not fall below the grid's Nyquist wavelength of $2 Delta x$. However, numerical and discretisation error can introduce unwanted, unphysical high-frequency modes. Non-linearities then transform these contributions beyond the grid's spatial cutoff, which alias into lower frequencies and can cause the simulation to diverge. To resolve this, we need to dampen or _dissipate_ unphysical high-frequency modes, while leaving the physical modes unchanged, and without harming the integration accuracy of the solver. 

So, let us analyse how to do this. We assume that we have a $1+1$-dimensional PDE that is in a first-order in time formulation,
$
  diff_t phi.alt = F(t,phi.alt).
$
We want to add a term that dissipates unphysical high-frequency modes while leaving low-frequency modes untouched. A first guess might be to add a diffusion term, turning it into
$
  diff_t phi.alt = F(t,phi.alt) + diff_x^2 phi.alt.
$
We can analyse what this does to different frequency components by switching to a spatial frequency representation, $phi.alt(t,x)->tilde(phi.alt)(t,k)$, where the PDE turns into
$
  diff_t tilde(phi.alt) = tilde(F)(t,tilde(phi.alt)) - k^2 tilde(phi.alt).
$
For modes where $k^2 tilde(phi.alt) >> tilde(F)(t,tilde(phi.alt))$, the diffusion term dominates, and the mode behaves as
$
  tilde(phi.alt) sim e^(-k^2 t).
$
That is, the mode is suppressed exponentially, with falloff $k^2$; the higher the frequency, the stronger the dissipation. However, this modification has two issues:

+ Even though the dissipation is weaker for low spatial frequencies $k$, it still modifies the behaviour of all modes; we would like to differentiate between low- and high-frequency modes more strongly. 

+ Since the dissipation term is not present in the physical PDE, it essentially acts as an error term. If introduced as above, without appropriate normalisation, it is an $fO(1)$ error term---completely dwarfing any numerical or discretisation errors, and causing the numerical solution to deviate strongly from the physical continuum solution.

Luckily, these issues are straightforward to resolve. To address the first, we can increase the power of $k$ in the exponential decay factor. Moving from $k^2$ to $k^(2 r)$ gives us a "knob" to adjust in the form of the integer $r$, which turns the falloff function from a Gaussian in the frequency domain into a steeper step-like threshold around $k=0$, leaving modes close to $0$ virtually unaffected while rapidly decaying those at larger $|k|$.

Since obtaining a $-k^(2 r)$ factor in the frequency domain requires an operator proportional to $diff_x^(2 r)$ in the spatial domain (recalling that $diff_x -> -i k$ and thus $diff_x^2 -> -k^2$), we must track the sign,
$
  diff_x^(2 r) -> (-1)^(r)k^(2 r).
$
Thus, the continuum dissipation term takes the form
$
  diff_t phi.alt = F(t,phi.alt) + (-1)^(r+1)diff_x^(2 r) phi.alt.
$
To resolve issue (ii) and prevent the dissipation from degrading the accuract of a $(2r-2)$-th order spatial discretisation scheme, we multiply the operator by $sigma Delta x^(2r-1)$, where $sigma > 0$ is a dimensionless, tunable parameter. This yields the modified continuous PDE
$
  diff_t phi.alt = F(t,phi.alt) + (-1)^(r+1) sigma Delta x^(2r - 1) diff_x^(2 r) phi.alt.
$
For resolved physical modes, where $k<<k_"Nyq"$, the damping term is of order $fO(Delta x^(2r-1))$, acting purely as a high-order truncation error that vanishes rapidly in the continuum limit $Delta x -> 0$. However, near the grid cutoff, where $k approx k_"Nyq" = pi\/Delta x$, the $k^(2 r)$ factor yields a dissipation rate scaling as $fO(1/Delta x)$, which rapidly suppresses unphysical grid-scale modes.

The final step is to discretise this operator on the spatial grid. Because the dissipation operator is already multiplied explicitly by $Delta x^(2r-1)$, any discretisation error introduced by the operator itself is pushed to even higher powers of $Delta x$. We are therefore free to use the simplest centered finite-difference stencil. Defining the forward and backward difference operators as
$
  D_+ f(x) = (f(x+Delta x) - f(x))/(Delta x), quad D_- f(x) = (f(x)-f(x-Delta x))/(Delta x),
$
we discretise the evolution equations as
$
  diff_t phi.alt = F(t,phi.alt) + (-1)^(r+1) Delta x^(2r-1) (D_+ D_-)^r phi.alt.
$
The term 
$
  fD_"KO" phi.alt = (-1)^(r+1) sigma Delta x^(2r-1) (D_+D_-)^r phi.alt
$
is the _Kreiss-Oliger dissipation operator of order $2r-1$_. For reasonable choices of $sigma << 1$, the $fD_"KO"$ term does not alter the CFL condition limiting the timestep for the dissipation-free PDE.

To get an explicit expression for the stencil $(D_+ D_-)^r$, it is useful to introduce the shift operator
$
  E^k f(x) := f(x + k Delta x), quad k in RR. 
$
We can write
$
  D_+ D_- f(x) = f(x+Delta x) - 2 f(x) + f(x-Delta x) = (E^(1\/2)-E^(-1\/2))^2 f(x),
$
so that according to the binomial theorem,
$
  (D_+ D_-)^r &= (E^(1\/2) - E^(-1\/2))^(2 r) = sum_(k=0)^(2r) binom(2r,k) (E^(1\/2))^(k) (-E^(-1\/2))^(2r-k)\
  &= sum_(k=0)^(2r)(-1)^(2r-k) binom(2r,k)E^(k-r) = sum_(k=-r)^r (-1)^(r-k) binom(2r,r+k)E^k
$
Thus,
$
  fD_"KO" phi.alt(x) = - sigma Delta x^(2r-1)sum_(k=-r)^r (-1)^k binom(2r,r+k) phi.alt(x + k Delta x).
$

== Implicit ODE and PDE Solvers
In this section, we consider numerical integration schemes for differential equations of the form
$
  diff_t phi.alt = F(t,phi.alt),
$<eqGeneralDE2>
where $F$ is either a function (for ODEs) or an operator/functional (for PDEs) acting on $phi.alt$. 

Formally, integrating @eqGeneralDE2 over a time interval $[t, t + Delta t]$ yields the exact integral equation
$
  phi.alt(t + Delta t) = phi.alt(t) + integral_(t)^(t+Delta t) F(t',phi.alt(t')) dt'.
$<eqExactStep>
Numerical time-stepping algorithms replace this integral with a discrete quadrature rule. Let $phi.alt^n approx phi.alt(n Delta t)$ denote the discrete numerical solution at step $n$. A general single-step time-integration scheme approximates the integral via a function $G$,
$
  phi.alt^(n+1) approx phi.alt^n + Delta t G(t_n, Delta t, phi.alt^n, phi.alt^(n+1)),
$
up to some local truncation error $fO(Delta t^(k+1))$, yielding a scheme of order $k$.

We classify the integration scheme based on how $G$ depends on the state vector:
+ _Explicit_: $G$ depends solely on known state values from the current or previous time steps ($phi.alt^n, phi.alt^(n-1), dots$). The update step for $phi.alt^(n+1)$ is a direct assignment operation.

+ _Implicit_: $G$ depends on the unknown state $phi.alt^(n+1)$ at the new time level (or intermediate stage evaluations requiring $phi.alt^(n+1)$). 

Consequently, for an implicit scheme, the update formula cannot be evaluated directly; instead, $phi.alt^(n+1)$ must be obtained by solving a system of algebraic equations (linear or non-linear) at every time step.

As the name would suggest, explicit/forward Euler, 
$
  phi.alt(t+Delta t) = phi.alt(t) + Delta t F(t, phi.alt(t)) + fO(Delta t)
$
is an example of an explicit integration scheme. Further examples are the Runge-Kutta integrators we considered in @sectRK. Another class of examples, which incorporates past values $phi.alt(t-Delta t), phi.alt(t- 2 Delta t),...$ exist, and are known as _Adams-Bashforth schemes_. 

=== Example: Backwards Euler
In this section, we introduce the most basic of implicit time integration schemes, the _backwards/implicit Euler method_. Although it has a rather large error, it will serve well in introducing a number of concepts, and provide a reference to compare more elaborate methods against.

The arguably simplest way to approximate the integral
$
  integral_(t)^(t+Delta t) F(t',phi.alt(t')) dt'
$
is to take the value of $phi.alt$ at one of the boundary times, $t$, $t+Delta t$, and multiply it by the width $Delta t$ of the integration interval. Since we are interested in an implicit scheme, we take the boundary value at $t+Delta t$, that is,
$
  integral_(t)^(t+Delta t) F(t',phi.alt(t')) dt' = Delta t thin F(t+Delta t, phi.alt(t+Delta t)) + fO(Delta t^2)
$
Neglecting to explicitly write the error term, this leads to the time-stepping scheme
$
  phi.alt(t+Delta t ) = phi.alt(t) + Delta t thin F(t+Delta t, phi.alt(t+Delta t)).
$
Clearly, this is an equation to be solved for $phi.alt(t+Delta t)$. In the case of $F$ being a function, it is purely algebraic; if $F$ however depends on derivatives of $phi.alt$, it is a differential equation. 

Let us go over some examples below.

+ *Exponential Decay:* For $phi.alt:RR->RR$, we consider the ODE
  $
    diff_t phi.alt(t) = -phi.alt(t).
  $ 
  Clearly, the solution space is given by $phi.alt(t) = phi.alt_0 e^(-t)$, $phi.alt_0 in RR$. The implicit time-stepping scheme prescribes the update equation
  $
    phi.alt(t+Delta t) = phi.alt(t) + Delta t thin (-phi.alt (t+Delta t)).
  $
  In this case---and many other ODE cases---it is possible to explicitly solve this for $phi.alt(t+Delta t)$:
  $
    phi.alt(t+Delta t) = phi.alt(t)/(1+Delta t).
  $
  This is a recursion, which can be made closed-form as
  $
    phi.alt(n Delta t) = phi.alt_0/(1+Delta t)^n
  $<eqImplicitEulerSolution>
  As expected from the analytical solution $phi.alt(t) = phi.alt_0 e^(-t)$, this is an exponential decay; however, each step truncates the exponential series to linear order. One timestep of the analytical solution scales it by a factor 
  $
    e^(-Delta t) = 1/e^(Delta t) = 1/(1+Delta t + 1/2 Delta t^2 + ...).
  $
  The implicit Euler scheme hence neglects the $fO(Delta t^2)$ terms in the denominator. 
  
  Lastly, we note that the stepping scheme is unconditionally stable---for any positive $Delta t$, @eqImplicitEulerSolution defines a stable exponential decay. The explicit/forward Euler scheme, which prescribes
  $
    phi.alt(t+Delta t) = phi.alt(t) + Delta t (-phi.alt(t)) = (1-Delta t) phi.alt(t),
  $
  has the closed-form solution
  $
    phi.alt(n Delta t) = (1-Delta t)^n phi.alt_0,
  $
  which for $Delta t > 1$ becomes oscillatory, and for $Delta t > 2$ starts growing exponentially in magnitude. Further, it is the truncation of the implicit solution:
  $
    phi.alt_0/(1+Delta t)^n = (1-Delta t + Delta t^2 + ...)^n phi.alt_0  approx  (1-Delta t)^n phi.alt_0 
  $
  What we can take away from this example is the following: implicit Euler handles exponential decay much more gracefully than explicit Euler. 

+ *1d Heat Equation, Discretised*:
  We now proceed to the one-dimensional heat equation,
  $
    diff_t phi.alt = diff_x^2 phi.alt.
  $
  Applying the implicit Euler method to this PDE leads to the time-step
  $
    phi.alt(t+Delta t,x) = phi.alt(t,x) + Delta t (diff_x^2 phi.alt)(t + Delta t,x).
  $
  Here, separating quantities evaluated at $t+Delta t$ to the left-hand side, we obtain the equation
  $
    (phi.alt - Delta t diff_x^2 phi.alt) (t+Delta t,x) = phi.alt(t,x).
  $
  This is a second-order linear ODE with constant coefficients, where the previous field state $phi.alt(t,x)$ acts as a source. Although technically, depending on the initial field configuration $phi.alt_0 (x)$, we can sometimes integrate this equation for each step, let us entertain a different (and more versatile) approach.

  Concretely, we discretise $phi.alt$ not just in time but also in space, with a grid spacing of $Delta x$. Since otherwise, notation will get heavy, let us introduce the shorthand
  $
    phi.alt^n_i = phi.alt(n Delta t, i Delta x).
  $
  Employing the @eq2ndDerivStencil3pt[stencil] for the second derivative, the discretised step equation reads
  $
    phi.alt_i^(n+1) - C (phi.alt_(i-1)^(n+1) - 2 phi.alt_i^(n+1) + phi.alt_(i-1)^(n+1)) = phi.alt^n_i
  $
  where $C = Delta t\/Delta x^2$. This is a linear system of equations---let us write it in matrix-vector form. Denoting $bold(phi.alt)^n = (phi.alt_i^n)#h(0em)_(i=1)^N$ for some grid size $N$, and additionally imposing Dirichlet boundary conditions, $phi.alt_(-1) = phi.alt_(N+1) = 0$, we can turn the step equation into
  $
    vM bold(phi.alt)^(n+1) = bold(phi.alt)^n
  $<eqLinearStepEqn>
  with
  $
    vM = mat(
      1 - 2C,     -C,       ,         ,       ;
          -C, 1 - 2C,     -C,         ,       ;
            ,     -C, 1 - 2C,   dots.down,    ;
            ,       , dots.down, dots.down,    -C;
            ,       ,       ,     -C,   1 - 2C;
      gap: #0.8em
    )
  $
  The first instinct might be to make the step equation explicit by writing
  $
    bold(phi.alt)^(n+1) = vM^(-1) bold(phi.alt)^n,
  $
  and precomputing the inverse $vM^(-1)$ once at simulation start. Although this works, this has a significant computational drawback: the matrix $vM$ is sparse, but $vM^(-1)$ is not. this makes the multiplication between $vM^(-1)$ and the state vector $bold(phi.alt)^n$ computationally expensive, especially for larger grids. Because of this, one usually opts to solve the linear system using other methods, such as an LDU decomposition or Gauss-Seidel.
+ *Heat Equation, Spectral*: 
  We now approach the heat equation in $d$ dimensions,
  $
   diff_t phi.alt = Delta phi.alt, quad Delta = sum_(i=1)^d diff_i^2,
  $
  on the domain $Omega = [0,pi]^d$ with Dirichlet boundary conditions, $phi.alt|_(diff Omega) = 0$. Its implicit Euler step equation is clearly
  $
    phi.alt(t+Delta t,vx) = phi.alt(t,vx) + Delta t (Delta phi.alt) (t+Delta t,vx)
  $
  This rearranges into an elliptic problem at each step,
  $
    (phi.alt-Delta t Delta phi.alt) (t+Delta t,vx) = phi.alt(t,vx).
  $<eqImplicitEulerStepHeatEqn>
  Just as in the $d=1$ case above, we could now discretise the Laplacian operator and proceed solving the linear system using Gauss-Seidel or any other linaer solver. Writing this down explicitly is more tedious, since the index gymnastics are far more involved, and there is not much gained conceptually from doing so. 

  Instead, we entertain another approach of approximating the solution $phi.alt$ and by that the action of $Delta$ on it---by means of a spectral expansion, in the special case where $d=3$. Concretely, we define the functions
  $
    psi_(i j k) (vx) &= sin(i x) sin(j y) sin(k z), quad i,j,k in NN^3.
  $
  Clearly, these functions form an orthogonal basis with respect to the standard $L^2$-inner product on the space of $L^2$-functions on $Omega$ satisfying the boundary conditions. We can hence express $phi.alt(t,vx)$ as a sum with time-variable coefficients $c_(i j k) (t)$,
  $
    phi.alt(t,vx) = sum_(i,j,k in NN) c_(i j k)(t) psi_(i j k)(vx).
  $
  Since computers cannot deal with infinities, we must truncate the sum at some highest-frequency mode $(N,N,N)$, turning it into
  $
    phi.alt(t,vx) = sum_(i, j, k = 1)^N c_(i j k)(t) psi_(i j k)(vx).
  $
  One might argue that this loses information just as discretising $phi.alt$ does. However, the kind of information that is lost is different; While discretising the $phi.alt$ leaves the function values at the grid points exact but makes derivatives introduce errors, truncating the basis expansion at some finite point (typically) allows derivatives to remain exact at the cost of introducing some error in raw function values. Let us examine this. The action of the Laplacian on a basis function $psi_(i j k)$ is
  $
    Delta psi_(i j k) = -(i^2 + j^2 + k^2) psi_(i j k).
  $
  Since the right-hand side involves no modes with larger indices, it can still be resolved fully by a series truncated at the mode $(N,N,N)$. Thus, the Laplacian applied to the truncated basis expansion of $phi.alt(t,vx)$ is exact, reading
  $
    Delta phi.alt = -sum_(i,j,k = 1)^N (i^2 + j^2 + k^2) c_(i j k) psi_(i j k).
  $
  We can insert this expansion into @eqImplicitEulerStepHeatEqn to obtain
  $
    sum_(i,j,k=1)^N (c_(i j k)^(n+1) + Delta t (i^2 + j^2 + k^2) c_(i j k)^(n + 1)) psi_(i j k) = sum_(i,j,k=1)^N c_(i j k)^n psi_(i j k).
  $<eqImplicitEulerStepHeatEqnExpansionInserted>
  Since basis expansions are unique, identification of coefficients implies the update equation
  $
    c_(i j k)^(n+1) = c^n_(i j k)/(1 + Delta t (i^2 + j^2 + k^2)).
  $
  This results in a geometric decay, with higher-frequency modes being suppressed more rapidly---just as expected for the heat equation.

  #remark[The reason we can invert @eqImplicitEulerStepHeatEqnExpansionInserted nicely for the explicit iteration step above is that the operator acting on the basis coefficient vector on the left-hand side is diagonal. In other words, the chosen basis diagonalises the Laplacian. Unfortunately, we cannot always do this; given a differential operator $L$, it is not always possible to find an explicit basis that diagonalises it so that the associated system
  $
    diff_t phi.alt = L phi.alt
  $
  can be time-integrated spectrally as cleanly as above. However, there are certain valuable properties of a basis ${psi_n}$ (with $n$ some generalised index) that typically can be achieved, which we now outline. 
  - The arguably most important advantage that spectral methods have over finite differences is that derivatives can be exact. For this to be the case, however, the truncation subpace $V_N = span{psi_1,...,psi_N}$ must be closed under differentiation. Just as we can adjust the resolution of a discretisation grid, we want to be able to adjust $N$ freely (up to stability and error considerations), which makes closure under differentiation require
    $
      diff psi_n in span{psi_1,...,psi_n} quad forall n in NN.
    $
    Here, $diff$ represents any of the relevant first derivatives.

  - Further, we should choose our basis such that $L$ becomes lower triangular with respect to it on $V_N$. This makes solving the equivalent of @eqImplicitEulerStepHeatEqn more straightforward to solve by employing a single Gauss elimination. Lower triangularity however is not a strict necessity, there are other forms where similarly efficient linear solvers exist

    Polynomial bases such as Legendre, Laguerre or Chebyshev polynomials work well for this application, since differentiation reduces their degree and hence ensures closure. For periodic domains, plane waves ($e^(i vk dot vx)$) are highly versatile, and for infinite domains, one can either map to a finite one with a conformal transformation and use a polynomial basis or employ Hermite functions if features are highly localised around the origin. 

    For linear equations with constant coefficients, the exact spatial differentiation of spectral methods is a clear gain over discrete finite differences. The only errors emerge from the approximation of initial conditions when projecting onto the basis and truncating, as well as the error introduced by the time integrator. However, as soon as terms which are non-linear in the basis functions emerge, we start to introduce additional truncation error. To illustrate this, consider terms of the form
    $
      f phi.alt quad "or" quad phi.alt^2.
    $
    When inserting basis expansions for the field $phi.alt$ and auxiliary functions such as $f$, products $psi_n psi_m$ of basis functions emerge. Although these can be expressed as a linear combination of $psi_n$'s again, this may introduce contributions from modes that are above the truncation cutoff $N$. Concretely, consider for example the real Fourier basis,
    $
      psi_n (x) = cos(n x), quad chi_n (x) = sin(n x).
    $
    The product of two such basis functions yields linear combinations like
    $
      psi_n (x)psi_m (x) = cos(n x) cos(m x) &= 1/2 cos((n-m)x) + 1/2 cos((n+m)x)\
      &= 1/2 psi_(n - m) (x) + 1/2 psi_(n + m)(x),
      \ \
      chi_n (x)psi_m (x) = sin(n x) cos(m x)&= 1/2 sin((n-m)x) + 1/2 sin((n+m)x)\
      &= 1/2 chi_(n - m)(x) + 1/2 chi_(n+m)(x),
      \ \
      chi_n (x)chi_m (x) = sin(n x) sin(m x)&= 1/2 cos((n-m)x) - 1/2 cos ((n+m) x)\
      &= 1/2 psi_(n-m)(x) - 1/2 psi_(n+m)(x).
    $
    This means that the product of two modes will alias as a low-frequency mode $psi_(n-m)$ or $chi_(n-m)$ _as well as_ a high-frequency mode $psi_(n+m)$ or $chi_(n+m)$. If $n+m > N$, the high-frequency alias is truncated, introducing an error. For this error not to be devastating, $N$ has to be chosen large enough so that modes close to it are low in amplitude.]<remarkSpectralMethods>

=== Example: Crank-Nicolson
A better way of approximating an integral than by taking one of its endpoint values multiplied by the interval width is to approximate the integrand as the linear polynomial passing through both endpoints. This is the so-called _trapezoidal_ integration rule,
$
  integral_a^b f(t) dt = (Delta t)/2 (f(b) + f(a)) + fO(Delta t^3),
$
which, in particular, improves the error from $fO(Delta t^2)$ to $fO(Delta t^3)$.

Applying this approximation to the exact @eqExactStep[step equation] leads to 
$
  phi.alt(t+Delta t) = phi.alt(t) + (Delta t)/2 (F(t,phi.alt(t)) + F(t+Delta t, phi.alt(t+Delta t))) + fO(Delta t^3).
$
This still leaves us with an equation to solve at each step, but improves the per-step error from $fO(Delta t^2)$ to $fO(Delta t^3)$.

Even though this is still implicit, the example of exponential decay is no longer unconditionally stable in any useful way. Let us consider this more precisely, starting from the differential equation
$
  diff_t phi.alt(t) = -phi.alt(t). 
$
The step equation (dropping error terms) then becomes
$
  phi.alt(t+Delta t) = phi.alt(t) - (Delta t)/2 phi.alt(t) - (Delta t)/2 phi.alt(t+Delta t),
$
which is easily solved for $phi.alt(t+Delta t)$ as
$
  phi.alt(t+Delta t) = (1-(Delta t)/2)/(1+(Delta t)/2) phi.alt(t)
$
Just as with the forward and backward Euler approaches for this equation, $phi.alt$ is multiplied by an approximation of $e^(-Delta t)$ in each timestep. 
$
#canvas({
  import draw: *

  plot.plot(
    size: (12, 8),
    x-label: $Delta t$,
    y-label: $g(Delta t)$,
    x-min: 0,
    x-max: 3,
    y-min: 0,
    y-max: 1.05,
    x-tick-step: 0.5,
    y-tick-step: 0.5,
    axis-style: "school-book",
    legend: "inner-north-east",
    legend-style: (
      stroke: none,
      fill: none,
      item-spacing: 1em, // Controls vertical gap between entries
    ),
    {
      // Exact exponential decay
      plot.add(
        dt => calc.exp(-dt),
        domain: (0, 3),
        label: $e^(-Delta t)$,
        style: (stroke: (paint: blue, thickness: 1.5pt))
      )

      // Forward Euler
      plot.add(
        dt => 1 - dt,
        domain: (0, 3),
        label: $1 - Delta t$,
        style: (stroke: (paint: red, thickness: 1.2pt, dash: "dashed"))
      )

      // Backward Euler
      plot.add(
        dt => 1 / (1 + dt),
        domain: (0, 3),
        label: $1 \/ (1 + Delta t)$,
        style: (stroke: (paint: green.darken(20%), thickness: 1.2pt, dash: "dotted"))
      )

      // Crank-Nicolson / Padé (1,1)
      plot.add(
        dt => (1 - dt / 2) / (1 + dt / 2),
        domain: (0, 3),
        label: $(1 - (Delta t) / 2) \/ (1 + (Delta t) / 2)$,
        style: (stroke: (paint: purple, thickness: 1.2pt, dash: "dash-dotted"))
      )
    }
  )
})
$
As can be seen in the plot above, the factor 
$
  g(Delta t) = (1-(Delta t)/2)/(1+(Delta t)/2)
$
is much closer to the exact solution's step factor of $e^(-Delta t)$ than the explicit Euler (red) and implicit Euler (green) factors. However, even though $|g(Delta t)|<=1$ for all $Delta t>=0$---implying unconditional stability---the factor becomes negative for $Delta t > 2$, producing unwanted oscillatory behaviour. However, the Crank-Nicolson procedure is still superior to explicit Euler, since the timestep may have a larger value while retaining the same integration accuracy.
=== Example: IMEX Schemes (Implicit-Explicit)
What we have seen so far is that, roughly speaking, implicit Euler is more or even unconditionally stable than explicit Euler---allowing for a timestep size $Delta t$ limited only by truncation error---but explicit Euler steps are much less costly computationally speaking, as they do not require solving algebraic or differential equations at each step. Hence, there is always a tradeoff between the two; stiff systems, like the heat equation
$
  diff_t phi.alt = Delta phi.alt,
$
are better approached using implicit Euler. This works well since one trades CFL-limited timesteps for having to solve a linear system of equations at each step---something for which plenty of efficient algorithms exist. However, if we introduce a nonlinear term, and consider e.g.
$
  diff_t phi.alt = Delta phi.alt - phi.alt^3,
$<eqNonLinExampleIMEX>
the implicit timesteps become much more costly---the $phi.alt^3$ term destroys the linearity of the system of equations. 

It would be nice if we could keep the advantage of large timesteps for the $Delta phi.alt$ smoothing term, while not having to invert the cubic term. It turns out, such methods exist---they fall under the category of _IMEX_ or _Implicit-Explicit_ integrators. Their starting point is a differential equation
$
  diff_t phi.alt = F(t,phi.alt) + G(t,phi.alt),
$<eqIMEXDE>
where on the right-hand side $F$ is a function or functional containing all stiff terms for which we want implicit-like timestepping, and $G(t,phi.alt)$ contains all other terms like nonlinear contributions, for which inversion is too computationally expensive. In the @eqNonLinExampleIMEX[example], we would pick
$
  F(t,phi.alt) = Delta phi.alt, quad G(t,phi.alt) = -phi.alt^3.
$
Clearly, the integral equation associated to @eqIMEXDE can be rearranged into the step equation
$
  phi.alt(t+Delta t) = phi.alt(t) + integral_t^(t+Delta t) F(t',phi.alt(t'))dt'+ integral_t^(t+Delta t) G(t',phi.alt(t'))dt'
$
where we have used the linearity of the integral to split it up into two separate integrations. The advantage of this is that we are free to pick approximation schemes for both of the terms. 

The simplest IMEX scheme consists of using the implicit Euler approximation
$
  integral_t^(t+Delta t) F(t',phi.alt(t'))dt' = F(t + Delta t,phi.alt(t+Delta t)) Delta t + fO(Delta t^2)
$
for the stiff term involving $F$ and the explicit Euler approximation
$
  integral_t^(t+Delta t) G(t',phi.alt(t'))dt' = G(t,phi.alt(t)) Delta t + fO(Delta t^2)
$
for the remainder involving $G$. This produces a step equation (neglecting error) reading
$
  phi.alt(t+Delta t) = phi.alt(t) + Delta t (F(t+Delta t, phi.alt(t+Delta t)) + G(t,phi.alt(t))).
$
Abusing notation slightly, we may rearrange this into
$
  (I - Delta t F(t + Delta t))(phi.alt(t + Delta t)) = phi.alt(t) + Delta t G(t,phi.alt(t)).
$
In this form all explicit evaluations at $t$ remain on the right-hand side, while the unknown value of $phi.alt(t+Delta t)$ is moved to the left-hand side. This leaves us to invert only the (typically but not necessarily linear) operator $I-Delta t F$, alleviating the need to solve non-linear contributions from $G$.

Applying this 1st-order IMEX scheme to our non-linear heat equation example @eqNonLinExampleIMEX with spatial finite differences yields:
$
  vM bold(phi.alt)^(n+1) = bold(phi.alt)^n - Delta t (bold(phi.alt)^n)^3,
$
where $vM$ is the same linear system matrix derived in @eqLinearStepEqn, and the power $(bold(phi.alt)^n)^3$ is evaluated component-wise. The non-linear term is computed purely as an explicit source update on the right-hand side, leaving $vM$ untouched.
== Newton-Raphson for Elliptic Equations
=== Recap of N-R for Root Finding in $d=1$
In this section, we discuss the Newton-Raphson method of finding a root of a function $f:RR->RR$, in preparation of the more general, higher dimensional case where we study maps $RR^m->RR^n$. 

The starting point is an initial guess $x_0 in RR$ for where a root of $f$ might be. To refine this guess, we approximate $f$ to linear order around $x_0$,
$
  f(x) approx f(x_0) + f'(x_0) (x-x_0).
$
This linear approximation is sure to have a root, which we use as the refined guess $x_1$ of the root of $f$ itself. That is, we set
$
  0 = f(x_0) + f'(x_0)(x_1-x_0).
$
Rearranging for $x_1$, we obtain
$
  x_1 = x_0 -f(x_0)/(f'(x_0)).
$
Repeating this step yields the _Newton-Raphson iteration_
$
  x_(n+1) = x_n - f(x_n)/(f'(x_n)).
$<eqNRiteration>
Clearly, roots of $f$ are fixed points of this iteration, since if for $x_*$ such that $f(x_*)=0$,
$
  x_* = x_* - f(x_*)/(f'(x_*)).
$
However, for this iteration to actually converge requires $f'(x_*) != 0$, and that the fixed point be an attractor. That is, there must be an interval $I = (x_*-epsilon,x_*+epsilon) subset.eq RR$ such that the iteration map
$
  g(x) = x-f(x)/(f'(x))
$
is a contraction on $I$. That is, we require $|g'(x_*)|<1$. Evaluating this, we get
$
  1 > lr(|1-underbrace((f'(x_*))/(f'(x_*)),=1) + (overbrace(f(x_*),=0)f''(x_*))/(f'(x_*))^2|,size:#55%) = 0,
$
which is unconditionally satisfied. Hence, such a contraction neighbourhood $I$ of the root $x_*$ exists, and by the Banach fixed point theorem, for any $x_0 in I$, the @eqNRiteration[Newton-Raphson iteration] converges to $x_*$.

Let us examine the nature of this convergence, provided that $x_*$ is a simple root ($f'(x_*) != 0$). To this end, we define the error $e_n$ as
$
  e_n = |x_* - x_n|,
$
the distance of $x_n$ to the true value of the root. We recall that $g(x_*) = x_*$, $g'(x_*) = 0$, so that
$
  x_(n+1) = g(x_n) &= underbrace(g(x_*),=x_*) + underbrace(g'(x_*),=0) (x_n-x_*) + 1/2 g''(x_*)(x_n-x_*)^2 + fO((x_n-x_*)^3) \
  &=x_* + 1/2 g''(x_*)(x_n-x_*)^2 + fO((x_n-x_*)^3).
$
Subtracting $x_*$ from both sides and taking the modulus, we obtain
$
  e_(n+1) = |x_(n+1) - x_*| <= 1/2 |g''(x_*)|e_n^2 + fO(e_n^3)
$
Since $e_n -> 0$ as $n->infty$, this implies
$
  lim_(n -> infty) e_(n+1)/e_n^2 <= 1/2|g''(x_*)| = 1/2 lr(|(f''(x_*))/(f'(x_*))|)
$
showing that the convergence order is quadratic if $f''(x_*)!= 0$, and at least cubic if $f''(x_*) = 0$.

=== Higher-Dimensional Generalisation
We now move to the multi-dimensional setting, where $f:RR^m->RR^n$. To derive a Newton-Raphson iteration for finding roots of $f$, we pick an initial guess $x_0 in RR^m$, and approximate $f$ to linear order around $x_0$, which yields
$
  f(x) approx f(x_0) + J(f)(x_0)(x-x_0),
$
where $[J(f)]_(i j) = diff_j f_i$ is the Jacobian matrix of $f$. As our refined guess $x_1 in RR^m$, we pick a root of $g$, that is, 
$
  J(f)(x_0)(x_1-x_0) = -f(x_0). 
$
This is now a linear system. If $n>m$, this generally has no solutions; if $n<m$, the system is underdetermined an $x_1$ is non-unique---we thus focus on the case where $n=m$. In that case, the system can be solved iff $det J(f)(x_0) != 0$, which is equivalent to the condition that $f'(x_0)!= 0$ in the one-dimensional case. 

These considerations lead us to define the _Newton-Raphson iteration_ for $f:RR^n->RR^n$ as
$
  x_(n+1) = x_n - J(f)(x_n)^(-1) f(x_n).   
$
Though this is nice to work with analytically, inverting the Jacobian explicitly is computationally expensive, whence one typically opts to solve
$
  J(f)(x_n) Delta x = -f(x_n)
$
and to then update $x_(n+1) = x_n + Delta x$.

== Adaptive Mesh Refinement
=== Refinement Conditions
The goal of adaptive mesh refinement (AMR) is to increase the resolution of a simulation wherever there are features in the field configuration which cannot be resolved adequately at the current resolution. Hence, we need a predicate to decide whether a grid cell should be refined or not; for this, we need to be able to detect features.

A first intuition for how to detect features might be something gradient-based---something like
$
  "refine if" quad |nabla phi.alt| > C
$
for some threshold value $C$. Unfortunately, this has multiple issues:

+ It is entirely decoupled from the current grid resolution. If large enough, the same gradient value will always lead to refinement, regardless of what resolution we already have. 

+ Further, the gradient is scale-dependent. A rescaling $phi.alt->lambda phi.alt$ for some real constant $lambda > 0$ does not affect the accuracy of a finite-difference stencil, since it simply gets scaled by the same factor. However, since $|nabla phi.alt| ->lambda|nabla phi.alt|$, the condition depends upon the scale of $phi.alt$.

+ For affine linear functions $phi.alt(vx) = va dot vx + vb$, finite difference approximations of its derivatives are exact. However, since $nabla phi.alt = va$, so we refine if $|va|>C$. This is bad: we do not need to increase the resolution for field configurations where the error of finite difference approximations vanishes already. 

Let us address these issues in order. The first is rather simple to fix: we simply multiply by the grid spacing $Delta x$, and turn the refinement condition into
$
  Delta x|nabla phi.alt| > C.
$
Effectively, this removes the division by the grid spacing in the finite difference approximation. Concretely, in the one-dimensional case, we have
$
  C < Delta x|nabla phi.alt| approx Delta x (phi.alt(x+Delta x)-phi.alt(x-Delta x))/(2 Delta x) = 1/2 (phi.alt(x+Delta x) - phi.alt(x-Delta x)).
$
We hence respect the grid spacing now; we refine if the difference between neighbouring cells' values becomes too large. This is already a much more sensible condition---we no longer refine indefinitely, but rather, until the change from cell to cell is small enough.

To address the second issue, (ii), we might think to rescale the gradient by $|phi.alt|$, turning it into
$
  (Delta x|nabla phi.alt|)/(|phi.alt|) > C.
$
Although this is now invariant under $phi.alt-> lambda phi.alt$, we have introduced two new issues: a potential division by zero, and a sensitivity to shift, $phi.alt->phi.alt + alpha$. The division by zero can be remedied by adding a typical scale/noise floor $phi.alt_0$ of the field $phi.alt$ to the denominator, turning it into
$
  (Delta x|nabla phi.alt|)/(|phi.alt| + phi.alt_0) > C.
$
This is resolution- and scale-invariant, so we have solved (i) and (ii), but (iii) is still unresolved---it got even worse---and we are now sensitive to shifts, $phi.alt -> phi.alt + alpha$. 

So, let us try to address (iii), together with this new issue. The shift-covariance issue is simple enough to deal with; any derivative of $phi.alt$ is shift-invariant, so we only use its derivatives in our condition. We are left to address (iii). Since we would like our condition's left-hand side to ideally evaluate to zero for affine-linear functions, we are bound to make it proportional to _second_ derivatives of $phi.alt$. The likely most infamous scalar second derivative of a function is its Laplacian, $Delta phi.alt$. Being a second derivative, it maps affine-linear functions $va dot vx + vb$ to zero, but unfortunately, it also annihilates the entire class of harmonic functions as well. For example, the function $phi.alt(x,y) = e^x sin y$ is harmonic, but its finite difference approximations of derivatives are not trivially exact. Thus, using $Delta phi.alt$ is not viable. 

The most general second derivative of $phi.alt$ is its Hessian, which we denote by
$
  H = [(diff^2 phi.alt)/(diff x^i diff x^j)]_(i,j=1)^d.
$
We would like to build a scalar from it that is sensitive to any second-order features, and which is isotropic. Under an orthogonal transformation $vx->R vx$, with $R^top R = I$, the Hessian transforms as $H->R H R^top$. Unfortunately, we have already ruled out $tr H = Delta phi.alt$ as an option, so the next-best isotropic scalar is the Frobenius norm
$
  ||H||_F^2 = tr(H^top H).
$
It is indeed isotropic, since
$
  tr(H^top H) -> tr((R H R^top)^top R H R^top) = tr(R H^top R^top R H R^top) = tr(H^top H).
$
In e.g. $d=3$, we can write it out explicitly as
$
  ||H||_F^2 = phi.alt_(x x)^2 + phi.alt_(y y)^2 + phi.alt_(z z)^2 + 2 phi.alt_(x y)^2 + 2 phi.alt_(y z)^2 + 2 phi.alt_(z x)^2.
$
Since all terms are non-negative, any non-zero second-order derivative will be picked up by $||H||_F$. Moreover, due to isotropy, it is agnostic to the feature's orientation---we have thus found an ingredient for a refinement condition solving (iii). 

Taking into account everything we have established in the above, a good refinement condition is
$
  "refine if" quad (Delta x^4||H||_F^2)/(Delta x^2|nabla phi.alt|^2 + phi.alt_0^2) > C.
$
This is a generalisation of the Löhner error estimate. At this point it is worth considering its behaviour in different regimes. In a region of small gradients---where $Delta x^2|nabla phi.alt|^2 << phi.alt_0^2$, such as at the top of a bell curve---the condition collapses to
$
  (Delta x^4)/(phi.alt_0^2)||H||_F^2 > C.
$
This means that refinement is proportional to the curvature $||H||_F^2$. In the other limiting case---$Delta x^2|nabla phi.alt|^2 >> phi.alt_0^2$, such as near shocks or wavefronts---the condition turns into
$
  Delta x^2 (||H||_F^2)/(|nabla phi.alt|^2) > C.
$
This is a measure for the relative rate of change of the gradient, loosely to be interpreted as $nabla log nabla phi.alt$. 

=== #text(fill: red)[Subcycling]

== #text(fill:red)[Discontinuous Galerkin Methods]