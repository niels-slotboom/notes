//tensor notation definition. T = tensor name, sink is indices, +a for ^a, -a for _a.
#let tensor(T, ..sink) = {
  let args = sink.pos()

  let (uppers, lowers) = ((), ())  // array, array
  let hphantom(s) = { hide($#s$) }  // Like Latex's \hphantom

  for i in range(args.len()) {
    let arg = args.at(i)
    let tuple = if arg.has("children") {
      arg.at("children")
    } else {
      ([+], sym.square)
    }
    assert(type(tuple) == array, message: "shall be array")

    let pos = tuple.at(0)
    let symbol = if tuple.len() >= 2 {
      tuple.slice(1).join()
    } else {
      sym.square
    }
    if pos == $+$.body  {
      let rendering = $#symbol$
      let space = $#symbol#h(-0.04cm)$
      uppers.push(rendering)
      lowers.push(hphantom(space))
    } else {  // Curiously, equality with [-] is always false, so we don't do it
      let rendering = $#symbol$
      let space = $#symbol#h(-0.04cm)$
      uppers.push(hphantom(space))
      lowers.push(rendering)
    }
  }
  (T,math.attach(hphantom(sym.zwj), t: uppers.join(), b: lowers.join())).join()
}

//cancel with a red stroke
#let cancelr(T) = math.cancel(T,stroke:(paint:color.red,thickness:0.8pt));
#let comment(T) = text(color.red)[*Comment:* #T];

#let inprod(T,S) = {$lr(angle.l #T , #S angle.r)$}
#let mhighlight(T) = {box([$#T$], fill:rgb("#fdff9c"), radius: 4pt, outset: 0.13em)}
#let rank = math.op("rank")

#let Ric = math.op("Ric")
#let mRic = overline(math.op("Ric"))
#let Div = math.op("div")

//various groups and algebras
#let SO = math.op("SO")
#let so = math.frak("so")

#let SU = math.op("SU")
#let su = math.frak("su")

#let Sp = math.op("Sp")
#let sp = math.frak("sp")

#let SL = math.op("SL")
#let sl = math.frak("sl")

#let GL = math.op("GL")
#let gl = math.frak("gl")

#let bottom-number(eq) = {
  set math.equation(number-align: bottom + right)
  eq
}

//miscellaneous math notation
#let const = math.op("const.")
#let diag = math.op("diag")
#let Box = sym.ballot
#let sim = sym.tilde
#let iff = math.arrow.l.r.double.long
#let ddot(a) = {
  math.accent(a,math.dot.double)
}
#let pm = math.plus.minus
#let mp = math.minus.plus
#let infty = math.infinity
#let span = math.op("span")
#let wedge = math.and
#let otimes = math.times.o
#let oplus = math.plus.o
#let star = (math.star.op,h(1pt) ).join()
#let flat = [\u{266D}]
#let sharp = [\u{266F}]
#let codim = math.op("codim")
#let dim = math.op("dim")
#let mnabla = math.macron(math.nabla)
#let tnabla = math.tilde(math.nabla)

#let diff = math.partial

//lowercase bold vectors
#let va = $bold(upright(a))$
#let vb = $bold(upright(b))$
#let vc = $bold(upright(c))$
#let vd = $bold(upright(d))$
#let ve = $bold(upright(e))$
#let vf = $bold(upright(f))$
#let vg = $bold(upright(g))$
#let vh = $bold(upright(h))$
#let vi = $bold(upright(i))$
#let vj = $bold(upright(j))$
#let vk = $bold(upright(k))$
#let vl = $bold(upright(l))$
#let vm = $bold(upright(m))$
#let vn = $bold(upright(n))$
#let vo = $bold(upright(o))$
#let vp = $bold(upright(p))$
#let vq = $bold(upright(q))$
#let vr = $bold(upright(r))$
#let vs = $bold(upright(s))$
#let vt = $bold(upright(t))$
#let vu = $bold(upright(u))$
#let vv = $bold(upright(v))$
#let vw = $bold(upright(w))$
#let vx = $bold(upright(x))$
#let vy = $bold(upright(y))$
#let vz = $bold(upright(z))$

//uppercase bold vectors
#let vA = $bold(upright(A))$
#let vB = $bold(upright(B))$
#let vC = $bold(upright(C))$
#let vD = $bold(upright(D))$
#let vE = $bold(upright(E))$
#let vF = $bold(upright(F))$
#let vG = $bold(upright(G))$
#let vH = $bold(upright(H))$
#let vI = $bold(upright(I))$
#let vJ = $bold(upright(J))$
#let vK = $bold(upright(K))$
#let vL = $bold(upright(L))$
#let vM = $bold(upright(M))$
#let vN = $bold(upright(N))$
#let vO = $bold(upright(O))$
#let vP = $bold(upright(P))$
#let vQ = $bold(upright(Q))$
#let vR = $bold(upright(R))$
#let vS = $bold(upright(S))$
#let vT = $bold(upright(T))$
#let vU = $bold(upright(U))$
#let vV = $bold(upright(V))$
#let vW = $bold(upright(W))$
#let vX = $bold(upright(X))$
#let vY = $bold(upright(Y))$
#let vZ = $bold(upright(Z))$

//lowercase math cal characters
#let fa = math.cal("a")
#let fb = math.cal("b")
#let fc = math.cal("c")
#let fd = math.cal("d")
#let fe = math.cal("e")
#let ff = math.cal("f")
#let fg = math.cal("g")
#let fh = math.cal("h")
#let fi = math.cal("i")
#let fj = math.cal("j")
#let fk = math.cal("k")
#let fl = math.cal("l")
#let fm = math.cal("m")
#let fn = math.cal("n")
#let fo = math.cal("o")
#let fp = math.cal("p")
#let fq = math.cal("q")
#let fr = math.cal("r")
#let fs = math.cal("s")
#let ft = math.cal("t")
#let fu = math.cal("u")
#let fv = math.cal("v")
#let fw = math.cal("w")
#let fx = math.cal("x")
#let fy = math.cal("y")
#let fz = math.cal("z")

//uppercase math cal characters
#let fA = math.cal("A")
#let fB = math.cal("B")
#let fC = math.cal("C")
#let fD = math.cal("D")
#let fE = math.cal("E")
#let fF = math.cal("F")
#let fG = math.cal("G")
#let fH = math.cal("H")
#let fI = math.cal("I")
#let fJ = math.cal("J")
#let fK = math.cal("K")
#let fL = math.cal("L")
#let fM = math.cal("M")
#let fN = math.cal("N")
#let fO = math.cal("O")
#let fP = math.cal("P")
#let fQ = math.cal("Q")
#let fR = math.cal("R")
#let fS = math.cal("S")
#let fT = math.cal("T")
#let fU = math.cal("U")
#let fV = math.cal("V")
#let fW = math.cal("W")
#let fX = math.cal("X")
#let fY = math.cal("Y")
#let fZ = math.cal("Z")

#let TM = $T cal(M)$
#let TN = $T cal(N)$
#let TS = $T cal(S)$

//highlithing in green for math mode
#let mhilite(T) = {
  (text(color.green)[$#T$])
}

#let differentialSymbol = $#h(0.1em)d$

//DIFFERENTIALS
//invariant measure with dimension D and coordinate V
#let meas(D,V) = {
  (math.sqrt("g"), math.thin, math.attach("d",t:D), V, math.thin).join()
}

//lowercase
#let da = $differentialSymbol a$
#let db = $differentialSymbol b$
#let dc = $differentialSymbol c$
#let dd = $differentialSymbol d$
#let de = $differentialSymbol e$
#let df = $differentialSymbol f$
#let dg = $differentialSymbol g$
#let dh = $differentialSymbol h$
#let di = $differentialSymbol i$
#let dj = $differentialSymbol j$
#let dk = $differentialSymbol k$
#let dl = $differentialSymbol l$
#let dm = $differentialSymbol m$
#let dn = $differentialSymbol n$
#let do = $differentialSymbol o$
#let dp = $differentialSymbol p$
#let dq = $differentialSymbol q$
#let dr = $differentialSymbol r$
#let ds = $differentialSymbol s$
#let dt = $differentialSymbol t$
#let du = $differentialSymbol u$
#let dv = $differentialSymbol v$
#let dw = $differentialSymbol w$
#let dx = $differentialSymbol x$
#let dy = $differentialSymbol y$
#let dz = $differentialSymbol z$

#let dms = $differentialSymbol macron(s)$

//uppercase
#let dA = $differentialSymbol A$
#let dB = $differentialSymbol B$
#let dC = $differentialSymbol C$
#let dD = $differentialSymbol D$
#let dE = $differentialSymbol E$
#let dF = $differentialSymbol F$
#let dG = $differentialSymbol G$
#let dH = $differentialSymbol H$
#let dI = $differentialSymbol I$
#let dJ = $differentialSymbol J$
#let dK = $differentialSymbol K$
#let dL = $differentialSymbol L$
#let dM = $differentialSymbol M$
#let dN = $differentialSymbol N$
#let dO = $differentialSymbol O$
#let dP = $differentialSymbol P$
#let dQ = $differentialSymbol Q$
#let dR = $differentialSymbol R$
#let dS = $differentialSymbol S$
#let dT = $differentialSymbol T$
#let dU = $differentialSymbol U$
#let dV = $differentialSymbol V$
#let dW = $differentialSymbol W$
#let dX = $differentialSymbol X$
#let dY = $differentialSymbol Y$
#let dZ = $differentialSymbol Z$

//lowercase vector
#let dva = $differentialSymbol va$
#let dvb = $differentialSymbol vb$
#let dvc = $differentialSymbol vc$
#let dvd = $differentialSymbol vd$
#let dve = $differentialSymbol ve$
#let dvf = $differentialSymbol vf$
#let dvg = $differentialSymbol vg$
#let dvh = $differentialSymbol vh$
#let dvi = $differentialSymbol vi$
#let dvj = $differentialSymbol vj$
#let dvk = $differentialSymbol vk$
#let dvl = $differentialSymbol vl$
#let dvm = $differentialSymbol vm$
#let dvn = $differentialSymbol vn$
#let dvo = $differentialSymbol vo$
#let dvp = $differentialSymbol vp$
#let dvq = $differentialSymbol vq$
#let dvr = $differentialSymbol vr$
#let dvs = $differentialSymbol vs$
#let dvt = $differentialSymbol vt$
#let dvu = $differentialSymbol vu$
#let dvv = $differentialSymbol vv$
#let dvw = $differentialSymbol vw$
#let dvx = $differentialSymbol vx$
#let dvy = $differentialSymbol vy$
#let dvz = $differentialSymbol vz$

//uppercase vector
#let dvA = $differentialSymbol vA$
#let dvB = $differentialSymbol vB$
#let dvC = $differentialSymbol vC$
#let dvD = $differentialSymbol vD$
#let dvE = $differentialSymbol vE$
#let dvF = $differentialSymbol vF$
#let dvG = $differentialSymbol vG$
#let dvH = $differentialSymbol vH$
#let dvI = $differentialSymbol vI$
#let dvJ = $differentialSymbol vJ$
#let dvK = $differentialSymbol vK$
#let dvL = $differentialSymbol vL$
#let dvM = $differentialSymbol vM$
#let dvN = $differentialSymbol vN$
#let dvO = $differentialSymbol vO$
#let dvP = $differentialSymbol vP$
#let dvQ = $differentialSymbol vQ$
#let dvR = $differentialSymbol vR$
#let dvS = $differentialSymbol vS$
#let dvT = $differentialSymbol vT$
#let dvU = $differentialSymbol vU$
#let dvV = $differentialSymbol vV$
#let dvW = $differentialSymbol vW$
#let dvX = $differentialSymbol vX$
#let dvY = $differentialSymbol vY$
#let dvZ = $differentialSymbol vZ$

//lowercase greek
#let dalpha = $differentialSymbol alpha$
#let dbeta = $differentialSymbol beta$
#let dgamma = $differentialSymbol gamma$
#let ddelta = $differentialSymbol delta$
#let depsilon = $differentialSymbol epsilon$
#let dzeta = $differentialSymbol zeta$
#let deta = $differentialSymbol eta$
#let dtheta = $differentialSymbol theta$
#let diota = $differentialSymbol iota$
#let dkappa = $differentialSymbol kappa$
#let dlambda = $differentialSymbol lambda$
#let dmu = $differentialSymbol mu$
#let dnu = $differentialSymbol nu$
#let dxi = $differentialSymbol xi$
#let domicron = $differentialSymbol omicron$
#let dpi = $differentialSymbol pi$
#let drho = $differentialSymbol rho$
#let dsigma = $differentialSymbol sigma$
#let dtau = $differentialSymbol tau$
#let dupsilon = $differentialSymbol upsilon$
#let dphi = $differentialSymbol phi$
#let dchi = $differentialSymbol chi$
#let dpsi = $differentialSymbol psi$
#let domega = $differentialSymbol omega$

//uppercase greek
#let dGamma = $differentialSymbol Gamma$
#let dDelta = $differentialSymbol Delta$
#let dTheta = $differentialSymbol Theta$
#let dLambda = $differentialSymbol Lambda$
#let dXi = $differentialSymbol Xi$
#let dPi = $differentialSymbol Pi$
#let dSigma = $differentialSymbol Sigma$
#let dUpsilon = $differentialSymbol Upsilon$
#let dPhi = $differentialSymbol Phi$
#let dPsi = $differentialSymbol Psi$
#let dOmega = $differentialSymbol Omega$

