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
#let rank = math.op("rank")

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
#let otimes = math.times.circle
#let oplus = math.plus.circle
#let star = (math.star.op,h(1pt) ).join()
#let flat = [\u{266D}]
#let sharp = [\u{266F}]
#let codim = math.op("codim")
#let dim = math.op("dim")
#let mnabla = math.macron(math.nabla)


//lowercase bold vectors
#let va = math.bold("a")
#let vb = math.bold("b")
#let vc = math.bold("c")
#let vd = math.bold("d")
#let ve = math.bold("e")
#let vf = math.bold("f")
#let vg = math.bold("g")
#let vh = math.bold("h")
#let vi = math.bold("i")
#let vj = math.bold("j")
#let vk = math.bold("k")
#let vl = math.bold("l")
#let vm = math.bold("m")
#let vn = math.bold("n")
#let vo = math.bold("o")
#let vp = math.bold("p")
#let vq = math.bold("q")
#let vr = math.bold("r")
#let vs = math.bold("s")
#let vt = math.bold("t")
#let vu = math.bold("u")
#let vv = math.bold("v")
#let vw = math.bold("w")
#let vx = math.bold("x")
#let vy = math.bold("y")
#let vz = math.bold("z")

//uppercase bold vectors
#let vA = math.bold("A")
#let vB = math.bold("B")
#let vC = math.bold("C")
#let vD = math.bold("D")
#let vE = math.bold("E")
#let vF = math.bold("F")
#let vG = math.bold("G")
#let vH = math.bold("H")
#let vI = math.bold("I")
#let vJ = math.bold("J")
#let vK = math.bold("K")
#let vL = math.bold("L")
#let vM = math.bold("M")
#let vN = math.bold("N")
#let vO = math.bold("O")
#let vP = math.bold("P")
#let vQ = math.bold("Q")
#let vR = math.bold("R")
#let vS = math.bold("S")
#let vT = math.bold("T")
#let vU = math.bold("U")
#let vV = math.bold("V")
#let vW = math.bold("W")
#let vX = math.bold("X")
#let vY = math.bold("Y")
#let vZ = math.bold("Z")

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

