#import "template.typ": *
#import "macros.typ": *
#import "@preview/xarrow:0.3.1": xarrow
#import "@preview/fletcher:0.4.5" as fletcher: diagram, node, edge
#import "@preview/cetz:0.5.1": canvas, draw
#import "@preview/cetz-plot:0.1.4": plot

#show: project.with(
  //title: "Dirac-Bergmann and Hamiltonian Field Theory",
  title: "Stuff I've derived",
  authors: (
    (name: "Niels Slotboom", email: "slotboom.n@gmail.com"),
  ),
)

#outline(indent:auto)
#pagebreak()
#include("num-methods.typ")
#pagebreak()
#include("num-relativity.typ")