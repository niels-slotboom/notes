#import "template.typ": *
#import "macros.typ": *
#import "@preview/xarrow:0.3.1": xarrow
#import "@preview/fletcher:0.4.5" as fletcher: diagram, node, edge

#show: project.with(
  //title: "Dirac-Bergmann and Hamiltonian Field Theory",
  title: "The ADM Action and Differential Geometry on Foliations",
  authors: (
    (name: "Niels Slotboom", email: "slotboom.n@gmail.com"),
  ),
)

#outline(indent:auto)

//#include("HamiltonianFormalism.typ")
#include("Submanifolds.typ") 
