// The project function defines how your document looks.
// It takes your content and some metadata and formats it.
// Go ahead and customize it to your liking!
#let project(title: "", authors: (), body) = {
  // Set the document's basic properties.
  set document(author: authors.map(a => a.name), title: title)
  set page(numbering: "1", number-align: center)
  set text(font: "New Computer Modern", lang: "en")
  show math.equation: set text(weight: 400)
  set heading(numbering: "1.1")

  // Title row.
  align(center)[
    #block(text(weight: 700, 1.75em, title))
  ]

  // Author information.
  pad(
    top: 0.5em,
    bottom: 0.5em,
    x: 2em,
    grid(
      columns: (1fr,) * calc.min(3, authors.len()),
      gutter: 1em,
      ..authors.map(author => align(center)[
        *#author.name* \
        #author.email
      ]),
    ),
  )

  // Main body.
  set par(justify: true)
  set page(margin: 0.9in)
  set par(leading: 0.55em, first-line-indent: 1.8em, justify: true)
  set text(font: "Computer Modern")
  //#set math.equation(numbering: "(1.1)")
  show raw: set text(font: "Cascadia Code")
  show par: set block(spacing: 0.55em)
  show heading: set block(above: 1.4em, below: 1em)
  set text(font: "New Computer Modern", lang: "en")
  // set math.equation(number-align: bottom + right)
  
  set heading(numbering: "1.1")
  show heading.where(level: 1): it => counter(math.equation).update(0) + it
  show heading.where(level: 2): it => counter(math.equation).update(0) + it
  
  // show math.equation.where(block:true): set math.lr(size: 115%)

  set math.equation(numbering: num => numbering("(1.1)", counter(heading).get().first(), num))

  body
}

