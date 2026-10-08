// Mini-project report — Uppsala University template (ported from the LaTeX title page)

// ---------------------------------------------------------------------------
//  Report metadata — fill in
// ---------------------------------------------------------------------------
#let course-name = "Numerical Linear Algebra"
#let course-code = "1TD452"
#let title = "Signal and Image Processing"
#let subtitle = "Mini-project 2: Iterative methods"
#let authors = ("Pontus Ahlberg", "Ivar Hammarberg", "Isac Persson")

// ---------------------------------------------------------------------------
//  Document settings
// ---------------------------------------------------------------------------
#set document(title: title, author: authors)
#set page(paper: "a4", margin: 2.5cm)
#set text(size: 12pt, lang: "en")
#set par(justify: true)
#set math.mat(delim: "[")
#set math.equation(numbering: "(1)")
// Number only labelled (i.e. referenced) equations.
#show math.equation: it => {
  if it.block and not it.has("label") [
    #counter(math.equation).update(n => n - 1)
    #math.equation(it.body, block: true, numbering: none)<unnumbered>
  ] else { it }
}
#set figure(gap: 0.8em)
#show figure.caption: set text(size: 0.9em)

// ---------------------------------------------------------------------------
//  Title page
// ---------------------------------------------------------------------------
#let hrule = line(length: 100%, stroke: 0.5mm)

#page(numbering: none)[
  #set align(center)
  #set par(justify: false)
  #set text(hyphenate: false)

  #v(1cm)
  #image("assets/UU_logo_CMYK_flat.pdf", width: 10cm)
  #v(1fr)

  #text(size: 14pt, smallcaps(course-name + " · " + course-code))
  #v(0.6cm)
  #hrule
  #v(0.5cm)
  #text(size: 24pt, weight: "bold", title)
  #v(0.3cm)
  #text(size: 15pt, style: "italic", subtitle)
  #v(0.5cm)
  #hrule
  #v(1.2cm)

  #text(size: 14pt, authors.join(linebreak()))

  #v(1fr)
  #text(size: 12pt, datetime.today().display("[day] [month repr:long] [year]"))
  #v(1cm)
]

// ---------------------------------------------------------------------------
//  Main text
// ---------------------------------------------------------------------------
// No page number on the table of contents; the main text starts at 1.
#set page(numbering: none)
#set heading(numbering: "1.1")
#import "@preview/pavemat:0.2.0": pavemat
#import "@preview/mannot:0.4.0": *

// One colour per role, reused in every annotated equation and in the Φ figure.
#let c-orth = rgb("#d97706")  // products that collapse to I
#let c-filt = rgb("#2563eb")  // filter factors φ_i
#let c-ls = rgb("#16a34a")    // plain least-squares coefficients
#let c-cut = rgb("#6b7280")   // discarded terms
// A soft highlight: `hl(x, #c-filt)` or `hl(x, #c-filt, #<tag>)` in math.
#let hl(body, c, ..tag) = markhl(body, c, ..tag, fill: c.transparentize(85%), radius: 2pt)
#let note = annot.with(annot-text-props: (size: 0.75em), leader-tip: none, leader-toe: none)
#set math.mat(delim: "[")


#outline()
#pagebreak()
#set page(numbering: "1")
#counter(page).update(1)

= Introduction


= Task 7


// ---------------------------------------------------------------------------
//  References
// ---------------------------------------------------------------------------
#set heading(numbering: none)
#bibliography("references.bib", title: "References")
