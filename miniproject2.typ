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
While direct methods with regularization are useful for computing the solution to least squares problems, they require dense matrix factorizations and are therefore not suitable for large-scale problems. In this project, we explore iterative methods for solving problems. Eller nåt sånt.

= Method
== Iterative Krylov subspace methods




= Task 7: Prove that the matrix is SPD
Let $A in RR^(m times n)$ and $L in RR^(p times n)$.
Assume that $lambda != 0$ and

$
  "rank"( mat(A; L) ) = n.
$

Define the matrix

$
  M = A^T A + lambda^2 L^T L.
$

We want to prove that $M$ is symmetric positive definite.

*Is $M$ symmetric?*

A matrix is symmetric if $M^T = M$. Taking the transpose gives

$
  M^T & = (A^T A + lambda^2 L^T L)^T \
      & =(A^T A)^T + lambda^2 (L^T L)^T \
      & = A^T A + lambda^2 L^T L = M.
$

Therefore, $M$ is symmetric.

*Is $M$ positive definite?*

A symmetric matrix $M$ is positive definite if

$
  z^T M z > 0
  quad "for all" z != 0.
$



$
  => z^T M z & = z^T (A^T A + lambda^2 L^T L) z \
             & = z^T A^T A z + lambda^2 z^T L^T L z \
             & = norm(A z)_2^2 + lambda^2 norm(L z)_2^2
$


Since squared Euclidean norms are always non-negative and
$lambda != 0$ implies $lambda^2 > 0$ we get:


$
  z^T M z
  = norm(A z)_2^2 + lambda^2 norm(L z)_2^2
  >= 0.
$

To prove positive definiteness, we must show that
the inequality is strict whenever $z != 0$. The assumption $"rank"( mat(A; L) ) = n$ means that the matrix obtained by stacking $A$ and $L$ has full column rank. Therefore, its null space contains only the zero vector:

$
  mat(A; L) z = 0
  quad & "Only if" quad z = 0 \
  => A z = 0 "and" L z = 0
  quad & "Only if" quad z = 0
$

Now let $z != 0$. It follows that $A z$ and $L z$
cannot both be zero at the same time. In other words at least on is greater than zero. Which means $M$ is symmetric positive definite.


If M is symmetric positive definite, we can assume that in the absence of rounding errors, the CG algorithm will converge exactly in at most n iterations. (Theorem 3.9, NLA_iterative_methods)

= Task 8: Implement the CG algorithm

#figure(
  stack(
    spacing: 0.5em,
    image("8.1.png", width: 100%),
    image("8.1.graph.png", width: 100%),
  ),
  caption: [Influence of $lambda$ on the CG reconstruction and convergence.],
) <fig:task8-lambda>

#figure(
  stack(
    spacing: 0.5em,
    image("8.2.png", width: 100%),
    image("8.2.graph.png", width: 100%),
  ),
  caption: [Influence of the tolerance on the CG reconstruction and convergence ($lambda = 0.05$).],
) <fig:task8-tol>

#figure(
  stack(
    spacing: 0.5em,
    image("8.3.png", width: 100%),
    image("8.3.graph.png", width: 100%),
  ),
  caption: [Influence of the maximum number of iterations on the CG reconstruction and convergence ($lambda = 0.01$).],
) <fig:task8-maxiter>

= Task 10
Preconditioner, om använda cg bevara spd
$L^{-1}A L^{-T}y=L^{-1}y=L^{-1}b, y=L^{T}x$
$M=LL^{T}$ är preconditioner då cholesky bevarar spd.
// ---------------------------------------------------------------------------
//  References
// ---------------------------------------------------------------------------
#set heading(numbering: none)
#bibliography("references.bib", title: "References")

= Task 11
SUPERAISVAR:LSQR can solve the regularized least-squares problem 
\(\min_x \|Ax-b\|_2^2+\lambda^2\|Lx\|_2^2\) by rewriting it as an ordinary least-squares problem with the augmented matrix \(\widetilde A=\begin{bmatrix}A\\\lambda L\end{bmatrix}\) and vector \(\widetilde b=\begin{bmatrix}b\\0\end{bmatrix}\). It requires only products with \(\widetilde A\) and \(\widetilde A^T\), which can be computed using \(A,A^T,L,L^T\). Therefore, LSQR avoids explicitly forming \(A^TA\) and \(L^TL\), reducing memory requirements and avoiding the additional rounding errors and potential sparsity loss associated with forming the normal equations.
CG can also avoid explicitly forming the normal equations by applying matrix-vector products through \(A^T(Av)+\lambda^2L^T(Lv)\). However, it solves the resulting normal-equation system, whose condition number is the square of that of the augmented matrix: \(\kappa_2(\widetilde A^T\widetilde A)=\kappa_2(\widetilde A)^2\), assuming full column rank. LSQR works directly with the augmented least-squares operator and generally has better numerical behavior for ill-conditioned problems. Both methods have similar per-iteration costs, typically requiring one forward and one transpose matrix-vector product, and can be effective for large sparse systems.
Regularization can improve conditioning by reducing the influence of poorly determined solution components, although it introduces a trade-off between fitting the data and limiting the solution norm. Neither LSQR nor CG guarantees an accurate solution for severely ill-conditioned problems, and both can exhibit semiconvergence, where later iterations begin fitting noise. In practice, LSQR is generally the preferred starting point for large, ill-conditioned least-squares problems, while CG may be appropriate when an effective preconditioner or an existing normal-equation solver is available. Early stopping and suitable regularization are important for obtaining a stable, meaningful solution.