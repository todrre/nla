// Mini-project report — Uppsala University template (ported from the LaTeX title page)

// ---------------------------------------------------------------------------
//  Report metadata — fill in
// ---------------------------------------------------------------------------
#let course-name = "Numerical Linear Algebra"
#let course-code = "1TD452"
#let title = "Signal and Image Processing"
#let subtitle = "Mini-project 1: Direct Methods"
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
// Short description of the application and the issue(s) investigated, in your
// own words. It should make sense without the mini-project description.
Many problems in signal and image processing have the same basic form. We want to recover an unknown signal $x$, but we can only observe a distorted and noisy version of it,
$ b = A x + e. $
Here $b$ is the observed data, $A$ is the _forward operator_ that describes how the measurement system (for example a blurring camera) distorts the signal, and $e$ is unknown noise.

A natural first attempt is to solve the least squares problem $limits(min)_x ||A x - b||_2^2$. But in  many applications the matrix $A$ is ill-conditioned. This means that small changes in $b$ cause large changes on the solution.

What we often do is apply _regularization_. We add a penalty term that rewards solutions we consider reasonable. This gives us the regularized least squares problem
$ min_x ||A x - b||_2^2 + lambda^2 ||L x||_2^2. $ <eq-general>
The parameter $lambda > 0$ controls the strength of the regularization. The matrix $L$ encodes what we know about the signal in advance and depends on the application. For example, $L = I$ favours solutions of small size, while a difference operator such as
$
  L = mat(
    1, -1, , , ;
    , 1, -1, , ;
    , , dots.down, dots.down, ;
    , , , 1, -1
  ) in RR^((n-1) times n)
$
favours smooth solutions, since $(L x)_i = x_i - x_(i+1)$ penalizes large jumps between neighbouring entries.

= Method
// Summarize the approach and explain the chosen algorithms and why they are
// appropriate. Short code snippets or pseudocode are fine; full code goes in the appendix.

== Tikhonov regularization and SVD
Singular value decomposition (SVD) is a useful tool for understanding least-squares problems. We split $A$ into simple building blocks.
$ A = U Sigma V^T, $
$U$ and $V$ are orthogonal matrices and $Sigma$ is a diagonal matrix holding the singular values $sigma_1 >= sigma_2 >= dots >= 0$ of $A$. We denote the columns of $U$ and $V$ by $u_i$ and $v_i$.

Using SVD, the ordinary least squares solution of $limits(min)_x ||A x - b||_2^2$ can be written as
$ x_"LS" = sum_(i=1)^n (u_i^T b) / sigma_i v_i. $
Here we can see the problem. When $A$ is ill-conditioned, some $sigma_i$ are tiny, and dividing by them blows up the noise in $b$. To prevent this, Tikhonov regularization can be used which adds a penalty on the size of $x$
$ min_x f(x) := ||A x - b||_2^2 + lambda^2 ||x||_2^2, $
which is @eq-general with $L = I$. Expanding $f$ gives
$
  f(x) = (A x - b)^T (A x - b) + lambda^2 x^T x
  = x^T A^T A x - 2 x^T A^T b + b^T b + lambda^2 x^T x.
$

Since $f$ is strictly convex, $x$ is the minimizer exactly when the gradient vanishes, $nabla f(x) = 2 (A^T A + lambda^2 I) x - 2 A^T b = 0$. This is the same as the _normal equations_
$ (A^T A + lambda^2 I) x = A^T b. $
The matrix $A^T A + lambda^2 I$ is positive definite and hence invertible, so the normal equations have exactly one solution, $x_"RLS" = (A^T A + lambda^2 I)^(-1) A^T b$.

#let c(x) = box(
  height: 1.6em,
  baseline: 50% - 0.3em,
  align(horizon, x),
)

To see what the regularization does, we insert the SVD into the normal equations and simplify using $U^T U = I$ and $V^T V = V V^T = I$:
$
                                        (A^T A + lambda^2 I) x_"RLS" & = A^T b \
  V Sigma^T cancel(U^T U) Sigma V^T + lambda^2 cancel(V V^T) x_"RLS" & = V Sigma^T U^T b \
                          V (Sigma^T Sigma + lambda^2 I) V^T x_"RLS" & = V Sigma^T U^T b \
              cancel(V^T V) (Sigma^T Sigma + lambda^2 I) V^T x_"RLS" & = cancel(V^T V) Sigma^T U^T b \
                             underbrace(
                               (Sigma^T Sigma + lambda^2 I), B= mat(
                                 sigma_1^2 + lambda^2, , , ;
                                 , sigma_2^2 + lambda^2, , ;
                                 , , dots.down, ;
                                 , , , sigma_n^2 + lambda^2;
                               )
                             ) underbrace(V^T x_"RLS", y)            & = Sigma^T U^T b \
                                                                => y & = B^(-1) Sigma^T U^T b \
                                                                     & = mat(
                                                                         #c($1/(sigma_1^2 + lambda^2)$), , , ;
                                                                         , #c($1/(sigma_2^2 + lambda^2)$), , ;
                                                                         , , #c($dots.down$), ;
                                                                         , , , #c($1/(sigma_n^2 + lambda^2)$)
                                                                       )
                                                                       mat(
                                                                         #c($sigma_1$), , , ;
                                                                         , #c($sigma_2$), , ;
                                                                         , , #c($dots.down$), ;
                                                                         , , , #c($sigma_n$)
                                                                       )
                                                                       mat(
                                                                         #c($u_1^T$) ;
                                                                         #c($u_2^T$) ;
                                                                         #c($dots.v$) ;
                                                                         #c($u_n^T$) ;
                                                                       ) b
$

Because $B^(-1)$ and $Sigma$ are both diagonal, each component of $y$ can be read off directly:
$ y_i = sigma_i / (sigma_i^2 + lambda^2) u_i^T b, quad i = 1, dots, n, $
and since $V$ is orthogonal we recover the solution as $x_"RLS" = V y$.

We assume $A$ has full rank, so every $sigma_i > 0$.
We can then split the coefficient and introduce the _filter factor_ $phi.alt_i (lambda)$:
$
  sigma_i / (sigma_i^2 + lambda^2) = 1 / sigma_i dot sigma_i^2 / (sigma_i^2 + lambda^2) = (phi.alt_i (lambda)) / sigma_i
  quad => quad
  y_i = (phi.alt_i (lambda)) / sigma_i u_i^T b
$

Writing $x_"RLS" = V y = sum_i y_i v_i$ as a sum then gives
$
  x_"RLS" = V y = sum_(i=1)^n sigma_i / (sigma_i^2 + lambda^2) (u_i^T b) v_i
  = sum_(i=1)^n phi.alt_i (lambda)(u_i^T b) / sigma_i v_i
$

Compared with $x_"LS"$, Tikhonov adds the filter factor $phi.alt_i (lambda)$, which keeps terms with large singular values and suppresses those with small ones. The parameter $lambda$ sets the transition.
- if $sigma_i >> lambda$, then $phi.alt_i (lambda) approx 1$ and the term is kept almost unchanged;
- if $sigma_i = lambda$, then $phi.alt_i (lambda) = 1/2$;
- if $sigma_i << lambda$, then $phi.alt_i (lambda) approx 0$ and the term is damped away.

== Truncated SVD as a filter
The filter-factor form is not specific to Tikhonov regularization. Any method that damps each term of $x_"LS"$ can be written as
$ x = sum_(i=1)^n phi.alt_i (u_i^T b) / sigma_i v_i, $
and methods differ only in their choice of $phi.alt_i$. Where Tikhonov damps the terms smoothly, the truncated SVD simply cuts them off. We pick a threshold $delta > 0$ and let $k$ be the number of singular values larger than $delta$. The filter factors are then
$ phi.alt_i = cases(1 & "if" i <= k quad (sigma_i > delta), 0 & "if" i > k quad (sigma_i <= delta)). $
Substituting this in, only the first $k$ terms survive:
$
  x = sum_(i=1)^k (u_i^T b) / sigma_i v_i
  + sum_(i=k+1)^n 0 dot (u_i^T b) / sigma_i v_i
  = x_"TLS".
$
The truncated least squares solution is thus a filter method with a sharp cut-off at $sigma_i = delta$, instead of the smooth transition of Tikhonov regularization.























== Generalized SVD
So far we have used $L = I$, which only penalizes the size of $x$. For a general $L$ the SVD of $A$ is no longer enough. Instead we use the generalized SVD (GSVD), which decomposes both matrices at once:
$ A = U Sigma_A W^(-1), quad L = V Sigma_L W^(-1). $ <eq-gsvd>
Here $U$ and $V$ are orthogonal, $W$ is invertible, and $Sigma_A$ and $Sigma_L$ are diagonal. The ratios $gamma_i = alpha_i \/ beta_i$ are the _generalized singular values_ of $(A, L)$.

== General-form regularization and GSVD
We first characterize the solution of the general problem in @eq-general. Expanding its objective $f(x)$ as in the case $L = I$ but now with a general $L$ gives
$
  f(x) & = (A x - b)^T (A x - b) + lambda^2 (L x)^T (L x) \
       & = x^T A^T A x - 2 x^T A^T b + b^T b + lambda^2 x^T L^T L x.
$
Since $f$ is convex, $x$ is a minimizer exactly when $nabla f(x) = 2 (A^T A + lambda^2 L^T L) x - 2 A^T b = 0$.  This is the normal equations for the general problem
$
  (A^T A + lambda^2 L^T L) x = A^T b.
$ <eq-normal-general>
The solution is unique when $A^T A + lambda^2 L^T L$ is invertible. This is the case if $A^T A + lambda^2 L^T L$ is positive definite, meaning, for $x!=0$,  $ x^T (A^T A + lambda^2 L^T L) x = ||A x||^2 + lambda^2 ||L x||^2 > 0 $
This implies that, when there is no unique solution, $A x = L x = 0$, which is true when $x$ is in the nullspace of both $A$ and $L$. So we require that $A$ and $L$ have no common nullspace, which is done by the choice of $L$.

Instead of solving @eq-normal-general directly, we use the GSVD @eq-gsvd, which shows what the regularization does to each component. Here
$
  Sigma_A = "diag"(alpha_1, alpha_2, dots, alpha_r, 0), quad
  Sigma_L = "diag"(beta_1, beta_2, dots, beta_r, 0),
$
and we assume all $alpha_i > 0$. Substituting the GSVD into the expanded $f(x)$ gives
$
  f(x) & = x^T W^(-T) Sigma_A^T cancel(U^T U) Sigma_A W^(-1) x - 2 x^T W^(-T) Sigma_A^T U^T b + b^T b \
       & quad + lambda^2 x^T W^(-T) Sigma_L^T cancel(V^T V) Sigma_L W^(-1) x.
$
With the change of variables $y = W^(-1) x$ this becomes
$
  f = y^T (Sigma_A^T Sigma_A + lambda^2 Sigma_L^T Sigma_L) y - 2 y^T Sigma_A^T U^T b + b^T b.
$
Since $W$ is invertible, minimizing over $x$ is the same as minimizing over $y$. Setting the gradient with respect to $y$ to zero gives
$
  nabla f(y) = 2 (Sigma_A^T Sigma_A + lambda^2 Sigma_L^T Sigma_L) y - 2 Sigma_A^T U^T b & = 0 \
              => quad underbrace((Sigma_A^T Sigma_A + lambda^2 Sigma_L^T Sigma_L), B) y & = Sigma_A^T U^T b.
$

$
  B := Sigma_A^T Sigma_A + lambda^2 Sigma_L^T Sigma_L & =
                                                        mat(
                                                          alpha_1^2, , , , , ;
                                                          , alpha_2^2, , , ;
                                                          , , dots.down, , , ;
                                                          , , , alpha_r^2, ;
                                                          , , , , , 0
                                                        ) + lambda^2 mat(
                                                          beta_1^2, , , , ;
                                                          , beta_2^2, , , ;
                                                          , , dots.down, , ;
                                                          , , , beta_r^2, ;
                                                          , , , , 0
                                                        ) \
                                                      & =
                                                        mat(
                                                          alpha_1^2 + lambda^2 beta_1^2, , , , ;
                                                          , alpha_2^2 + lambda^2 beta_2^2, , , ;
                                                          , , dots.down, , ;
                                                          , , , alpha_r^2 + lambda^2 beta_r^2, ;
                                                          , , , , 0
                                                        )
$
#pagebreak()
Solving for $y$ gives
$
  y & = B^(-1) Sigma_A^T U^T b \
    & = mat(
        #c($1/(alpha_1^2 + lambda^2 beta_1^2)$), , , ;
        , #c($1/(alpha_2^2 + lambda^2 beta_2^2)$), , ;
        , , #c($dots.down$), ;
        , , , #c($1/(alpha_n^2 + lambda^2 beta_n^2)$)
      )
      mat(
        #c($alpha_1$), , , ;
        , #c($alpha_2$), , ;
        , , #c($dots.down$), ;
        , , , #c($alpha_n$)
      )
      mat(
        #c($u_1^T$) ;
        #c($u_2^T$) ;
        #c($dots.v$) ;
        #c($u_n^T$) ;
      ) b
$
Because $B^(-1)$ and $Sigma_A$ are both diagonal, each component of $y$ can be read off directly:
$
  y_i = alpha_i / (alpha_i^2 + lambda^2 beta_i^2) u_i^T b, quad i = 1, dots, n,
$
and we recover the solution as $x = W y$.

As for Tikhonov, we split the coefficient and introduce the filter factor $phi.alt_i (lambda)$
$
  alpha_i / (alpha_i^2 + lambda^2 beta_i^2) = 1 / alpha_i dot alpha_i^2 / (alpha_i^2 + lambda^2 beta_i^2) = (phi.alt_i (lambda)) / alpha_i,
  quad
  phi.alt_i (lambda) = alpha_i^2 / (alpha_i^2 + lambda^2 beta_i^2) = gamma_i^2 / (gamma_i^2 + lambda^2).
$
Writing $x = W y = sum_i y_i w_i$, where $w_i$ are the columns of $W$, as a sum then gives
$
  x = sum_(i=1)^n alpha_i / (alpha_i^2 + lambda^2 beta_i^2) (u_i^T b) w_i
  = sum_(i=1)^n phi.alt_i (lambda) (u_i^T b) / alpha_i w_i.
$

Here $gamma_i = alpha_i \/ beta_i$ compares how well $A$ measures the component $w_i$ with how much $L$ penalizes it:
- $gamma_i >> lambda$: $phi.alt_i (lambda) approx 1$, the term is kept;
- $gamma_i << lambda$: $phi.alt_i (lambda) approx 0$, the term is damped.
Since $alpha_i^2 + beta_i^2 = 1$, the damped terms have small $alpha_i$ and large $beta_i$. When $L = I$, only $A$ decides what is damped. For a general $L$, $alpha_i$ and $beta_i$ depend on both $A$ and $L$. So by choosing $L$ we decide which components are treated as unstable.

== One-dimensional deblurring
We now apply the method to a concrete problem: removing blur from a one-dimensional signal. A blurred signal is modelled as
$
  b(s) = integral_(-1)^1 K(s, t) x(t) dif t + e(s), quad K(s, t) = c e^(-beta (s - t)^2),
$
where the Gaussian kernel $K$ spreads each point of $x$ over its neighbours, and a larger $beta$ gives a narrower blur.
As true signal we use $x(t) = e^t sin(pi t)$. We use the box quadrature rule to turn the integral into a sum. We split $[-1, 1]$ into $n$ boxes of width $h = 2 \/ n$ with midpoints $t_i$. The integral is then approximated by the sum of the box areas,
$
  b(t_i) approx sum_(j=1)^n h K(t_i, t_j) x(t_j) + e(t_i).
$
This is $b = A x + e$ with $A_(i j) = h K(t_i, t_j)$, and we take the noise as $e = delta dot cal(N)(0, 1)$. We choose $c = sqrt(beta \/ pi)$ so that the blur does not change the size of the signal, and use the first-difference matrix from the introduction as $L$.

The unregularized solution is computed with a standard least-squares solver. For the regularized solution we compute the GSVD of $(A, L)$ once and use $x = W y$ with $y_i = alpha_i \/ (alpha_i^2 + lambda^2 beta_i^2) u_i^T b$.
To choose $lambda$ without knowing the true signal, we use the _discrepancy principle_ @KARL2005183. The idea is that the solution should fit the data as well as the noise allows, but not better. So we pick the $lambda$ where the residual is as large as the noise, $||A x_lambda - b||_2 = ||e||_2 approx delta sqrt(n)$.


= Results
// Results for each part. Show that the program works: figures, (trimmed)
// output, and so on.
// Use the same true signal and parameters as in the Method section throughout.

== Regularized vs. unregularized solution
Since $"cond"(A) approx 10^19$, the unregularized solution is only amplified noise, with a relative error of about $10^10$. The regularized solution follows the true signal well. Its error is $0.055$, compared with $0.225$ for the data $b$, see @fig-task5-1.

#figure(
  image("assets/task5_1.svg", width: 100%),
  caption: [Regularized (left) and unregularized (right) solution.],
) <fig-task5-1>

#pagebreak()
== Filter factors
Components with $gamma_i >> lambda$ are kept and those with $gamma_i << lambda$ are suppressed. Increasing $lambda$ moves the cut-off to the right, so more components are removed, see @fig-task5-2.

#figure(
  image("assets/task5_2.svg", width: 65%),
  caption: [Filter factors for different $lambda$. The dotted vertical lines show the point where $gamma_i = lambda$.],
) <fig-task5-2>

== Effect of $lambda$
Too small $lambda$ lets the noise through (under-regularization), too large $lambda$ smooths away the signal (over-regularization). The smallest error, $0.054$, is at $lambda approx 0.093$, close to our hand-picked $lambda = 0.1$, see @fig-task5-3.

#figure(
  image("assets/task5_3_error.svg", width: 100%),
  caption: [Relative error against $lambda$.],
) <fig-task5-3>

== Effect of the noise level
The unregularized error grows in proportion to $delta$, while the regularized error grows only slowly and the best $lambda$ increases with $delta$, see @tab-noise and @fig-task5-4.

#figure(
  table(
    columns: 4,
    align: (right, right, right, right),
    stroke: none,
    table.hline(),
    table.header([$delta$], [error, unregularized], [best $lambda$], [error, regularized]),
    table.hline(stroke: 0.5pt),
    [$10^(-3)$], [$7.1 dot 10^8$], [$0.027$], [$0.040$],
    [$10^(-2)$], [$7.1 dot 10^9$], [$0.093$], [$0.054$],
    [$10^(-1)$], [$7.1 dot 10^10$], [$0.66$], [$0.116$],
    table.hline(),
  ),
  caption: [Errors for different noise levels.],
) <tab-noise>

#figure(
  image("assets/task5_4.svg", width: 100%),
  caption: [Regularized solutions for different noise levels.],
) <fig-task5-4>

== Automatic choice of $lambda$
The discrepancy principle gives a somewhat larger $lambda$ and a slightly larger error than the hand-picked and the best value, see @tab-lambda.

#figure(
  table(
    columns: 3,
    align: (left, right, right),
    stroke: none,
    table.hline(),
    table.header([method], [$lambda$], [relative error]),
    table.hline(stroke: 0.5pt),
    [discrepancy principle], [$0.28$], [$0.084$],
    [chosen by hand], [$0.1$], [$0.055$],
    [best (from @fig-task5-3)], [$0.093$], [$0.054$],
    table.hline(),
  ),
  caption: [Choice of $lambda$.],
) <tab-lambda>

= Discussion
// Further comments: observations, open questions, possible extensions, and
// any discussion questions from the mini-project description.
Since $"cond"(A) approx 10^19$, even small noise blows up without regularization. Regularization removes the components with small generalized singular values, where this happens. More noise needs a larger $lambda$, so regularization becomes more important when the noise level grows.

In practice the true signal is unknown, so the best $lambda$ cannot be found directly. The discrepancy principle only needs the noise level and still gives a good reconstruction, although it slightly over-regularizes, which is typical for the method @KARL2005183.

// ---------------------------------------------------------------------------
//  References
// ---------------------------------------------------------------------------
#set heading(numbering: none)
#bibliography("references.bib", title: "References")
