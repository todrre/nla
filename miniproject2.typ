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

// ---------------------------------------------------------------------------
//  Numbers from task6.ipynb (written by save_results, updated on every rerun)
// ---------------------------------------------------------------------------
#let t6 = json("assets/task6.json")
#let t6A = t6.memory.A
#let t6L = t6.memory.L
#let mth(s) = math.equation(eval(s, mode: "math"))
// x with a fixed number of decimals, as a string
#let fixed(x, d) = {
  let r = calc.round(x * calc.pow(10.0, d))
  let digits = str(int(r))
  if d == 0 { return digits }
  while digits.len() <= d { digits = "0" + digits }
  digits.slice(0, digits.len() - d) + "." + digits.slice(digits.len() - d)
}
// x with a number of significant digits, as a string (x > 0)
#let sigs(x, digits: 3) = fixed(x, calc.max(0, digits - 1 - calc.floor(calc.log(x))))
// integer with thin spaces between the thousands: 90 000
#let thousands(x) = {
  let d = str(int(calc.round(x)))
  let out = ""
  for (i, ch) in d.clusters().enumerate() {
    if i > 0 and calc.rem(d.len() - i, 3) == 0 { out += " thin " }
    out += ch
  }
  mth(out)
}
#let num(x, digits: 3) = mth(sigs(x, digits: digits))
// a dot 10^b
#let sci(x, digits: 3) = {
  let e = calc.floor(calc.log(x))
  let m = sigs(x / calc.pow(10.0, e), digits: digits)
  if m.starts-with("10") {
    e += 1
    m = sigs(1.0, digits: digits)
  }
  if float(m) == 1.0 and digits == 1 { mth("10^(" + str(e) + ")") } else { mth(m + " dot 10^(" + str(e) + ")") }
}
#let pow10(x) = sci(x, digits: 1)
#let mem(bytes) = if bytes >= 1e9 [#num(bytes / 1e9) GB] else [#num(bytes / 1e6) MB]
#let percent(M) = 100 * M.nnz / (M.rows * M.cols)
#let t6beta = thousands(t6.params.beta)
#let t6eps = pow10(t6.params.eps)
#set math.mat(delim: "[")


#outline()
#pagebreak()
#set page(numbering: "1")
#counter(page).update(1)

= Introduction
In this project we study iterative methods for solving linear systems. We apply these methods to recover a sharp image from a blurred one. The task boils down to a very large least squares problem, and with a dense matrix of that size the direct methods introduced earlier become too expensive. We therefore make the matrix artificially sparse and use iterative methods that are built for sparse problems. The deblurring problem is, however, ill-posed, so small amounts of noise can ruin the solution, and we show how regularization turns it into a solvable one.


= Method
Direct methods like GSVD need too much storage and computation for large image problems. We instead use iterative methods, which refine the solution step by step using mainly matrix-vector products and a few stored vectors.

== Images as vectors
#grid(
  columns: (2fr, 1fr),
  gutter: 1em,
  align: top,

  [
    A grayscale image is a matrix $X in RR^(n_1 times n_2)$ of pixel brightnesses from $0$ (black) to $1$ (white). Stacking its columns gives

    $x_"true" = "vec"(X) in RR^n, quad n = n_1 n_2$

    We resize the image to $300 times 300$ pixels, so $n = 90 000$. The goal is to recover the sharp image $x_"true"$ from the blurred and noisy observation $b = A x_"true" + e$ <eq:model> where $A$ is the blurring matrix and $e$ is the noise.
  ],

  [
    #v(-1em)
    #figure(
      image("true_image.png", width: 100%),
      caption: [Original $X_"true"$ image],
    )
  ],
)

== Least squares

Only $b$ is known, not $x_"true"$ or $e$. Without noise we could solve $A x = b$, but noise makes an exact solution unreliable, so we solve

$
  min_x norm(A x - b)_2^2
$

which minimizes the residual $A x - b$, the difference between the predicted blurred image and $b$. Least squares alone also fits the noise in $b$, which motivates regularization.

== Regularization

To suppress noise, we penalize non-smooth reconstructions. $L$ takes differences between neighbouring pixels, so $norm(L x)_2^2$ is large when neighbours differ strongly. We minimize

$
  f(x) = underbrace(norm(A x - b)_2^2, "fit the observed data") + underbrace(lambda^2 norm(L x)_2^2, "encourage smoothness")
$ <eq:tik>

where $lambda > 0$ balances the two terms: too small leaves noise, too large loses details and edges.

#pagebreak()
== The normal equations
We write out $f$ and set its gradient to zero:
$
  f(x) & = norm(A x - b)_2^2 + lambda^2 norm(L x)_2^2 = (A x - b)^T (A x - b) + lambda^2 (L x)^T (L x) \
  & = x^T (A^T A + lambda^2 L^T L) x - 2 x^T A^T b + b^T b \
  mid(arrow.r.double) nabla f(x) & = 2 (A^T A + lambda^2 L^T L) x - 2 A^T b = 0
$ <eq:f-expanded>
This gives the linear system
$
  M x = c, quad M = A^T A + lambda^2 L^T L, quad c = A^T b
$ <eq:normal>
This is the system we solve in the rest of the report.

== Symmetric positive definite matrix
Let $A in RR^(n times n)$ and $L in RR^(p times n)$.
Assume that $lambda != 0$ and

$
  "rank"( mat(A; L) ) = n
$

Using
$M = A^T A + lambda^2 L^T L$ we want to prove that $M$ is symmetric positive definite.
A matrix is symmetric if $M^T = M$. Taking the transpose gives

$
  M^T & = (A^T A + lambda^2 L^T L)^T && =(A^T A)^T + lambda^2 (L^T L)^T \
      & = A^T A + lambda^2 L^T L     && = M
$

Therefore, $M$ is symmetric. A symmetric matrix is positive definite if

$
  z^T M z > 0
  quad "for all" z != 0
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
  >= 0
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


== Iterative Krylov subspace methods
An iterative method starts from a guess $x_0$ (we use $x_0 = 0$, a black image) and improves it step by step, $x_1, x_2, dots$. To see how close $x_k$ is to solving $M x = c$ we use the residual.
$
  r_k = c - M x_k
$
which is zero exactly at the solution. We stop when the relative residual $norm(r_k)_2 \/ norm(c)_2$ is below a tolerance or after a maximum number of iterations.

In the first step the only direction we know is the residual $r_0$, so $x_1 = x_0 + alpha_0 r_0$ for some number $alpha_0$. In Krylov subspace methods a new direction is always found for each iteration, and the step length in that direction is chosen so $x_k$ is the best possible approximation that the method can find in $k$ steps. The subspace is built from $r_0$ and its products with $M$: first $M r_0$, then $M (M r_0) = M^2 r_0$, and so on. $x_k - x_0$ is therefore a linear combination of
$
  cal(K)_k (M, r_0) = "span"{r_0, M r_0, dots, M^(k-1) r_0}
$

The two methods used in this project are Conjugate Gradient (CG) and Minimum Residual (MINRES).

== Conjugate Gradient
Since $M$ is SPD, solving $M x = c$ is the same as minimizing
$
  phi.alt(x) = 1/2 x^T M x - c^T x
$

We can use the update rule
$
  x_(k+1) = x_k + alpha_k p_k
$
where $p_k$ is the search direction and $alpha_k$ is the step length. We choose $p_k$ to be $M$-orthogonal to all previous search directions. $p_0$ is chosen to be the initial residual $r_0 = c - M x_0$. The search direction is updated according to
$ p_(k+1) = r_(k+1) + beta_k p_k $ where
$
  beta_k = (r_(k+1)^T r_(k+1)) / (r_k^T r_k)
$
and
$
  alpha_k = (r_k^T r_k) / (p_k^T M p_k)
$

When M is SPD, CG converges in at most n iterations, and the error decreases monotonically. The convergence rate depends on the condition number of M, which is improved by regularization.

In this project, CG is implemented in Python, adapted from the algorithm in @nla_iter.


== MINRES
MINRES (minimal residual) solves $M x = c$ for symmetric $M$, which unlike in CG does not need to be positive definite. In iteration $k$ it chooses
$
  x_k = op("arg min", limits: #true)_(x in x_0 + cal(K)_k (M, r_0)) norm(c - M x)_2
$
the iterate with the smallest residual in the Krylov subspace

Since $M$ is symmetric, the Lanczos process reduces this to a small least squares problem with a tridiagonal matrix. This problem is solved with Givens rotations, which gives the short update
$
  x_(k+1) = x_k + tau_k p_k
$

Like CG, MINRES needs one matrix-vector product per iteration and a fixed number of stored vectors. CG minimizes the error in the $M$-norm, MINRES the residual in the 2-norm.

== LSQR
There are other iterative methods that can be used to solve the problem. LSQR is a third option, it is not implemented or tested in this project. This section can be skipped without missing anything in the results.
LSQR @paige1982lsqr @lsqr_se is an iterative method for least squares problems $min_x norm(hat(A) x - hat(b))_2$. Unlike CG and MINRES it does not use the normal equations, but works on the least squares problem directly.

The regularized problem @eq:tik becomes an ordinary least squares problem if we stack $A$ on top of $lambda L$:
$
  norm(mat(A; lambda L) x - mat(b; 0))_2^2 = norm(A x - b)_2^2 + lambda^2 norm(L x)_2^2
$
So we apply LSQR with $hat(A) = mat(A; lambda L)$ and $hat(b) = mat(b; 0)$. In each iteration LSQR only multiplies a vector by $hat(A)$ and by $hat(A)^T$:
$
  hat(A) v = mat(A v; lambda L v), quad hat(A)^T mat(u_1; u_2) = A^T u_1 + lambda L^T u_2
$
These are sparse matrix-vector products so no matrix-matrix product is needed.

LSQR and CG solve the same problem, since $M = hat(A)^T hat(A)$ and $c = hat(A)^T hat(b)$, and in exact arithmetic they give the same iterates. Both use four sparse products per iteration and store only a few vectors. The difference is the rounding errors. CG works with $M$, and $kappa(M) = kappa(hat(A))^2$, so the errors are amplified more. LSQR works with $hat(A)$, which makes it more accurate for ill-conditioned problems, that is for small $lambda$. The disadvantage is that LSQR is more complicated and does slightly more work per iteration, so for a well-conditioned problem CG is just as good and simpler @paige1982lsqr.

#pagebreak()
== Blur matrix $A$
To blur the image, the vector $x_"true"$ is multiplied by a blurring matrix $A$.
$A$ comes from the continuous blur model (1.8)
$
  b(s) = integral K(s, t) x(t) dif t + e(s)
$
where $b(s)$ is the blurred image and $e(s)$ is the noise, which in this case is Gaussian.
$K(s,t)$ is called the _point spread function_ and determines how much each pixel is affected by every other pixel in the image. We use the Gaussian kernel (1.9) with $B = beta I$, where $beta$ controls how strongly pixels far away contribute:
$
  K(s, t) = c e^(-beta norm(s - t)_2^2)
$

To make sure that $K(s,t)$ does not change the overall brightness, $c$ is chosen as $c = beta \/ pi$.

Since the image is discrete, the integral is approximated by a sum, where $h$ is the width of each pixel:
$
  (A x)_i = sum_j h^2 K(s_i, t_j) x_j
$

Finally, since $norm(s - t)_2^2 = (s_1 - t_1)^2 + (s_2 - t_2)^2$, the Gaussian kernel factors into one kernel for the $x$-direction and one for the $y$-direction. We can therefore build one-dimensional blur matrices $A_1$ and $A_2$ for each direction, and the Kronecker product gives $A$ for the whole image:
$
  A = A_2 times.o A_1
$
The order of the factors follows from $x = "vec"(X)$ stacking the columns of $X$: $A_1 in RR^(n_1 times n_1)$ blurs every column (vertical direction) and $A_2 in RR^(n_2 times n_2)$ every row (horizontal direction), so $A x = "vec"(A_1 X A_2^T)$.

$A$ is now a dense #thousands(t6A.rows) $times$ #thousands(t6A.cols) matrix, which requires #mem(t6A.dense_bytes) of memory. This is unreasonably large, and pixels far apart barely affect each other anyway. So instead of storing $A$ as a dense matrix, entries below a threshold ($epsilon$) can be set to zero and the matrix can be stored in CSR (compressed sparse row) format.

#pagebreak()
== The regularization matrix $L$
To penalize images that are not smooth, we need a matrix $L$ that measures how much neighbouring pixels differ. We start in one dimension, where the first-order difference matrix
$
  L_1 = mat(1, -1; , 1, -1; , , dots.down, dots.down; , , , 1, -1) in RR^((n_1 - 1) times n_1), quad (L_1 v)_i = v_i - v_(i+1)
$
computes the difference between each pair of neighbouring entries of a vector $v$.

Since $x_"true"$ stacks the columns of $X$, applying $L_1$ to every column gives all vertical differences, likewise $L_2$ applied to every row gives all horizontal differences. Adding them together gives L:
$
  L = L_G = mat(I_(n_2) times.o L_1; L_2 times.o I_(n_1)) in RR^(#thousands(t6L.rows) times #thousands(t6L.cols))
$
Given that $L x = 0$ is only attained when all pixels have the same value, and $A x = 0$ when we have the original image exactly we can conclude that $"rank" mat(A; L) = n$

== Sparsity pattern of $A$
To see how much we gain from storing $A$ as a sparse matrix, we look at its sparsity pattern, that is, where its nonzero entries are (@fig:task6-spy-A). With $beta = #t6beta$ and the threshold $epsilon = #t6eps$, $A$ has #sci(t6A.nnz) nonzeros, which is only $#num(percent(t6A)) %$ of all entries. Stored in CSR format, $A$ therefore needs #mem(t6A.sparse_bytes) instead of #mem(t6A.dense_bytes) as a dense matrix.

#figure(
  image("assets/task6_1_spy.svg", width: 85%),
  caption: [Sparsity pattern of $A$ ($beta = #t6beta$, $epsilon = #t6eps$). Left: all of $A$. Right: the first $3 n_1$ rows and columns, showing the banded blocks of $A_2 times.o A_1$.],
) <fig:task6-spy-A>

#pagebreak()
== Sparsity pattern of $L$
$L$ is even sparser than $A$ (@fig:task6-spy-L). Since every row has exactly two nonzeros, $L$ has #thousands(t6L.nnz) nonzeros, which is only $#num(percent(t6L), digits: 2) %$ of all entries. In CSR format, $L$ needs #mem(t6L.sparse_bytes) instead of #mem(t6L.dense_bytes) as a dense matrix.

#figure(
  image("assets/task6_2_spy.svg", width: 100%),
  caption: [Sparsity pattern of $L$. Left: all of $L$. Middle: the first rows of the top half. Right: the first $3 n_1$ rows of the bottom half.],
) <fig:task6-spy-L>

== Influence of $beta$ and $epsilon$
To see how the choice of parameters affects the sparsity, we vary $beta$ and $epsilon$ and count the nonzeros of $A$ (@tab:sparsity and @fig:task6-params). Since $beta$ controls how strongly pixels far away contribute, a larger $beta$ gives a narrower Gaussian. Fewer entries are then above $epsilon$, so $A$ becomes sparser.

A larger $epsilon$ also makes $A$ sparser, since more entries are set to zero. The price is that $A$ becomes a worse approximation of the real blur, because more of the Gaussian is cut away.

#grid(
  columns: (auto, 1fr),
  gutter: 1.5em,
  align: horizon,
  [#figure(
    table(
      columns: 4,
      align: (center, center, right, right),
      stroke: none,
      table.hline(),
      [$beta$], [$epsilon$], [nnz$(A)$], [nonzero %],
      table.hline(stroke: 0.5pt),
      ..t6
        .sweep
        .map(r => (
          thousands(r.beta),
          pow10(r.eps),
          sci(r.nnz),
          num(100 * r.nnz / (t6A.rows * t6A.cols)),
        ))
        .flatten(),
      table.hline(),
    ),
    caption: [Number of nonzeros of $A$ for different $beta$ and $epsilon$.],
  ) <tab:sparsity>],

  [#figure(
    image("assets/task6_1_params.svg", width: 100%),
    caption: [Sparsity pattern of $A$ for different $beta$ (rows) and $epsilon$ (columns).],
  ) <fig:task6-params>],
)

#pagebreak()
= Results
== Influence of the regularization parameter $lambda$ on CG
We use the tolerance $10^(-8)$, at most $300$ iterations, and $lambda in {0, 0.01, 0.05, 0.3, 1}$. The results are shown in @fig:task8-lambda.
With small $lambda$, the relative residual does not shrink much, but the reconstruction error is smaller, except for $lambda$ = 0, where the error grows. For larger $lambda$ the residual shrinks faster, but the error is larger.

#figure(
  stack(
    spacing: 0.5em,
    image("assets/task8_1_images.svg", width: 100%),
    image("assets/task8_1_convergence.svg", width: 100%),
  ),
  caption: [Influence of $lambda$ on the CG reconstruction and convergence.],
) <fig:task8-lambda>

== Influence of the stopping tolerance on CG
We fix $lambda = 0.05$ and at most $300$ iterations, and vary the tolerance. The results are shown in @fig:task8-tol. We see that the error and residual decrease with the tolerance, but the error levels off at for tolerances below $0.001$. This means it is unnecessary to use a lower tolerance.

#figure(
  stack(
    spacing: 0.5em,
    image("assets/task8_2_images.svg", width: 100%),
    image("assets/task8_2_convergence.svg", width: 100%),
  ),
  caption: [Influence of the tolerance on the CG reconstruction and convergence ($lambda = 0.05$).],
) <fig:task8-tol>

== The maximum number of iterations for CG
We use a weak regularization, $lambda = 0.01$, and tolerance $0$, so CG always runs exactly `max_iter` iterations The results are shown in @fig:task8-maxiter. We see a similar pattern as for the tolerance: the error decreases with the number of iterations, but levels off for `max_iter` higher than 20.

#figure(
  stack(
    spacing: 0.5em,
    image("assets/task8_3_images.svg", width: 100%),
    image("assets/task8_3_convergence.svg", width: 100%),
  ),
  caption: [Influence of the maximum number of iterations on the CG reconstruction and convergence ($lambda = 0.01$).],
) <fig:task8-maxiter>

== Comparison of CG and MINRES
CG and MINRES are compared in @fig:task9-conv and @fig:task9-relres.
#figure(
  stack(
    spacing: 0.5em,
    image("assets/task9_1_images.svg", width: 100%),
    image("assets/task9_1_convergence.svg", width: 100%),
  ),
  caption: [CG and MINRES for $lambda = 0.05$, tolerance $10^(-8)$, at most $300$ iterations.],
) <fig:task9-conv>

#figure(
  image("assets/task9_2_relres.svg", width: 100%),
  caption: [Relative residual of CG and MINRES for three values of $lambda$.],
) <fig:task9-relres>

#figure(
  table(
    columns: 8,
    align: (center, left, right, right, right, right, right, right),
    stroke: none,
    table.hline(),
    [$lambda$], [*Method*], [*Iter.*], [*$M v$*], [*Sparse mat-vec*], [*Time [s]*], [*Rel. res.*], [*Error*],
    table.hline(stroke: 0.5pt),
    [$0.01$], [CG], [$247$], [$248$], [$992$], [$2.68$], [$9.6 dot 10^(-9)$], [$0.032$],
    [], [MINRES], [$217$], [$218$], [$872$], [$2.21$], [$9.9 dot 10^(-9)$], [$0.032$],
    [$0.05$], [CG], [$66$], [$67$], [$268$], [$0.66$], [$9.8 dot 10^(-9)$], [$0.036$],
    [], [MINRES], [$64$], [$65$], [$260$], [$0.7$], [$8.3 dot 10^(-9)$], [$0.036$],
    [$0.3$], [CG], [$18$], [$19$], [$76$], [$0.18$], [$3.5 dot 10^(-9)$], [$0.07$],
    [], [MINRES], [$17$], [$18$], [$72$], [$0.16$], [$9.5 dot 10^(-9)$], [$0.07$],
    table.hline(),
  ),
  caption: [CG and MINRES with tolerance $10^(-8)$ and at most $300$ iterations: iterations, products with $M$, sparse matrix-vector products, time, final relative residual and relative reconstruction error.],
) <tab:task9>

== Why CG is expected to outperform MINRES
$M$ is symmetric positive definite, which is exactly the case CG is designed for. MINRES is designed for symmetric indefinite systems, which is a more general case. In every iteration CG gives the smallest possible value of the regularized least squares objective $f$ in the Krylov space. MINRES instead minimizes $norm(c - M x)_2$. That makes CG a better fit here.

== Preconditioned CG
The idea of preconditioning a problem is to turn it into a similar problem with a smaller condition number. To precondition for CG, the matrix has to stay SPD. A way of achieving this is with an incomplete Cholesky factorization $P = C C^T approx M$, where $C$ is sparse and lower triangular, we get the SPD system
$
  C^(-1) M C^(-T) tilde(x) = C^(-1) c, quad tilde(x) = C^T x.
$
Applying CG to this system with $tilde(M) = C^(-1) M C^(-T)$ gives
$
  alpha_k = (tilde(r)_k^T tilde(r)_k) / (tilde(p)_k^T tilde(M) tilde(p)_k), quad
  tilde(x)_(k+1) = tilde(x)_k + alpha_k tilde(p)_k, quad
  tilde(r)_(k+1) = tilde(r)_k - alpha_k tilde(M) tilde(p)_k, \
  beta_k = (tilde(r)_(k+1)^T tilde(r)_(k+1)) / (tilde(r)_k^T tilde(r)_k), quad
  tilde(p)_(k+1) = tilde(r)_(k+1) + beta_k tilde(p)_k,
$
and the solution is recovered from $x = C^(-T) tilde(x)$. In practice, $tilde(M)$ is never formed. With $r = C tilde(r)$ and $p = C^(-T) tilde(p)$, the same iteration can be carried out in the original variables, where every iteration only needs one product with $M$ and one solve with $P$ (two triangular solves).

We compare CG and PCG, both starting from $x_0 = 0$ and stopping when the relative residual is below $10^(-6)$ or after $500$ iterations (@fig:task10, @tab:task10).

#figure(
  image("assets/task10_cg_vs_pcg.svg", width: 100%),
  caption: [Relative residual of CG and PCG (incomplete Cholesky) for $lambda = 0.01, 0.05, 0.1$.],
) <fig:task10>

#figure(
  table(
    columns: 6,
    align: (center, left, right, right, right, right),
    stroke: none,
    table.hline(),
    [$lambda$], [*Method*], [*Iter.*], [*Time [s]*], [*Rel. res.*], [*Error*],
    table.hline(stroke: 0.5pt),
    [$0.01$], [CG], [$138$], [$1.66$], [$9.9 dot 10^(-7)$], [$0.0316$],
    [], [PCG], [$13$], [$1.99 + 15.6$], [$6.8 dot 10^(-7)$], [$0.0316$],
    [$0.05$], [CG], [$43$], [$0.58$], [$9.7 dot 10^(-7)$], [$0.0365$],
    [], [PCG], [$4$], [$0.57 + 12.6$], [$1.2 dot 10^(-7)$], [$0.0365$],
    [$0.1$], [CG], [$26$], [$0.35$], [$9.8 dot 10^(-7)$], [$0.0461$],
    [], [PCG], [$3$], [$0.47 + 12.7$], [$2.3 dot 10^(-7)$], [$0.0461$],
    table.hline(),
  ),
  caption: [CG and PCG with tolerance $10^(-6)$: iterations, time (for PCG: solve + building the incomplete Cholesky factor), final relative residual and relative reconstruction error.],
) <tab:task10>

@tab:task10 shows that PCG needs about $10$ times fewer iterations than CG for all three values of $lambda$. Despite this, PCG is not faster: the solve itself takes about as long as CG, since each iteration is much more expensive, and building the incomplete Cholesky factor adds time. Both methods reach the tolerance and give the same reconstruction error, since they solve the same system to the same tolerance.

= Summary
For solving inverse problems, direct methods and iterative Krylov methods have significant differences in terms of computational complexity, memory requirements, numerical stability, scalability, and suitability for small and large problems. Direct methods require significant memory and computational resources, using dense matrix operations, especially for large problems. However, they are often more numerically stable and can provide better solutions. Krylov subspace methods are more memory-efficient and scalable, making them suitable for large problems, but they may be less stable and require careful tuning of parameters. Direct methods are advantageous for small problems due to their accuracy, while iterative methods are preferred in large-scale applications.


// ---------------------------------------------------------------------------
//  References
// ---------------------------------------------------------------------------
#pagebreak()
#set heading(numbering: none)
#bibliography("references.bib", title: "References")


