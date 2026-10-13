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
  if m.starts-with("10") { e += 1; m = sigs(1.0, digits: digits) }
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
In this project we study iterative methods for solving linear systems. We apply these methods to to recover a sharp image from a blurred one. The task boils down to a very large least squares problem, and with a dense matrix of that size the direct methods introduced earlier become too expensive. We therefore make the matrix artificially sparse and use iterative methods that are built for sparse problems. The deblurring problem is, however, ill-posed, so small amounts of noise can ruin the solution, and we show how regularization turns it into a solvable one.


#pagebreak()
= Method
Direct methods like GSVD need too much storage and computation for large image problems. We instead use iterative methods, which refine the solution step by step using mainly matrix-vector products and a few stored vectors.

== Images as vectors
#grid(
  columns: (2fr, 1fr),
  gutter: 1em,
  align: top,

  [
    A grayscale image is a matrix $X in RR^(n_1 times n_2)$ of pixel brightnesses from $0$ (black) to $1$ (white). Stacking its columns gives

    $x_"true" = "vec"(X) in RR^n, quad n = n_1 n_2.$

    We resize the image to $300 times 300$ pixels, so $n = 90 000$. The goal is to recover the sharp image $x_"true"$ from the blurred and noisy observation

    $b = A x_"true" + e,$ <eq:model>

    where $A$ is the blurring matrix and $e$ is the noise.
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
  min_x norm(A x - b)_2^2,
$

which minimizes the residual $A x - b$, the difference between the predicted blurred image and $b$. Least squares alone also fits the noise in $b$, which motivates regularization.

== Regularization

To suppress noise, we penalize non-smooth reconstructions. $L$ takes differences between neighbouring pixels, so $norm(L x)_2^2$ is large when neighbours differ strongly. We minimize

$
  f(x) = underbrace(norm(A x - b)_2^2, "fit the observed data") + underbrace(lambda^2 norm(L x)_2^2, "encourage smoothness"),
$ <eq:tik>

where $lambda > 0$ balances the two terms: too small leaves noise, too large loses details and edges.

== The normal equations
We write out $f$ and set its gradient to zero:
$
    f(x) & = norm(A x - b)_2^2 + lambda^2 norm(L x)_2^2 = (A x - b)^T (A x - b) + lambda^2 (L x)^T (L x) \
             & = x^T (A^T A + lambda^2 L^T L) x - 2 x^T A^T b + b^T b \
  mid(arrow.r.double) nabla f(x) & = 2 (A^T A + lambda^2 L^T L) x - 2 A^T b = 0.
$ <eq:f-expanded>
This gives the linear system
$
  M x = c, quad M = A^T A + lambda^2 L^T L, quad c = A^T b.
$ <eq:normal>
This is the system we solve in the rest of the report.

== Symmetric positive definite matrix
Let $A in RR^(m times n)$ and $L in RR^(p times n)$.
Assume that $lambda != 0$ and

$
  "rank"( mat(A; L) ) = n.
$

Using the matrix

$
  M = A^T A + lambda^2 L^T L.
$

we want to prove that $M$ is symmetric positive definite.

A matrix is symmetric if $M^T = M$. Taking the transpose gives

$
  M^T & = (A^T A + lambda^2 L^T L)^T \
      & =(A^T A)^T + lambda^2 (L^T L)^T \
      & = A^T A + lambda^2 L^T L = M.
$

Therefore, $M$ is symmetric.

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


== Iterative Krylov subspace methods
An iterative method starts from a guess $x_0$ (we use $x_0 = 0$, a black image) and improves it step by step, $x_1, x_2, dots$. To see how close $x_k$ is to solving $M x = c$ we use the residual.
$
  r_k = c - M x_k,
$
which is zero exactly at the solution. We stop when the relative residual $norm(r_k)_2 \/ norm(c)_2$ is below a tolerance $epsilon$ or after a maximum number of iterations.

In the first step the only direction we know is the residual $r_0$, so $x_1 = x_0 + alpha_0 r_0$ for some number $alpha_0$. In Krylov subspace methods a new direction is always found for each iteration, and the step length in that direction is chosen so $x_m$ is the best possible approximation that the method can find in $m$ steps. The subspace is built from $r_0$ and its products with $M$: first $M r_0$, then $M (M r_0) = M^2 r_0$, and so on. $x_m - x_0$ is therefore a linear combination of
$
  cal(K)_m (M, r_0) = "span"{r_0, M r_0, dots, M^(m-1) r_0},
$

The two methods used in this project are Conjugate Gradient (CG) and Minimum Residual (MINRES).

== Conjugate Gradient
Since $M$ is SPD, solving $M x = c$ is the same as minimizing
$
  phi.alt(x) = 1/2 x^T M x - c^T x.
$

We can use the update rule
$
  x_(m+1) = x_m + alpha_m p_m
$
where $p_m$ is the search direction and $alpha_m$ is the step length. We choose $p_m$ to be A-orthogonal to all previous search directions. $p_0$ is chosen to be the initial residual $r_0 = c - M x_0$. The search direction is updated according to
$ p_(m+1) = r_(m+1) + beta_m p_m $ where
$
  beta_m = (r_(m+1)^T r_(m+1)) / (r_m^T r_m)
$
and
$
  alpha_m = (r_m^T r_m) / (p_m^T M p_m).
$

When M is SPD, CG converges in at most n iterations, and the error decreases monotonically. The convergence rate depends on the condition number of M, which is improved by regularization.

In this project, CG is implemented in Python, adapted from the algorithm in @nla_iter.


== MINRES
MINRES (minimal residual) solves $M x = c$ for symmetric $M$, which unlike in CG does not need to be positive definite. In iteration $m$ it chooses
$
  x_m = op("arg min", limits: #true)_(x in x_0 + cal(K)_m (M, r_0)) norm(c - M x)_2,
$
the iterate with the smallest residual in the Krylov subspace

Since $M$ is symmetric, the Lanczos process reduces this to a small least squares problem with a tridiagonal matrix. This problem is solved with Givens rotations, which gives the short update
$
  x_(m+1) = x_m + tau_m p_m.
$

Like CG, MINRES needs one matrix-vector product per iteration and a fixed number of stored vectors. CG minimizes the error in the $M$-norm, MINRES the residual in the 2-norm.

== Blur matrix $A$
To blur the image the vector $x_"true"$ is multiplied by a blurring matrix $A$. 
$A$ is given by the continuous blur model (1.8)
$
  b(s) = integral K(s, t) x(t) dif t + e(s),
$
where $b(s)$ is the final blurred image, $e(s)$ represents the noise which is gaussian noise in this case.
$K(s,t)$ is called _point spread function_ and weighs how each pixel is affected by all other pixels in the image, we use the Gaussian kernel (1.9), with $B = beta I$, where $beta$ controls how much pixels far away affect:
$
  K(s, t) = c e^(-beta norm(s - t)_2^2).
$

To make sure $K(s,t)$ does not change the overall brightness, c is chosen as $c = beta \/ pi$

As the problem is discrete the integral can be approximated by the sum, where $h$ represents the width of each pixel

$
  (A x)_i = sum_j h^2 K(s_i, t_j) x_j
$

Finnaly as the problem is 2 dimensional and our A now works in one dimension, a kronecker product between the x direction and y direction gives us a A for the whole image.$
  A = A_2 times.o A_1.
$

A is now a dense $90 thin 000 times 90 thin 000$

*Thresholding.* The Gaussian is never exactly zero, so $A$ is a dense  matrix with $8.1 dot 10^9$ entries, about $65$ GB. But the weights of pixels far apart are tiny: for $beta = #t6beta$, pixels more than about #num(t6.params.radius_px, digits: 1) pixels apart get a weight below #t6eps. We therefore set all entries $|a_(i j)| < epsilon$ to zero, so every pixel only depends on its neighbours and $A$ becomes sparse.


*Sparse storage.* We store $A$ in the CSR (compressed sparse row) format, which only stores the nonzero entries: for each nonzero its value ($8$ bytes) and its column index ($4$ bytes), plus one pointer per row ($4$ bytes) telling where the row starts.

== The data $b$
The data is generated by !!!!eq:model, $b = A x_"true" + e$, where every entry of $e$ is an independent standard normal number times $sigma = #num(t6.params.sigma, digits: 1)$. Since the pixel values lie in $[0, 1]$, this is noise of about $#num(100 * t6.params.sigma, digits: 1) %$ of the brightness range. The true image and $b$ are shown in @fig:task8-lambda.

== The regularization matrix $L$
In one dimension, the first-order difference matrix
$
  L_1 = mat(1, -1; , 1, -1; , , dots.down, dots.down; , , , 1, -1) in RR^((n_1 - 1) times n_1), quad (L_1 v)_i = v_i - v_(i+1),
$
computes the differences between neighbouring entries of a vector. Applied to every column of $X$, $L_1 X$ contains all vertical differences, and $"vec"(L_1 X) = (I_(n_2) times.o L_1) x$. In the same way, $X L_2^T$ contains all horizontal differences, and $"vec"(X L_2^T) = (L_2 times.o I_(n_1)) x$. Stacking the two gives the discrete gradient
$
  L = L_G = mat(I_(n_2) times.o L_1; L_2 times.o I_(n_1)) in RR^(179 thin 400 times 90 thin 000).
$
So $norm(L x)_2^2$ is the sum of the squared differences over all pairs of vertically and horizontally neighbouring pixels. This is the smoothness penalty from the introduction. Every row of $L$ is one such difference and has exactly two nonzeros, $1$ and $-1$.

Note that $L x = 0$ only when all pixels have the same value (a constant image), while $A$ does not map a constant image to zero, since all entries of $A$ are non-negative and every row has nonzeros. Hence $A z = 0$ and $L z = 0$ can only both hold for $z = 0$, which means that $"rank" mat(A; L) = n$. This is the assumption we need in Task 7.

== Sparsity pattern of $A$
@fig:task6-spy-A shows the sparsity pattern of $A$, that is, where its nonzero entries are. $A$ has #sci(t6A.nnz) nonzeros, which is $#num(percent(t6A)) %$ of all entries. On average each row has #num(t6A.nnz / t6A.rows, digits: 2) nonzeros out of #thousands(t6A.cols): every pixel of the blurred image depends on about #num(t6A.nnz / t6A.rows, digits: 2) pixels of the sharp image, roughly a disk with a radius of #num(t6.params.radius_px, digits: 1) pixels. The pattern is a band of blocks: each block is a band (from $A_1$), and the blocks themselves form a band (from $A_2$).

#figure(
  image("assets/task6_1_spy.svg", width: 85%),
  caption: [Sparsity pattern of $A$ ($beta = #t6beta$, $epsilon = #t6eps$). Left: all of $A$. Right: the first $3 n_1$ rows and columns, showing the banded blocks of $A_2 times.o A_1$.],
) <fig:task6-spy-A>

== Sparsity pattern of $L$
$L$ has #thousands(t6L.nnz) nonzeros ($2$ per row), which is only $#num(percent(t6L), digits: 2) %$ of all entries (@fig:task6-spy-L). In the top half ($I_(n_2) times.o L_1$, vertical neighbours) the two nonzeros of a row are in neighbouring columns $i$ and $i + 1$. In the bottom half ($L_2 times.o I_(n_1)$, horizontal neighbours) they are $n_1 = 300$ columns apart, since horizontal neighbours are $n_1$ positions apart in $x$.

#figure(
  image("assets/task6_2_spy.svg", width: 100%),
  caption: [Sparsity pattern of $L$. Left: all of $L$. Middle: the first rows of the top half. Right: the first $3 n_1$ rows of the bottom half.],
) <fig:task6-spy-L>

== Memory
A dense $m times n$ matrix stores every entry as an $8$ byte number, so it needs $8 m n$ bytes. In CSR format a matrix needs $12 thin "nnz" + 4 (m + 1)$ bytes, where nnz is the number of nonzeros. @tab:memory shows that dense storage is impossible on a normal computer, while the sparse matrices fit easily.

#figure(
  table(
    columns: 4,
    align: (left, right, right, right),
    stroke: none,
    table.hline(),
    [*Matrix*], [*Size*], [*Dense*], [*Sparse (CSR)*],
    table.hline(stroke: 0.5pt),
    ..(([$A$], t6A), ([$L$], t6L)).map(((name, M)) => (
      name, [#thousands(M.rows) $times$ #thousands(M.cols)], mem(M.dense_bytes), mem(M.sparse_bytes),
    )).flatten(),
    table.hline(),
  ),
  caption: [Memory needed for dense and sparse storage of $A$ and $L$ ($beta = #t6beta$, $epsilon = #t6eps$).],
) <tab:memory>

== Influence of $beta$ and $epsilon$
@tab:sparsity and @fig:task6-params show how the number of nonzeros of $A$ depends on $beta$ and $epsilon$.
- A larger $beta$ gives a narrower Gaussian, so fewer entries are above $epsilon$ and $A$ becomes sparser.
- A larger $epsilon$ removes more entries, so $A$ also becomes sparser. The price is that $A$ is a worse approximation of the real blur, since we cut away more of the tails of the Gaussian.

The memory and the cost of a matrix-vector product ($approx 2 thin "nnz"$ operations) both grow linearly with nnz. For our parameters, one product with $A$ costs about #sci(2 * t6A.nnz, digits: 2) operations instead of #sci(2 * t6A.rows * t6A.cols, digits: 2) for a dense $A$, about #thousands(t6A.rows * t6A.cols / t6A.nnz) times fewer. Since every iteration in Tasks 8–9 uses four sparse products, this is what makes the iterative methods fast. $L$ does not depend on $beta$ or $epsilon$, and it is so sparse that its cost is negligible compared to $A$.

#figure(
  table(
    columns: 4,
    align: (center, center, right, right),
    stroke: none,
    table.hline(),
    [$beta$], [$epsilon$], [nnz$(A)$], [nonzero %],
    table.hline(stroke: 0.5pt),
    ..t6.sweep.map(r => (
      thousands(r.beta), pow10(r.eps), sci(r.nnz), num(100 * r.nnz / (t6A.rows * t6A.cols)),
    )).flatten(),
    table.hline(),
  ),
  caption: [Number of nonzeros of $A$ for different $beta$ and $epsilon$.],
) <tab:sparsity>

#figure(
  image("assets/task6_1_params.svg", width: 80%),
  caption: [Sparsity pattern of $A$ for different $beta$ (rows) and $epsilon$ (columns).],
) <fig:task6-params>

== CG results (från AI)
== The regularization parameter $lambda$
We use the tolerance $10^(-8)$, at most $300$ iterations, and $lambda in {0, 0.01, 0.05, 0.3, 1}$ (@fig:task8-lambda).

#figure(
  stack(
    spacing: 0.5em,
    image("8.1.png", width: 100%),
    image("8.1.graph.png", width: 100%),
  ),
  caption: [Influence of $lambda$ on the CG reconstruction and convergence.],
) <fig:task8-lambda>

- *$lambda = 0$ (no regularization).* This is plain least squares. The residual stalls around $10^(-4)$ and the error grows to $1.42$. CG starts to fit the noise, exactly as predicted in the introduction, and the image is destroyed.
- *$lambda = 0.01$.* The image is better, but there is still ringing noise, and CG needs all $300$ iterations because $M$ is still badly conditioned.
- *$lambda = 0.05$* gives the smallest error ($0.126$) of the values we tried. This is the compromise between noise and smoothing.
- *$lambda = 0.3$ and $lambda = 1$.* The image is too smooth and the edges are blurred, so the error grows again (over-regularization).
- *Speed.* A larger $lambda$ makes $M$ better conditioned, so CG converges much faster ($39$–$47$ iterations instead of more than $300$), as the bound above predicts.
- The final objective is larger for larger $lambda$, simply because the penalty term $lambda^2 norm(L x)_2^2$ is weighted more.

== The stopping tolerance
We fix $lambda = 0.05$ and at most $300$ iterations, and vary the tolerance (@fig:task8-tol).

#figure(
  stack(
    spacing: 0.5em,
    image("8.2.png", width: 100%),
    image("8.2.graph.png", width: 100%),
  ),
  caption: [Influence of the tolerance on the CG reconstruction and convergence ($lambda = 0.05$).],
) <fig:task8-tol>

CG is deterministic and always starts from $x_0 = 0$, so all runs follow the same iterates. The tolerance only decides where a run stops (the dots in the plots).
- `tol` $= 10^(-2)$ stops after $4$ iterations with a blurry image (error $0.155$).
- From `tol` $= 10^(-4)$ ($38$ iterations) the error stays at $0.126$. Smaller tolerances cost more iterations ($89$ and $140$) without making the image better.
- The error and the objective level off after about $30$ iterations, so a moderate tolerance is enough. Solving the normal equations to high accuracy is pointless when the equations themselves only describe the image up to the noise.

== The maximum number of iterations
We use a weak regularization, $lambda = 0.01$, and `tol` $= 0$, so CG always runs exactly `max_iter` iterations (@fig:task8-maxiter).

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
With weak regularization the number of iterations itself acts as a regularizer. This is called _semi-convergence_:
- The error first decreases and is smallest around $20$–$25$ iterations ($0.128$, about as good as $lambda = 0.05$).
- After that the error grows to $0.21$ and ringing artifacts appear, even though the residual and the objective keep decreasing.

The reason is the SVD picture from the introduction. The first Krylov vectors are dominated by the directions with large singular values, which are smooth and carry the real image. The directions with small singular values, where the noise is amplified, only enter the Krylov space after many iterations. Stopping early therefore keeps the noisy directions out of the solution. A small residual does not mean a good reconstruction: with a small $lambda$, stopping early gives a better image.

== Summary of the convergence behaviour
- The objective decreases in every iteration, as it must, since CG minimizes it over a growing space.
- The residual decreases (not always monotonically), and the rate is set by $kappa(M)$: a larger $lambda$ gives faster convergence.
- The reconstruction error is what we actually care about. It levels off (good $lambda$) or even grows again (small $lambda$), so the best reconstruction is obtained by a good $lambda$ combined with a moderate tolerance.


== Experiment setup
Both methods start from $x_0 = 0$, solve $M x = c$ with the same operator, and stop when the relative residual $norm(c - M x_k)_2 \/ norm(c)_2$ is below the tolerance, so the stopping criteria are identical. We count the products with $M$ by wrapping the operator in a counter. Every product with $M$ is four sparse products ($A$, $A^T$, $L$, $L^T$). The timed runs record nothing during the iterations (computing the error and objective would dominate the time), and each run is timed three times, keeping the fastest.

== Results
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
    [$lambda$], [*Method*], [*Iter.*], [*$M v$*], [*Sparse $m v$*], [*Time [s]*], [*Rel. res.*], [*Error*],
    table.hline(stroke: 0.5pt),
    [$0.01$], [CG], [$300$], [$301$], [$1204$], [$15.9$], [$1.2 dot 10^(-6)$], [$0.212$],
    [], [MINRES], [$300$], [$301$], [$1204$], [$15.1$], [$2.5 dot 10^(-7)$], [$0.211$],
    [$0.05$], [CG], [$140$], [$141$], [$564$], [$11.7$], [$9.4 dot 10^(-9)$], [$0.126$],
    [], [MINRES], [$130$], [$131$], [$524$], [$9.1$], [$9.2 dot 10^(-9)$], [$0.126$],
    [$0.3$], [CG], [$39$], [$40$], [$160$], [$2.5$], [$6.9 dot 10^(-9)$], [$0.140$],
    [], [MINRES], [$38$], [$39$], [$156$], [$2.6$], [$7.7 dot 10^(-9)$], [$0.140$],
    table.hline(),
  ),
  caption: [CG and MINRES with tolerance $10^(-8)$ and at most $300$ iterations: iterations, products with $M$, sparse matrix-vector products, time, final relative residual and relative reconstruction error.],
) <tab:task9>

*Convergence rate.* MINRES reaches the tolerance in slightly fewer iterations than CG ($130$ instead of $140$ for $lambda = 0.05$, $38$ instead of $39$ for $lambda = 0.3$), and for $lambda = 0.01$ it ends with a $5$ times smaller residual after $300$ iterations (@fig:task9-relres, @tab:task9). This is expected: MINRES minimizes exactly the residual we use to stop, while the CG residual is not minimized and decreases less smoothly. The reconstruction error and the objective, however, are almost the same for both methods (@fig:task9-conv). In the first $20$–$30$ iterations CG even has a slightly smaller error, and both end at the same error. The extra residual reduction of MINRES does not give a better image.

*Matrix-vector products.* Both methods use one product with $M$ per iteration plus one for the initial residual, so $k + 1$ products with $M$, or $4(k + 1)$ sparse products, for $k$ iterations. The difference only comes from the number of iterations.

*Computational time.* The time per iteration is about the same for both methods, since the four sparse products with $A$, $A^T$ ($3 dot 10^7$ nonzeros) and $L$, $L^T$ dominate the cost, not the vector operations. The measured times varied between runs by up to about $30 %$ (the computer was not idle), so the small differences in @tab:task9 are not significant.

== Why CG is expected to outperform MINRES here
- *CG is made for this kind of matrix.* $M$ is symmetric positive definite, which is exactly the case CG is designed for. MINRES is designed for symmetric indefinite systems, and we do not need that generality.
- *CG minimizes what we actually want to minimize.* As shown in Task 8, $norm(x - x_*)_M^2 = f(x) - f(x_*)$, so in every iteration CG gives the smallest possible value of the regularized least squares objective $f$ in the Krylov space. MINRES instead minimizes $norm(c - M x)_2$, which is only the size of the gradient of $f$. This is why MINRES wins in the residual plot, while CG is at least as good in the error and the objective.
- *CG is cheaper per iteration.* MINRES does more work per iteration (the Lanczos recurrence and the Givens rotations) and stores more vectors ($v_(k-1), v_k, w, p_(k-1), p_k, x$ instead of $x, r, p, M p$ for CG). For our problem the matrix-vector products dominate, so the difference is small, but CG is never more expensive.
- *Neither method fixes the conditioning.* Both search in the same Krylov space, so both are limited by the same $kappa(M)$, which for small $lambda$ is close to $kappa(A)^2$. This is why a preconditioner is needed to speed up the convergence (Task 10).

= Summary
We wanted to recover a sharp image $x_"true"$ from a blurred, noisy image $b = A x_"true" + e$. Solving $A x = b$ or the plain least squares problem amplifies the noise, because the blur matrix $A$ has many tiny singular values. Adding the regularization term $lambda^2 norm(L x)_2^2$, which penalizes differences between neighbouring pixels, damps the noise, and the minimizer satisfies the normal equations $(A^T A + lambda^2 L^T L) x = A^T b$. For a $300 times 300$ image this is a $90 thin 000 times 90 thin 000$ system, far too large for direct methods, but $A$ and $L$ are sparse (Task 6). The matrix is symmetric positive definite (Task 7), so CG can be used, and it only needs sparse matrix-vector products (Task 8). The best reconstruction is obtained with a good $lambda$ ($approx 0.05$) and a moderate tolerance. With a small $lambda$, stopping early acts as a regularizer. MINRES (Task 9) uses the same number of products per iteration and reaches a small residual slightly faster, but gives the same image, and CG is the natural method for this SPD problem.



= Discussion
// Task 13. Summarize the differences between direct methods (SVD/GSVD) and
// iterative Krylov methods for solving inverse problems. Your discussion
// should include
//   - computational complexity,
//   - memory requirements,
//   - numerical stability,
//   - scalability,
//   - suitability for small and large problems,
//   - advantages and disadvantages of each approach.


= Task 11
SUPERAISVAR:LSQR can solve the regularized least-squares problem
\(\min_x \|Ax-b\|_2^2+\lambda^2\|Lx\|_2^2\) by rewriting it as an ordinary least-squares problem with the augmented matrix \(\widetilde A=\begin{bmatrix}A\\\lambda L\end{bmatrix}\) and vector \(\widetilde b=\begin{bmatrix}b\\0\end{bmatrix}\). It requires only products with \(\widetilde A\) and \(\widetilde A^T\), which can be computed using \(A,A^T,L,L^T\). Therefore, LSQR avoids explicitly forming \(A^TA\) and \(L^TL\), reducing memory requirements and avoiding the additional rounding errors and potential sparsity loss associated with forming the normal equations.
CG can also avoid explicitly forming the normal equations by applying matrix-vector products through \(A^T(Av)+\lambda^2L^T(Lv)\). However, it solves the resulting normal-equation system, whose condition number is the square of that of the augmented matrix: \(\kappa_2(\widetilde A^T\widetilde A)=\kappa_2(\widetilde A)^2\), assuming full column rank. LSQR works directly with the augmented least-squares operator and generally has better numerical behavior for ill-conditioned problems. Both methods have similar per-iteration costs, typically requiring one forward and one transpose matrix-vector product, and can be effective for large sparse systems.
Regularization can improve conditioning by reducing the influence of poorly determined solution components, although it introduces a trade-off between fitting the data and limiting the solution norm. Neither LSQR nor CG guarantees an accurate solution for severely ill-conditioned problems, and both can exhibit semiconvergence, where later iterations begin fitting noise. In practice, LSQR is generally the preferred starting point for large, ill-conditioned least-squares problems, while CG may be appropriate when an effective preconditioner or an existing normal-equation solver is available. Early stopping and suitable regularization are important for obtaining a stable, meaningful solution.


// ---------------------------------------------------------------------------
//  References
// ---------------------------------------------------------------------------
#set heading(numbering: none)
#bibliography("references.bib", title: "References")


