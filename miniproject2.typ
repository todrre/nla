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
== Goal
The goal of part two is to recover a sharp image from a blurred one. The task boils down to a very large least squares problem, and with a dense matrix of that size the direct methods introduced earlier become too expensive. We therefore make the matrix artificially sparse and use iterative methods that are built for sparse problems. The deblurring problem is however ill-posed, so small amounts of noise can ruin the solution, and we show how regularization turns it into a solvable one.


#pagebreak()
= Method

== Iterative methods
An iterative method compared to a direct method (GSVD) solves the problem by iterativly guessing and improving a solution. The main reason to use iterative methods instead of direct, is because sparse-matrix products are much cheaper then dense ones.

== Least squares
Since $b$ contains noise, we do not ask for $A x = b$ exactly. We ask for the $x$ that fits the data best:
$
  min_x norm(A x - b)_2^2.
$
But the best fit to noisy data is again $x = A^(-1) b$, so least squares also fits the noise.

== Regularization
Real images are mostly smooth, noise is not. So we also punish rough images. The matrix $L$ takes differences between neighbouring pixels (Task 6), so $norm(L x)_2$ is large when the image is rough. We minimize
$
  f(x) = underbrace(norm(A x - b)_2^2, "fit the data") + underbrace(lambda^2 norm(L x)_2^2, "be smooth").
$ <eq:tik>
The parameter $lambda > 0$ decides how much smoothness counts: too small gives noise, too large gives blurry edges.

== Images as vectors

A grayscale image is a matrix $X in RR^(n_1 times n_2)$, where each entry represents a pixel brightness between $0$ (black) and $1$ (white). We represent the image as a single vector by stacking its columns:

$
  x_"true" = "vec"(X) in RR^n, quad n = n_1 n_2.
$

In our project, the image is resized to $300 times 300$ pixels, so $n = 90 000$. The goal is to recover this original sharp image from its blurred and noisy version $b$. We model the observation as

$
  b = A x_"true" + e,
$ <eq:model>

where $A$ is the blurring matrix and $e$ is the noise.

== Least squares

Since we observe $b$, but do not know the original image $x_"true"$ or the noise $e$, we want to estimate the sharp image. If there were no noise, we would have $A x = b$. In practice, noise makes an exact solution unreliable, so instead we find the image $x$ whose blurred version $A x$ is as close as possible to the observed image $b$. This gives the least squares problem:

$
  min_x norm(A x - b)_2^2.
$

The quantity $A x - b$ is the residual, and its squared norm measures the difference between the predicted blurred image and the observed image. However, since $b$ contains noise, minimizing the residual can also make the reconstruction fit the noise instead of the true image. Least squares alone does not prevent this, which is why we introduce regularization.

== Regularization

Least squares tries to fit the observed image $b$, but since $b$ contains noise, the reconstructed image $x$ may also contain noise. To reduce this effect, we add a penalty that encourages the reconstructed image to be smooth.

The matrix $L$ measures differences between neighbouring pixels (Task 6), so $norm(L x)_2^2$ is large when neighbouring pixels differ a lot. We therefore minimize

$
  f(x) = underbrace(norm(A x - b)_2^2, "fit the observed data") + underbrace(lambda^2 norm(L x)_2^2, "encourage smoothness").
$ <eq:tik>

The first term is the least squares objective, which makes the blurred reconstruction $A x$ fit the observed image $b$. The second term penalizes roughness in the reconstructed image. The parameter $lambda > 0$ controls the balance between these goals: if $lambda$ is too small, the reconstruction may contain noise; if it is too large, important details and edges may be lost.

== Normal equations
We write out $f$ and set its gradient to zero:
$
        f(x) & = norm(A x - b)_2^2 + lambda^2 norm(L x)_2^2 \
             & = (A x - b)^T (A x - b) + lambda^2 (L x)^T (L x) \
             & = x^T (A^T A + lambda^2 L^T L) x - 2 x^T A^T b + b^T b, \
  nabla f(x) & = 2 (A^T A + lambda^2 L^T L) x - 2 A^T b = 0.
$ <eq:f-expanded>
This gives the linear system
$
  M x = c, quad M = A^T A + lambda^2 L^T L, quad c = A^T b.
$ <eq:normal>
This is the system we solve in the rest of the report.

== Iterative Krylov subspace methods
An iterative method starts from a guess $x_0$ (we use $x_0 = 0$, a black image) and improves it step by step, $x_1, x_2, dots$. To see how close $x_k$ is to solving $M x = c$ we use the _residual_
$
  r_k = c - M x_k,
$
which is zero exactly at the solution. We stop when the relative residual $norm(r_k)_2 \/ norm(c)_2$ is below a tolerance `tol`, or after `max_iter` steps.

In the first step the only direction we know is the residual $r_0$, so $x_1 = x_0 + alpha_0 r_0$ for some number $alpha_0$. In krylov methods a new direction is always found for each iteration. And the length we walk in that direction is chosen optimally, so $x_k$ is the best possible approximation that the method can find in $k$ steps.

For every new step the space in wich we pick the direction grows by one dimension. The standard krylov sub space is build from the residual $r_0$ and its products with $M$: first $M r_0$, then $M (M r_0) = M^2 r_0$, and so on. After $k$ steps, $x_k - x_0$ is therefore a linear combination of
$
  cal(K)_k (M, r_0) = "span"{r_0, M r_0, dots, M^(k-1) r_0},
$

Since every space contains the previous one, the optimal $x_k$ can only get better from one iteration to the next. After at most $n$ iterations the space contains the exact solution, but in practice far fewer are enough. Krylov methods therefore build a better basis of the same space, one vector per iteration, where every new vector is orthogonal to the earlier ones. 

The two methods used in this project are Conjugate Gradient (CG) and Minimum Residual (MINRES).

== Conjugate Gradient
Since $M$ is SPD, solving $M x = c$ is the same as minimizing
$
  phi.alt(x) = 1/2 x^T M x - c^T x.
$

We can use the update rule
$
  x_(m+1) = x_m + alpha_m p_m,
$
where 

In this project, CG is implemented in Python, adapted from the algorithm in @nla_iter.

For CG there is the classical bound
$
  norm(x_k - x_*)_M <= 2 ((sqrt(kappa) - 1) / (sqrt(kappa) + 1))^k norm(x_0 - x_*)_M, quad kappa = kappa(M) = lambda_max (M) / lambda_min (M).
$
If $kappa$ is close to $1$ the error decreases quickly. If $kappa$ is large the factor is close to $1$ and convergence is slow. For $lambda = 0$ we have $kappa(M) = kappa(A^T A) = kappa(A)^2$, which is enormous. The regularization helps here too: $A^T A$ is tiny for rough images (the blur erases them), while $L^T L$ is large exactly for rough images (they have large differences between neighbours). Adding $lambda^2 L^T L$ therefore lifts the smallest eigenvalues of $M$, so $kappa(M)$ decreases and CG converges faster when $lambda$ grows. The two terms complement each other: $L^T L$ is zero for constant images, but $A^T A$ is not small there.

*What "optimal" means* depends on the method:
- *CG* (Task 8) minimizes the error $norm(x - x_*)_M$, which is the same as minimizing $f$. It needs $M$ to be symmetric positive definite.
- *MINRES* (Task 9) minimizes the residual $norm(c - M x)_2$. It only needs $M$ to be symmetric.

= Task 6: Building the problem
In this task we build the concrete $x_"true"$, $A$, $b$ and $L$ for our image and store $A$ and $L$ as sparse matrices. We use the parameters in @tab:params.

#figure(
  table(
    columns: 3,
    align: (left, center, left),
    stroke: none,
    table.hline(),
    [*Parameter*], [*Value*], [*Meaning*],
    table.hline(stroke: 0.5pt),
    [$n_1 = n_2$], [$300$], [image size in pixels],
    [$beta$], [$1000$], [width of the blur (larger $beta$ = less blur)],
    [$epsilon$], [$10^(-4)$], [entries $|a_(i j)| < epsilon$ are set to zero],
    [$sigma$], [$10^(-2)$], [noise level],
    table.hline(),
  ),
  caption: [Parameters used in Tasks 6–9 unless stated otherwise.],
) <tab:params>

== The image
The image from the course page is resized to $300 times 300$ pixels, converted to grayscale and scaled to values in $[0, 1]$. This gives $X_"true"$, and $x_"true" = "vec"(X_"true") in RR^(90 thin 000)$.

== The blur matrix $A$
The continuous blur model (1.8) in the project description is
$
  b(s) = integral K(s, t) x(t) dif t + e(s),
$
which says that the brightness at the point $s$ of the blurred image is a weighted integral of the sharp image. The weight $K(s, t)$, called the _point spread function_ (PSF), tells how much light from the point $t$ ends up at the point $s$. We use the Gaussian kernel (1.9) with $B = beta I$:
$
  K(s, t) = c e^(-beta norm(s - t)_2^2).
$
Points close to each other get a large weight and points far apart almost zero. The parameter $beta$ controls the width: the Gaussian has standard deviation $1 \/ sqrt(2 beta) approx 0.022$ for $beta = 1000$, which is about $3.4$ pixels. The constant $c = beta \/ pi$ makes the weights sum to one, so that blurring does not change the overall brightness.

*Discretization.* We place the image on the square $[-1, 1]^2$. Each pixel is a small square with side $h = 2 \/ 300$, and its position $t_j$ is its midpoint. Replacing the integral by a sum over the pixels (each pixel contributes its value times its area $h^2$, the box rule from Task 5) gives
$
  (A x)_i = sum_j h^2 K(s_i, t_j) x_j, quad "so" quad a_(i j) = h^2 c e^(-beta norm(s_i - t_j)_2^2).
$
Pixels outside the image are not part of the sum, which is the same as assuming that everything outside the image is black (zero boundary conditions). This is why the blurred image gets a dark frame (@fig:task8-lambda).

*Kronecker structure.* The Gaussian splits into a vertical and a horizontal part,
$
  e^(-beta norm(s - t)_2^2) = e^(-beta (s_1 - t_1)^2) e^(-beta (s_2 - t_2)^2),
$
so a 2D blur is the same as first blurring every column of the image vertically and then every row horizontally. With the one-dimensional blur matrices (built exactly like $A$ in Task 5)
$
  (A_1)_(i k) = h sqrt(beta / pi) e^(-beta (t_i - t_k)^2) in RR^(n_1 times n_1), quad A_2 in RR^(n_2 times n_2) "built the same way",
$
blurring the columns of $X$ is $A_1 X$ and blurring its rows is $X A_2^T$, so the blurred image is $A_1 X A_2^T$. Using the standard rule $"vec"(A_1 X A_2^T) = (A_2 times.o A_1) "vec"(X)$ we get
$
  A = A_2 times.o A_1.
$
Here $times.o$ is the _Kronecker product_: $A_2 times.o A_1$ is the block matrix whose block $(k, l)$ is $(A_2)_(k l) A_1$. So $A$ is an $n_2 times n_2$ grid of blocks, each of size $n_1 times n_1$.

*Thresholding.* The Gaussian is never exactly zero, so $A$ is a dense $90 thin 000 times 90 thin 000$ matrix with $8.1 dot 10^9$ entries, about $65$ GB. But the weights of pixels far apart are tiny: for $beta = 1000$, pixels more than about $10$ pixels apart get a weight below $10^(-4)$. We therefore set all entries $|a_(i j)| < epsilon$ to zero, so every pixel only depends on its neighbours and $A$ becomes sparse.

We cannot form the dense $A$ and then threshold it, since it does not fit in memory. Instead we use that $a_(i j) = (A_1)_(i_1 k_1) (A_2)_(i_2 k_2)$ is a product of an entry of $A_1$ and an entry of $A_2$. An entry of $A_1$ that is so small that even multiplied by the largest entry of $A_2$ it is below $epsilon$ can never give an entry of $A$ above $epsilon$, so it is removed already in $A_1$ (and the same for $A_2$). We then form the sparse Kronecker product and remove the remaining entries below $epsilon$. The result is the same as thresholding $A$ directly.

*Sparse storage.* We store $A$ in the CSR (compressed sparse row) format, which only stores the nonzero entries: for each nonzero its value ($8$ bytes) and its column index ($4$ bytes), plus one pointer per row ($4$ bytes) telling where the row starts.

== The data $b$
The data is generated by !!!!eq:model, $b = A x_"true" + e$, where every entry of $e$ is an independent standard normal number times $sigma = 0.01$. Since the pixel values lie in $[0, 1]$, this is noise of about $1 %$ of the brightness range. The true image and $b$ are shown in @fig:task8-lambda.

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
@fig:task6-spy-A shows the sparsity pattern of $A$, that is, where its nonzero entries are. $A$ has $3.05 dot 10^7$ nonzeros, which is $0.376 %$ of all entries. On average each row has $339$ nonzeros out of $90 thin 000$: every pixel of the blurred image depends on about $339$ pixels of the sharp image, roughly a disk with a radius of $10$ pixels. The pattern is a band of blocks: each block is a band (from $A_1$), and the blocks themselves form a band (from $A_2$).

#figure(
  image("assets/task6_1_spy.svg", width: 85%),
  caption: [Sparsity pattern of $A$ ($beta = 1000$, $epsilon = 10^(-4)$). Left: all of $A$. Right: the first $3 n_1$ rows and columns, showing the banded blocks of $A_2 times.o A_1$.],
) <fig:task6-spy-A>

== Sparsity pattern of $L$
$L$ has $358 thin 800$ nonzeros ($2$ per row), which is only $0.0022 %$ of all entries (@fig:task6-spy-L). In the top half ($I_(n_2) times.o L_1$, vertical neighbours) the two nonzeros of a row are in neighbouring columns $i$ and $i + 1$. In the bottom half ($L_2 times.o I_(n_1)$, horizontal neighbours) they are $n_1 = 300$ columns apart, since horizontal neighbours are $n_1$ positions apart in $x$.

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
    [$A$], [$90 thin 000 times 90 thin 000$], [$64.8$ GB], [$0.37$ GB],
    [$L$], [$179 thin 400 times 90 thin 000$], [$129$ GB], [$5.0$ MB],
    table.hline(),
  ),
  caption: [Memory needed for dense and sparse storage of $A$ and $L$ ($beta = 1000$, $epsilon = 10^(-4)$).],
) <tab:memory>

== Influence of $beta$ and $epsilon$
@tab:sparsity and @fig:task6-params show how the number of nonzeros of $A$ depends on $beta$ and $epsilon$.
- A larger $beta$ gives a narrower Gaussian, so fewer entries are above $epsilon$ and $A$ becomes sparser.
- A larger $epsilon$ removes more entries, so $A$ also becomes sparser. The price is that $A$ is a worse approximation of the real blur, since we cut away more of the tails of the Gaussian.

The memory and the cost of a matrix-vector product ($approx 2 thin "nnz"$ operations) both grow linearly with nnz. For our parameters, one product with $A$ costs about $6 dot 10^7$ operations instead of $1.6 dot 10^10$ for a dense $A$, about $265$ times fewer. Since every iteration in Tasks 8–9 uses four sparse products, this is what makes the iterative methods fast. $L$ does not depend on $beta$ or $epsilon$, and it is so sparse that its cost is negligible compared to $A$.

#figure(
  table(
    columns: 4,
    align: (center, center, right, right),
    stroke: none,
    table.hline(),
    [$beta$], [$epsilon$], [nnz$(A)$], [nonzero %],
    table.hline(stroke: 0.5pt),
    [$500$], [$10^(-3)$], [$2.43 dot 10^7$], [$0.300$],
    [$500$], [$10^(-4)$], [$5.13 dot 10^7$], [$0.633$],
    [$500$], [$10^(-5)$], [$7.99 dot 10^7$], [$0.987$],
    [$1000$], [$10^(-3)$], [$1.63 dot 10^7$], [$0.201$],
    [$1000$], [$10^(-4)$], [$3.05 dot 10^7$], [$0.376$],
    [$1000$], [$10^(-5)$], [$4.42 dot 10^7$], [$0.545$],
    [$2000$], [$10^(-3)$], [$1.07 dot 10^7$], [$0.132$],
    [$2000$], [$10^(-4)$], [$1.70 dot 10^7$], [$0.210$],
    [$2000$], [$10^(-5)$], [$2.50 dot 10^7$], [$0.308$],
    table.hline(),
  ),
  caption: [Number of nonzeros of $A$ for different $beta$ and $epsilon$.],
) <tab:sparsity>

#figure(
  image("assets/task6_1_params.svg", width: 80%),
  caption: [Sparsity pattern of $A$ for different $beta$ (rows) and $epsilon$ (columns).],
) <fig:task6-params>

= Task 7: Prove that the matrix is SPD
We now know what we have to solve: the normal equations @eq:normal, $M x = c$ with $M = A^T A + lambda^2 L^T L$. Before choosing a method, we need to know what kind of matrix $M$ is, because the CG method only works for _symmetric positive definite_ (SPD) matrices. A matrix is symmetric if $M^T = M$, and positive definite if $z^T M z > 0$ for every vector $z != 0$. Geometrically, SPD means that the quadratic function $phi(x) = 1/2 x^T M x - c^T x$ is a "bowl" with exactly one lowest point, and that point is the solution of $M x = c$ (Task 8 explains why this matters for CG).

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

*Why this makes CG applicable.* CG finds the solution of $M x = c$ by walking downhill on the bowl $phi(x) = 1/2 x^T M x - c^T x$ (Task 8). This only works because $M$ is SPD:
- the bowl has exactly one lowest point, and it is the solution $x_*$ (since $nabla phi(x) = M x - c$),
- the step length in CG divides by $p^T M p$, which is positive for every direction $p != 0$,
- $norm(z)_M = sqrt(z^T M z)$ is a norm, so "the smallest error in the $M$-norm", which CG minimizes, makes sense.

For our image problem the assumption $"rank" mat(A; L) = n$ holds (see the end of the description of $L$ in Task 6), so $M$ is SPD for every $lambda != 0$ and CG can be used. For $lambda = 0$ we get $M = A^T A$, which in exact arithmetic is SPD when $A$ is invertible, but numerically it is singular because of the tiny singular values (introduction). We will see the effect of this in Task 8.1.

= Task 8: Solving the normal equations with CG
== The idea of CG

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

= Task 9: Solving the normal equations with MINRES
== The idea of MINRES
MINRES searches in the same Krylov space $x_0 + cal(K)_k (M, r_0)$ as CG, but chooses the $x_k$ with the smallest _residual_,
$
  x_k = arg min_(x in x_0 + cal(K)_k (M, r_0)) norm(c - M x)_2.
$
Since $r = -1/2 nabla f$, MINRES makes the gradient of the objective as small as possible, while CG makes the objective itself as small as possible. MINRES only needs $M$ to be symmetric, so it also works for symmetric _indefinite_ matrices, where CG can break down. Our $M$ is SPD, so both can be used.

== The algorithm
*Lanczos.* MINRES first builds an orthonormal basis $v_1, dots, v_k$ of the Krylov space. Start with $v_1 = r_0 \/ rho$, where $rho = norm(r_0)_2$. To get the next basis vector we compute $M v_k$ and remove its components along the earlier basis vectors. Because $M$ is symmetric, $M v_k$ is automatically orthogonal to all but the last two basis vectors, so only a three-term recurrence is needed:
$
  beta_k v_(k+1) = M v_k - alpha_k v_k - beta_(k-1) v_(k-1), quad alpha_k = v_k^T M v_k, quad beta_k = norm(M v_k - alpha_k v_k - beta_(k-1) v_(k-1))_2.
$
(Here $alpha_k, beta_k$ are the Lanczos coefficients, not the CG step length or the blur parameter.) In matrix form, with $V_k = [v_1, dots, v_k]$,
$
  M V_k = V_(k+1) underline(T)_k, quad underline(T)_k = mat(alpha_1, beta_1; beta_1, alpha_2, dots.down; , dots.down, dots.down, beta_(k-1); , , beta_(k-1), alpha_k; , , , beta_k) in RR^((k+1) times k).
$

*The small least squares problem.* Every $x$ in the search space can be written $x = x_0 + V_k y$ with $y in RR^k$. Then
$
  c - M x = r_0 - M V_k y = rho v_1 - V_(k+1) underline(T)_k y = V_(k+1) (rho e_1 - underline(T)_k y),
$
and since $V_(k+1)$ has orthonormal columns, $norm(c - M x)_2 = norm(rho e_1 - underline(T)_k y)_2$. The huge minimization over $RR^n$ has become a small $(k+1) times k$ least squares problem.

*Givens rotations.* The small problem is solved with a QR factorization of $underline(T)_k$, built with _Givens rotations_. A Givens rotation $mat(c, s; -s, c)$ with $c^2 + s^2 = 1$ rotates two rows so that one chosen entry becomes zero. Since $underline(T)_k$ is tridiagonal, each new column $(beta_(k-1), alpha_k, beta_k)$ only needs the two previous rotations plus one new rotation that zeroes $beta_k$. Applying the same rotations to the right-hand side $rho e_1$, its last entry $tau$ satisfies $|tau| = norm(r_k)_2$, so the residual norm is known in every iteration without computing $c - M x_k$. Finally, $x_k$ is updated with a multiple of a new direction $p_k$, which is computed from $v_k$ and the two previous directions. Like CG, MINRES therefore only stores a few vectors of length $n$ and uses one product with $M$ per iteration. Our implementation is `minres` in the notebook, and it uses the same `LinearOperator` for $M$ as CG.

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
