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

= Task 1
We want to solve the regularized least-squares problem:
$min_x f(x) := ||A x - b||_2^2 + lambda^2 ||x||_2^2.$
Expanding and differentiating it gives
$
  f(x) = (A x - b)^T (A x - b) + lambda^2 x^T x
  = x^T A^T A x - 2 x^T A^T b + b^T b + lambda^2 x^T x.
$

$f$ is strictly convex, so $x$ minimizes $f$ if and only if $nabla f(x) = 2 (A^T A + lambda^2 I) x - 2 A^T b = 0$, that is, if and only if
$ (A^T A + lambda^2 I) x = A^T b, $
which is the normal equation. Since $A^T A + lambda^2 I$ is positive definite, it is invertible. The normal equations therefore have exactly one solution $x_"RLS" = (A^T A + lambda^2 I)^(-1) A^T b$. Which was to be shown.

#line(length: 100%, stroke: (paint: gray, thickness: 1pt, cap: "round"))

#let c(x) = box(
  height: 1.6em,
  baseline: 50% - 0.3em,
  align(horizon, x),
)

Now let $A = U Sigma V^T$ be the SVD of $A$. Substituting into the normal equations, using that $U$ and $V$ are orthogonal $=> U^(-1) = U^T$ and $V^(-1) = V^T$.
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

Since $B^(-1)$ and $Sigma$ are diagonal, the $i$-th component of $y$ is
$ y_i = sigma_i / (sigma_i^2 + lambda^2) u_i^T b, quad i = 1, dots, n. $
Finally, since $y = V^T x_"RLS"$ and $V$ is orthogonal, multiplying by $V$ yields $x_"RLS" = V y$, which was to be shown.

#line(length: 100%, stroke: (paint: gray, thickness: 1pt, cap: "round"))

Since $A$ has full rank, all $sigma_i > 0$, so we may write $sigma_i / (sigma_i^2 + lambda^2) = 1 / sigma_i dot sigma_i^2 / (sigma_i^2 + lambda^2) = phi_i (lambda) / sigma_i$, i.e. $y_i = phi_i (lambda) / sigma_i u_i^T b$. Expanding $x_"RLS" = V y = sum_i y_i v_i$ then gives
$
  x_"RLS" = V y = sum_(i=1)^n sigma_i / (sigma_i^2 + lambda^2) (u_i^T b) v_i
  = sum_(i=1)^n hl(phi_i (lambda), #c-filt, #<phi>) hl((u_i^T b) / sigma_i, #c-ls, #<ls>) v_i,
  quad phi_i (lambda) := sigma_i^2 / (sigma_i^2 + lambda^2)
  #note(<phi>, pos: (bottom + right, top + right), dy: 1.3em, dx: -0.6em)[filter factor]
  #note(<ls>, pos: (bottom + left, top + left), dy: 0.8em, dx: 0.6em)[coefficient of $x_"LS"$]
$
#v(1.6em)

The filter factor $phi_i (lambda)$ acts as a low-pass filter on the terms of the least-squares solution: for $lambda << sigma_i$ it approaches $1$ (term kept), at $lambda = sigma_i$ it equals $1/2$, and for $lambda >> sigma_i$ it approaches $0$ (term damped).



#pagebreak()
= Task 2

Starting from the (unregularized) least-squares solution expressed via the SVD,
$ x_"LS" = sum_(i=1)^r (u_i^T b) / sigma_i v_i, $
small singular values $sigma_i$ can amplify noise in $b$. As in Task 1, we introduce a filter factor $phi.alt_i (lambda)$ that damps the contribution of each term:
$
  y_i = (phi_i (lambda)) / sigma_i u_i^T b, quad
  x = V y = sum_(i=1)^r phi_i (lambda) (u_i^T b) / sigma_i v_i.
$

For the truncated SVD, the filter factor is chosen as
$ phi_i = cases(1 & "if" i <= k quad (sigma_i > delta), 0 & "if" i > k quad (sigma_i <= delta)). $

Substituting this in, only the first $k$ terms survive:
$
  x = sum_(i=1)^k (u_i^T b) / sigma_i v_i
  + sum_(i=k+1)^r 0 dot (u_i^T b) / sigma_i v_i
  = x_"TLS"
$

#pagebreak()
= Task 3
We consider the opimization problem
$
  min_x f(x):= ||A x - b||_2^2 + lambda^2 ||L x||_2^2
$
By expanding the terms, we get
$
  f(x) = (A x - b)^T (A x - b) + lambda^2 (L x)^T (L x)
  = x^T A^T A x - 2 x^T A^T b + b^T b + lambda^2 x^T L^T L x.
$

$x$ is minimized when $nabla f(x) = 2 (A^T A + lambda^2 L^T L) x - 2 A^T b = 0$, that is, when
$
  (A^T A + lambda^2 L^T L) x = A^T b
$

This solution is unique when $A^T A + lambda^2 L^T L$ is invertible.


#pagebreak()
= Task 4
We consider the opimization problem
$ min_x f(x):= ||A x - b||_2^2 + lambda^2 ||L x||_2^2 $
Using GSVD we devide the matrices $A$ and $L$ into the form
$
  A = U Sigma_A W^(-1), quad quad
  L = V Sigma_L W^(-1)
$
where
$
  Sigma_A = mat(
    alpha_1, , , , 0;
    , alpha_2, , , ;
    , , dots, , dots.v;
    , , , alpha_r, ;
    0, , dots, , 0
  ), quad
  Sigma_L = mat(
    beta_1, , , , 0;
    , beta_2, , , ;
    , , dots, , dots.v;
    , , , beta_r, ;
    0, , dots, , 0
  )
$

We then expand $f(x)$ and substitute the GSVD of $A$ and $L$ into the normal equations:
$
  f(x) &= (A x - b)^T (A x - b) + lambda^2 (L x)^T (L x)
  = x^T A^T A x - 2 x^T A^T b + b^T b + lambda^2 x^T L^T L x
  \ &=
  x^T (U Sigma_A W^(-1))^T U Sigma_A W^(-1) x - 2x^T (U Sigma_A W^(-1))^T b + b^T b + lambda^2 x^T (V Sigma_L W^(-1))^T V Sigma_L W^(-1) x
  \ &=
  x^T W^(-1^(T)) Sigma_A^T cancel(U^T U) Sigma_A W^(-1) x -2x^T W^(-1^(T)) Sigma_A^T U^T b + b^T b + lambda^2 x^T W^(-1^(T)) Sigma_L^T cancel(V^T V) Sigma_L W^(-1) x
  \ {y = W^(-1) x}
  \ &=
  y^T Sigma_A^T Sigma_A y - 2 y^T Sigma_A^T U^T b + b^T b + lambda^2 y^T Sigma_L^T Sigma_L y
  \ &=
  y^T (Sigma_A^T Sigma_A + lambda^2 Sigma_L^T Sigma_L) y - 2 y^T Sigma_A^T U^T b + b^T b
$

differentiating $f(y)$ gives
$
  nabla f(y) = 2 (Sigma_A^T Sigma_A + lambda^2 Sigma_L^T Sigma_L) y - 2 Sigma_A^T U^T b = 0
$

$
  (Sigma_A^T Sigma_A + lambda^2 Sigma_L^T Sigma_L) y = Sigma_A^T U^T b
$

$
  B := Sigma_A^T Sigma_A + lambda^2 Sigma_L^T Sigma_L =
  mat(
    alpha_1^2, , , , 0;
    , alpha_2^2, , , ;
    , , dots, , dots.v;
    , , , alpha_r^2, ;
    0, , dots, , 0
  ) + lambda^2 mat(
    beta_1^2, , , , 0;
    , beta_2^2, , , ;
    , , dots, , dots.v;
    , , , beta_r^2, ;
    0, , dots, , 0
  ) =
  mat(
    alpha_1^2 + lambda^2 beta_1^2, , , , 0;
    , alpha_2^2 + lambda^2 beta_2^2, , , ;
    , , dots, , dots.v;
    , , , alpha_r^2 + lambda^2 beta_r^2, ;
    0, , dots, , 0
  )
$

this gives us the final relation
$
  y = B^(-1) Sigma_A^T U^T b
  = mat(
    1/(alpha_1^2 + lambda^2 beta_1^2), , , , 0;
    , 1/(alpha_2^2 + lambda^2 beta_2^2), , , ;
    , , dots, , dots.v;
    , , , 1/(alpha_r^2 + lambda^2 beta_r^2), ;
    0, , dots, , 0
  ) mat(
    alpha_1, , , , 0;
    , alpha_2, , , ;
    , , dots, , dots.v;
    , , , alpha_r, ;
    0, , dots, , 0
  ) vec(u_1, dots.v, u_r)^T b
$

component wise, this gives
$
  y_i = alpha_i / (alpha_i^2 + lambda^2 beta_i^2) u_i^T b
  = hl(phi_i (lambda), #c-filt, #<phi>) hl((u_i^T b) / alpha_i, #c-ls, #<ls>) v_i
  quad phi_i (lambda) := alpha_i^2 / (alpha_i^2 + lambda^2 beta_i^2)
  #note(<phi>, pos: (bottom + right, top + right), dy: 1.3em, dx: -0.6em)[filter factor]
  #note(<ls>, pos: (bottom + left, top + left), dy: 0.8em, dx: 0.6em)[coefficient of $x_"LS"$]
$


#pagebreak()
= Task 5
A samlar in data
L vad vi redan vet eller antar om bilden/signalen.

Debluring = göra en suddig bild skarp igen.

$
  integral_(RR^d) K(s,t) x(t) d t + e(s) = b(s)
$

t en punkt i skarpa bilden. s en punkt i suddiga bilden.

K ( s,t), point spread function (PSF), det som gör bilden suddig.
b(s) den suddiga bilden.
e(s) brus
x(t) den skarpa bilden.

så formeln beskriver från skarp bild till suddig bild. Vi vill göra tvärtom.


A matrisen ärbluring operatorn

=== val av regulariserings matris $L$
Valet beror på specifika problemet.
enklaste valet är $L = I$ (identity matrix) ....
elle

$
  L = mat(
    1, -1, , , ;
    , 1, -1, , ;
    , , dots.down, dots.down, ;
    , , , 1, -1;
  )
$

finite discretization of the first derivative. T


=== Sammanfattning: från suddig till skarp signal
Vi har en sann signal $x(t)$ på $[a, b]$ och samplar den i $n$ punkter $t_1, dots, t_n$. Mätningen suddas ut av Gausskärnan $K(s,t) = c exp(-beta (s - t)^2)$, som ersätter varje punkt med ett viktat medelvärde av sina grannar. Integralen approximeras med boxregeln ($h = (b - a)\/n$):
$
  b_i approx sum_(j=1)^n h K(t_i, t_j) x_j quad => quad A_(i j) = h c exp(-beta (t_i - t_j)^2), quad b = A x + e.
$

Suddningen tar bort detaljer, så singulärvärdena för $A$ avtar snabbt mot noll. Direkt invertering ger $A^(-1) b = x + A^(-1) e$, där bruset förstärks kraftigt av de små singulärvärdena. Vi löser därför det regulariserade problemet
$
  min_x ||A x - b||_2^2 + lambda^2 ||L x||_2^2,
$
där $L$ är första differensen. $||L x||$ straffar stora hopp mellan grannvärden, så lösningen blir jämn, och $lambda$ balanserar datapassning mot jämnhet.

GSVD:n diagonaliserar $A$ och $L$ samtidigt, och enligt Task 4 blir lösningen komponentvis
$
  x_lambda = W y, quad y_i = phi_i (lambda) (u_i^T b) / alpha_i, quad phi_i (lambda) = alpha_i^2 / (alpha_i^2 + lambda^2 beta_i^2).
$
Filterfaktorn behåller komponenter där $alpha_i >> lambda beta_i$ och dämpar de brusdominerade där $alpha_i << lambda beta_i$. GSVD:n räknas bara en gång, så det blir billigt att testa många $lambda$ och välja det som ger bäst rekonstruktion.
