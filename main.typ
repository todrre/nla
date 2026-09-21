= Task 1
We want to solve the regularized least-squares problem
$ min_x quad f(x) := ||A x - b||_2^2 + lambda^2 ||x||_2^2. $

Expanding the objective gives
$ f(x) = (A x - b)^T (A x - b) + lambda^2 x^T x
       = x^T A^T A x - 2 x^T A^T b + b^T b + lambda^2 x^T x. $

Since $f$ is quadratic (and convex), $nabla f(x) = 0$ is a necessary and sufficient condition for optimality:
$ nabla f(x) = 2 A^T A x - 2 A^T b + 2 lambda^2 x = 0. $

This simplifies to the normal equations
$ (A^T A + lambda^2 I) x = A^T b. $



Let $A = U E V^T$ be the SVD of $A$. Substituting into the normal equations and using $U^T U = I$:
$ (A^T A + lambda^2 I) x_"RLS" &= A^T b \
  (V E^T U^T U E V^T + lambda^2 I) x_"RLS" &= V E^T U^T b \
  V (E^2 + lambda^2 I) V^T x_"RLS" &= V E^T U^T b \
  (E^2 + lambda^2 I) x_"RLS" &= E^T U^T b $

Let $B := E^2 + lambda^2 I$, which is diagonal. Then
$ x_"RLS" = B^(-1) E^T U^T b. $

Since $B$ is diagonal, this holds component-wise: writing $y_i := v_i^T x_"RLS"$,
$ y_i = sigma_i / (sigma_i^2 + lambda^2) u_i^T b. $

This gives us the final relation
$ x_"RLS" = sum_(i=1)^n v_i sigma_i / (sigma_i^2 + lambda^2) u_i^T b = V y. $

The filter factor $sigma_i / (sigma_i^2 + lambda^2)$ can be interpreted as a low-pass filter: as $lambda << sigma_i$ it approaches $1$, at $lambda = sigma_i$ it equals $1/2$, and as $lambda >> sigma_i$ it approaches $0$.
= Task 2

Starting from the (unregularized) least-squares solution expressed via the SVD,
$ x_"LS" = sum_(i=1)^n (u_i^T b) / sigma_i v_i, $
small singular values $sigma_i$ can amplify noise in $b$. As in Task 1, we introduce a filter factor $phi_i (lambda)$ that damps the contribution of each term:
$ y_i = phi_i (lambda) / sigma_i u_i^T b, quad
  y = sum_(i=1)^n y_i v_i = sum_(i=1)^n phi_i (lambda) / sigma_i (u_i^T b) v_i. $

For the truncated SVD, the filter factor is chosen as
$ phi_i (lambda) = cases(1 & "if" i <= k, 0 & "if" i > k). $

Substituting this in, only the first $k$ terms survive:
$ V y = sum_(i=1)^k (u_i^T b) / sigma_i v_i + sum_(i=k+1)^n 0 dot (u_i^T b) / sigma_i v_i = x_"TLS". $