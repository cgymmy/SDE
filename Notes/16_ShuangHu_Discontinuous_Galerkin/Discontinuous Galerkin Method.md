> [!remark] Remark. 
> This note concentrates on the discontinuous Galerkin method for hyperbolic equation or conservation system.

## Brief Introduction

Discontinuous Galerkin (DG) methods are a class of finite element methods that employ **completely discontinuous** basis functions. In essence, these methods seek an approximate **weak solution** to a partial differential equation (PDE) of the form $\mathcal{L}u = f$ in a **broken Sobolev space**, which allows for discontinuities at the boundaries of mesh elements. Compared to classical finite element methods, DG methods exhibit the following characteristics:
- **High parallel efficiency**: The discontinuous nature of the basis functions enables efficient parallel computation.
- **High flexibility**: DG methods can handle complex geometries and adapt to varying solution behaviors.
- **Complete freedom in polynomial degree adaptation (p-adaptivity)**: The polynomial degree in each element can be adjusted independently of neighboring elements.
- **Extremely local data structure**: The methods rely on localized data, reducing communication overhead in parallel implementations.
However, DG methods also have some drawbacks:
- **Increased degrees of freedom**: The discontinuous nature leads to a larger number of unknowns compared to continuous methods.
- **Condition number of the stiffness matrix**: The stiffness matrix in DG methods often has a higher condition number, which can affect solver performance.
- **Complex error analysis**: The analysis of errors in DG methods is more intricate due to the discontinuities and the need for specialized techniques.
## DG method for conservation laws
First, we discuss general DG methods for conservation laws, then derive its error and stability analysis.
### One-Dimensional Time-Dependent Conservation Laws
#### Spatial Discretization

We begin our discussion with the following system:
$$
\left\{
\begin{aligned}
u_{t} + (f(u))_{x} &= 0, \\
u(x, 0) &= u_{0}(x),
\end{aligned}
\right.
$$
where $x \in \Omega := [0, 1]$, and $u$ is equipped with **periodic boundary conditions**. First, we define a mesh on $\Omega$: let $h := \frac{1}{N}$ , $x_{j+\frac{1}{2}} := jh$, and $I_{j} := \left[ x_{j-\frac{1}{2}}, x_{j+\frac{1}{2}} \right]$. We then define the finite element space:
$$
V_{h}^{k} := \{ u \in L^{2}(\Omega) : u \in \mathcal{P}^{k}(I_{j}), \, 1 \leq j \leq N \}.
$$
Next, we derive the weak form on $I_{j}$: for $v \in V_{h}^{k}$,
$$
\int_{I_{j}} u_{t} v \, \mathrm{d}x + \int_{I_{j}} (f(u))_{x} v \, \mathrm{d}x = 0.
$$
Applying integration by parts, we obtain:
$$
\int_{I_{j}} u_{t} v \, \mathrm{d}x - \int_{I_{j}} f(u) v_{x} \, \mathrm{d}x + \left. (f(u) v) \right|_{x_{j+\frac{1}{2}}^{-}} - \left. (f(u) v) \right|_{x_{j-\frac{1}{2}}^{+}} = 0.
$$
We seek $u \in L^{2}(V_{h}^{k}; [0, T])$ that satisfies the above weak form. Since $u$ and $v$ are both discontinuous, we must determine the values of $f(u)$ and $v$ at the mesh boundaries. 

Determining $v\left( x_{j+\frac{1}{2}}^{-} \right)$ and $v\left( x_{j-\frac{1}{2}}^{+} \right)$ is straightforward because $v$ is a test function; we simply choose its left or right limit. However, the situation is different for $\left. f(u) \right|_{x_{j+\frac{1}{2}}^{-}}$. 

To solve the weak form, **boundary conditions** for $u|_{I_{j}}$ are required. Therefore, it is inappropriate to simply set $\left. f(u) \right|_{x_{j+\frac{1}{2}}^{-}} = f\left( u_{j+\frac{1}{2}}^{-} \right)$. For example, if we consider $V_{h}^{0}$, the equation reduces to $\int_{I_{j}} u_{t} v \, \mathrm{d}x = 0$, which does not allow us to determine the value of $u$. Thus, we must derive $\left. f(u) \right|_{x_{j+\frac{1}{2}}^{-}}$ and $\left. f(u) \right|_{x_{j-\frac{1}{2}}^{+}}$ using appropriate approximations. This leads to the concept of the **numerical flux**.

Typically, the numerical flux for a function $f$ is a function of the left and right limits of $u$ denoted as $\hat{f}(u^{-}, u^{+})$. In general, $\hat{f}$ should satisfy the following properties:
- **Consistency**: $\hat{f}(u, u) = f(u)$.
- **Continuity**: $\hat{f}(u^{-}, u^{+})$ is at least **Lipschitz continuous** with respect to both arguments $u^{-}$ and $u^{+}$.
- **Monotonicity**: $\hat{f}(u^{-}, u^{+})$ is a non-decreasing function of $u^{-}$ and a non-increasing function of $u^{+}$.

Some well-known monotone numerical fluxes include:
- **The Lax-Friedrichs flux**:
  $$
  \hat{f}^{LF}(u^{-}, u^{+}) = \frac{1}{2} \left( f(u^{-}) + f(u^{+}) - \alpha (u^{+} - u^{-}) \right), \quad \alpha = \max_{u} |f'(u)|;
  $$
- **The Godunov flux**:
  $$
  \hat{f}^{God}(u^{-}, u^{+}) = 
  \left\{
  \begin{aligned}
  \min_{u^{-} \leq u \leq u^{+}} f(u), & \quad \text{if } u^{-} < u^{+}, \\
  \max_{u^{+} \leq u \leq u^{-}} f(u), & \quad \text{if } u^{-} \geq u^{+};
  \end{aligned}
  \right.
  $$
- **The Engquist-Osher flux**:
  $$
  \hat{f}^{EO} = \int_{0}^{u^{-}} \max(f'(u), 0) \, \mathrm{d}u + \int_{0}^{u^{+}} \min(f'(u), 0) \, \mathrm{d}u + f(0).
  $$

Finally, the spatial discretization scheme for the one-dimensional conservation law is given by:
$$
\int_{I_{j}} u_{t} v \, \mathrm{d}x - \int_{I_{j}} f(u) v_{x} \, \mathrm{d}x + \hat{f}_{j+\frac{1}{2}} \left. v \right|_{x_{j+\frac{1}{2}}^{-}} - \hat{f}_{j-\frac{1}{2}} \left. v \right|_{x_{j-\frac{1}{2}}^{+}} = 0.
$$
#### Time discretization
For hyperbolic problems, we often use a class of **high order nonlinearly stable Runge-Kutta** time discretizations. The most popular scheme in this class is the following third order RK method for solving 
$$
u_{t}=\mathcal{L}(u,t):
$$
$$
\begin{align}
u^{(1)}&= u^{n}+\Delta t \mathcal{L}(u^{n},t^{n}), \\
u^{(2)}&= \frac{3}{4}u^{n}+\frac{1}{4}u^{(1)}+\frac{1}{4}\Delta t \mathcal{L}(u^{(1)},t^{n}+\Delta t), \\
u^{n+1}&=\frac{1}{3}u^{n}+\frac{2}{3}u^{(2)}+\frac{2}{3}\Delta t \mathcal{L}\left( u^{(2)},t^{n}+\frac{1}{2}\Delta t \right).
\end{align}
$$
By combining the expressions for time discretization and spatial discretization, we derive the DG scheme for the one-dimensional conservation law.
#### $L^{2}$ stability
Since the weak solution of conservation law may not be unique, we seek the **entropy weak solution**. The entropy weak solution $u$ satisfies the following **entropy inequality**
$$
U(u)_{t}+F(u)_{x}\le 0
$$
in distribution sense, for any *convex* function $U(u)$ and the corresponding entropy flux $F(u)=\int U'(u)f'(u)\mathrm{d}u$. 

The solution $u_{h}$ to the semi-discretize DG scheme satisfies the following **cell entropy inequality**:
$$
\frac{\mathrm{d}}{\mathrm{d}t}\int_{I_{i}}U(u_{h})\mathrm{d}x+\hat{F}_{i+\frac{1}{2}}-\hat{F}_{i-\frac{1}{2}}\le 0
$$
for the square entropy $U(u):=\frac{1}{2}u^{2}$ and some consistent entropy flux $\hat{F}$.
**Proof:** Recall the semi-discrete scheme:
$$
\int_{I_{j}} u_{t} v \, \mathrm{d}x - \int_{I_{j}} f(u) v_{x} \, \mathrm{d}x + \hat{f}_{j+\frac{1}{2}} \left. v \right|_{x_{j+\frac{1}{2}}^{-}} - \hat{f}_{j-\frac{1}{2}} \left. v \right|_{x_{j-\frac{1}{2}}^{+}} = 0.
$$
Mark the notation 
$$
B_{j}(u_{h};v_{h}):=\int_{I_{j}} (u_{h})_{t} v_{h} \, \mathrm{d}x - \int_{I_{j}} f(u_{h}) (v_{h})_{x} \, \mathrm{d}x + \hat{f}_{j+\frac{1}{2}} \left. v_{h} \right|_{x_{j+\frac{1}{2}}^{-}} - \hat{f}_{j-\frac{1}{2}} \left. v_{h} \right|_{x_{j-\frac{1}{2}}^{+}} = 0,
$$
since $B_{j}(u_{h};u_{h})=0$, we have:
$$
\frac{\mathrm{d}}{\mathrm{d}t}\int_{I_{j}}U(u_{h})\mathrm{d}x-\int_{I_{j}}f(u_{h})\mathrm{d}u_{h}+\hat{f}_{j+\frac{1}{2}}u_{h}\left( x_{j+\frac{1}{2}}^{-} \right)-\hat{f}_{j-\frac{1}{2}}u_{h}\left( x_{j-\frac{1}{2}}^{+} \right),
$$
then mark $\tilde{F}(u):=\int_{}^{u}f(u)\mathrm{d}u$ and $\hat{F}_{j+\frac{1}{2}}:=-\tilde{F}\left( u_{h}\left( x_{j+\frac{1}{2}}^{-} \right) \right)+\hat{f}_{j+\frac{1}{2}}u_{h}\left( x_{j+\frac{1}{2}}^{-} \right)$, we have:
$$
\frac{\mathrm{d}}{\mathrm{d}t}\int_{I_{j}}U(u_{h})\mathrm{d}x+\hat{F}_{j+\frac{1}{2}}-\hat{F}_{j-\frac{1}{2}}+\Theta_{j-\frac{1}{2}}=0,
$$
where 
$$
\Theta_{j-\frac{1}{2}}=-\tilde{F}\left( u_{h}\left( x_{j-\frac{1}{2}}^{-} \right) \right)+\hat{f}_{j-\frac{1}{2}}u_{h}\left( x_{j-\frac{1}{2}}^{-} \right)+\tilde{F}\left( u_{h}\left( x_{j-\frac{1}{2}}^{+} \right) \right)-\hat{f}_{j-\frac{1}{2}}u_{h}\left( x_{j-\frac{1}{2}}^{+} \right).
$$
By Lagrangian M.V.T:
$$
\begin{align}
\Theta_{j-\frac{1}{2}}&=\hat{f}_{j-\frac{1}{2}}\left( u_{h}\left( x_{j-\frac{1}{2}}^{-} \right)-u_{h}\left( x_{j-\frac{1}{2}}^{+} \right) \right)+f(\xi )\left( u_{h}\left( x_{j-\frac{1}{2}}^{+} \right)-u_{h}\left( x_{j-\frac{1}{2}}^{-} \right) \right) \\
&=(\hat{f}(\xi ,\xi )-\hat{f}(u^{-},u^{+}))\left( u_{h}\left( x_{j-\frac{1}{2}}^{+} \right)-u_{h}\left( x_{j-\frac{1}{2}}^{-} \right) \right)\ge0.
\end{align}
$$
By monotonic property, the final inequality holds. It means the cell entropy inequality holds.

Sum up the cell entropy inequality w.r.t $j$, we have the $L^{2}$ stability:
$$
\frac{\mathrm{d}}{\mathrm{d}t} \int_{0}^{1}(u_{h})^{2}\mathrm{d}x\le0.
$$
#### Limiters and total variation
When we apply high order DG scheme, the numerical solution may appear high oscillation, specially on the points near discontinuities. To handle this issue, we often need to apply **nonlinear limiters** to control these oscillations.

We first consider the forward Euler (FE) time discretization for the DG scheme. By the discription of FE scheme, the "pre-process" of the time step $n+1$ shows:
$$
\int_{I_{j}}\frac{u_{h}^{n+1,pre}-u_{h}^{n}}{\Delta t}v_{h}\mathrm{d}x-\int_{I_{j}}f(u_{h}^{n})(v_{h})_{x}\mathrm{d}x+\hat{f}_{j+\frac{1}{2}}^{n}v_{h}\left( x_{i+\frac{1}{2}}^{-} \right)-\hat{f}_{j-\frac{1}{2}}^{n}v_{h}\left( x_{i-\frac{1}{2}}^{+} \right)=0,
$$
after this procedure, we should act a **limiting procedure** from $u_{h}^{n,pre}$ to $u_{h}^{n}$. This procedure should satisfy the following two conditions:
* It should not change the **cell averaged** of $u_{h}^{n,pre}$, this is for the **conservation property** of the DG method.
* It should not affect the **accuracy** of the scheme in **smooth regions**, i.e. in the smooth region, $u_{h}^{n}=u_{h}^{n,pre}$.

In this note, we present one of limiters as an example.

We denote 
$$
\bar{u}_{i}:=\frac{1}{\Delta x_{i}}\int_{I_{i}}u_{h}\mathrm{d}x
$$
and 
$$
\tilde{u}_{i}:=u_{h}\left( x_{i+\frac{1}{2}}^{-} \right)-\bar{u}_{i},\quad \tilde{\tilde{u}}_{i}:=\bar{u}_{i}-u_{h}\left( x_{i-\frac{1}{2}}^{+} \right).
$$
The limiter should not change $\bar{u}_{i}$, but it can change the range $\tilde{u}_{i}$ and $\tilde{\tilde{u}}_{i}$. The **minmod limiter** changes $\tilde{u}_{i}$ and $\tilde{\tilde{u}}_{i}$ into:
$$
\tilde{u}_{i}^{(\text{mod})}=m(\tilde{u}_{i},\Delta_{+}\bar{u}_{i},\Delta_{-}\bar{u}_{i}),\quad \tilde{\tilde{u}}_{i}^{\text{mod}}=m(\tilde{\tilde{u}}_{i},\Delta_{+}\bar{u}_{i},\Delta_{-}\bar{u}_{i}),
$$
where 
$$
\Delta_{+} \bar{u}_{i}=\bar{u}_{i+1}-\bar{u}_{i},\quad \Delta_{-}\bar{u}_{i}=\bar{u}_{i}-\bar{u}_{i-1},
$$
and the minmod function $m$ is:
$$
m(a_{1},\ldots,a_{l})=
\left\{
\begin{aligned}
&s\min(|a_{1}|,\ldots,|a_{l}|),&\quad \text{if}\;s=\text{sign}(a_{1})=\ldots=\text{sign}(a_{l});\\
&0,& \text{otherwise}.
\end{aligned}
\right.
$$
Then we can modify $u_{h}^{n,pre}$. The new point values are given by:
$$
u_{h}^{(\text{mod})}\left( x_{i+\frac{1}{2}}^{-} \right)=\bar{u}_{i}+\tilde{u}_{i}^{(\text{mod})},\quad u_{h}^{(\text{mod})}\left( x_{i-\frac{1}{2}}^{+} \right)=\bar{u}_{i}-\tilde{\tilde{u}}_{i}^{\text{mod}}.
$$
This recovery is unique for $P^{k}$ polynomials with $k\le 2$. By these modified point value, we give the modified numerical solutions.

The DG scheme with modification has total variation stability. First, we present *Harten Lemma* .

> [!lemma] Lemma.[Harten]
If a scheme can be written in the form 
$$
u_{i}^{n+1}=u_{i}^{n}+C_{i+\frac{1}{2}}\Delta_{+}u_{i}^{n}-D_{i-\frac{1}{2}}\Delta_{-}u_{i}^{n}
$$
with periodic or compactly supported boundary conditions, where $C_{i+\frac{1}{2}}$ and $D_{i-\frac{1}{2}}$ may be nonlinear functions of the grid values $u_{j}^{n}$ for $j=i-p,\ldots,i+q$ with some $p,q\ge0$, satisfying:
$$
C_{i+\frac{1}{2}}\ge0,\quad D_{i+\frac{1}{2}}\ge 0,\quad C_{i+\frac{1}{2}}+D_{i+\frac{1}{2}}\le 1,\quad \forall i
$$
then the scheme is TVD:
$$
TV(u^{n+1})\le TV(u^{n}),
$$
where the total variation seminorm is defined by 
$$
TV(u):=\sum\limits_{i}|\Delta_{+}u_{i}|.
$$

**Proof:** Since
$$
u_{i}^{n+1}=u_{i}^{n}+C_{i+\frac{1}{2}}\Delta_{+}u_{i}^{n}-D_{i-\frac{1}{2}}\Delta_{-}u_{i}^{n}
$$
and
$$
u_{i+1}^{n+1}=u_{i+1}^{n}+C_{i+\frac{3}{2}}\Delta_{+}u_{i
+1}^{n}-D_{i+\frac{1}{2}}\Delta_{-}u_{i+1}^{n},
$$
we have:
$$
\Delta_{+} u_{i}^{n+1}=\left( 1-C_{i+\frac{1}{2}}-D_{i+\frac{1}{2}} \right)\Delta_{+}u_{i}^{n}+C_{i+\frac{3}{2}}\Delta_{+}u_{i+1}^{n}+D_{i-\frac{1}{2}}\Delta_{+}u_{i-1}^{n},
$$
then sum up, we have:
$$
\sum\limits_{i}\Delta_{+}u_{i}^{n+1}=\sum\limits_{i}\left(\left( 1-C_{i+\frac{1}{2}}-D_{i+\frac{1}{2}} \right)\Delta_{+}u_{i}^{n}+C_{i+\frac{3}{2}}\Delta_{+}u_{i+1}^{n}+D_{i-\frac{1}{2}}\Delta_{+}u_{i-1}^{n}\right),
$$
it means
$$
\sum\limits_{i}| \Delta _{+} u_{i}^{n+1} |\le \sum\limits_{i}\left( 1-C_{i+\frac{1}{2}}-D_{i+\frac{1}{2}} \right)| \Delta _{+}u_{i}^{n} |+\sum\limits_{i}C_{i+\frac{1}{2}}| \Delta _{+}u_{i}^{n} |+\sum\limits_{i}D_{i+\frac{1}{2}}| \Delta _{+}u_{i}^{n} |=\sum\limits_{i}| \Delta _{+}u_{i}^{n} |.
$$

For DG scheme, we define the *total variation in the means* semi-norm, or **TVM**, as:
$$
TVM(u_{h}):=\sum\limits_{i}| \Delta _{+}\bar{u}_{i} |,
$$
then we have the following proposition:
> [!proposition] Proposition. 
For periodic or compactly supported boundary conditions, the solution $u_{h}^{n}$ of DG scheme with forward Euler (FE) time discretization and the "pre-processing"  by the limiter, is total variation diminishing in the means **(TVDM)**, i.e.
$$
TVM(u_{h}^{n+1})\le TVM(u_{h}^{n}).
$$

**Proof:** Choose $v_{h}=1$ on the scheme, we have:
$$
\int_{I_{i}}\frac{u_{h}^{n+1,pre}-u_{h}^{n}}{\Delta t}\mathrm{d}x+\hat{f}^{n}_{i+\frac{1}{2}}-\hat{f}_{i-\frac{1}{2}}^{n}=0,
$$
mark $\lambda_{i}:=\frac{\Delta t}{\Delta x_{i}}$, it means:
$$
\bar{u}_{i}^{n+1,pre}=\bar{u}_{i}-\lambda_{i}(\hat{f}(\bar{u}_{i}+\tilde{u}_{i},\bar{u}_{i+1}-\tilde{\tilde{u}}_{i+1})-\hat{f}(\bar{u}_{i-1}+\tilde{u}_{i-1},\bar{u}_{i}-\tilde{\tilde{u}}_{i})),
$$
where the right hand omits the time parameter $n$. 
Then, we try to write this scheme as the Harten form. We write:
$$
\begin{align}
C_{i+\frac{1}{2}}&:=-\lambda_{i}\frac{\hat{f}(\bar{u}_{i}+\tilde{u}_{i},\bar{u}_{i+1}-\tilde{\tilde{u}}_{i+1})-\hat{f}(\bar{u}_{i}+\tilde{u}_{i},\bar{u}_{i}-\tilde{\tilde{u}}_{i})}{\Delta_{+}\bar{u}_{i}}, \\
D_{i-\frac{1}{2}}&:=\lambda _{i}\frac{\hat{f}(\bar{u}_{i}+\tilde{u}_{i},\bar{u}_{i}-\tilde{\tilde{u}}_{i})-\hat{f}(\bar{u}_{i-1}+\tilde{u}_{i-1},\bar{u}_{i}-\tilde{\tilde{u}}_{i})}{\Delta _{-}\bar{u}_{i}}.
\end{align}
$$
Now, we should estimate $C_{i+\frac{1}{2}}$ and $D_{i-\frac{1}{2}}$. We denote:
$$
C_{i+\frac{1}{2}}=-\lambda_{i}\frac{\hat{f}(\bar{u}_{i}+\tilde{u}_{i},\bar{u}_{i+1}-\tilde{\tilde{u}}_{i+1})-\hat{f}(\bar{u}_{i}+\tilde{u}_{i},\bar{u}_{i}-\tilde{\tilde{u}}_{i})}{\bar{u}_{i+1}-\tilde{\tilde{u}}_{i+1}-\bar{u}_{i}+\tilde{\tilde{u}}_{i}}\frac{\Delta _{+}\bar{u}_{i}+\tilde{\tilde{u}}_{i}-\tilde{\tilde{u}}_{i+1}}{\Delta _{+}\bar{u}_{i}}:=-\lambda _{i}A.B
$$
Since $\hat{f}$ is monotony, we have $A\le0$ and $A\ge -L_{2}$. Then, from the limiter operator, we have $| \tilde{\tilde{u}}_{i} |\le | \Delta_{+}u_{i} |$ and $|\tilde{\tilde{u}}_{i+1}|\le | \Delta_{-}u_{i+1} |=| \Delta_{+}u_{i} |$, and 
$$
\text{sign}(\tilde{\tilde{u}}_{i})=\text{sign}(\tilde{\tilde{u}}_{i+1})=\text{sign}(\Delta _{+}u_{i}),
$$
then $0\le B\le 2$. Therefore, $0\le C_{i+\frac{1}{2}}\le 2L_{2}\lambda_{i}$. Similarly, $0\le D_{i-\frac{1}{2}}\le 2L_{1}\lambda_{i}$. 
Therefore, if we choose $\Delta t$ such that $\forall i$, $\lambda_{i}\le \frac{1}{2(L_{1}+L_{2})}$, we have $0\le C_{i+\frac{1}{2}}+D_{i+\frac{1}{2}}\le 1$. By Harten's Lemma, the scheme is TVDM.

In the same way, for any *strong stability preserving (SSP)* Runge-Kutta (RK) time discretizations, the DG schemes are all TVDM.

Then, we need to verify that the limiter **doesn't affact accuracy** in smooth regions. By Taylor expansion:
$$
\tilde{u}_{i}=\frac{1}{2}u_{x}(x_{i})\Delta x_{i}+O(h^{2}),\quad \tilde{\tilde{u}}_{i}=\frac{1}{2}u_{x}(x_{i})\Delta x_{i}+O(h^{2}).
$$
> [!remark] Remark. 
$$
\begin{align} 
 u\left( x_{i}+\frac{h}{2} \right)-\bar{u}_{i}&=\frac{1}{h}\int_{-\frac{h}{2}}^{\frac{h}{2}}\left[ u\left( x_{i}+\frac{h}{2} \right)-u(x_{i}+t) \right]\mathrm{d}t \\
&=\frac{1}{h}\int_{-\frac{h}{2}}^{\frac{h}{2}}\left[ u(x_{i})-u(x_{i}+t)+\frac{h}{2}u_{x}(x_{i})+O(h^{2}) \right]\mathrm{d}t \\
&=\frac{h}{2}u_{x}(x_{i})+O(h^{2}).
\end{align}
$$

While:
$$
\Delta _{+}\bar{u}_{i}=\frac{1}{2}u_{x}(x_{i})(\Delta x_i+\Delta x_{i+1})+O(h^{2}),\quad \Delta _{-}\bar{u}_{i}=\frac{1}{2}u_{x}(x_{i})(\Delta x_{i}+\Delta x_{i-1})+O(h^{2}).
$$
When we are in a smooth and monotone region, the first argument in the minmod function is of the same sign as the second and third arguments and is smaller in magnitude when $h$ is small. Therefore, in this case, the modified values will take the unmodified $\tilde{u}_{i}$ and $\tilde{\tilde{u}}_{i}$, which doesn't affact the accuracy.

On the other hand, the TVD limiter does kill accuracy at smooth extrema. They are **at most second order accurate** for smooth but non-monotone solutions. Therefore, in practice we often use a total variation bounded (TVB) corrected limiter:
$$
\tilde{m}(a_{1},\ldots,a_{l})=\left\{\begin{aligned}
&a_{1},\quad & \text{if }| a_{1} |\le Mh^{2};\\
&m(a_{1},\ldots,a_{l}),& \text{otherwise}.
\end{aligned}
\right.
$$
to correct the result near smooth extrema.

#### Error estimates
If we assume the exact solution of conservation law is smooth, we can obtain optimal $L^{2}$ error estimates. For simplicity, in this note, we will give the proof for semi-discrete DG scheme and the linear conservation law:
$$
u_{t}+u_{x}=0
$$
with the monotone flux $\hat{f}(u^{-},u^{+}):=u^{-}.$ 
> [!theorem] Theorem.
The solution $u_{h}$ of the semidiscrete DG scheme for the linear conservation law with a smooth solution $u$ satisfies:
$$
\left\|u-u_{h}\right\|\le Ch^{k+1},
$$
where $C$ depends on $u$ and its derivatives but is independent of $h$.

**Proof:** The DG scheme shows 
$$
B_{i}(u_{h};v_{h})=0
$$
for all $v_{h}\in V_{h}$ and for all $i$. The exact solution of the PDE satisfies:
$$
B_{i}(u;v_{h})=0
$$
for all $i$. Since $B_{i}$ is a bilinear form in this case, 
$$
B_{i}(u-u_{h};v_{h})=0.
$$
Now, we take a special projection $P$ into $V_{h}$. For a given smooth function $w$, we set 
$$
\int_{I_{i}}(Pw(x)-w(x))v_{h}(x)\mathrm{d}x=0\quad \forall v_{h}\in P^{k-1}(I_{i});\quad Pw\left( x_{i+\frac{1}{2}}^{-} \right)=w\left( x_{i+\frac{1}{2}} \right).
$$
Standard approximation theory implies, for a smooth function $w$, 
$$
\left\|Pw(x)-w(x)\right\|\le C(w)h^{k+1}.
$$
Now take 
$$
v_{h}:=Pu-u_{h},
$$
and denote 
$$
e_{h}:=Pu-u_{h},\quad \varepsilon _{h}:=u-Pu.
$$
Since $B_{i}(u-u_{h};v_{h})=0$, we have 
$$
B_{i}(e_{h};e_{h})=-B_{i}(\varepsilon _{h};e_{h}).
$$
First, consider $B_{i}(e_{h};e_{h})$. From the proof for $L^{2}$ stability, we have:
$$
B_{i}(e_{h};e_{h})=
\frac{1}{2}\frac{\mathrm{d}}{\mathrm{d}t}\int_{I_{i}}e_{h}^{2}\mathrm{d}x+\hat{F}_{i+\frac{1}{2}}-\hat{F}_{i-\frac{1}{2}}+\Theta_{i-\frac{1}{2}},
$$
with $\Theta_{i-\frac{1}{2}}\ge 0$. Then, consider $B_{i}(\varepsilon_{h};e_{h})$, we have:
$$
B_{i}(\varepsilon _{h};e_{h})=\int_{I_{i}}(\varepsilon _{h})_{t}e_{h}\mathrm{d}x-\int_{I_{i}}\varepsilon_{h}(e_{h})_{x}\mathrm{d}x+\varepsilon_{h}\left( x_{i+\frac{1}{2}}^{-} \right)e_{h}\left( x_{i+\frac{1}{2}}^{-} \right)-\varepsilon _{h}\left( x_{i-\frac{1}{2}}^{-} \right)e_{h}\left( x_{i-\frac{1}{2}}^{+} \right).
$$
By the definition, $\varepsilon_{h}\left( x_{i+\frac{1}{2}}^{-} \right)=0$, it means 
$$
-B_{i}(\varepsilon _{h};e_{h})=-\int_{I_{i}}(\varepsilon _{h})_{t}e_{h}\mathrm{d}x\le \frac{1}{2}\left( \int_{I_{i}}((\varepsilon _{h})_{t})^{2}\mathrm{d}x +\int_{I_{i}}e_{h}^{2}\mathrm{d}x\right).
$$
Therefore:
$$
\frac{\mathrm{d}}{\mathrm{d}t}\int_{0}^{1}e_{h}^{2}\mathrm{d}x\le \int_{0}^{1}e_{h}^{2}\mathrm{d}x
+Ch^{2k+2}.
$$
By the Gronwall inequality and the initial error, we have $\left\|u-u_{h}\right\|\le Ch^{k+1}$.

### Multi-Dimensional Conservation Laws
In this section, we consider the multi-dimensional cases with arbitrary triangulations. The two dimentional time dependent conservation law:
$$
u_{t}+(f(u))_{x}+(g(u))_{y}=0.
$$
One of the polynomial basis on a triangle cell $\Delta_{i}$ shows 
$\{x_{i}y_{j}\}_{i+j\le k}$, so there are $K=\frac{(k+1)(k+2)}{2}$ degrees of freedom (DoFs) per cell freedom.

Multiply the equation by a test function $v(x,y)$, then integrate over the cell $\Delta_{j}$, we have:
$$
\frac{\mathrm{d}}{\mathrm{d}t}\int_{\Delta _{j}}u(x,y,t)v(x,y)\mathrm{d}x \mathrm{d}y-\int_{\Delta _{j}}F(u)\cdot \nabla v\mathrm{d}x \mathrm{d}y+\int_{\partial \Delta _{j}}F(u)\cdot n v \mathrm{d}s=0,
$$
where $F=(f,g)$, and $n$ is the outward unit normal of the cell boundary $\partial \Delta_{j}$. The line integral is typically discretized by a **Gaussian quadrature**:
$$
\int_{\partial \Delta _{j}}F\cdot nv \mathrm{d}s \approx | \partial \Delta _{j} | \sum\limits_{k=1}^{q}\omega_{k}F(u(G_{k},t))\cdot n v(G_{k}),
$$
where $F(u(G_{k},t))\cdot n$ is replaced by a **numerical flux**. We can choose the simple **Lax-Friedrichs flux**, which is given by:
$$
F(u(G_{k},t))\cdot n\approx \frac{1}{2}[(F(u^{-}(G_{k},t))+F(u^{+}(G_{k},t)))\cdot n- \alpha(u^{+}(G_{k},t)-u^{-}(G_{k},t))].
$$
where $\alpha$ is taken as an upper bound for the eigenvalues of the Jacobian in the $n$ direction, and $u^{-}$, $u^{+}$ are the values of $u$ inside $\Delta_{j}$ and outside $\Delta_{j}$, respectively.

The cell entropy inequality holds for arbitrary triangulation, and the limiter can also be defined for arbitrary triangulation. For multi-dimensional cases, we can prove the **maximum norm stability** of the limited scheme.

For nonlinear hyperbolic equations including symmetrizable systems, if the solution of the PDE is smooth, $L^{2}$ error estimates of $O(h^{k+1/2}+\Delta t^{2})$ where $\Delta t$ is the time step. For upwind fluxes, the optimal $O(h^{k+1}+\Delta t^{2})$ error estimate can be obtained.

## DG method for convection diffusion models

 In this section, we discuss the time dependent convection diffusion equations:
 $$
u_{t}+\sum\limits_{i=1}^{d}f_{i}(u)_{x_{i}}-\sum\limits_{i=1}^{d}\sum\limits_{j=1}^{d}(a_{ij}(u)u_{x_{j}})_{x_{i}}=0,
$$
where $(a_{ij}(u))$ is a **symmetric** and **semi-positive definite** matrix. We will discuss the local discontinuous Galerkin (LDG) method.

Since the solution space is not regular enough to handle higher derivatives, DG methods cannot be directly applied. The idea of LDG methods for time dependent PDEs with higher derivatives is to rewrite the equation into a first order system, then apply the DG method o the system. A key ingredient for the success of such methods is the correct design of **interface numerical fluxes.** These fluxes must be designed to guarantee **stability** and **local solvability** to approximate the derivatives. 

Next, we will focus on the one dimensional convection diffusion equations and its LDG scheme.

### LDG scheme
Now, we consider the one dimensional **convection diffusion equation**:
$$
u_{t}+(f(u))_{x}=(a(u)u_{x})_{x},
$$
with $a(u)\ge 0$. First, we should split this equation to a first order system.

We should decompose $a(u)u_{x}$ into two factors, marked as $a(u)u_{x}:=A(u)B(u)$, where $A(u)$ is a parameter and $B(u)$ is derived by an auxiliary PDE related to $A(u)$ *FOR UNIFORM NUMERICAL FLUX? I GUESS*. Since $a(u)\ge 0$, we choose 
$$
A(u):=\sqrt{a(u)},\quad B(u):=\int_{}^{u}A(u)\mathrm{d}u,
$$
the system shows:
$$
\begin{align}
u_{t}+(f(u))_{x}&=(A(u)q)_{x}, \\
q&=(B(u))_{x}.
\end{align}
$$
Then apply mixed FEM scheme to this system, for all test functions $v_{h},p_{h}\in V_{h}^{k}$, the semi-discretize scheme shows:
$$
\begin{align}
\int_{I_{i}}(u_{h})_{t}v_{h}\mathrm{d}x-\int_{I_{i}}(f(u)-A(u)q)(v_{h})_{x}\mathrm{d}x+(\hat{f}-\hat{A}\hat{q})|_{i+\frac{1}{2}}(v_{h})^{-}_{i+\frac{1}{2}}-(\hat{f}-\hat{A}\hat{q})_{i-\frac{1}{2}}(v_{h})^{+}_{i-\frac{1}{2}}&=0, \\
\int_{I_{i}}q_{h}p_{h}\mathrm{d}x+\int_{I_{i}}B(u_{h})(p_{h})_{x}\mathrm{d}x-\hat{B}_{i+\frac{1}{2}}(p_{h})_{i+\frac{1}{2}}^{-}+\hat{B}_{i-\frac{1}{2}}(p_{h})_{i-\frac{1}{2}}^{+}&=0.
\end{align}
$$
Here, all the "hat" terms are the numerical fluxes. Now, we discuss *alternating fluxes*, defined as:
$$
\hat{q}:=q_{h}^{+},\quad \hat{B}:=B(u_{h}^{-}), \quad \hat{A}:= \frac{B(u_{h}^{+})-B(u_{h}^{-})}{u_{h}^{+}-u_{h}^{-}};
$$
or
$$
\hat{q}:=q_{h}^{-},\quad \hat{B}:=B(u_{h}^{+}),\quad \hat{A}:=\frac{B(u_{h}^{+})-B(u_{h}^{-})}{u_{h}^{+}-u_{h}^{-}}
$$
is also available. 

From the second equation, we can solve $q_{h}$ *explicitly* and *locally* in terms of $u_{h}$. So this method is referred to as the **local** DG method.
### Stability
In this section, we derive the **cell entropy inequality** for the semi-discrete LDG method.
> [!proposition] Proposition.
The solution $u_{h}$, $q_{h}$ to the semi-discrete LDG scheme satisfies:
$$
\frac{1}{2}\frac{\mathrm{d}}{\mathrm{d}t}\int_{I_{i}}(u_{h})^{2}\mathrm{d}x+\int_{I_{i}}(q_{h})^{2}\mathrm{d}x+\hat{F}_{i+\frac{1}{2}}-\hat{F}_{i-\frac{1}{2}}\le 0,
$$
for some consistent entropy flux
$$
\hat{F}_{i+\frac{1}{2}}:=\hat{F}\left( u_{h}\left( x_{i+\frac{1}{2}}^{-},t \right),q_{h}\left( x_{i+\frac{1}{2}}^{-},t \right);u_{h}\left( x_{i+\frac{1}{2}}^{+},t \right),q_{h}\left( x_{i+\frac{1}{2}}^{+},t \right) \right)
$$
satisfying $\hat{F}(u,u)=F(u)-ub(u)q$, where $F(u):=\int_{}^{u}uf'(u)\mathrm{d}u$.

**Proof:** First, we introduce a short-hand notation
$$
\begin{align}
B_{i}(u_{h},q_{h};v_{h},p_{h})=& \int_{I_{i}}(u_{h})_{t}v_{h}\mathrm{d}x -\int_{I_{i}}(f(u_{h})-A(u_{h})q_{h})(v_{h})_{x}\mathrm{d}x \\
&+(\hat{f}-\hat{A}\hat{q})|_{i+\frac{1}{2}}(v_{h})_{i+\frac{1}{2}}^{-}-(\hat{f}-\hat{A}\hat{q})|_{i-\frac{1}{2}}(v_{h})_{i-\frac{1}{2}}^{+} \\
 \\
&+\int_{I_{i}}q_{h}p_{h}\mathrm{d}x+\int_{I_{i}}B(u_{h})(p_{h})_{x}\mathrm{d}x-\hat{B}_{i+\frac{1}{2}}(q_{h})_{i+\frac{1}{2}}^{-}+\hat{B}_{i-\frac{1}{2}}(q_{h})_{i-\frac{1}{2}}^{+} \\
=&0.
\end{align}
$$
If we take $v_{h}:=u_{h}$, $p_{h}:=q_{h}$, we obtain $B_{i}(u_{h},q_{h};u_{h},q_{h})=0$. Then, mark $\tilde{F}(u):=\int_{}^{u}f(u)\mathrm{d} u$,  
$$
\hat{F}:=-\tilde{F}(u_{h}^{-})+\hat{f}u_{h}^{-}-\hat{b}q_{h}^{+}u_{h}^{-}
$$
and 
$$
\Theta :=-\tilde{F}(u_{h}^{-})+\hat{f}u_{h}^{-}+\tilde{F}(u_{h}^{+})-\hat{f}u_{h}^{+},
$$
we have:
$$
B_{i}(u_{h},q_{h};u_{h},q_{h})=\frac{1}{2}\frac{\mathrm{d}}{\mathrm{d}t}\int_{I_{i}}(u_{h})^{2}\mathrm{d}x+\int_{I_{i}}(q_{h})^{2}\mathrm{d}x+\hat{F}_{i+\frac{1}{2}}-\hat{F}_{i-\frac{1}{2}}+\Theta_{i-\frac{1}{2}}=0.
$$
We readily have $\Theta\ge 0$, which completes the proof.

> [!corollary] Corollary. 
For periodic or compactly supported boundary conditions, the solution $u_{h},q_{h}$ to the semi-discrete LDG scheme satisfies the $L^{2}$ stability:
$$
\frac{\mathrm{d}}{\mathrm{d}t}\int_{0}^{1}(u_{h})^{2}\mathrm{d}x+2 \int_{0}^{1}(q_{h})^{2}\mathrm{d}x\le 0.
$$

### Error Estimates
If we assume the exact solution is smooth, we can obtain optimal $L^{2}$ error estimates. For simplicity, we will give the proof for the heat equation 
$$
u_{t}=u_{xx}
$$
defined on $[0,1]$ with periodic boundary conditions.

> [!proposition] Proposition.
The solution $u_{h}$ and $q_{h}$ to the semi-discrete DG scheme with a smooth solution $u$ satisfies:
$$
\int_{0}^{1}(u-u_{h})^{2}\mathrm{d}x+\int_{0}^{t}\int_{0}^{1}(u_{x}(x,\tau )-q_{h}(x,\tau ))^{2}\mathrm{d}x \mathrm{d}\tau \le Ch^{2(k+1)},
$$
where $C$ depends on $u$ and its derivatives but is independent of $h$.

**Proof:** similar to the convection equation.

