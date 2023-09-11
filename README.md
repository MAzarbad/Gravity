# Gravity
Suppose $`N`$ objects with masses $`m_1,m_2,\ldots,m_N`$ and locations(at time $`t`$) $`P_1(t),P_2(t),\ldots,P_N(t)`$.
```math
P_i(t) = \big(x_i(t),y_i(t),z_i(t)\big) \quad \quad i=1,2,\ldots,N
```
the $`j`$th object feels a force (due to gravity) equal to:

```math
F_j(t) = Gm_j \sum_{\overset{i=1}{i \neq j}}^N m_i \frac{P_i(t)-P_j(t)}{|P_i(t)-P_j(t)|^3}
```
Where $`G=6.674 \times 10^{-11}`$. on the other hand we know that
```math
F_j(t) = m_j a_j(t) =  m_j P_j{''}(t)
```
therefore
```math
P_j{''}(t) = G\sum_{\overset{i=1}{i \neq j}}^N m_i \frac{P_i(t)-P_j(t)}{|P_i(t)-P_j(t)|^3} \quad \quad j=1,2,\ldots,N
```

so for finding $`P_1(t), P_2(t), \ldots,P_N(t)`$ we need to solve a System of $`N`$differential equations.
for solving the equation let assume that we have locations and speeds at time zero:
```math
P_1(0), P_2(0), \ldots,P_N(0) \quad \text{and} \quad P_1'(0), P_2'(0), \ldots,P_N'(0)
```
for a small number $`\Delta>0`$
```math
P_j'((k+1)\Delta) \approx \Delta P_j{''}(k \Delta) + P_j'(k\Delta)= \Delta  G \sum_{i \neq j} m_i \frac{P_i(k \Delta)-P_j(k \Delta)}{|P_i(k \Delta)-P_j(k \Delta)|^3}+ P_j'(k\Delta)
```
```math
P_j((k+1)\Delta) \approx \Delta P_j'(k \Delta) + P_j(k\Delta)
```
