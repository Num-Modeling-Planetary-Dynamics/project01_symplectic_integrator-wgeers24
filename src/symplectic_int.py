import numpy as np
import matplotlib.pyplot as plt
import
#Hamiltonian
\mathcal{H} =\mathcal{H}{Kepler} + \mathcal{H}{Interaction} +\mathcal{H}_{Sun}
#component hamiltonians
\begin{aligned} \mathcal{H}{Kepler} &= \sum{i=i}^{n}\left(\frac{\left|\mathbf{P}i\right|i^2}{2m_i}-\frac{\mathcal{G}m_0m_i}{\mathbf{Q}i}\right)\ \mathcal{H}{Interaction} &= -\sum{i=1}^{n-1}\sum{j=i+1}^{n}\frac{\mathcal{G}m_im_j}{\left|\mathbf{Q}i-\mathbf{Q}j\right|}\ \mathcal{H}{Sun} &= \frac{1}{2m_0}\left|\sum{i=1}^n\mathbf{P}_i\right|^2 \end{aligned}
#standard coordinate variable
\mathbf{q}=\mathbf{r}
#momentum variable
\mathbf{p}=m\mathbf{v}
#algorithim, kick-drift-kick leapfrog integration
E_{Sun}\left(\frac{\Delta t}{2}\right)E_{int}\left(\frac{\Delta t}{2}\right)E_{Kep}\left(\Delta t\right)E_{int}\left(\frac{\Delta t}{2}\right)E_{Sun}\left(\frac{\Delta t}{2}\right)

