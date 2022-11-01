%matplotlib inline
from orbitgeometry import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
import xarray as xr
rc('text', usetex=True)
rc('text.latex', preamble=r'\usepackage{amsmath}')

a = 1.0
mu = 1.0
e = 0.3
P = 2 * np.pi * np.sqrt(a**3/mu)
rp = a * (1.0 - e)
vp = np.sqrt(mu * (2.0 / rp - 1.0 / a))
y0 = np.array([rp,0.0,0.0,vp])
t0 = 0.0
nperiods = 300
tend = nperiods * P
def step_one_leapfrog(yold,h):
    xold = yold[0:2]
    vold = yold[2:4]
    xhalf = xold + 0.5 * h * vold
    vdot = - mu * xhalf / np.sqrt((xhalf[0] ** 2 + xhalf[1] ** 2) ** 3)
    vnew = vold + h * vdot
    xnew = xhalf + 0.5 * h * vnew
    return np.concatenate((xnew, vnew))

step_vals = np.array([0.01,0.1/3.,0.1,1./3.])
hvals = step_vals * P
sim_results_all = []
conserved_results_all = []
for h in hvals:
    y = y0
    tvals = np.arange(start=t0, stop=tend+h, step=h)
    itheta = np.ceil(P/h).astype(np.int64)
    Energy_vals = []
    delta_theta_vals = []
    h_vals = []
    yhist = []
    thist = []

for i,t in enumerate(tvals):
    yhist.append(y)
    thist.append(t)
    if i % itheta == 0:
        rmag = np.sqrt(y[0]**2 + y[1]**2)
        vmag2 = y[2]**2 + y[3]**2
        Energy_vals.append(vmag2 / 2.0 - 1.0 / rmag)
        hmag = np.cross(y[0:2],y[2:4])
        h_vals.append(np.sqrt(np.dot(hmag,hmag)))
        delta_theta_vals.append(np.arctan2(y[1],y[0]))
        y = step_one_leapfrog(y,h)

sim_results = xr.DataArray(data=yhist,dims=["t","vec"],
            coords={"t": ("t", thist),"vec": ("vec", ["r_x","r_y","v_x","v_y"])})
sim_results = sim_results.to_dataset(dim="vec")
sim_results_all.append(sim_results)

conserved_results = xr.Dataset(
    {
        "E": (["nP"], Energy_vals),
        "Eerr": (["nP"],np.abs((Energy_vals - Energy_vals[0])/Energy_vals[0])),
        "hmag": (["nP"], h_vals),
        "hmagerr": (["nP"],np.abs((h_vals - h_vals[0])/h_vals[0])),
        "DeltaTheta": (["nP"], delta_theta_vals),
    },
    coords={
        "nP": ("nP", np.arange(start=0,stop=nperiods+1,dtype=int))
    },
    )
    conserved_results_all.append(conserved_results)
fig, axs = plt.subplots(ncols=2,nrows=2,figsize=(12,12))
axes_fontsize = 16
legend_fontsize = 12

def rorbit(theta,a,e):
    p = a * (1 - e**2)
    return p / (1 + e * np.cos(theta))

theta = np.linspace(start=0.0,stop=2*np.pi,num=300,endpoint=True)
r = rorbit(theta,a,e)
xorbit = r * np.cos(theta)
yorbit = r * np.sin(theta)

axs[0][0].plot(xorbit,yorbit,label="True",color='black',linewidth=3,linestyle="--")
for i,ds in enumerate(sim_results_all):
    r_x = ds['r_x'].values
    r_y = ds['r_y'].values
    axs[0][0].plot(r_x,r_y,label=f"h = {step_vals[i]:.2f}P")
axs[0][0].legend(fontsize=legend_fontsize,loc='lower right')
axs[0][0].set_xlabel("$r_x$",fontsize=axes_fontsize)
axs[0][0].set_ylabel("$r_y$",fontsize=axes_fontsize)
axs[0][0].set_xlim((-2,2))
axs[0][0].set_ylim((-2,2))

for i,ds in enumerate(conserved_results_all):
    ds['Eerr'].plot.line(ax=axs[0][1])
axs[0][1].set_xlabel("$nP$",fontsize=axes_fontsize)
axs[0][1].set_ylabel("$\\left|(\mathcal{E} - \mathcal{E}_0)/\mathcal{E}_0\\right|$",fontsize=axes_fontsi
axs[0][1].set_yscale('log')

for i,ds in enumerate(conserved_results_all):
    ds['hmagerr'].plot.line(ax=axs[1][0])
axs[1][0].set_xlabel("$nP$",fontsize=axes_fontsize)
axs[1][0].set_ylabel("$\\left|(\mathbf{h} - \mathbf{h}_0)/\mathbf{h}_0\\right|$",fontsize=axes_fontsize)
axs[1][0].set_yscale('log')

for i,ds in enumerate(conserved_results_all):
    ds['DeltaTheta'].plot.line(ax=axs[1][1])
axs[1][1].set_xlabel("$nP$",fontsize=axes_fontsize)
axs[1][1].set_ylabel("$\Delta \\theta $",fontsize=axes_fontsize)
for axr in axs:
    for ax in axr:
        ax.tick_params(labelsize=axes_fontsize)
plt.show()
