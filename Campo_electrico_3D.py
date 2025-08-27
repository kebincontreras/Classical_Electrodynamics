import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

q = 1.0
k = 1.0
eps = 1e-6

n = 7
x = np.linspace(-1.5, 1.5, n)
y = np.linspace(-1.5, 1.5, n)
z = np.linspace(-1.5, 1.5, n)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

def campo_electrico_3d(X, Y, Z, q=1.0, k=1.0, eps=1e-6):
    Rx = X
    Ry = Y
    Rz = Z
    r2 = Rx**2 + Ry**2 + Rz**2
    r2_safe = r2 + (r2 == 0) * eps
    R3 = r2_safe**1.5
    Ex = k * q * Rx / R3
    Ey = k * q * Ry / R3
    Ez = k * q * Rz / R3
    Ex = np.where(r2 == 0, 0.0, Ex)
    Ey = np.where(r2 == 0, 0.0, Ey)
    Ez = np.where(r2 == 0, 0.0, Ez)
    return Ex, Ey, Ez

Ux, Uy, Uz = campo_electrico_3d(X, Y, Z, q=q, k=k, eps=eps)

dx = x[1] - x[0]
dy = y[1] - y[0]
dz = z[1] - z[0]

dUx_dx = np.gradient(Ux, dx, axis=0)
dUy_dy = np.gradient(Uy, dy, axis=1)
dUz_dz = np.gradient(Uz, dz, axis=2)
divE = dUx_dx + dUy_dy + dUz_dz

eps0 = 1.0
rho = eps0 * divE
rho = np.nan_to_num(rho, nan=0.0, posinf=0.0, neginf=0.0)

Xf = X.flatten()
Yf = Y.flatten()
Zf = Z.flatten()
Uf = Ux.flatten()
Vf = Uy.flatten()
Wf = Uz.flatten()
Rhof = rho.flatten()
N = Xf.size

vmin = float(np.min(Rhof))
vmax = float(np.max(Rhof))
if vmin == vmax:
    vmin -= 1e-6
    vmax += 1e-6
cmap = plt.cm.seismic

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(x.min(), x.max())
ax.set_ylim(y.min(), y.max())
ax.set_zlim(z.min(), z.max())
ax.set_box_aspect([1,1,1])
ax.set_title('Campo eléctrico 3D: flechas coloreadas por \rho (num.)')

quiver = None
scat = ax.scatter([], [], [], s=0)

from matplotlib.cm import ScalarMappable
sm = ScalarMappable(cmap=cmap)
sm.set_clim(vmin, vmax)
cb = fig.colorbar(sm, ax=ax, shrink=0.6, pad=0.1)
cb.set_ticks([vmin, vmax])
cb.set_ticklabels(['ρ<0 (sumidero)', 'ρ>0 (fuente)'])

fig.text(0.5, 0.94, r'$\mathbf{E}(\mathbf{r}) = k q \; \mathbf{r} / |\mathbf{r}|^3$', ha='center', fontsize=10)

def update(frame):
    global quiver
    if quiver is not None:
        try:
            quiver.remove()
        except Exception:
            if hasattr(quiver, 'collections'):
                for coll in quiver.collections:
                    try:
                        coll.remove()
                    except Exception:
                        pass
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(y.min(), y.max())
    ax.set_zlim(z.min(), z.max())
    ax.set_box_aspect([1,1,1])
    ax.set_title('Campo eléctrico 3D: flechas coloreadas por \\rho (num.)')

    end = min(frame, N-1)
    Xs = Xf[:end+1]
    Ys = Yf[:end+1]
    Zs = Zf[:end+1]
    Us = Uf[:end+1]
    Vs = Vf[:end+1]
    Ws = Wf[:end+1]
    Rs = Rhof[:end+1]

    normed = (Rs - vmin) / (vmax - vmin)
    normed = np.clip(normed, 0, 1)
    colors = cmap(normed)

    quiver = ax.quiver(Xs, Ys, Zs, Us, Vs, Ws, length=0.2, colors=colors, normalize=True)

    return quiver,

total_frames = N
ani = animation.FuncAnimation(fig, update, frames=total_frames, interval=100, blit=False)

plt.show()

SAVE_GIF = True
if SAVE_GIF:
    try:
        ani.save('campo_3d.gif', writer=animation.PillowWriter(fps=8))
        print("GIF guardado como 'campo_3d.gif'")
    except Exception as e:
        print('No se pudo guardar GIF:', e)
