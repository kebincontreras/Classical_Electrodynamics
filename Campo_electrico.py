import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

x = np.linspace(-2, 2, 20)
y = np.linspace(-2, 2, 20)
X, Y = np.meshgrid(x, y)

def campo_electrico(X, Y, q=1.0, k=1.0, eps=1e-6):
    Rx = X
    Ry = Y
    r2 = Rx**2 + Ry**2
    r2_safe = r2 + (r2 == 0) * eps
    R3 = r2_safe**1.5
    Ex = k * q * Rx / R3
    Ey = k * q * Ry / R3
    Ex = np.where(r2 == 0, 0.0, Ex)
    Ey = np.where(r2 == 0, 0.0, Ey)
    return Ex, Ey

U_full, V_full = campo_electrico(X, Y)

fig, ax = plt.subplots(figsize=(6,6))
ax.set_xlim(-2.2, 2.2)
ax.set_ylim(-2.2, 2.2)
ax.set_aspect('equal')

positions = np.column_stack((X.flatten(), Y.flatten()))
N = positions.shape[0]

dx = x[1] - x[0]
dy = y[1] - y[0]
dEx_dx = np.gradient(U_full, dx, axis=1)
dEy_dy = np.gradient(V_full, dy, axis=0)
divE = dEx_dx + dEy_dy
eps0 = 1.0
rho = eps0 * divE
rho = np.nan_to_num(rho, nan=0.0, posinf=0.0, neginf=0.0)

vmin = float(np.min(rho))
vmax = float(np.max(rho))
im = ax.pcolormesh(X, Y, rho, cmap='seismic', shading='auto', alpha=0.6, vmin=vmin, vmax=vmax)
cb = fig.colorbar(im, ax=ax, label=r'$\rho$ (arb.)')
cb.set_ticks([vmin, vmax])
cb.set_ticklabels(['ρ<0 (sumidero)', 'ρ>0 (fuente)'])
cont = ax.contour(X, Y, rho, levels=[0], colors='k', linewidths=1)

fig.text(0.5, 0.98, r'$\mathbf{E}(\mathbf{r}) = k\,q\;\dfrac{\mathbf{r}}{|\mathbf{r}|^3}$', ha='center', va='top', fontsize=11)
fig.text(0.5, 0.945, r'$E_x = k\,q\dfrac{x}{(x^2+y^2)^{3/2}}$', ha='center', va='top', fontsize=10)
fig.tight_layout(rect=[0,0,1,0.90])

rho_flat = rho.flatten()
quiver = ax.quiver(positions[:,0], positions[:,1], np.zeros(N), np.zeros(N), rho_flat, cmap='seismic', angles='xy', scale_units='xy', scale=5)
quiver.set_array(rho_flat)
quiver.set_clim(vmin, vmax)

info_text = ax.text(-2.0, 2.0, '', fontsize=10, va='top', bbox=dict(facecolor='white', alpha=0.7))

def init():
    quiver.set_offsets(positions)
    quiver.set_UVC(np.zeros(N), np.zeros(N))
    info_text.set_text('')
    return quiver, info_text

def update(frame):
    U_arr = np.zeros(N)
    V_arr = np.zeros(N)
    end = min(frame, N-1)
    U_flat = U_full.flatten()
    V_flat = V_full.flatten()
    U_arr[:end+1] = U_flat[:end+1]
    V_arr[:end+1] = V_flat[:end+1]
    quiver.set_UVC(U_arr, V_arr)
    rho_flat = rho.flatten()
    rho_val = rho_flat[end]
    info_text.set_text(f'index={end}\n rho={rho_val:.3e}')
    return quiver, info_text

total_vectors = X.size

ani = animation.FuncAnimation(
    fig,
    update,
    frames=total_vectors,
    init_func=init,
    blit=True,
    interval=100,
    repeat=False
)

plt.show()

SAVE_GIF = True
if SAVE_GIF:
    try:
        ani.save('campo.gif', writer=animation.PillowWriter(fps=10))
        print("GIF guardado como 'campo.gif'")
    except Exception as e:
        print('No se pudo guardar GIF:', e)

