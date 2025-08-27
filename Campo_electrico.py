import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Definición de la malla de puntos
x = np.linspace(-2, 2, 20)
y = np.linspace(-2, 2, 20)
X, Y = np.meshgrid(x, y)

# Función que calcula el campo eléctrico de una carga puntual en el origen
def campo_electrico(X, Y, q=1.0, k=1.0, eps=1e-6):
    Rx = X
    Ry = Y
    # Evitar división por cero en el origen
    r2 = Rx**2 + Ry**2
    r2_safe = r2 + (r2 == 0) * eps
    R3 = r2_safe**1.5
    Ex = k * q * Rx / R3
    Ey = k * q * Ry / R3
    # Forzamos campo cero exactamente en el origen para evitar inf/nan
    Ex = np.where(r2 == 0, 0.0, Ex)
    Ey = np.where(r2 == 0, 0.0, Ey)
    return Ex, Ey

# Calculamos el campo completo una sola vez
U_full, V_full = campo_electrico(X, Y)

# Configuración de la figura y el quiver
fig, ax = plt.subplots(figsize=(6,6))
ax.set_xlim(-2.2, 2.2)
ax.set_ylim(-2.2, 2.2)
ax.set_aspect('equal')
#ax.set_title("Construcción paso a paso del campo eléctrico")

# Preparamos posiciones fijas para todos los puntos; los componentes se actualizarán paso a paso
positions = np.column_stack((X.flatten(), Y.flatten()))
N = positions.shape[0]

# Calculamos la divergencia numérica: div(E) = dEx/dx + dEy/dy -> rho = eps0 * div(E)
dx = x[1] - x[0]
dy = y[1] - y[0]
# dEx/dx: derivada de U_full respecto a x (axis=1), dEy/dy: derivada de V_full respecto a y (axis=0)
dEx_dx = np.gradient(U_full, dx, axis=1)
dEy_dy = np.gradient(V_full, dy, axis=0)
divE = dEx_dx + dEy_dy
eps0 = 1.0  # constante arbitraria para visualización; si quieres valores físicos usa 8.854e-12
rho = eps0 * divE
# Reemplazar valores no finitos para evitar que el colormap pinte 'bad' como negro
rho = np.nan_to_num(rho, nan=0.0, posinf=0.0, neginf=0.0)

# Mostrar rho como mapa de color y contorno rho=0 (líneas rectas donde corresponde)
vmin = float(np.min(rho))
vmax = float(np.max(rho))
im = ax.pcolormesh(X, Y, rho, cmap='seismic', shading='auto', alpha=0.6, vmin=vmin, vmax=vmax)
cb = fig.colorbar(im, ax=ax, label=r'$\rho$ (arb.)')
# Etiquetas de la colorbar indicando sumidero/fuente (sin valores numéricos)
cb.set_ticks([vmin, vmax])
cb.set_ticklabels(['ρ<0 (sumidero)', 'ρ>0 (fuente)'])
cont = ax.contour(X, Y, rho, levels=[0], colors='k', linewidths=1)

# Mostrar la ecuación del campo en dos líneas usando mathtext simple (sin \begin{array})
fig.text(0.5, 0.98, r'$\mathbf{E}(\mathbf{r}) = k\,q\;\dfrac{\mathbf{r}}{|\mathbf{r}|^3}$', ha='center', va='top', fontsize=11)
fig.text(0.5, 0.945, r'$E_x = k\,q\dfrac{x}{(x^2+y^2)^{3/2}}$', ha='center', va='top', fontsize=10)
# Ajustar layout para evitar solapamientos y dejar espacio para el texto superior
fig.tight_layout(rect=[0,0,1,0.90])

# Inicializa el quiver con todos los offsets fijos y componentes inicialmente en cero
# Los colores de las flechas estarán dados por los valores de rho en cada posición
rho_flat = rho.flatten()
# Create quiver without vmin/vmax (those are not valid kwargs); set array and clim after
quiver = ax.quiver(positions[:,0], positions[:,1], np.zeros(N), np.zeros(N), rho_flat, cmap='seismic', angles='xy', scale_units='xy', scale=5)
# Ensure the quiver uses the same color range as the rho map
quiver.set_array(rho_flat)
quiver.set_clim(vmin, vmax)

# Texto informativo para mostrar rho en el punto actual
info_text = ax.text(-2.0, 2.0, '', fontsize=10, va='top', bbox=dict(facecolor='white', alpha=0.7))

# Función de inicialización: componentes en cero
def init():
    quiver.set_offsets(positions)
    quiver.set_UVC(np.zeros(N), np.zeros(N))
    info_text.set_text('')
    return quiver, info_text

# Función de actualización: va activando vectores hasta el frame actual y muestra rho en el punto
def update(frame):
    # Creamos arrays de componentes con ceros y activamos hasta el índice "frame"
    U_arr = np.zeros(N)
    V_arr = np.zeros(N)
    end = min(frame, N-1)
    U_flat = U_full.flatten()
    V_flat = V_full.flatten()
    U_arr[:end+1] = U_flat[:end+1]
    V_arr[:end+1] = V_flat[:end+1]
    quiver.set_UVC(U_arr, V_arr)
    # Mostrar rho en el punto actual
    rho_flat = rho.flatten()
    rho_val = rho_flat[end]
    info_text.set_text(f'index={end}\n rho={rho_val:.3e}')
    return quiver, info_text

# Número total de vectores
total_vectors = X.size

# Creamos la animación
ani = animation.FuncAnimation(
    fig,
    update,
    frames=total_vectors,
    init_func=init,
    blit=True,
    interval=100,   # milisegundos entre frames
    repeat=False
)

# Mostrar la figura interactiva primero
plt.show()

# Después de cerrar la ventana, opcionalmente guardar la animación como GIF
SAVE_GIF = True
if SAVE_GIF:
    try:
        ani.save('campo.gif', writer=animation.PillowWriter(fps=10))
        print("GIF guardado como 'campo.gif'")
    except Exception as e:
        print('No se pudo guardar GIF:', e)
