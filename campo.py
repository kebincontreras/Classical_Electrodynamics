import numpy as np
import matplotlib.pyplot as plt

# Parámetros
a = 1.0         # radio del hemisferio
E0 = 1.0        # magnitud del campo externo

# Crear figura única
plt.figure(figsize=(12, 10))

# Mallado según especificaciones del problema
y = np.linspace(-3*a, 3*a, 60)
z = np.linspace(0, 6*a, 60)
Y, Z = np.meshgrid(y, z)
R = np.sqrt(Y**2 + Z**2)

# Evitar división por cero
R[R == 0] = 1e-9

# Campo eléctrico según la derivación completa del potencial
# 
# DERIVACIÓN DEL CAMPO ELÉCTRICO:
# ================================
# El campo eléctrico se obtiene como: E(r) = −∇Φ(r)
# 
# Partiendo del potencial: Φ(r) = Φ₀ - (1 - a³/r³) E₀ · r
# 
# Definimos T(r) = (1 - a³/r³) E₀ · r y calculamos gradientes:
# ∇(E₀ · r) = E₀
# ∇(1/r³) = ∇(r⁻³) = -3r⁻⁴∇r = -3r/r⁵
# ∇(a³/r³) = a³∇(1/r³) = -3a³ r/r⁵
# 
# Aplicando regla del producto para ∇[(a³/r³) E₀ · r]:
# ∇[(a³/r³) E₀ · r] = (a³/r³) ∇(E₀ · r) + (E₀ · r)∇(a³/r³)
#                   = (a³/r³) E₀ + (E₀ · r)(-3a³ r/r⁵)
#                   = (a³/r³) E₀ - 3a³ (E₀ · r)r/r⁵
# 
# Por tanto, el gradiente completo es:
# ∇T(r) = ∇(E₀ · r) - ∇[(a³/r³) E₀ · r]
#       = E₀ - (a³/r³) E₀ + 3a³ (E₀ · r)r/r⁵
# 
# El campo eléctrico resulta:
# E(r) = −∇Φ(r) = ∇T(r) = E₀ - (a³/r³) E₀ + 3a³ (E₀ · r)r/r⁵
# 
# Factorizando: E(r) = E₀ + (a³/r³)[3r̂(r̂ · E₀) − E₀]
# 
# ANÁLISIS FÍSICO:
# - Lejos del hemisferio (r ≫ a): a³/r³ → 0, entonces E → E₀
# - Cerca del hemisferio (r ∼ a): término dipolar ∝ a³/r³ domina
# - En el polo (r̂ ∥ E₀): E = E₀ + [3E₀ - E₀] = 3E₀ (efecto punta)
# - El término a³/r³[3r̂(r̂ · E₀) − E₀] es la perturbación dipolar
#
Ey = (3 * a**3 * Y * Z / R**5) * E0
Ez = E0 * (1 - a**3 / R**3 + 3 * a**3 * Z**2 / R**5)

# Magnitud del campo
Emag = np.sqrt(Ey**2 + Ez**2)

# Normalizar vectores para mostrar dirección
Ey_norm = np.divide(Ey, Emag, out=np.zeros_like(Ey), where=(Emag != 0) & np.isfinite(Emag))
Ez_norm = np.divide(Ez, Emag, out=np.zeros_like(Ez), where=(Emag != 0) & np.isfinite(Emag))

# Contornos de magnitud del campo (colores de fondo)
contour = plt.contourf(Y, Z, Emag, levels=25, cmap='plasma', alpha=0.9)
plt.contour(Y, Z, Emag, levels=12, colors='white', alpha=0.4, linewidths=0.8)

# Vectores del campo eléctrico (dirección)
skip = 4  # Saltar puntos para claridad visual
plt.quiver(Y[::skip, ::skip], Z[::skip, ::skip], 
           Ey_norm[::skip, ::skip], Ez_norm[::skip, ::skip], 
           color='white', scale=28, width=0.004, alpha=0.95)

# Hemisferio conductor
theta = np.linspace(0, np.pi, 100)
y_hemisphere = a * np.sin(theta)
z_hemisphere = a * np.cos(theta)

plt.plot(y_hemisphere, z_hemisphere, 'k-', linewidth=5)
plt.plot(-y_hemisphere, z_hemisphere, 'k-', linewidth=5)
plt.fill_between(y_hemisphere, z_hemisphere, a, alpha=1.0, color='black')
plt.fill_between(-y_hemisphere, z_hemisphere, a, alpha=1.0, color='black')

# Plano conductor
plt.plot([-3*a, -a], [0, 0], 'k-', linewidth=7)
plt.plot([a, 3*a], [0, 0], 'k-', linewidth=7)

# Campo externo E₀ (flechas indicativas)
plt.arrow(-2.7*a, 5.2*a, 0, 0.6*a, head_width=0.25*a, head_length=0.2*a, 
          fc='red', ec='red', linewidth=3, alpha=0.9)
plt.arrow(2.7*a, 5.2*a, 0, 0.6*a, head_width=0.25*a, head_length=0.2*a, 
          fc='red', ec='red', linewidth=3, alpha=0.9)
plt.text(-2.5*a, 6.0*a, r'$\vec{E_0}$', fontsize=16, color='red', fontweight='bold')
plt.text(2.3*a, 6.0*a, r'$\vec{E_0}$', fontsize=16, color='red', fontweight='bold')

# Configuración del gráfico
plt.xlim(-3*a, 3*a)
plt.ylim(0, 6*a)
plt.xlabel('y', fontsize=16)
plt.ylabel('z', fontsize=16)
plt.title('CÁLCULO DEL CAMPO ELÉCTRICO ALREDEDOR DEL HEMISFERIO CONDUCTOR\n' + 
          r'$\mathbf{E}(\mathbf{r}) = -\nabla\Phi(\mathbf{r}) = \mathbf{E}_0 + \frac{a^3}{r^3}[3\hat{\mathbf{r}}(\hat{\mathbf{r}} \cdot \mathbf{E}_0) - \mathbf{E}_0]$',
          fontsize=14, fontweight='bold', pad=20)

plt.grid(True, alpha=0.4, linewidth=1)
plt.gca().set_aspect('equal')

plt.tight_layout()

# Guardar imagen
plt.savefig('campo_electrico_hemisferio.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')

plt.show()

print("✓ Imagen guardada como 'campo_electrico_hemisferio.png'")
print("\nEsta figura muestra la derivación completa del campo eléctrico:")
print("• Fórmula final: E(r) = E₀ + (a³/r³)[3r̂(r̂ · E₀) − E₀]")
print("• Magnitud del campo: colores (plasma colormap)")
print("• Dirección del campo: vectores blancos normalizados")
print("• Hemisferio conductor: región negra sólida")
print("• Campo externo E₀: flechas rojas")
print("• Comportamiento físico:")
print("  - Lejos: r ≫ a → E ≈ E₀ (campo uniforme)")
print("  - Cerca: r ∼ a → perturbación dipolar significativa")
print("  - Polo: E = 3E₀ (intensificación por efecto punta)")
print("• Rango: -3a ≤ y ≤ 3a, 0 ≤ z ≤ 6a")
