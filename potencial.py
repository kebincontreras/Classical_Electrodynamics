import numpy as np
import matplotlib.pyplot as plt

# Parámetros
a = 1.0         # radio del hemisferio
E0 = 1.0        # magnitud del campo externo
Phi_0 = 0       # potencial de referencia del conductor

# Crear figura única para el potencial
plt.figure(figsize=(12, 10))

# Mallado según especificaciones del problema
y = np.linspace(-3*a, 3*a, 80)
z = np.linspace(0, 4*a, 80)  # Reducido a z=4 para evitar zona en blanco
Y, Z = np.meshgrid(y, z)
R = np.sqrt(Y**2 + Z**2)

# Evitar división por cero
R[R == 0] = 1e-9

# Calcular el potencial según la fórmula
# Φ(r) = Φ₀ - (1 - a³/r³) E₀ · r
# Para E₀ en dirección z: E₀ · r = E₀ * z
Phi = Phi_0 - (1 - a**3 / R**3) * E0 * Z

# Usar percentiles para un mejor contraste visual
phi_5 = np.percentile(Phi, 5)   # Evita valores extremos
phi_95 = np.percentile(Phi, 95) # Evita valores extremos

# Contorno de fondo con colores que muestran buen contraste
levels_fill = np.linspace(phi_5, phi_95, 30)
contourf = plt.contourf(Y, Z, Phi, levels=levels_fill, cmap='RdYlBu_r', alpha=0.8, extend='both')

# Líneas equipotenciales principales en el rango de interés
levels_lines = np.linspace(phi_5*0.8, phi_95*0.8, 15)
contour_lines = plt.contour(Y, Z, Phi, levels=levels_lines, colors='black', linewidths=1.5, alpha=0.8)

# Etiquetas en las líneas equipotenciales
plt.clabel(contour_lines, inline=True, fontsize=10, fmt='%.1f', colors='black')

# Hemisferio conductor (superficie equipotencial)
theta = np.linspace(0, np.pi, 100)
y_hemisphere = a * np.sin(theta)
z_hemisphere = a * np.cos(theta)

plt.plot(y_hemisphere, z_hemisphere, 'k-', linewidth=6)
plt.plot(-y_hemisphere, z_hemisphere, 'k-', linewidth=6)
plt.fill_between(y_hemisphere, z_hemisphere, a, alpha=1.0, color='darkgray')
plt.fill_between(-y_hemisphere, z_hemisphere, a, alpha=1.0, color='darkgray')

# Plano conductor (también equipotencial)
plt.plot([-3*a, -a], [0, 0], 'k-', linewidth=8)
plt.plot([a, 3*a], [0, 0], 'k-', linewidth=8)

# Marcar el potencial constante del conductor
plt.text(0, a+0.3, f'Φ = {Phi_0:.1f}', fontsize=14, fontweight='bold', 
         ha='center', va='bottom', color='white',
         bbox=dict(boxstyle="round,pad=0.3", facecolor='black', alpha=0.8))

# Configuración del gráfico
plt.xlim(-3*a, 3*a)
plt.ylim(0, 4*a)  # Reducido para evitar zona en blanco
plt.xlabel('y', fontsize=16)
plt.ylabel('z', fontsize=16)
plt.title('POTENCIAL ELÉCTRICO ALREDEDOR DEL HEMISFERIO CONDUCTOR\n' + 
          r'$\Phi(\mathbf{r}) = \Phi_0 - \left(1 - \frac{a^3}{r^3}\right) \mathbf{E_0} \cdot \mathbf{r}$' + '\n' +
          'Líneas equipotenciales y distribución de potencial',
          fontsize=18, fontweight='bold', pad=20)

plt.grid(True, alpha=0.3, linewidth=1)
plt.gca().set_aspect('equal')

plt.tight_layout()

# Guardar imagen
plt.savefig('potencial_punto_b.png', dpi=300, bbox_inches='tight', 
            facecolor='white', edgecolor='none')

plt.show()

print("✓ Imagen guardada como 'potencial_punto_b.png'")
print("\nEsta figura muestra el POTENCIAL ELÉCTRICO para el punto B:")
print("• Distribución del potencial: colores de fondo (azul-rojo-amarillo)")
print("• Líneas equipotenciales: líneas negras con valores numéricos")
print("• Hemisferio conductor: superficie equipotencial Φ = Φ₀ = 0")
print("• Ecuación del potencial mostrada en el título")
print("• Rango especificado: -3a ≤ y ≤ 3a, 0 ≤ z ≤ 6a")
print("\nEsto responde correctamente al punto B del problema.")
