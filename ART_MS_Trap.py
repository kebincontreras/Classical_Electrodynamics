#!/usr/bin/env python3
"""
TRAMPA ELECTROSTÁTICA ART-MS - Réplica exacta de la Figura 2

Este programa replica exactamente la Figura 2 del artículo:
"Numerical Simulation of Autoresonant Ion Oscillations in an Anharmonic Electrostatic Trap"

Geometría: Trampa cilíndrica con electrodo washer central a -500V y electrodos cilíndricos a tierra
Coordenadas: (r, z) con simetría axial, donde y=0 corresponde al plano r-z
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def main():
    # Parámetros para replicar exactamente la Figura 2 del artículo
    nr = 150                    # Puntos radiales (alta resolución como el artículo)
    nz = 400                    # Puntos axiales (alta resolución como el artículo)
    itmax = 4000                # Iteraciones máximas
    
    # Geometría EXACTA según el artículo (Figura 1b)
    R_max = 6e-3                # Radio máximo: 6 mm
    L_total = 16e-3             # Longitud total: 16 mm
    
    # Washer electrode: ubicado en el centro, con ancho de 2mm (como en la figura)
    washer_width = 2e-3         # Ancho de 2mm (1mm a cada lado del centro)
    washer_r_max = 4e-3         # Radio máximo de la washer: 4mm (visible en figura)
    
    V_washer = -500.0           # Voltaje del washer: -500V (como en el artículo)
    V_cylinder = 0.0            # Electrodos cilíndricos: 0V (tierra)
    
    # Voltajes
    V_washer = -500.0           # Electrodo washer central: -500V
    V_cylinder = 0.0            # Electrodos cilíndricos: 0V (tierra)
    
    # Discretización
    hr = R_max / (nr - 1)       # Paso radial
    hz = L_total / (nz - 1)     # Paso axial
    
    tol = 1e-6                  # Tolerancia
    
    # Arrays coordenadas
    r = np.linspace(0, R_max, nr)
    z = np.linspace(0, L_total, nz)
    
    # Matriz de potencial
    phi = np.zeros((nr, nz))
    phi_old = np.zeros((nr, nz))
    
    # CONDICIONES DE FRONTERA según geometría EXACTA del artículo (Figura 2)
    # Electrodo washer central: ubicado exactamente en el centro (z = 8mm)
    z_center_idx = int(0.5 * nz)                      # Centro exacto en z
    washer_half_width_idx = int(washer_width/(2*hz))  # Mitad del ancho en índices
    
    z_washer_start = z_center_idx - washer_half_width_idx  # Inicio del washer
    z_washer_end = z_center_idx + washer_half_width_idx    # Fin del washer  
    r_washer_max = int(washer_r_max / hr)                  # Radio máximo en índices
    
    # Inicializar condiciones de frontera
    apply_boundary_conditions(phi, nr, nz, V_washer, V_cylinder, 
                            z_washer_start, z_washer_end, r_washer_max)
    
    print("Resolviendo ecuación de Laplace cilíndrica para trampa ART-MS...")
    print("Geometría: R_max = 6mm, L = 16mm, V_washer = -500V")
    
    # Resolver ecuación de Laplace cilíndrica
    phi = solve_laplace_cylindrical(phi, nr, nz, hr, hz, tol, itmax,
                                   z_washer_start, z_washer_end, r_washer_max,
                                   V_washer, V_cylinder)
    
    # Crear figura exactamente como en el artículo
    create_figure_2_replica(r, z, phi, R_max, L_total, hr, hz, 
                          z_washer_start, z_washer_end, r_washer_max)
    
    print("✓ Figura 2 del artículo ART-MS replicada exitosamente!")

def apply_boundary_conditions(phi, nr, nz, V_washer, V_cylinder, 
                            z_washer_start, z_washer_end, r_washer_max):
    """Aplicar condiciones de frontera de la trampa ART-MS con geometría rectangular"""
    
    # Electrodo washer central rectangular (como en la imagen)
    # Forma rectangular sólida desde r=0 hasta r_washer_max
    for j in range(z_washer_start, z_washer_end + 1):
        for i in range(r_washer_max + 1):
            phi[i, j] = V_washer
    
    # Extensiones verticales del washer (paredes laterales)
    for i in range(r_washer_max + 1):
        phi[i, z_washer_start] = V_washer
        phi[i, z_washer_end] = V_washer
    
    # Electrodos cilíndricos externos (r = R_max, todos los z)
    for j in range(nz):
        phi[nr-1, j] = V_cylinder
    
    # Extremos axiales (z = 0 y z = L) - conectados a tierra
    for i in range(nr):
        phi[i, 0] = V_cylinder
        phi[i, nz-1] = V_cylinder
    
    # Eje de simetría (r = 0): condición natural, se maneja en la iteración

def solve_laplace_cylindrical(phi, nr, nz, hr, hz, tol, itmax,
                            z_washer_start, z_washer_end, r_washer_max,
                            V_washer, V_cylinder):
    """Resolver usando SOR adaptativo para alta resolución y estabilidad"""
    
    phi_old = np.zeros_like(phi)
    
    # SOR adaptativo: empezar conservador, acelerar gradualmente
    omega_start = 0.5   # Factor inicial conservador
    omega_max = 1.4     # Factor máximo (clásico para SOR)
    
    for iteration in range(itmax):
        phi_old[:] = phi[:]
        max_error = 0.0
        
        # Calcular omega adaptativo
        progress = min(iteration / (itmax * 0.3), 1.0)  # 30% del camino para omega máximo
        omega = omega_start + progress * (omega_max - omega_start)
        
        hr_hz_ratio2 = (hr / hz)**2
        hz_hr_ratio2 = (hz / hr)**2
        
        # Para i ≠ 0 usando el método SOR mejorado
        for i in range(1, nr-1):
            for j in range(1, nz-1):
                if not is_electrode_point(i, j, z_washer_start, z_washer_end, r_washer_max):
                    
                    # Coeficientes exactos del enunciado
                    alpha = 1.0 / (2.0 * i)
                    
                    # Términos de la ecuación diferencial cilíndrica
                    A = (1.0 + alpha) / (hr**2)  # Coeficiente de phi[i+1,j]
                    B = (1.0 - alpha) / (hr**2)  # Coeficiente de phi[i-1,j]  
                    C = 1.0 / (hz**2)            # Coeficiente de phi[i,j±1]
                    
                    # Forma estándar: (A+B+2C)φ[i,j] = A·φ[i+1,j] + B·φ[i-1,j] + C·(φ[i,j+1]+φ[i,j-1])
                    denominator = A + B + 2.0*C
                    phi_new = (A * phi[i+1, j] + B * phi[i-1, j] + C * (phi[i, j+1] + phi[i, j-1])) / denominator
                    
                    # SOR update
                    phi[i, j] = (1.0 - omega) * phi[i, j] + omega * phi_new
                    
                    error = abs(phi[i, j] - phi_old[i, j])
                    max_error = max(max_error, error)
        
        # Para i = 0 (eje de simetría) con SOR
        for j in range(1, nz-1):
            if not is_electrode_point(0, j, z_washer_start, z_washer_end, r_washer_max):
                
                # En r=0, usar L'Hôpital: ∂²φ/∂r² + (1/r)∂φ/∂r → 2∂²φ/∂r²
                # Por lo tanto: 2∂²φ/∂r² + ∂²φ/∂z² = 0
                A = 2.0 / (hr**2)    # Coeficiente doble por L'Hôpital
                C = 1.0 / (hz**2)    # Coeficiente normal en z
                
                denominator = A + 2.0*C
                phi_new = (A * phi[1, j] + C * (phi[0, j+1] + phi[0, j-1])) / denominator
                
                # SOR update
                phi[0, j] = (1.0 - omega) * phi[0, j] + omega * phi_new
        
        # Condiciones de frontera después de cada iteración
        apply_boundary_conditions(phi, nr, nz, V_washer, V_cylinder,
                                z_washer_start, z_washer_end, r_washer_max)
        
        # Reportar progreso
        if iteration % 200 == 0:
            print(f"Iteración {iteration}, Error máximo: {max_error:.2e}, ω = {omega:.3f}")
        
        if max_error < tol:
            print(f"Convergencia alcanzada en {iteration} iteraciones con ω = {omega:.3f}")
            break
            
        # Control de estabilidad mejorado
        if max_error > 1e6:
            print("¡Inestabilidad detectada! Reduciendo factor de relajación...")
            omega = max(0.3, omega * 0.8)  # Reducir omega dinámicamente
            if omega < 0.4:
                print("Reiniciando con condiciones conservadoras...")
                phi.fill(0.0)
                apply_boundary_conditions(phi, nr, nz, V_washer, V_cylinder,
                                        z_washer_start, z_washer_end, r_washer_max)
                omega = 0.3  # Restart conservativo
    
    return phi

def is_electrode_point(i, j, z_washer_start, z_washer_end, r_washer_max):
    """Determinar si un punto está en un electrodo (condición de frontera fija)"""
    
    # Electrodo washer rectangular (toda la región)
    if (z_washer_start <= j <= z_washer_end and i <= r_washer_max):
        return True
    
    # Paredes laterales del washer
    if ((j == z_washer_start or j == z_washer_end) and i <= r_washer_max):
        return True
        
    return False

def create_figure_2_replica(r, z, phi, R_max, L_total, hr, hz, 
                          z_washer_start, z_washer_end, r_washer_max):
    """Crear réplica exacta de la Figura 2 del artículo"""
    
    # Convertir a mm para la visualización
    r_mm = r * 1000
    z_mm = z * 1000
    R_mm, Z_mm = np.meshgrid(z_mm, r_mm)  # Nota: invertido para tener z en x, r en y
    
    # Crear figura con espaciado correcto para evitar que se corte
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    fig.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.1, wspace=0.25)
    
    # Panel (a): Potencial eléctrico EXACTO como en la Figura 2a del artículo
    
    # Contornos de potencial con niveles específicos (como en el artículo)
    # Niveles visibles en la figura: -500, -400, -300, -275, -250, -225, -200, -175, -150, -125, -100, -75, -50, -25, 0
    main_levels = [-500, -400, -300, -275, -250, -225, -200, -175, -150, -125, -100, -75, -50, -25, 0]
    
    # Contornos principales (líneas negras como en el artículo)
    cs_main = ax1.contour(R_mm, Z_mm, phi, levels=main_levels, colors='black', linewidths=0.8)
    ax1.clabel(cs_main, inline=True, fontsize=8, fmt='%g')
    
    # Contornos adicionales para más detalle (más líneas finas)
    detailed_levels = np.linspace(-500, 0, 50)  # 50 niveles para suavidad
    cs_detailed = ax1.contour(R_mm, Z_mm, phi, levels=detailed_levels, colors='black', linewidths=0.4, alpha=0.6)
    
    # Fondo de colores (plasma como en el artículo)  
    cs_filled = ax1.contourf(R_mm, Z_mm, phi, levels=detailed_levels, cmap='plasma', alpha=0.8)
    
    # Calcular campo eléctrico E = -∇φ EXACTO
    dz_mm = (z[1] - z[0]) * 1000  # mm
    dr_mm = (r[1] - r[0]) * 1000  # mm
    
    # Gradientes usando diferencias finitas de segundo orden
    dphi_dr = np.zeros_like(phi)
    dphi_dz = np.zeros_like(phi)
    
    # Calcular ∂φ/∂r (gradiente radial) con diferencias centradas
    for i in range(1, phi.shape[0]-1):
        for j in range(phi.shape[1]):
            dphi_dr[i,j] = (phi[i+1,j] - phi[i-1,j]) / (2*dr_mm)
    
    # Calcular ∂φ/∂z (gradiente axial) con diferencias centradas  
    for i in range(phi.shape[0]):
        for j in range(1, phi.shape[1]-1):
            dphi_dz[i,j] = (phi[i,j+1] - phi[i,j-1]) / (2*dz_mm)
    
    # Campo eléctrico: E = -∇φ
    Er = -dphi_dr  # Er = -∂φ/∂r
    Ez = -dphi_dz  # Ez = -∂φ/∂z
    
    # Vectores de campo eléctrico (flechas ROJAS como en el artículo)
    skip_r = 8  # Saltar puntos para claridad visual
    skip_z = 12 # Más espaciado en z
    
    magnitude = np.sqrt(Er**2 + Ez**2)
    
    # Filtrar vectores en regiones de interés (fuera del electrodo)
    mask_field = np.zeros_like(magnitude, dtype=bool)
    for i in range(0, phi.shape[0], skip_r):
        for j in range(0, phi.shape[1], skip_z):
            # Solo mostrar vectores fuera del electrodo y con magnitud significativa
            if not is_electrode_point(i, j, z_washer_start, z_washer_end, r_washer_max):
                if magnitude[i,j] > 0.05 * np.max(magnitude):
                    mask_field[i,j] = True
    
    # Normalizar vectores para visualización consistente
    Er_vis = np.divide(Er, magnitude, out=np.zeros_like(Er), where=(magnitude > 0))
    Ez_vis = np.divide(Ez, magnitude, out=np.zeros_like(Ez), where=(magnitude > 0))
    
    # Dibujar vectores ROJOS como en el artículo original
    for i in range(0, phi.shape[0], skip_r):
        for j in range(0, phi.shape[1], skip_z):
            if mask_field[i,j]:
                # Posición en mm
                z_pos = j * dz_mm / 1000 * 1000  # Convertir a mm para plot
                r_pos = i * dr_mm / 1000 * 1000  # Convertir a mm para plot
                
                # Componentes del vector (escalados para visualización)
                scale_factor = 2.0  # Factor de escala para tamaño de flechas
                ez_scaled = Ez_vis[i,j] * scale_factor
                er_scaled = Er_vis[i,j] * scale_factor
                
                # Dibujar flecha roja
                ax1.arrow(z_pos, r_pos, ez_scaled, er_scaled, 
                         head_width=0.15, head_length=0.1, 
                         fc='red', ec='red', alpha=0.8, width=0.02)
    
    # Configuración del panel (a) EXACTA como el artículo
    ax1.set_xlabel('z (mm)', fontsize=14, fontweight='bold')
    ax1.set_ylabel('x (mm)', fontsize=14, fontweight='bold')  # Nota: artículo usa x, no r
    ax1.set_title('a)', fontsize=16, fontweight='bold', loc='left')
    ax1.set_xlim(0, L_total*1000)        # 0 a 16 mm
    ax1.set_ylim(-R_max*1000, R_max*1000) # -6 a +6 mm (simétrico como en el artículo)
    ax1.grid(False)  # Sin grid para coincidir con el artículo
    
    # Dibujar geometría de electrodos EXACTA como en el artículo
    # Electrodo washer rectangular NEGRO (como en la Figura 2a)
    washer_z_start_mm = (z_washer_start * hz) * 1000
    washer_z_end_mm = (z_washer_end * hz) * 1000  
    washer_r_mm = (r_washer_max * hr) * 1000
    
    # Rectángulo del washer SUPERIOR
    rect_upper = plt.Rectangle((washer_z_start_mm, 0), 
                              washer_z_end_mm - washer_z_start_mm, 
                              washer_r_mm, 
                              facecolor='black', edgecolor='black', linewidth=0, zorder=10)
    ax1.add_patch(rect_upper)
    
    # Rectángulo del washer INFERIOR (simétrico)
    rect_lower = plt.Rectangle((washer_z_start_mm, -washer_r_mm), 
                              washer_z_end_mm - washer_z_start_mm, 
                              washer_r_mm, 
                              facecolor='black', edgecolor='black', linewidth=0, zorder=10)
    ax1.add_patch(rect_lower)
    
    # Electrodos cilíndricos externos (líneas horizontales NEGRAS)
    ax1.axhline(y=R_max*1000, color='black', linewidth=4, alpha=1.0, zorder=10)
    ax1.axhline(y=-R_max*1000, color='black', linewidth=4, alpha=1.0, zorder=10)
    
    # Eje central (línea en r=0)
    ax1.axhline(y=0, color='black', linewidth=1, alpha=0.5)
    
    # Barra de colores con etiquetas exactas como el artículo
    cbar1 = plt.colorbar(cs_filled, ax=ax1, shrink=0.8, pad=0.02)
    cbar1.set_label('V', fontsize=14, rotation=0, labelpad=15)
    
    # Configurar ticks de la barra de colores
    cbar1.set_ticks([-500, -400, -300, -200, -100, 0])
    cbar1.set_ticklabels(['-500', '-400', '-300', '-200', '-100', '0'])
    
    # Panel (b): Perfil de potencial a lo largo del eje z (x=0) EXACTO como el artículo
    
    phi_axis = phi[0, :]  # Potencial en r=0 (eje central)
    
    # Línea azul sólida (potencial anarmónico obtenido) - EXACTO como artículo
    ax2.plot(z_mm, phi_axis, 'b-', linewidth=2.5, label='Anharmonic potential')
    
    # Línea verde punteada (aproximación armónica) - EXACTO como artículo
    # Encontrar el mínimo del potencial
    z_min_idx = np.argmin(phi_axis)
    z_min = z_mm[z_min_idx]
    phi_min = phi_axis[z_min_idx]
    
    # Crear aproximación armónica más realista
    # Usar la curvatura en el fondo del pozo para definir la parábola
    if z_min_idx > 10 and z_min_idx < len(z_mm) - 10:
        # Estimar curvatura usando puntos alrededor del mínimo
        window = 5  # Ventana para calcular curvatura
        z_fit = z_mm[z_min_idx-window:z_min_idx+window+1]
        phi_fit = phi_axis[z_min_idx-window:z_min_idx+window+1]
        
        # Ajuste parabólico: φ = a(z-z0)² + φ0
        from numpy.polynomial import Polynomial
        z_centered = z_fit - z_min
        poly = Polynomial.fit(z_centered, phi_fit, deg=2)
        coeffs = poly.convert().coef
        
        # Crear la aproximación armónica completa
        phi_harmonic = coeffs[0] + coeffs[2] * (z_mm - z_min)**2
        
        # Solo mostrar la parte relevante (como en el artículo)
        mask_harmonic = (z_mm >= 2) & (z_mm <= 14)  # Rango visible en el artículo
        ax2.plot(z_mm[mask_harmonic], phi_harmonic[mask_harmonic], 'g--', linewidth=2, alpha=0.8, label='Harmonic approximation')
    
    # Configuración del panel (b) EXACTA como el artículo
    ax2.set_xlabel('z (mm)', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Electrostatic potential (V)', fontsize=14, fontweight='bold')
    ax2.set_title('b)', fontsize=16, fontweight='bold', loc='left')
    ax2.set_xlim(0, 16)  # 0 a 16 mm como el artículo
    ax2.set_ylim(-400, 0)  # Rango del artículo
    ax2.grid(True, alpha=0.3, linestyle='--')
    
    # Líneas verticales verdes punteadas (como en el artículo original)
    # Estas marcan los límites del washer electrode
    washer_z_start_mm = (z_washer_start * hz) * 1000
    washer_z_end_mm = (z_washer_end * hz) * 1000
    ax2.axvline(x=washer_z_start_mm, color='green', linestyle='--', linewidth=1.5, alpha=0.7)
    ax2.axvline(x=washer_z_end_mm, color='green', linestyle='--', linewidth=1.5, alpha=0.7)
    
    # Marcar el mínimo con punto rojo
    ax2.plot(z_min, phi_min, 'ro', markersize=8, label=f'Minimum: {phi_min:.0f} V', zorder=5)
    
    # Anotación de oscilaciones mq (como en el artículo)
    ax2.text(0.05, 0.95, 'mq\noscillations', transform=ax2.transAxes, 
             fontsize=12, ha='left', va='top', fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.4', facecolor='yellow', alpha=0.8, edgecolor='black'))
    
    # Leyenda
    ax2.legend(loc='lower right', fontsize=11, framealpha=0.9)
    
    # Título general como en el artículo
    fig.suptitle('Figure 2. (a) Electric potential in the y = 0 plane, with the corresponding electrostatic field vectors E (red arrows). ' +
                 '(b) Electric potential profile along the z-axis (x = 0), comparing the obtained...',
                 fontsize=11, y=0.02, wrap=True)
    
    # Ajustes finales para evitar recortes
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.12, top=0.95)
    
    # Guardar y mostrar
    plt.savefig('Figure_2_ART_MS_Replica.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    
    print("✓ Figura 2 del artículo ART-MS replicada exitosamente!")
    return fig

if __name__ == "__main__":
    main()
