#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ART-MS Trap - Versión mejorada adaptada del código C++
Simulación del potencial eléctrico en trampa electrostática ART-MS
Basado en el código C++ proporcionado con mejores condiciones de frontera
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.gridspec as gridspec

def solve_art_ms_cylindrical(quick_test=False):
    """
    Resuelve la ecuación de Laplace en coordenadas cilíndricas para el ART-MS
    Ecuación: ∂²φ/∂r² + (1/r)∂φ/∂r + ∂²φ/∂z² = 0
    Dominio: r ∈ [0, R_max], z ∈ [0, L_total]
    
    Parameters:
    -----------
    quick_test : bool
        Si True, usa parámetros reducidos para pruebas rápidas
        Si False, usa parámetros de alta resolución
    """
    
    # Parámetros del problema
    V_ground = 0.0    # Potencial de tierra (electrodos cilíndricos y extremos)
    V_washer = -500.0 # Potencial NEGATIVO del washer electrode
    
    # ¿Por qué -500V y no +500V?
    # - Los iones positivos son atraídos hacia potenciales NEGATIVOS
    # - El washer crea un "pozo de potencial" (mínimo energético)
    # - -500V es típico en trampas iónicas ART-MS experimentales
    
    # Dimensiones del dominio cilíndrico
    R_max = 0.006   # 6 mm (radio máximo)
    L_total = 0.016 # 16 mm (longitud total en z)
    
    # Parámetros ajustables según el modo
    if quick_test:
        print("🏃‍♂️ MODO RÁPIDO: Resolución reducida para pruebas")
        nr = 60    # Puntos radiales (reducido de 120)
        nz = 120   # Puntos axiales (reducido de 240)  
        max_iter = 500   # Iteraciones SOR (reducido de 2000)
        tolerance = 1e-4  # Tolerancia relajada
        omega = 1.8      # Factor SOR optimizado para convergencia rápida
    else:
        print("🎯 MODO ALTA RESOLUCIÓN: Cálculo completo")
        nr = 120   # Puntos radiales 
        nz = 240   # Puntos axiales
        max_iter = 2000  # Iteraciones SOR
        tolerance = 1e-6  # Tolerancia estricta
        omega = 1.9      # Factor SOR para precisión
    
    # Crear arrays de coordenadas cilíndricas
    r = np.linspace(0, R_max, nr)  # r ∈ [0, 6mm]
    z = np.linspace(0, L_total, nz)  # z ∈ [0, 16mm]
    
    # Paso espacial
    hr = R_max / (nr - 1)
    hz = L_total / (nz - 1)
    
    print(f"COORDENADAS CILÍNDRICAS (r, z)")
    print(f"Malla: {nr}×{nz} puntos")
    print(f"Paso espacial: hr = {hr*1000:.3f} mm, hz = {hz*1000:.3f} mm")
    print(f"Dominio: r ∈ [0, {R_max*1000:.1f}] mm, z ∈ [0, {L_total*1000:.1f}] mm")
    
    # Inicializar matriz de potencial
    phi = np.zeros((nr, nz))
    phi_old = np.zeros((nr, nz))
    
    # Definir región del washer electrode en coordenadas cilíndricas
    # Washer centrado en z = 8mm (como en el artículo)
    z_center = 0.008  # 8 mm (centro del washer)
    washer_width = 0.0004  # 0.4 mm de ancho (más delgado)
    washer_r_inner = 0.0008  # 0.8 mm (gap interno para iones)
    washer_r_outer = 0.0055   # 5.5 mm (casi todo el radio, como en la figura)
    
    # Convertir a índices de malla
    z_washer_start_idx = int((z_center - washer_width/2) / hz)
    z_washer_end_idx = int((z_center + washer_width/2) / hz)
    r_washer_inner_idx = int(washer_r_inner / hr)
    r_washer_outer_idx = min(int(washer_r_outer / hr), nr-2)  # No tocar la frontera
    
    print(f"Washer electrode: z ∈ [{(z_center-washer_width/2)*1000:.1f}, {(z_center+washer_width/2)*1000:.1f}] mm")
    print(f"Washer electrode: r ∈ [{washer_r_inner*1000:.1f}, {washer_r_outer*1000:.1f}] mm")
    print(f"Índices washer: z[{z_washer_start_idx}:{z_washer_end_idx}], r[{r_washer_inner_idx}:{r_washer_outer_idx}]")
    
    # Aplicar condiciones de frontera cilíndricas
    # Extremos axiales (z = 0 y z = L_total)
    phi[:, 0] = V_ground      # z = 0
    phi[:, -1] = V_ground     # z = L_total
    
    # Pared cilíndrica externa (r = R_max)
    phi[-1, :] = V_ground
    
    # Washer electrode: SOLO el anillo (NO el centro)
    # El centro (r < r_inner) queda LIBRE para crear el pozo de potencial
    for j in range(z_washer_start_idx, z_washer_end_idx + 1):
        for i in range(r_washer_inner_idx, r_washer_outer_idx + 1):
            phi[i, j] = V_washer
    
    # IMPORTANTE: El centro (r < r_inner) NO se fija, permitiendo el pozo
    print(f"Washer aplicado: z[{z_washer_start_idx}:{z_washer_end_idx}], r[{r_washer_inner_idx}:{r_washer_outer_idx}]")
    print(f"Centro libre: r[0:{r_washer_inner_idx}] para crear pozo de potencial")
    
    print(f"\nResolviendo ecuación de Laplace cilíndrica...")
    print(f"Ecuación: ∂²φ/∂r² + (1/r)∂φ/∂r + ∂²φ/∂z² = 0")
    print(f"Tolerancia: {tolerance:.1e}, Iteraciones máximas: {max_iter}")
    print(f"Factor SOR: {omega}")
    
    # Proceso iterativo usando SOR
    for it in range(1, max_iter + 1):
        max_error = 0.0
        phi_old[:, :] = phi[:, :]
        
        # Iterar sobre puntos internos
        for i in range(1, nr - 1):
            for j in range(1, nz - 1):
                
                # Verificar si estamos en el washer electrode
                washer_condition = (
                    (z_washer_start_idx <= j <= z_washer_end_idx) and
                    (r_washer_inner_idx <= i <= r_washer_outer_idx)
                )
                
                if washer_condition:
                    continue  # Potencial fijo en el washer
                
                # Ecuación de diferencias finitas en coordenadas cilíndricas
                # ∂²φ/∂r² + (1/r)∂φ/∂r + ∂²φ/∂z² = 0
                
                r_i = r[i]  # Radio en el punto i
                
                if r_i > 1e-10:  # Para r > 0
                    # Coeficientes de la ecuación discretizada
                    # ∂²φ/∂r² ≈ (φ[i+1] - 2φ[i] + φ[i-1])/hr²
                    # (1/r)∂φ/∂r ≈ (φ[i+1] - φ[i-1])/(2*r*hr)
                    # ∂²φ/∂z² ≈ (φ[j+1] - 2φ[j] + φ[j-1])/hz²
                    
                    A = 1.0 + hr/(2.0*r_i)  # Coeficiente de φ[i+1,j]
                    B = 1.0 - hr/(2.0*r_i)  # Coeficiente de φ[i-1,j]
                    C = (hr/hz)**2           # Coeficiente de φ[i,j±1]
                    
                    phi_new = (A*phi[i+1,j] + B*phi[i-1,j] + C*(phi[i,j+1] + phi[i,j-1])) / (2.0 + 2.0*C)
                    
                else:  # Para r = 0 (eje de simetría)
                    # Aplicar regla de L'Hôpital: lim(r→0) (1/r)∂φ/∂r = ∂²φ/∂r²
                    # Ecuación se convierte en: 2∂²φ/∂r² + ∂²φ/∂z² = 0
                    C = (hr/hz)**2
                    phi_new = (2.0*phi[i+1,j] + C*(phi[i,j+1] + phi[i,j-1])) / (2.0 + 2.0*C)
                
                # Aplicar SOR
                phi[i,j] = (1.0 - omega) * phi[i,j] + omega * phi_new
        
        # Reaplicar condiciones de frontera después de cada iteración
        phi[:, 0] = V_ground
        phi[:, -1] = V_ground
        phi[-1, :] = V_ground
        
        # Reaplicar washer electrode
        for j in range(z_washer_start_idx, z_washer_end_idx + 1):
            for i in range(r_washer_inner_idx, r_washer_outer_idx + 1):
                phi[i, j] = V_washer
        
        # Calcular error máximo
        for i in range(1, nr - 1):
            for j in range(1, nz - 1):
                er_abs = abs(phi[i, j] - phi_old[i, j])
                if er_abs > max_error:
                    max_error = er_abs
        
        # Mostrar progreso
        if it % 2000 == 0 or max_error <= tolerance:
            print(f"Iteración {it}: Error = {max_error:.2e}")
        
        # Verificar convergencia
        if max_error <= tolerance:
            print(f"¡Convergencia alcanzada en {it} iteraciones!")
            break
    
    if it >= max_iter:
        print(f"¡Máximo de iteraciones ({max_iter}) alcanzado!")
    
    return phi, r, z, nr, nz, z_washer_start_idx, z_washer_end_idx, r_washer_inner_idx, r_washer_outer_idx

def plot_art_ms_cylindrical(quick_test=False):
    """
    Crea las gráficas de la simulación del ART-MS en coordenadas cilíndricas
    
    Parameters:
    -----------
    quick_test : bool
        Si True, usa parámetros reducidos para visualización rápida
        Si False, usa parámetros de alta resolución
    """
    
    # Resolver el problema en coordenadas cilíndricas
    phi, r, z, nr, nz, z_washer_start_idx, z_washer_end_idx, r_washer_inner_idx, r_washer_outer_idx = solve_art_ms_cylindrical(quick_test=quick_test)
    
    # Convertir a mm para visualización
    r_mm = r * 1000  # r en mm
    z_mm = z * 1000  # z en mm
    
    # Para mostrar simetría completa, crear dominio r ∈ [-6, 6] mm
    r_full_mm = np.concatenate([-r_mm[::-1][:-1], r_mm])
    nr_full = len(r_full_mm)
    
    # Extender matriz de potencial con simetría
    phi_full = np.zeros((nr_full, nz))
    phi_full[:nr-1, :] = phi[::-1][:-1, :]  # Parte negativa (espejada)
    phi_full[nr-1:, :] = phi                # Parte positiva
    
    # Crear meshgrids para plotting
    Z_mm, R_full_mm = np.meshgrid(z_mm, r_full_mm)
    
    # Crear figura con dos subplots
    fig = plt.figure(figsize=(16, 7))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.3, 1])
    
    # ==================== PANEL (a): Mapa 2D del potencial ====================
    ax1 = fig.add_subplot(gs[0])
    
    # Niveles de contorno de -500V a 0V
    levels = np.linspace(-500, 0, 25)
    
    # Mapa de colores INVERTIDO: Rojo (bajo/negativo) -> Azul (alto/positivo)
    # Como lo pediste: Rojo = -500V, Azul = 0V
    from matplotlib.colors import LinearSegmentedColormap
    colors = ['#FF0000', '#FF8000', '#FFFF00', '#00FF00', '#00FFFF', '#0000FF', '#000080']
    cmap_custom = LinearSegmentedColormap.from_list('art_ms', colors, N=100)
    
    # Contornos con relleno
    contourf = ax1.contourf(Z_mm, R_full_mm, phi_full, levels=levels, cmap=cmap_custom, extend='both')
    
    # Líneas equipotenciales
    contour_lines = ax1.contour(Z_mm, R_full_mm, phi_full, levels=levels[::3], 
                               colors='black', linewidths=0.8, alpha=0.7)
    ax1.clabel(contour_lines, inline=True, fontsize=8, fmt='%g V')
    
    # Calcular y dibujar campo eléctrico en coordenadas cilíndricas
    Er = np.zeros_like(phi)  # Campo radial
    Ez = np.zeros_like(phi)  # Campo axial
    
    # Calcular gradientes (E = -∇φ)
    for i in range(1, nr - 1):
        for j in range(1, nz - 1):
            Er[i, j] = -(phi[i+1, j] - phi[i-1, j]) / (2 * (r[1] - r[0]))
            Ez[i, j] = -(phi[i, j+1] - phi[i, j-1]) / (2 * (z[1] - z[0]))
    
    # Extender campos con simetría
    Er_full = np.zeros((nr_full, nz))
    Ez_full = np.zeros((nr_full, nz))
    Er_full[:nr-1, :] = -Er[::-1][:-1, :]  # Campo radial cambia signo
    Er_full[nr-1:, :] = Er
    Ez_full[:nr-1, :] = Ez[::-1][:-1, :]   # Campo axial mantiene signo
    Ez_full[nr-1:, :] = Ez
    
    # Vectores del campo eléctrico
    skip = 10
    Z_sub = Z_mm[::skip, ::skip]
    R_sub = R_full_mm[::skip, ::skip] 
    Ez_sub = Ez_full[::skip, ::skip]
    Er_sub = Er_full[::skip, ::skip]
    
    # Normalizar para mejor visualización
    magnitude = np.sqrt(Ez_sub**2 + Er_sub**2)
    Ez_norm = np.where(magnitude > 0, Ez_sub / magnitude, 0)
    Er_norm = np.where(magnitude > 0, Er_sub / magnitude, 0)
    
    ax1.quiver(Z_sub, R_sub, Ez_norm, Er_norm, 
              color='red', alpha=0.8, scale=25, width=0.003)
    
    # Dibujar geometría del trap
    z_washer_start_mm = z_mm[z_washer_start_idx]
    z_washer_end_mm = z_mm[z_washer_end_idx]
    r_washer_inner_mm = r_mm[r_washer_inner_idx]
    r_washer_outer_mm = r_mm[r_washer_outer_idx]
    
    # Washer electrode (rectángulos superior e inferior)
    washer_rect_upper = Rectangle((z_washer_start_mm, r_washer_inner_mm), 
                                 z_washer_end_mm - z_washer_start_mm, 
                                 r_washer_outer_mm - r_washer_inner_mm,
                                 facecolor='darkgray', edgecolor='black', linewidth=2, alpha=0.9)
    ax1.add_patch(washer_rect_upper)
    
    washer_rect_lower = Rectangle((z_washer_start_mm, -r_washer_outer_mm), 
                                 z_washer_end_mm - z_washer_start_mm, 
                                 r_washer_outer_mm - r_washer_inner_mm,
                                 facecolor='darkgray', edgecolor='black', linewidth=2, alpha=0.9)
    ax1.add_patch(washer_rect_lower)
    
    # Fronteras del dominio
    ax1.axhline(y=r_full_mm[-1], color='black', linewidth=4)   # r = +R_max
    ax1.axhline(y=r_full_mm[0], color='black', linewidth=4)    # r = -R_max
    ax1.axvline(x=0, color='black', linewidth=4)               # z = 0
    ax1.axvline(x=z_mm[-1], color='black', linewidth=4)        # z = L_total
    ax1.axhline(y=0, color='gray', linewidth=1, alpha=0.5)     # Eje de simetría
    
    # Configuración
    ax1.set_xlabel('z (mm)', fontsize=12)
    ax1.set_ylabel('r (mm)', fontsize=12)
    ax1.set_title('(a) Potencial eléctrico - Coordenadas Cilíndricas', fontsize=14, fontweight='bold')
    ax1.set_xlim(0, z_mm[-1])
    ax1.set_ylim(r_full_mm[0], r_full_mm[-1])
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal')
    
    # Colorbar
    cbar1 = fig.colorbar(contourf, ax=ax1, shrink=0.8)
    cbar1.set_label('Potencial eléctrico (V)', fontsize=11)
    cbar1.ax.tick_params(labelsize=10)
    
    # ==================== PANEL (b): Perfil axial del potencial ====================
    ax2 = fig.add_subplot(gs[1])
    
    # Potencial a lo largo del eje z (r=0)
    phi_axis = phi[0, :]  # Potencial en r=0 (eje de simetría)
    
    # Plot principal
    ax2.plot(z_mm, phi_axis, 'b-', linewidth=3, label='Potencial anarmónico obtenido')
    
    # Aproximación armónica simple
    phi_min = np.min(phi_axis)
    z_min_idx = np.argmin(phi_axis)
    z_min = z_mm[z_min_idx]
    
    # Parábola ajustada
    k_harmonic = phi_min / ((z_mm[-1]/2 - z_min)**2)
    phi_harmonic = phi_min + k_harmonic * (z_mm - z_min)**2
    
    ax2.plot(z_mm, phi_harmonic, 'g--', linewidth=2.5, label='Aproximación armónica')
    
    # Marcar región de oscilaciones con punto rojo y flechas
    # Encontrar una posición representativa para las oscilaciones (cerca del mínimo)
    z_osc = z_min + 1.0  # 1 mm del mínimo
    z_osc_idx = np.argmin(np.abs(z_mm - z_osc))
    phi_osc = phi_axis[z_osc_idx]
    
    # Punto rojo para marcar las oscilaciones
    ax2.plot(z_osc, phi_osc, 'ro', markersize=8, markerfacecolor='red', 
             markeredgecolor='darkred', markeredgewidth=1.5)
    
    # Flechas indicando oscilaciones
    arrow_length = 0.8
    ax2.annotate('', xy=(z_osc - arrow_length, phi_osc), xytext=(z_osc, phi_osc),
                arrowprops=dict(arrowstyle='<-', color='red', lw=2))
    ax2.annotate('', xy=(z_osc + arrow_length, phi_osc), xytext=(z_osc, phi_osc),
                arrowprops=dict(arrowstyle='->', color='red', lw=2))
    
    # Texto "m.g oscillations"
    ax2.text(z_osc, phi_osc + 30, 'm.g\noscillations', 
             ha='center', va='bottom', fontsize=10, color='red',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                      edgecolor='red', alpha=0.8))
    
    # Marcar región del washer
    ax2.axvspan(z_washer_start_mm, z_washer_end_mm, alpha=0.3, color='gray', label='Washer electrode')
    
    # Configuración
    ax2.set_xlabel('z (mm)', fontsize=12)
    ax2.set_ylabel('Potencial eléctrico (V)', fontsize=12) 
    ax2.set_title('(b) Perfil del potencial a lo largo del eje z (r=0)', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    ax2.set_xlim(0, z_mm[-1])
    ax2.set_ylim(-500, 50)
    
    # Ajustar layout
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.25)
    
    # Guardar
    filename = 'Figure_2_ART_MS_Replica.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight', 
               facecolor='white', edgecolor='none')
    print(f"\n✅ Figura guardada como '{filename}'")
    
    plt.show()
    
    return fig, phi, r, z

# ==================== FUNCIÓN PRINCIPAL ====================
def main(quick_test=False):
    """
    Función principal que ejecuta la simulación completa
    
    Parameters:
    -----------
    quick_test : bool
        Si True, ejecuta versión rápida para pruebas
        Si False, ejecuta versión de alta resolución
    """
    
    print("="*70)
    print("SIMULACIÓN ART-MS TRAP - COORDENADAS CILÍNDRICAS")
    print("Basada en código C++ con malla rectangular y condiciones precisas")
    if quick_test:
        print("🏃‍♂️ MODO RÁPIDO ACTIVADO - Para pruebas y ajustes")
    else:
        print("🎯 MODO ALTA RESOLUCIÓN - Cálculo completo")
    print("="*70)
    
    # Ejecutar simulación
    fig, phi, r, z = plot_art_ms_cylindrical(quick_test=quick_test)
    
    # Estadísticas finales
    print("\n" + "="*70)
    print("✅ SIMULACIÓN COMPLETADA EXITOSAMENTE")
    print(f"✅ Potencial mínimo: {np.min(phi):.1f} V")
    print(f"✅ Potencial máximo: {np.max(phi):.1f} V") 
    print(f"✅ Rango de cálculo: z ∈ [0, 16] mm, r ∈ [0, 6] mm")
    print(f"✅ Resolución de malla: {len(z)} × {len(r)} puntos")
    print("="*70)
    
    return fig, phi, r, z

def test_quick():
    """Función para pruebas rápidas durante desarrollo"""
    print("🏃‍♂️ EJECUTANDO PRUEBA RÁPIDA...")
    return main(quick_test=True)

def run_full():
    """Función para simulación completa de alta resolución"""
    print("🎯 EJECUTANDO SIMULACIÓN COMPLETA...")
    return main(quick_test=False)

# ==================== EJECUCIÓN ====================
if __name__ == "__main__":
    import sys
    
    # Permitir seleccionar modo desde línea de comandos
    if len(sys.argv) > 1 and sys.argv[1] == "--quick":
        test_quick()
    elif len(sys.argv) > 1 and sys.argv[1] == "--full":
        run_full()
    else:
        # Por defecto, modo rápido para desarrollo
        print("💡 TIPS:")
        print("  - Usa 'python ART_MS_Trap.py --quick' para pruebas rápidas")
        print("  - Usa 'python ART_MS_Trap.py --full' para alta resolución")
        print("  - Sin argumentos = modo rápido por defecto\n")
        test_quick()
