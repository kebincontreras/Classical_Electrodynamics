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

def solve_art_ms_improved():
    """
    Resuelve la ecuación de Laplace para el ART-MS usando el método mejorado
    Adaptado del código C++ con malla rectangular y condiciones de frontera precisas
    """
    
    # Parámetros del problema (como en el código C++)
    V1, V2, V3, V4 = 0.0, 0.0, 0.0, 0.0  # Condiciones de frontera externas (tierra)
    V5 = -500.0  # Potencial del washer electrode
    Vmax = abs(V5)
    
    # Dimensiones de la zona de cálculo
    Lx = 0.016  # 16 mm (longitud total)
    Ly = 0.012  # 12 mm (ancho total: -6 a +6 mm)
    
    # Puntos de malla (alta resolución como en C++)
    npx = 251  # puntos en x (dirección z en nuestro caso)
    npy = 201  # puntos en y (dirección r en nuestro caso)
    
    # Crear arrays de coordenadas
    x = np.linspace(0, Lx, npx)  # x ∈ [0, 16mm] (dirección z)
    y = np.linspace(-0.006, 0.006, npy)  # y ∈ [-6mm, 6mm] (dirección r)
    
    # Paso espacial
    hx = Lx / (npx - 1)
    hy = Ly / (npy - 1)
    
    print(f"Malla: {npx}x{npy} puntos")
    print(f"Paso espacial: hx = {hx*1000:.3f} mm, hy = {hy*1000:.3f} mm")
    print(f"Dominio: x ∈ [0, {Lx*1000:.1f}] mm, y ∈ [{y[0]*1000:.1f}, {y[-1]*1000:.1f}] mm")
    
    # Normalizar voltajes
    V1_norm = V1 / Vmax
    V2_norm = V2 / Vmax  
    V3_norm = V3 / Vmax
    V4_norm = V4 / Vmax
    V5_norm = V5 / Vmax
    
    # Inicializar matriz de potencial
    v = np.zeros((npx, npy))
    v_old = np.zeros((npx, npy))
    
    # Aplicar condiciones de frontera iniciales
    # Fronteras externas a tierra
    v[0, :] = V1_norm      # x = 0 (z = 0)
    v[-1, :] = V2_norm     # x = Lx (z = 16mm)  
    v[:, 0] = V3_norm      # y = -6mm (r = -6mm)
    v[:, -1] = V4_norm     # y = +6mm (r = +6mm)
    
    # Parámetros de convergencia
    itmax = 1000
    tol = 1e-6
    
    print(f"\nResolviendo ecuación de Laplace...")
    print(f"Tolerancia: {tol:.1e}, Iteraciones máximas: {itmax}")
    
    # Proceso iterativo (Gauss-Seidel mejorado)
    for it in range(1, itmax + 1):
        max_error = 0.0
        
        # Guardar iteración anterior
        v_old[:, :] = v[:, :]
        
        # Iterar sobre puntos internos
        for i in range(1, npx - 1):
            for j in range(1, npy - 1):
                
                # Definir región del washer electrode (como en C++)
                # Zona vertical delgada centrada en x = 8mm con ancho pequeño
                x_washer_center = 0.008  # 8 mm (centro en z)
                x_washer_width = 0.0004  # 0.4 mm de ancho (muy delgado)
                y_washer_inner = 0.0002  # 0.2 mm (gap central pequeño)
                y_washer_outer = 0.006   # 6 mm (radio máximo)
                
                # Condición del washer: zona vertical delgada excluyendo región central
                washer_condition = (
                    (x[i] >= x_washer_center - x_washer_width/2) and 
                    (x[i] <= x_washer_center + x_washer_width/2) and
                    (
                        (y[j] >= -y_washer_outer and y[j] <= -y_washer_inner) or
                        (y[j] >= y_washer_inner and y[j] <= y_washer_outer)
                    )
                )
                
                if washer_condition:
                    # Potencial fijo en el washer
                    v[i, j] = V5_norm
                else:
                    # Ecuación de diferencias finitas mejorada (como en C++)
                    v[i, j] = (
                        hy**2 * (v[i+1, j] + v[i-1, j]) +
                        hx**2 * (v[i, j+1] + v[i, j-1])
                    ) / (2.0 * (hx**2 + hy**2))
        
        # Reaplicar condiciones de frontera después de cada iteración
        v[0, :] = V1_norm
        v[-1, :] = V2_norm
        v[:, 0] = V3_norm
        v[:, -1] = V4_norm
        
        # Calcular error máximo
        for i in range(1, npx - 1):
            for j in range(1, npy - 1):
                er_abs = abs(v[i, j] - v_old[i, j])
                if er_abs > max_error:
                    max_error = er_abs
        
        # Mostrar progreso
        if it % 1000 == 0 or max_error <= tol:
            print(f"Iteración {it}: Error = {max_error:.2e}")
        
        # Verificar convergencia
        if max_error <= tol:
            print(f"¡Convergencia alcanzada en {it} iteraciones!")
            break
    
    if it >= itmax:
        print(f"¡Máximo de iteraciones ({itmax}) alcanzado!")
    
    # Desnormalizar el potencial
    v_real = v * Vmax
    
    return v_real, x, y, npx, npy

def plot_art_ms_improved():
    """
    Crea las gráficas de la simulación mejorada del ART-MS
    """
    
    # Resolver el problema
    v, x, y, npx, npy = solve_art_ms_improved()
    
    # Convertir a mm para visualización
    x_mm = x * 1000  # z en mm
    y_mm = y * 1000  # r en mm
    
    # Crear meshgrids para plotting
    X_mm, Y_mm = np.meshgrid(x_mm, y_mm)
    
    # Crear figura con dos subplots
    fig = plt.figure(figsize=(16, 7))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.3, 1])
    
    # ==================== PANEL (a): Mapa 2D del potencial ====================
    ax1 = fig.add_subplot(gs[0])
    
    # Niveles de contorno de 0 a -500V
    levels = np.linspace(-500, 0, 25)
    
    # Mapa de colores personalizado (azul -> verde -> rojo)
    from matplotlib.colors import LinearSegmentedColormap
    colors = ['#0000FF', '#00FFFF', '#00FF00', '#FFFF00', '#FF0000']
    cmap_custom = LinearSegmentedColormap.from_list('art_ms', colors, N=100)
    
    # Contornos con relleno
    contourf = ax1.contourf(X_mm, Y_mm, v.T, levels=levels, cmap=cmap_custom, extend='both')
    
    # Líneas equipotenciales
    contour_lines = ax1.contour(X_mm, Y_mm, v.T, levels=levels[::3], 
                               colors='black', linewidths=0.8, alpha=0.7)
    ax1.clabel(contour_lines, inline=True, fontsize=8, fmt='%g V')
    
    # Calcular y dibujar campo eléctrico
    Ex = np.zeros_like(v)  # Campo en dirección x (z)
    Ey = np.zeros_like(v)  # Campo en dirección y (r)
    
    # Calcular gradientes (E = -∇φ)
    for i in range(1, npx - 1):
        for j in range(1, npy - 1):
            Ex[i, j] = -(v[i+1, j] - v[i-1, j]) / (2 * (x[1] - x[0]))
            Ey[i, j] = -(v[i, j+1] - v[i, j-1]) / (2 * (y[1] - y[0]))
    
    # Vectores del campo eléctrico
    skip = 10
    X_sub = X_mm[::skip, ::skip]
    Y_sub = Y_mm[::skip, ::skip] 
    Ex_sub = Ex.T[::skip, ::skip]
    Ey_sub = Ey.T[::skip, ::skip]
    
    # Normalizar para mejor visualización
    magnitude = np.sqrt(Ex_sub**2 + Ey_sub**2)
    Ex_norm = np.where(magnitude > 0, Ex_sub / magnitude, 0)
    Ey_norm = np.where(magnitude > 0, Ey_sub / magnitude, 0)
    
    ax1.quiver(X_sub, Y_sub, Ex_norm, Ey_norm, 
              color='red', alpha=0.8, scale=25, width=0.003)
    
    # Dibujar geometría del trap
    # Washer electrode (líneas verticales)
    ax1.axvline(x=7.8, color='black', linewidth=6, alpha=0.8)
    ax1.axvline(x=8.2, color='black', linewidth=6, alpha=0.8)
    
    # Región central (gap)
    rect_gap = Rectangle((7.8, -0.2), 0.4, 0.4, 
                        facecolor='white', edgecolor='black', linewidth=2)
    ax1.add_patch(rect_gap)
    
    # Fronteras externas
    ax1.axhline(y=-6, color='black', linewidth=4)
    ax1.axhline(y=6, color='black', linewidth=4)
    ax1.axvline(x=0, color='black', linewidth=4)
    ax1.axvline(x=16, color='black', linewidth=4)
    
    # Configuración
    ax1.set_xlabel('z (mm)', fontsize=12)
    ax1.set_ylabel('r (mm)', fontsize=12)
    ax1.set_title('(a) Potencial eléctrico en el plano y=0', fontsize=14, fontweight='bold')
    ax1.set_xlim(0, 16)
    ax1.set_ylim(-6, 6)
    ax1.grid(True, alpha=0.3)
    ax1.set_aspect('equal')
    
    # Colorbar
    cbar1 = fig.colorbar(contourf, ax=ax1, shrink=0.8)
    cbar1.set_label('Potencial eléctrico (V)', fontsize=11)
    cbar1.ax.tick_params(labelsize=10)
    
    # ==================== PANEL (b): Perfil axial del potencial ====================
    ax2 = fig.add_subplot(gs[1])
    
    # Potencial a lo largo del eje z (r=0, que corresponde a y=0)
    j_center = npy // 2  # Índice del centro en y
    phi_axis = v[:, j_center]  # Potencial en r=0
    
    # Plot principal
    ax2.plot(x_mm, phi_axis, 'b-', linewidth=3, label='Potencial anarmónico obtenido')
    
    # Aproximación armónica simple
    phi_min = np.min(phi_axis)
    x_min_idx = np.argmin(phi_axis)
    x_min = x_mm[x_min_idx]
    
    # Parábola ajustada
    k_harmonic = phi_min / ((x_mm[-1]/2 - x_min)**2)
    phi_harmonic = phi_min + k_harmonic * (x_mm - x_min)**2
    
    ax2.plot(x_mm, phi_harmonic, 'g--', linewidth=2.5, label='Aproximación armónica')
    
    # Marcar región del washer
    ax2.axvspan(7.8, 8.2, alpha=0.3, color='gray', label='Washer electrode')
    
    # Añadir punto rojo en el mínimo del potencial obtenido
    ax2.scatter(x_min, phi_min, color='red', s=100, label='Oscilaciones (m q)')
    
    # Configuración
    ax2.set_xlabel('z (mm)', fontsize=12)
    ax2.set_ylabel('Potencial eléctrico (V)', fontsize=12) 
    ax2.set_title('(b) Perfil del potencial a lo largo del eje z', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=10)
    ax2.set_xlim(0, 16)
    ax2.set_ylim(-500, 50)
    
    # Ajustar layout
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.25)
    
    # Guardar
    filename = 'ART_MS_Trap_mejorado.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight', 
               facecolor='white', edgecolor='none')
    print(f"\n✅ Figura guardada como '{filename}'")
    
    plt.show()
    
    return fig, v, x, y

# ==================== FUNCIÓN PRINCIPAL ====================
def main():
    """Función principal que ejecuta la simulación completa"""
    
    print("="*70)
    print("SIMULACIÓN ART-MS TRAP - VERSIÓN MEJORADA")
    print("Basada en código C++ con malla rectangular y condiciones precisas")
    print("="*70)
    
    # Ejecutar simulación
    fig, v, x, y = plot_art_ms_improved()
    
    # Estadísticas finales
    print("\n" + "="*70)
    print("✅ SIMULACIÓN COMPLETADA EXITOSAMENTE")
    print(f"✅ Potencial mínimo: {np.min(v):.1f} V")
    print(f"✅ Potencial máximo: {np.max(v):.1f} V") 
    print(f"✅ Rango de cálculo: z ∈ [0, 16] mm, r ∈ [-6, 6] mm")
    print(f"✅ Resolución de malla: {len(x)} × {len(y)} puntos")
    print("="*70)
    
    return fig, v, x, y

# ==================== EJECUCIÓN ====================
if __name__ == "__main__":
    main()
