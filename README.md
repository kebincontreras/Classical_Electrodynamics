# Taller de Electrodinámica Clásica

Este repositorio contiene las simulaciones y visualizaciones para el taller de electrodinámica clásica, implementando soluciones numéricas para problemas de campos eléctricos y potenciales.

## Contenido

### 1. Potencial Eléctrico y Equipotenciales

![Potencial Eléctrico](potencial.png)

Para ejecutar la simulación del potencial eléctrico:
```bash
python potencial.py
```


### 2. Campo Eléctrico

![Campo Eléctrico Hemisferio](campo_electrico_hemisferio.png)

Para ejecutar la simulación del campo eléctrico:
```bash
python campo.py
```

### 3. Simulación ART-MS Trap

![Réplica Figure 2 ART-MS](Figure_2_ART_MS_Replica.png)

Para ejecutar la simulación del ART-MS Trap:
```bash
python ART_MS_Trap.py
```

### 4. Simulaciones Adicionales

![Campo Eléctrico 2D](campo.gif)

![Campo Eléctrico 3D](campo_3d.gif)

Para ejecutar simulaciones adicionales del campo eléctrico:
```bash
python Campo_electrico.py      # Para la visualización 2D
python Campo_electrico_3D.py   # Para la visualización 3D
```

## Ejecución Automática

Para ejecutar todas las simulaciones de forma automática, utiliza el archivo batch:
```bash
run_simulations.bat
```

Este archivo creará un entorno virtual, instalará las dependencias necesarias y ejecutará todas las simulaciones en secuencia.

## Requisitos

- Python 3.7+
- NumPy
- Matplotlib
- SciPy

## Instalación Manual

```bash
python -m venv venv
venv\Scripts\activate
pip install numpy matplotlib scipy
```
