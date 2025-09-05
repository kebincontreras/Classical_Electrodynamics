@echo off
echo ===============================================
echo     Taller de Electrodinamica Clasica
echo     Ejecutando simulaciones automaticamente
echo ===============================================
echo.

REM Crear entorno virtual si no existe
if not exist "venv" (
    echo Creando entorno virtual...
    python -m venv venv
    echo Entorno virtual creado.
    echo.
)

REM Activar entorno virtual
echo Activando entorno virtual...
call venv\Scripts\activate
echo.

REM Instalar dependencias
echo Instalando dependencias...
pip install -r requirements.txt
echo Dependencias instaladas.
echo.

REM Ejecutar simulacion de potencial electrico
echo ===============================================
echo 1. Ejecutando simulacion de potencial electrico...
echo ===============================================
python potencial.py
echo.
pause

REM Ejecutar simulacion de campo electrico
echo ===============================================
echo 2. Ejecutando simulacion de campo electrico...
echo ===============================================
python campo.py
echo.
pause

REM Ejecutar simulacion ART-MS Trap
echo ===============================================
echo 3. Ejecutando simulacion ART-MS Trap...
echo ===============================================
python ART_MS_Trap.py
echo.

echo ===============================================
echo    Todas las simulaciones han sido ejecutadas
echo ===============================================
echo.
pause

REM Desactivar entorno virtual
deactivate
