@echo off
setlocal
cd /d "%~dp0"

if not exist "build\Release\UMPS.exe" (
    echo UMPS has not been built. Run Build-UMPS.bat first.
    pause
    exit /b 1
)

"build\Release\UMPS.exe" "%~dp0models\beam3\beam3_tet.veg"

if errorlevel 1 pause