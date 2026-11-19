@echo off
setlocal

set "UMPS_VCPKG=C:\dev\vcpkg"

where git >nul 2>&1
if errorlevel 1 (
    echo Git is missing. Install Git for Windows, then run this again.
    goto :failed
)

if not exist "%UMPS_VCPKG%\." (
    git clone https://github.com/microsoft/vcpkg.git "%UMPS_VCPKG%"
    if errorlevel 1 goto :failed
)

if not exist "%UMPS_VCPKG%\bootstrap-vcpkg.bat" (
    echo The folder exists but does not contain a valid vcpkg installation.
    goto :failed
)

if not exist "%UMPS_VCPKG%\vcpkg.exe" (
    call "%UMPS_VCPKG%\bootstrap-vcpkg.bat" -disableMetrics
    if errorlevel 1 goto :failed
)

set "VCPKG_ROOT=%UMPS_VCPKG%"

"%UMPS_VCPKG%\vcpkg.exe" install ^
    opengl:x64-windows ^
    freeglut:x64-windows ^
    glew:x64-windows ^
    tbb:x64-windows ^
    arpack-ng:x64-windows

if errorlevel 1 goto :failed

echo.
echo SUCCESS: vcpkg and UMPS dependencies installed.
echo Location: %UMPS_VCPKG%
pause
exit /b 0

:failed
echo.
echo Installation stopped. Check the error above.
pause
exit /b 1