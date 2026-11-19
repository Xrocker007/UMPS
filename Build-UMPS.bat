@echo off
setlocal
cd /d "%~dp0"

where git >nul 2>&1
if errorlevel 1 (
    echo Git is required and must be available on PATH.
    goto :failed
)

where cmake >nul 2>&1
if errorlevel 1 (
    echo CMake is required and must be available on PATH.
    goto :failed
)

rem Optional first argument: path to vcpkg.
set "UMPS_VCPKG=%~1"
if not defined UMPS_VCPKG set "UMPS_VCPKG=%VCPKG_ROOT%"
if not defined UMPS_VCPKG set "UMPS_VCPKG=C:\dev\vcpkg"

if not exist "%UMPS_VCPKG%\scripts\buildsystems\vcpkg.cmake" (
    echo Cannot find vcpkg at "%UMPS_VCPKG%".
    echo Run Build-UMPS.bat "path-to-vcpkg"
    goto :failed
)

set "VCPKG_ROOT=%UMPS_VCPKG%"

if not exist "external\VegaFEM\CMakeLists.txt" (
    if not exist "external" mkdir external
    git clone https://github.com/FabienPean/VegaFEM.git external\VegaFEM
    if errorlevel 1 goto :failed
)

"%UMPS_VCPKG%\vcpkg.exe" install opengl:x64-windows freeglut:x64-windows glew:x64-windows tbb:x64-windows arpack-ng:x64-windows
if errorlevel 1 goto :failed

rem Optional second argument: folder containing MKLConfig.cmake.
set "UMPS_MKL=%~2"

if not defined UMPS_MKL (
    for /d %%D in ("%ProgramFiles(x86)%\Intel\oneAPI\mkl\*") do (
        if exist "%%~fD\lib\cmake\mkl\MKLConfig.cmake" (
            set "UMPS_MKL=%%~fD\lib\cmake\mkl"
        )
    )
)

if not exist "%UMPS_MKL%\MKLConfig.cmake" (
    echo Cannot find Intel oneMKL.
    echo Install it, or supply its CMake folder as the second argument:
    echo Build-UMPS.bat "path-to-vcpkg" "path-to-MKL-CMake-folder"
    goto :failed
)

cmake -S . -B build ^
    -G "Visual Studio 17 2022" -A x64 ^
    "-DCMAKE_TOOLCHAIN_FILE=%UMPS_VCPKG%\scripts\buildsystems\vcpkg.cmake" ^
    -DVCPKG_TARGET_TRIPLET=x64-windows ^
    "-DMKL_DIR=%UMPS_MKL%" ^
    "-DVEGA_SOURCE_DIR=%CD%\external\VegaFEM"

if errorlevel 1 goto :failed

cmake --build build --config Release --target UMPS --parallel
if errorlevel 1 goto :failed

echo.
echo UMPS built successfully. Use Run-UMPS.bat to launch it.
pause
exit /b 0

:failed
echo.
echo Build stopped. Review the message above.
pause
exit /b 1