@echo off
rem
rem This file compiles source using g77/gcc on Windows-NT
rem
rem   written by: Paul Adams and Ralf Grosse-Kunstleve
rem
rem   copyright Yale University
rem
rem ==========================================================================

echo compiling Fortran source
echo compiler: g77 
echo flags: -c -fno-globals -Wno-globals -O3 -malign-double -funroll-loops -ffast-math

for %%i in (*.f) do echo compiling file: %%i & g77 -c -fno-globals -Wno-globals -O3 -malign-double -funroll-loops -ffast-math %%i

echo compiling C source
echo compiler: gcc 
echo flags: -c -O -ffast-math

for %%i in (*.c) do echo compiling file: %%i & gcc -c -O -ffast-math %%i

echo linking cns_solve.exe
echo linker: g77

g77 -o cns_solve.exe *.o
