rem
rem run all tests under Windows-NT
rem

for %%i in (*.inp) DO cns < %%i > %%~ni.out

rem
rem remove any data files
rem

del *.dat*
del *.DAT*
del *.list*

rem
