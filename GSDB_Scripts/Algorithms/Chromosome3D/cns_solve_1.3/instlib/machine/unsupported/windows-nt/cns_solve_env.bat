@echo off
rem
rem This file sets up the appropriate environmental variables
rem for CNSsolve under Windows-NT
rem
rem   written by: Paul Adams and Ralf Grosse-Kunstleve
rem
rem   copyright Yale University
rem
rem ==========================================================================
rem
rem ****** Important: define the location of the CNSsolve directory ******
rem
rem CHANGE THE NEXT LINE TO POINT TO THE LOCATION OF THE CNSsolve DIRECTORY

	set CNS_SOLVE=U:\cns_solve

rem
rem ==========================================================================
rem

rem ****** Windows-98 ONLY ******
rem
rem if you want to enable convenient command-line editing,
rem remove the "rem" at the beginning of the next line.
rem C:\WINDOWS\COMMAND\DOSKEY /INSERT

rem ==========================================================================
rem
rem general environmental variables
rem

set CNS_LIB=%CNS_SOLVE%\libraries
set CNS_MODULE=%CNS_SOLVE%\modules
set CNS_TOPPAR=%CNS_LIB%\toppar
set CNS_CONFDB=%CNS_LIB%\confdb
set CNS_XTALLIB=%CNS_LIB%\xtal
set CNS_NMRLIB=%CNS_LIB%\nmr
set CNS_XRAYLIB=%CNS_LIB%\xray
set CNS_XTALMODULE=%CNS_MODULE%\xtal
set CNS_NMRMODULE=%CNS_MODULE%\nmr
set CNS_HELPLIB=%CNS_SOLVE%\helplib

rem
rem ==========================================================================
rem
rem the full path to the CNSsolve executable directory
rem   if there is an error about too many parameters for the set
rem   statement you should quote the PATH variable: "%PATH%"

set path=%PATH%;%CNS_SOLVE%\windows-nt

rem
rem ==========================================================================
rem
