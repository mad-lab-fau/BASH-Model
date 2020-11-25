@echo off

echo.Enter or Drag'n'Drop the input files here:
echo.

echo..osim file:
set /P osim=

echo.scale file:
set /P scale=

echo..mot file:
set /P mot=

echo..sto file:
set /P sto=

echo.limit frames to:
set /P frames=

echo.Baseline-Model directory:
set /P model=

IF NOT "%osim%" == "" (
	set osim=--osim %osim%
)
IF NOT "%scale%" == "" (
	set scale=--scale %scale%
)
IF NOT "%mot%" == "" (
	set mot=--mot %mot%
)
IF NOT "%sto%" == "" (
	set sto=--sto %sto%
)
IF NOT "%frames%" == "" (
	set frames=--frames %frames%
)
IF NOT "%model%" == "" (
	set model=--model %model%
)

echo.
echo.Starting Application: .\build\x64\Release\SCAPE.exe %osim% %scale% %mot% %sto% %frames% %model%
.\build\x64\Release\SCAPE.exe %osim% %scale% %mot% %sto% %frames% %model%

@pause