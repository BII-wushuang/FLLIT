@echo off
echo Executing FLLIT under Docker
FOR /F "tokens=4 delims= " %%i in ('route print ^| find " 0.0.0.0"') do set ipv4=%%i
docker run -v %cd%:/FLLIT -e DISPLAY=%ipv4%:0 feijianke/fllit