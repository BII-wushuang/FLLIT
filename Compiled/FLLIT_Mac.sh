#!/bin/sh
echo Executing FLLIT under Docker
ipv4=$(ifconfig | grep "inet " | grep -Fv 127.0.0.1 | awk '{print $2}')
fllitdir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd)
docker run -v $fllitdir:/FLLIT -e DISPLAY=$ipv4:0 feijianke/fllit
