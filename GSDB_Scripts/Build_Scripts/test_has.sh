#!/bin/sh
n=6
output=$(tr -cd '[:alnum:]' < /dev/urandom | head -c$n)
echo $output