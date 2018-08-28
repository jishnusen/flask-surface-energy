#!/bin/sh

gcc -O3 -c -fPIC numerical_integral_pkg/*.c 
gcc *.o -shared -o numerical_integral.so

rm *.o
