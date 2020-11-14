#! /bin/bash


for f in $(find src -iname "*.h" -o -iname "*.hpp" -o -iname "*.cpp" -o -iname "*.c"); do 
    expand -t 4 $f > /tmp/e
    mv /tmp/e $f
done
