#! /bin/bash


for f in $(find src -iname "*.h" -o -iname "*.hpp"); do 
    echo $f
done
