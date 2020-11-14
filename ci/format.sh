#! /bin/bash


for f in $(find src -iname "*.h" -o -iname "*.hpp" -o -iname "*.cpp" -o -iname "*.c"); do 
    expand -t 4 $f > /tmp/e
    mv -f /tmp/e $f
done

for f in $(find src -iname "*.h" -o -iname "*.hpp" -o -iname "*.cpp" -o -iname "*.c"); do 
    if ! grep -q License $f
    then
        cat lic.txt $f > /tmp/e
        mv -f /tmp/e $f
    else
        echo "License exists in $f"
    fi
done
