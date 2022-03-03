#!/bin/bash
gcc -ansi -Wall -Wextra -Werror -pedantic-errors test.c matrix.c -lm -o test
if [ "$?" -ne 0 ]; then
    echo "Compilation failed!"
    exit 1
fi
./test