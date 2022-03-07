#!/bin/bash

# assembling and linking
gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c matrix.c -lm -o matrix.o
gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c graph.c -lm -o graph.o
gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c eigen.c -lm -o eigen.o
gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c spkmeans.c -lm -o spkmeans.o
gcc -ansi -Wall -Wextra -Werror -pedantic-errors matrix.o graph.o eigen.o spkmeans.o -lm -o spkmeans


# clean up
rm *.o
