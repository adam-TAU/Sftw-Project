#!/bin/bash

# assembling and linking
gcc -g -ansi -Wall -Wextra -Werror -pedantic-errors -c matrix.c -lm -o matrix.o
gcc -g -ansi -Wall -Wextra -Werror -pedantic-errors -c graph.c -lm -o graph.o
gcc -g -ansi -Wall -Wextra -Werror -pedantic-errors -c eigen.c -lm -o eigen.o
gcc -g  -ansi -Wall -Wextra -Werror -pedantic-errors -c spkmeans.c -lm -o spkmeans.o
gcc -g -ansi -Wall -Wextra -Werror -pedantic-errors -c spkmeans_goals.c  -lm -o spkmeans_goals.o
gcc -g -ansi -Wall -Wextra -Werror -pedantic-errors matrix.o graph.o eigen.o spkmeans.o spkmeans_goals.o -lm -o spkmeans


# clean up
rm *.o
