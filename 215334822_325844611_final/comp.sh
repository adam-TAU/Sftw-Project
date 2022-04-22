#!/bin/bash

# assembling and linking
gcc -ansi -Wall -Wextra -Werror -pedantic-errors matrix.c graph.c eigen.c spkmeans.c spkmeans_goals.c -lm -o spkmeans
