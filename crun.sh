#!/bin/bash

find *.c | sed -e 's/\.c//g' | xargs -I % sh -c 'echo "----%----" && gcc -o ./out/%.o %.c && ./out/%.o'