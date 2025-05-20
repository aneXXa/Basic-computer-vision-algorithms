@echo off
:: Creating dynamic library
gcc filters.c -shared -fPIC -O2 -o filters.dll
:: Compiling working program
gcc main.c -o bcva.exe -L. -lfilters -lm