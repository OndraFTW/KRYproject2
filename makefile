#Autor: Ondřej Šlampa, xslamp01@stud.fit.vutbr.cz
#Projekt: KRY 2
#Popis: Program pro práci s RSA šifrováním.
#Obsah souboru: makefile

COMP=gcc
FLAGS=-std=c99 -Wall -Wextra -Ofast -lgmp

.PHONY: clean valg

all:kry.c
	$(COMP) -o kry kry.c $(FLAGS)
clean:
	rm -f kry
valg:
	valgrind --tool=memcheck --leak-check=yes --track-origins=yes ./kry -g 8
