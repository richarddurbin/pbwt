
pbwt: pbwtMain.o pbwtCore.o pbwtIO.o pbwtMatch.o pbwtImpute.o utils
	gcc -g -o pbwt pbwtMain.o pbwtCore.o pbwtIO.o pbwtMatch.o pbwtImpute.o hash.o dict.o array.o utils.o -lm

pbwtMain.o: pbwtMain.c pbwt.h utils.h
	gcc -g -c pbwtMain.c

pbwtCore.o: pbwtCore.c pbwt.h utils.h
	gcc -g -c pbwtCore.c

pbwtIO.o: pbwtIO.c pbwt.h utils.h
	gcc -g -c pbwtIO.c

pbwtMatch.o: pbwtMatch.c pbwt.h utils.h
	gcc -g -c pbwtMatch.c

pbwtImpute.o: pbwtImpute.c pbwt.h utils.h
	gcc -g -c pbwtImpute.c

#################################

utils: hash.o dict.o array.o utils.o
	touch utils

utils.h: hash.h dict.h array.h

hash.o: hash.c utils.h
	gcc -g -c hash.c

dict.o: dict.c utils.h
	gcc -g -c dict.c

array.o: array.c utils.h
	gcc -g -c array.c

utils.o: utils.c utils.h
	gcc -g -c utils.c
