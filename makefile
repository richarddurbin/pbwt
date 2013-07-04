
CFLAGS= -g

all: pbwt

.PHONY:test
test:
	./test/test.pl

pbwt: pbwtMain.o pbwtCore.o pbwtIO.o pbwtMatch.o pbwtImpute.o pbwtMerge.o utils
	gcc $(CFLAGS) -o pbwt pbwtMain.o pbwtCore.o pbwtIO.o pbwtMatch.o pbwtImpute.o pbwtMerge.o hash.o dict.o array.o utils.o -lm

pbwtMain.o: pbwtMain.c pbwt.h utils.h
	gcc $(CFLAGS) -c pbwtMain.c

pbwtCore.o: pbwtCore.c pbwt.h utils.h
	gcc $(CFLAGS) -c pbwtCore.c

pbwtIO.o: pbwtIO.c pbwt.h utils.h
	gcc $(CFLAGS) -c pbwtIO.c

pbwtMatch.o: pbwtMatch.c pbwt.h utils.h
	gcc $(CFLAGS) -c pbwtMatch.c

pbwtImpute.o: pbwtImpute.c pbwt.h utils.h
	gcc $(CFLAGS) -c pbwtImpute.c

pbwtMerge.o: pbwtMerge.c pbwt.h utils.h
	gcc $(CFLAGS) -c pbwtMerge.c

#################################

utils: hash.o dict.o array.o utils.o
	touch utils

utils.h: hash.h dict.h array.h

hash.o: hash.c utils.h
	gcc $(CFLAGS) -c hash.c

dict.o: dict.c utils.h
	gcc $(CFLAGS) -c dict.c

array.o: array.c utils.h
	gcc $(CFLAGS) -c array.c

utils.o: utils.c utils.h
	gcc $(CFLAGS) -c utils.c

clean:
	rm -f *.o pbwt
