
CFLAGS= -g -O3
CPPFLAGS=-I$(HTSDIR)
HTSDIR=../htslib
HTSLIB=$(HTSDIR)/libhts.a
LDLIBS=-lpthread -lz -lm $(HTSLIB)

all: pbwt

test:
	./test/test.pl

PBWT_OBJS=pbwtMain.o pbwtCore.o pbwtSample.o pbwtIO.o pbwtMatch.o pbwtImpute.o pbwtPaint.o pbwtLikelihood.o pbwtMerge.o pbwtGeneticMap.o pbwtHtslib.o
UTILS_OBJS=hash.o dict.o array.o utils.o
AUTOZYG_OBJS=autozygExtract.o

pbwt: $(PBWT_OBJS) $(UTILS_OBJS)
	$(LINK.c) $^ $(LDLIBS) -o $@

autozygExtract: $(AUTOZYG_OBJS)
	$(LINK.c) $^ $(LDLIBS) -o $@

utils.o: hash.o dict.o array.o
	$(COMPILE.c) -o $@ utils.c
	touch utils

install: all
	install -d $(PREFIX)
	install pbwt $(PREFIX)

clean:
	$(RM) *.o pbwt *~

.PHONY: all test clean

