
CFLAGS= -g -O3
CPPFLAGS=-I$(HTSDIR)
HTSDIR=../htslib
HTSLIB=$(HTSDIR)/libhts.a
LDLIBS=-lpthread -lz -lm $(HTSLIB)

all: pbwt

PBWT_COMMIT_HASH = ""
ifneq "$(wildcard .git)" ""
PBWT_COMMIT_HASH = $(shell git describe --always --long --dirty)
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PBWT_COMMIT_HASH)",$(shell cat version.h)),,force))
endif
version.h:
	echo '#define PBWT_COMMIT_HASH "$(PBWT_COMMIT_HASH)"' > $@


force:

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

.PHONY: all test clean force

