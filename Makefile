CFLAGS= -g -O3
CPPFLAGS=-I$(HTSDIR)
HTSDIR=../htslib
HTSLIB=$(HTSDIR)/libhts.a
LDLIBS=-lpthread $(HTSLIB) -lz -lm -lbz2 -llzma -lcurl

all: pbwt

PBWT_COMMIT_HASH = ""
ifneq "$(wildcard .git)" ""
PBWT_COMMIT_HASH = $(shell git describe --always --long --dirty)
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PBWT_COMMIT_HASH)",$(shell cat version.h)),,force))
endif

version.h:
	echo '#define PBWT_COMMIT_HASH "$(PBWT_COMMIT_HASH)"' > $@

force:

test: all
	./test/test.pl

PBWT_OBJS=pbwtMain.o pbwtCore.o pbwtSample.o pbwtIO.o pbwtMatch.o pbwtImpute.o pbwtPaint.o pbwtLikelihood.o pbwtMerge.o pbwtGeneticMap.o pbwtHtslib.o
UTILS_OBJS=hash.o dict.o array.o utils.o
UTILS_HEADERS=utils.h array.h dict.h hash.h
AUTOZYG_OBJS=autozygExtract.o

pbwt: $(PBWT_OBJS) $(UTILS_OBJS)
	$(LINK.c) $^ $(LDLIBS) -o $@

pbwtMain.o: version.h
$(PBWT_OBJS): pbwt.h $(UTILS_HEADERS)
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

autozygExtract: $(AUTOZYG_OBJS)
	$(LINK.c) $^ $(LDLIBS) -o $@

install: all
	install -d $(PREFIX)
	install pbwt $(PREFIX)

clean:
	$(RM) *.o pbwt *~ version.h

.PHONY: all test clean force install

