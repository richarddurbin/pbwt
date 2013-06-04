#include "pbwt.h"

#include <string.h>
#include <errno.h>

// Synced reading of multiple PBWTs. 
// Note: all operations in-memory for now, buffered input in future?
typedef struct
{
	// all operations in-memory for now, buffered input in future?
	int n;          // number of PBWTs
	PBWT **pbwt;
    int *cpos;      // current position of each PBWT (site index)
    int mpos;       // minimum position across all PBWT (position)
    Update **update;
    int *unpacked;
}
pbwt_reader_t;

pbwt_reader_t *pbwt_reader_init(const char **fnames, int nfiles)
{
	pbwt_reader_t *reader = calloc(1,sizeof(pbwt_reader_t));
	reader->n        = nfiles;
	reader->pbwt     = malloc(sizeof(PBWT*)*nfiles);
    reader->cpos     = malloc(sizeof(int)*nfiles);
    reader->update   = malloc(sizeof(Update*)*nfiles);
    reader->unpacked = malloc(sizeof(int)*nfiles);

	int i;
	for (i=0; i<nfiles; i++)
	{
		FILE *fp = fopen(fnames[i],"r");
		if ( !fp ) die("failed to open %s: %s\n", fnames[i], strerror(errno));
		reader->pbwt[i] = pbwtRead(fp);
		fclose(fp);

		int j = strlen(fnames[i]);
		char *fname = malloc(sizeof(char)*(j+2));
		memcpy(fname, fnames[i], j);
		memcpy(fname+j-4,"sites",6);
		fp = fopen(fname,"r");
		if ( !fp ) die("failed to open %s: %s\n", fname, strerror(errno));
		free(fname);

		pbwtReadSites(reader->pbwt[i], fp);
		fclose(fp);

        reader->update[i] = updateCreate(reader->pbwt[i]->M, 0);
        reader->unpacked[i] = 0;
        reader->cpos[i] = 0;
    }
    for (i=1; i<nfiles; i++)
    {
        if ( strcmp(reader->pbwt[0]->chrom,reader->pbwt[i]->chrom) ) 
            die("Different chromosomes: %s in %s vs %s in %s\n", reader->pbwt[0]->chrom,fnames[0],reader->pbwt[i]->chrom,fnames[i]);
    }
	return reader;
}

void pbwt_reader_destroy(pbwt_reader_t *reader)
{
	int i;
	for (i=0; i<reader->n; i++)
    {
		pbwtDestroy(reader->pbwt[i]);
        updateDestroy(reader->update[i]);
    }
	free(reader->pbwt);
	free(reader->update);
	free(reader->unpacked);
    free(reader->cpos);
	free(reader);
}

// Return value: 0 if all PBWTs finished or minimum current position 
int pbwt_reader_next(pbwt_reader_t *reader)
{
    int i, min_pos = INT_MAX;
    for (i=0; i<reader->n; i++)
    {  
        PBWT *p = reader->pbwt[i];
        int j   = reader->cpos[i];
        Site *site = arrp(p->sites, j, Site);

        // todo:
        //  - check sequence names, for now assuming single block
        //  - conflicting variations at the same position

        while ( j < p->N && site->x <= reader->mpos  )
            site = arrp(p->sites, ++j, Site);
        reader->cpos[i] = j;

        if ( j < p->N && site->x < min_pos ) min_pos = site->x;
    }
    reader->mpos = min_pos==INT_MAX ? 0 : min_pos;
    return reader->mpos;
}

PBWT *pbwtMerge(const char **fnames, int nfiles)
{
	pbwt_reader_t *reader = pbwt_reader_init(fnames, nfiles);

    int nhaps = 0, i;
    for (i=0; i<nfiles; i++) nhaps += reader->pbwt[i]->M;
    PBWT *out_pbwt  = pbwtCreate(nhaps);
    out_pbwt->yz    = arrayCreate(4096*32, uchar);
    uchar *yz       = myalloc(nhaps, uchar);
    Update *out_up  = updateCreate(nhaps, 0);
    uchar *yseq     = myalloc(nhaps, uchar);
    out_pbwt->sites = arrayReCreate(out_pbwt->sites, reader->pbwt[0]->N, Site);
    out_pbwt->variationDict = dictCreate(32);
    out_pbwt->chrom = strdup(reader->pbwt[0]->chrom);

    int pos, j;
	while ( (pos=pbwt_reader_next(reader)) )
	{
        char *cvar = NULL;

        // read and merge
        int ihap = 0;
        for (i=0; i<nfiles; i++)
        {
            Update *u  = reader->update[i];
            PBWT *p    = reader->pbwt[i];
            Site *site = arrp(p->sites, reader->cpos[i], Site);
            if ( site->x!=pos )
            {
                for (j=0; j<p->M; j++) yseq[ihap + j] = 0;  // missing site, using ref for now
            }
            else
            {
                reader->unpacked[i] += unpack3(arrp(p->yz,reader->unpacked[i],uchar), p->M, u->y, 0);
                for (j=0; j<p->M; j++) yseq[ihap + u->a[j]] = u->y[j];
                updateForwardsA(u);
                if ( !cvar ) cvar = dictName(p->variationDict, site->varD);
            }
            ihap += p->M;
        }

        // pack merged haplotypes
        for (j=0; j<nhaps; j++)
            out_up->y[j] = yseq[out_up->a[j]];
        int nyPack = pack3(out_up->y, out_pbwt->M, yz);
        for (j=0; j<nyPack; j++)
            array(out_pbwt->yz,arrayMax(out_pbwt->yz),uchar) = yz[j];
        updateForwardsA(out_up);

        // insert new site
        arrayExtend(out_pbwt->sites, out_pbwt->N+1);
        Site *site = arrayp(out_pbwt->sites, out_pbwt->N, Site);
        site->x = pos;
        dictAdd(out_pbwt->variationDict, cvar, &site->varD);

        out_pbwt->N++;
	}

    free(yseq);
    updateDestroy(out_up);
    pbwt_reader_destroy(reader);
	return out_pbwt;
}


