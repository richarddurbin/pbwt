/*  File: utils.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Modified by Daniel Lawson (dan.lawson@bristol.ac.uk) in December 2014, adding gzip output for paintSparse
 *  Copyright (C) Genome Research Limited, 2011
 * -------------------------------------------------------------------
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *-------------------------------------------------------------------
 * Description: includes standard system headers and own headers
 * Exported functions:
 * HISTORY:
 * Last edited: Jan 20 01:01 2015 (rd)
 * adding gzip output for paintSparse
 * Created: Wed Jan  5 16:13:48 2011 (rd)
 *-------------------------------------------------------------------
 */

#include <stdio.h>		/* FILE etc. */
#include <stdlib.h>		/* malloc(), free(), ... notation */
#include <string.h>		/* memset() */
#include <limits.h>		/* INT_MAX etc. */
#include <errno.h>
#include "zlib.h"

#ifndef BOOL_DEFINED
#define BOOL_DEFINED
typedef char BOOL ;
#define TRUE 1
#define FALSE 0
#endif

#include "array.h"
#include "dict.h"
#include "hash.h"

void die (char *format, ...) ;
void warn (char *format, ...) ;
#define myalloc(n,type)	(type*)_myalloc((n)*sizeof(type))
void *_myalloc (long size) ;
#define mycalloc(n,type) (type*)_mycalloc(n,sizeof(type))
void *_mycalloc (long number, int size) ;
FILE *fopenTag (char* root, char* tag, char* mode) ;
gzFile gzopenTag (char* root, char* tag, char* mode) ;
char *fgetword (FILE *f) ;	/* not threadsafe */
void timeUpdate (FILE *f) ;	/* report to stderr resources used since last called */

/************************/
