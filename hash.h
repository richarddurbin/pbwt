/*  File: hash.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2011
 *-------------------------------------------------------------------
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
 * Description: based on acedb code from Jean Thierry-Mieg and Richard Durbin 1999-2004
      - can hash integers, floats, pointers
      - returns an integer to be used as index into a local array
      - index of first item to be stored is 1, second is 2 etc.
      - support removal of objects, with free-list to reuse indices of removed objects
 * Exported functions:
 * HISTORY:
 * Last edited: May  5 11:05 2013 (rd)
 * Created: Thu Jan 13 10:57:20 2011 (rd)
 *-------------------------------------------------------------------
 */

#ifndef HASH_DEFINED
#define HASH_DEFINED

typedef void* HASH ;
typedef union {long int i ; float f ; void* p ;} HASHKEY ;
static  HASHKEY _hk ;

/* define keys so as to allow ourselves to hash value 0 */
#define HASH_INT(x) (_hk.i = (x)^INT_MAX, _hk)
#define HASH_FLOAT(x) (_hk.f = (x), _hk.i ^= INT_MAX, _hk)
#define HASH_PTR(x) (_hk.p = (x), _hk)

HASH hashCreate (int n) ;
void hashDestroy (HASH h) ;
void hashClear (HASH h) ;
int  hashFind (HASH h, HASHKEY k) ; /* if found, returns index, else returns 0 */
int  hashAdd  (HASH h, HASHKEY k) ; /* returns index for the key, creating if necessary */
BOOL hashRemove (HASH h, HASHKEY k) ; /* if found, remove and return TRUE, else FALSE */

/* iterator to get all key-value pairs, in arbitrary order */
/* note that if nothing is removed then values are incrementing from 1 to n_added */
/* note also that adding while iterating can invalidate the iterator */
void hashInitIterator (HASH h) ;
int  hashNextKeyValue (HASH h, HASHKEY *kp, int *ip) ; /* returns value, 0 if done */

void hashStats (void) ;		/* overall stats on package/performance */

#endif

/*********** end of file ************/
