/*  File: dict.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
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
 * Description: based on acedb code from Jean Thierry-Mieg and Richard Durbin 1999-2004
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 23 16:18 2014 (rd)
 * Created: July 2003 (rd)
 *-------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

/****************************************/

static void* remap (void *old, int oldSize, int newSize)
{
  void* new = _mycalloc (newSize, 1) ;
  memcpy (new, old, oldSize) ;
  free (old) ;
  return new ;
}

/****************************************/

static int hashString (char *cp, int n, int isDiff)
{
  int i ;
  unsigned int j, x = 0 ;
  int rotate = isDiff ? 21 : 13 ;
  int leftover = 8 * sizeof(int) - rotate ;

  while (*cp)
    x = (*cp++) ^ ((x >> leftover) | (x << rotate)) ;

  for (j = x, i = n ; i < sizeof(int) ; i += n)
    j ^= (x >> i) ;
  j &= (1 << n) - 1 ;

  if (isDiff)
    j |= 1 ;

  return j ;
}

/*****************************/

DICT *dictCreate (int size)
{
  DICT *dict = mycalloc (1,DICT) ;

  for (dict->dim = 10, dict->size = 1024 ; dict->size < size ; ++dict->dim, dict->size *= 2) ;
  dict->table = mycalloc (dict->size, int) ;
  dict->names = mycalloc (dict->size/2, char*) ;
  return dict ; 
}

/*****************************/

void dictDestroy (DICT *dict)
{
  int i ;
  for (i = 1 ; i <= dict->max ; ++i) free (dict->names[i]) ;
  free (dict->names) ;
  free (dict->table) ;
  free (dict) ;
}

/*****************************/

static int newPos ;		/* communication between dictFind() and dictAdd() */

int dictFind (DICT *dict, char *s, int *ip)
{
  int i, x, d ;

  if (!dict) die ("dictAdd received null dict\n") ;
  if (!s) die ("dictAdd received null string\n") ;

  x = hashString (s, dict->dim, 0) ;
  if (!(i = dict->table[x]))
    { newPos = x ; 
      return 0 ; 
    }
  else if (!strcmp (s, dict->names[i]))
    { if (ip) *ip = i-1 ; 
      return 1 ; 
    }
  else
    { d = hashString (s, dict->dim, 1) ;
      while (1)
	{ x = (x + d) & ((1 << dict->dim) - 1) ;
	  if (!(i = dict->table[x]))
	    { newPos = x ; 
	      return 0 ; 
	    }
	  else if (!strcmp (s, dict->names[i]))
	    { if (ip) *ip = i-1 ; 
	      return 1 ; 
	    }
	}
    }
}

/*****************************/

int dictAdd (DICT *dict, char *s, int *ip)
{
  int i, x ;

  if (dictFind (dict, s, ip)) return 0 ;

  i = ++dict->max ;
  dict->table[newPos] = i ;
  dict->names[i] = myalloc (strlen(s) + 1, char) ;
  strcpy (dict->names[i], s) ;
  if (ip) *ip = i-1 ;

  if (dict->max > 0.3 * dict->size) /* double table size and remap */
    { int *newTable ;
      ++dict->dim ; dict->size *= 2 ;
      dict->names = (char**) remap (dict->names, (dict->max+1)*sizeof(char*), (dict->size/2)*sizeof(char*)) ;
      newTable = mycalloc (dict->size, int) ;
      for (i = 1 ; i <= dict->max ; ++i)
	{ s = dict->names[i] ;
	  x = hashString (s, dict->dim, 0) ;
	  if (!newTable[x])
	    newTable[x] = i ;
	  else
	    { int d = hashString (s, dict->dim, 1) ;
	      while (1)
		{ x = (x + d) & ((1 << dict->dim) - 1) ;
		  if (!newTable[x])
		    { newTable[x] = i ; break ; }
		}
	    }
	}
      free (dict->table) ; dict->table = newTable ;
    }

  return 1 ;
}

/*****************************/

char* dictName (DICT *dict, int i)
{ return dict->names[i+1] ; }

/*********** end of file ***********/
