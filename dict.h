/*  File: dict.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2003-2008
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
 * Description: header file for DICT package, string hash tables
                developed from the corresponding functions in acedb
		Jean Thierry-Mieg and Richard Durbin 1989-
 * Exported functions:
 * HISTORY:
 * Last edited: Sep 23 16:18 2014 (rd)
 * Created: Sat Dec 20 09:34:14 2008 (rd)
 *-------------------------------------------------------------------
 */

#ifndef DICT_DEFINED
#define DICT_DEFINED

typedef struct {
  char* *names ;
  int *table ;
  int max ;			/* current number of entries */
  int dim ;
  int size ;			/* 2^dim = size of tables */
} DICT ;

DICT *dictCreate (int size) ;
void dictDestroy (DICT *dict) ;
int dictAdd (DICT *dict, char* string, int *index) ;
int dictFind (DICT *dict, char *string, int *index) ;
char* dictName (DICT *dict, int i) ;

#define dictMax(dict)  ((dict)->max)

#endif

/*********** end of file ***********/
