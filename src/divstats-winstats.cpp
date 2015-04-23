/* divstats -- a program to calculate window-based diversity statistics
   Copyright (C) 2015  Zachary A Szpiech

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
*/
#include "divstats-winstats.h"



double pi_window(HaplotypeData *hapData, pair_t* snpIndex) {
  //int startSnpIndex; int endSnpIndex;
  double pi = 0;
  double denominator = (hapData->nhaps) * (hapData->nhaps - 1) * 0.5;
  int length = snpIndex->end - snpIndex->start + 1;
  if (length == 0)
  {
    pi = 0;
  }
  else
  {
    for (int i = 0; i < hapData->nhaps; i++)
    {
      for (int j = i + 1; j < hapData->nhaps; j++)
      {
        pi += hamming_dist_ptr(hapData->data[i] + snpIndex->start, hapData->data[j] + snpIndex->start, length);
      }
    }
  }
  return (pi / denominator);
}

