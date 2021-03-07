#ifndef  _HAS_HOMOLGY_H
#define _HAS_HOMOLGY_H

int Treat_pair_Percent_homology(int *Datapair, int Lengthtar, int Lengthcan);
float Two_sequence_identity(int *seq1, int *seq2, int Len1, int Len2,  int Nbase, int *pam_tab,int *Score,int (*blast)[21]);
//Nbase=0 for normal, Nbase=-1 for normal and best macthed segment(Segment_Length,Score[0] = Frag100_beg;,Score[1] = Ncommon_match; Score[2] = Segment_Length;)
//Nbase=10 for normal and multiple alignment score()
#endif
