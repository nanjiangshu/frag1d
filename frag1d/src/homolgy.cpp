#include <cstdlib>
#include <cassert>
#include <algorithm>
#include "base.h"
#include "homolgy.h"
using namespace std;

int Treat_pair_Percent_homology(int *Datapair, int Lengthtar, int Lengthcan)/*{{{*/
{
    int Numbre_Pair, ik,ij, ip, Len2_beg, Len1_end, Len_cal_tar,Len_cal_can,Npos,nbeg;
    int Len_Segment, Npair_last_can,Npair_last_tar;
    BYTE byte0, byte1 , *Point_can,*Point_tar;
    bool Number_seg_less500;
    int Xtar_seg[500], Ycan_seg[500], XYLen_seg[500], Number_seg;//suppose the number of consecutive segments < 500
    //int Number_seg_tmp;
    //int NSquare;
    //int Low_end_tar;
    //int Low_end_can; 
    //int High_end_tar; 
    //int High_end_can, NBoundary;
    //int tar_NSquare, can_NSquare ;
    int Max_Numbre_Pair, Keep_Segment, Icode, Num_Data, Max_Point_Line,Point_Line , Max_line_threshhold;
    int Num_blank, BBLANK, Num_Check;
    //int  Max_data;



    //Datapair[]--; -1 for no point, 1000 for point, 0-500 for the number of consecutive segments
    byte0 = BYTE (0);
    byte1 = BYTE (1);
    Numbre_Pair = 2000;
    Keep_Segment = 3;
    BBLANK = 3;
    Max_Point_Line = 0;
    Point_Line = 0;
    Max_line_threshhold = 5;//percent
    Num_Data = Lengthtar*Lengthcan;
    //keep only those pairs that are at least 3 in consecutive
    Number_seg_less500 = true;
    while (   Number_seg_less500  )
    {
        Number_seg_less500 = false;//suppose one time
        Len2_beg = Lengthcan - Keep_Segment;//first position--Lengthcan - 1, at leats 3 segment
        Number_seg = 0;
        for (ik=Len2_beg; ik>=0; ik--) //from top-left to origin
        {
            Len_cal_tar = 0;
            Len_Segment = 0;
            Point_Line = 0;
            Num_blank = 0;
            Num_Check = 0;
            for (ij=ik; ij<Lengthcan; ij++)
            {
                Npos = Len_cal_tar*Lengthcan + ij;
                assert (  Npos < Num_Data  );
                if (  Datapair[Npos] >= 1000  )
                {
                    Len_Segment++;
                    Point_Line++;
                    Num_blank = 0;
                }
                else
                {
                    Num_blank++;
                    if  (   (Num_blank>=BBLANK) || (Len_cal_tar==(Lengthtar-1))   )
                    {
                        if (  Len_Segment < Keep_Segment  )//thrown away the points which are spread-- Len_Segment points
                        {
                            for (ip=0; ip<Num_Check; ip++)
                            {
                                Npos = (Len_cal_tar-ip)*Lengthcan + ij-ip;
                                assert (  Npos < Num_Data  );
                                if (  (Npos>=0) && (Npos<Num_Data)  )
                                {
                                    Datapair[Npos] = -1;
                                }
                            }
                        }
                        else if (  Len_Segment >= Keep_Segment  )
                        {
                            for (ip=0; ip<Num_Check; ip++)
                            {
                                Npos = (Len_cal_tar-ip-1)*Lengthcan + ij-ip-1;
                                assert (  Npos < Num_Data  );
                                if (  (Npos>=0) && (Npos<Num_Data)  )
                                {
                                    if (  Datapair[Npos] >= 0  )
                                    {
                                        Datapair[Npos] = Number_seg;
                                    }
                                }
                            }
                            Xtar_seg[Number_seg] = Len_cal_tar - Len_Segment;
                            Ycan_seg[Number_seg] = ij - Len_Segment;
                            XYLen_seg[Number_seg] = Len_Segment;
                            Number_seg++;
                            if (  Number_seg >= 500  )
                            {
                                //AfxMessageBox("Number_seg >= 500 in CCommonClass::Treat_pair_Percent_homology(BYTE *Datapair, int Lengthtar, int Lengthcan)");
                                //return Numbre_Pair ;
                                ik = -1;
                                Keep_Segment = Keep_Segment + 3;
                                Number_seg_less500 = true;
                                Len_Segment = 0;
                                for (ip=0; ip<Num_Data; ip++)
                                {
                                    if (  Datapair[ip] >= 0   )
                                    {
                                        Datapair[ip] = 1000;
                                    }
                                }
                                break;//ij
                            }
                        }
                        Num_Check = 0;
                        Num_blank = 0;
                        Len_Segment = 0;
                    }
                }
                Len_cal_tar++; 
                Num_Check++;
                if (  Len_cal_tar >= Lengthtar  )
                {
                    break;
                }
            }
            if (  Point_Line > Max_Point_Line )
            {
                Max_Point_Line = Point_Line ;
            }
        }



        Len1_end = Lengthcan - 2;
        for (ik=Len1_end; ik>=0; ik--)  //from origin to right domn
        {
            Len_cal_can = 0;
            Len_Segment = 0;
            Point_Line = 0;
            Num_blank = 0;
            Num_Check = 0;
            nbeg = Lengthcan - 1 - ik;
            for (ij=nbeg; ij<Lengthtar; ij++)
            {
                Npos = ij*Lengthcan + Len_cal_can;
                assert (  Npos < Num_Data  );
                if (  Datapair[Npos] >=1000  )
                {
                    Len_Segment++;
                    Point_Line++;
                    Num_blank = 0;
                }
                else
                {
                    Num_blank++;
                    if (   (Num_blank>=BBLANK) || (Len_cal_can==(Lengthcan-1))   )
                    {
                        if (  Len_Segment < Keep_Segment  )//thrown away the points which are spread-- Len_Segment points
                        {
                            for (ip=0; ip<Num_Check; ip++)
                            {
                                Npos = (ij-ip)*Lengthcan + Len_cal_can-ip;
                                assert (  Npos < Num_Data  );
                                if (  (Npos>=0) && (Npos<Num_Data)  )
                                {
                                    Datapair[Npos] = -1;
                                }
                            }
                        }
                        else if (  Len_Segment >= Keep_Segment  )
                        {
                            for (ip=0; ip<Num_Check; ip++)
                            {
                                Npos = (ij-ip-1)*Lengthcan + Len_cal_can-ip-1;
                                assert (  Npos < Num_Data  );
                                if (  (Npos>=0) && (Npos<Num_Data)  )
                                {
                                    if (  Datapair[Npos] >= 0  )
                                    {
                                        Datapair[Npos] = Number_seg;
                                    }
                                }
                            }
                            Xtar_seg[Number_seg] = ij - Len_Segment;
                            Ycan_seg[Number_seg] = Lengthcan - Len_Segment;
                            XYLen_seg[Number_seg] = Len_Segment;
                            Number_seg++;
                            if (  Number_seg >= 500  )
                            {
                                //AfxMessageBox("Number_seg >= 500 in CCommonClass::Treat_pair_Percent_homology(BYTE *Datapair, int Lengthtar, int Lengthcan)");
                                //return Numbre_Pair;
                                ik = -1;
                                Keep_Segment = Keep_Segment + 3;
                                Number_seg_less500 = true;
                                Len_Segment = 0;
                                for (ip=0; ip<Num_Data; ip++)
                                {
                                    if (  Datapair[ip] >= 0   )
                                    {
                                        Datapair[ip] = 1000;
                                    }
                                }
                                break;//ij
                            }
                        }
                    }
                    Len_Segment = 0;
                    Num_blank = 0;
                    Num_Check = 0;
                }
                Len_cal_can++;
                Num_Check++;
                if (  Len_cal_can >= Lengthcan  )
                {
                    break;
                }
            }
            if (  Point_Line > Max_Point_Line )
            {
                Max_Point_Line = Point_Line ;
            }
        }
    }




    Point_can = new BYTE [Lengthcan+1];
    Point_tar = new BYTE [Lengthtar+1];
    Max_Numbre_Pair = 0;
    Number_seg = 0;


    for (ij=0; ij<Lengthcan; ij++)
    {
        Point_can[ij] = byte0;
    }
    for (ij=0; ij<Lengthtar; ij++)
    {
        Point_tar[ij] = byte0;
    }


    Numbre_Pair = 0;
    for (ij=0; ij<Lengthtar; ij++)
    {
        for (ip=0; ip<Lengthcan; ip++)
        {
            Npos = ij*Lengthcan + ip;
            if  (  Datapair[Npos] >= 0  ) 
            {
                Point_can[ip] = byte1;
                Point_tar[ij] = byte1;
                Numbre_Pair++;
            }
        }
    }

    Npair_last_can = 0;
    for (ij=0; ij<Lengthcan; ij++)
    {
        if (  Point_can[ij] == byte1  )
        {
            Npair_last_can++;
        }
    }

    Npair_last_tar = 0;
    for (ij=0; ij<Lengthtar; ij++)
    {
        if (  Point_tar[ij] == byte1  )
        {
            Npair_last_tar++;
        }
    }

    //Numbre_Pair = __min(Npair_last_can,Npair_last_tar);
    Max_Numbre_Pair = Numbre_Pair;
    Len_Segment = min(Lengthcan,Lengthtar);//shorter sequence
    Icode = Max_Point_Line*100/Len_Segment;
    if (  Icode > Max_line_threshhold  )
    {
        //Max_Numbre_Pair = Numbre_Pair + (Icode-Max_line_threshhold)*Len_Segment;

    }




    delete [] Point_can;
    delete [] Point_tar;

    return Max_Numbre_Pair;
}/*}}}*/
float Two_sequence_identity(int *seq1, int *seq2, int Len1, int Len2, int Nbase, int *pam_tab,int *Score, int (*blast)[21])/*{{{*/
{

    float X_Identity = 0.0; 
    //float XPer100;
    //char *CResult1,*CResult2;
    //char Res1,Res2;
    //int Align_Len,ip,ij,Ncommon_match,Nscore[5];
    //int *Value_align, NUM_Match, Min_Length, Segment_Length, Frag100_beg, Single_step, Length1,Length2;

    /*
       Segment_Length = 100;
       XPer100 = 100;
       Min_Length = __min(Len1,Len2);
       Score[0] = 0;

       Align_Len = Aligment_Two_sequence_500(seq1,seq2,Len1,Len2, CResult, pam_tab, Nscore, "temp");


       Value_align = new int [Align_Len+1];
       NUM_Match = 0;
       for (ip=0; ip<Align_Len; ip++)
       {
       Res1 = CResult[0].GetAt(ip);
       Res2 = CResult[1].GetAt(ip);
       if (  Res1 == Res2  )
       {
       Value_align[ip] = 1;//match,  impossible for two blank, 
       NUM_Match++;
       }
       else
       {
       if (  Res1 == ' ' )//Res1==blank, Res2==aas
       {
       Value_align[ip] = -1;
       }
       else if (  Res2 == ' ' )//Res2==blank, Res1==aas
       {
       Value_align[ip] = -2;
       }
       else//not matched
       {
       Value_align[ip] = 0;
       }
       }
       }
       X_Identity = 0;

       if (  Nbase == 0  )//common sequence identity of alignment
       {
       X_Identity = NUM_Match*XPer100/Min_Length;
       Score[0] = int (NUM_Match*XPer100*10/Len1+0.5);
       Score[1] = int (NUM_Match*XPer100*10/Len2+0.5);
       }
       else if (  Nbase == -1  )//best matched segment that has at least Segment_Length long if two sequences are longer than 100
       {  

       X_Identity = NUM_Match*XPer100/Min_Length;
       if (  Min_Length <= Segment_Length )
       {
       Score[0] = 0;
       Score[1] = int (X_Identity+0.5);
       }
       else
       {
    //search the 100 long segment which has the highest identity
    for (ip=0; ip<Align_Len; ip++)
    {
    if (  Value_align[ip] < 0 )
    {
    Value_align[ip] = -1;
    }
    }
    Ncommon_match = 0;
    Frag100_beg = 0;
    for (ip=0; ip<Align_Len-Segment_Length; ip++)
    {
    Single_step = 0;
    for (ij=ip; ij<ip+Segment_Length; ij++)
    {
    Single_step = Single_step + Value_align[ip];
}
if (  Single_step > Ncommon_match  )
{
    Ncommon_match = Single_step;
    Frag100_beg = ip;
}
}
Score[0] = Frag100_beg;
Score[1] = Ncommon_match;
Score[2] = Segment_Length;
}
}
else if (  Nbase == 10  )//multiple alignment score
{
    X_Identity = NUM_Match*XPer100/Min_Length;
    Length1 = 0;
    Length2 = 0;
    Ncommon_match = 0;
    Single_step = 0;
    for (ip=0; ip<Align_Len; ip++)
    {
        if (  Value_align[ip] != -1  )//sequnece 1 "Seq1" is the aas
        {
            if (  Value_align[ip] >= 0  )//both aas
            {
                Single_step = Single_step + blast[Length1][seq2[Length2]];
                Ncommon_match++;
            }
            Length1++;
        }
        if (  Value_align[ip] != -2  )//sequnece 2 "Seq2" is the aas
        {
            Length2++;
        }
        if  (   (Length1>Len1) || (Length2>Len2)   )
        {
            AfxMessageBox("Ncommon_match > Len1 in Two_sequence_identity");
        }
    }
    Score[0] = Ncommon_match;
    Score[1] = Single_step/Ncommon_match;
}
else //take those away if there are consecutive three blanks
{
    Ncommon_match = 0;
    for (ip=0; ip<Align_Len; ip++)
    {
        if (  Value_align[ip] < 0 )
        {
            Ncommon_match++;
        }
        else
        {
            if (  Ncommon_match >= 3  )
            {
                Frag100_beg = Frag100_beg + Ncommon_match;

            }
            Ncommon_match = 0;
        }
    }
    X_Identity = NUM_Match*XPer100/(Align_Len-Frag100_beg);

}

delete  [] Value_align;
*/
return X_Identity;


}/*}}}*/
