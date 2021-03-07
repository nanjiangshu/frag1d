/*
 * =====================================================================================
 *       Filename:  reportFrag1D.cpp
 *    Description: report the prediction result of Frag1D 
 *
 *        Version:  1.0
 *        Created:  2009-11-08 22:29:24 Sunday Week 45 
 *       Revision:  none
 *       Compiler:  gcc
 *         Author:  Nanjiang Shu (Shu), nanjiang.shu@mmk.su.se
 *        Company:  Structural Chemistry, Stockholm Univesity
 * =====================================================================================
 */

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <set>
#include "array.h"
#include "mytemplate.h"
#include "myfunc.h"
using namespace std;

#if defined(_Windows) || defined(__WINDOWS__) || \
    defined(__WIN32__) || defined(WIN32) || \
defined(__WINNT__) || defined(__NT__)
#   ifndef WINDOWS
#       define WINDOWS
#   endif
#endif

#define NUM_HSR_STATE 3
#define NUM_SHAPE8_STATE 8
#define NUM_SHAPE3_STATE 3

int method_HSR  = 1 ;
int method_SHAPE  = 0 ;
int typeConfidence = 1; /*default using the normlized confidence*/

char HSRalphabet[] = "HSR-";
char S8Alphabet[] = "SRUVKATG-";
char S3Alphabet[] = "SHT-";

int maxSeqLength = 8000;
bool isProofReading = true;
int method_ProofReading = 1;

bool isHSRrestricted = false;
char restrictHSRList[MAX_PATH+1] = "";
string frag1d_version="Frag1D version 1.3";

string reference/*{{{*/="\
# Reference:\n\
#  Tuping Zhou*, Nanjiang Shu* and Sven Hovmoller. A Novel Method for\n\
#  Accurate One-dimensional Protein Structure Prediction based on Fragment\n\
#  Matching, Bioinformatics, 2010;26(4):470-477. (*Co-first author)\n\
";/*}}}*/

string usage=/*{{{*/"\n\
usage: reportFrag1D [-outname STR] [-outpath DIR] [-l LISTFILE] [-h|--help]\n\
                    -datapath DIR ID [ID ...]\n\
Report the result of Frag1D prediction\n\
\n\
  -datapath DIR    Set the path with Frag1D prediction\n\
  -outpath  DIR    Set output path\n\
  -outname  STR    Set the rootname of the output file, (default: query)\n\
                       the output file will be named \n\
                       $outpath/$outname.predfrag1d \n\
                       $outpath/$outname.predfrag1d.html\n\
  -l  LISTFILE     Set the file containing a list of IDs\n\
  -h|--help        Print this help message and exit\n\
\n\
Created on 2009-11-08, updated 2011-10-19, Nanjiang Shu\n\
\n\
Examples:\n\
";/*}}}*/

string headerAnnotation=/*{{{*/"\
# Explanation:\n\
# Num      Residue index in the sequence.\n\
# AA       One-letter amino acid code.\n\
# Sec      Predicted three state secondary structure,\n\
#          (H, S and R).\n\
# ConfSec  Confidence of the predicted secondary structure.\n\
# S8       Predicted 8 state Shape String,\n\
#          (R, S, U, V, A, K, G and T).\n\
# ConfS8   Confidence of the predicted 8 state Shape String.\n\
# S3       Predicted 3 state Shape String\n\
#          (R, S, U, V -> S; A, K -> H; G, T - > T). \n\
# ConfS3   Confidence of the predicted 3 state Shape String.\n\
#";
/*}}}*/
void PrintHelp() {
    cout << usage << endl;
}
void PrintVerboseHelp() { }


void WriteHeader(FILE *fpout)/*{{{*/
{
    fprintf(fpout,"# One-dimensional structure prediction by %s (c) Shu.\n", frag1d_version.c_str());
    fprintf(fpout,"%s\n", reference.c_str());
    fprintf(fpout,"%s\n", headerAnnotation.c_str());
}/*}}}*/
double GetRawHSRConfidence(int probH, int probS, int probR)/*{{{*/
/*Determine raw confidence of the prediction based on the probability on H, S
 * and R
 * 2009-07-21, Nanjiang*/
{
    double sumHSR =  double(probH+probS+probR) + 1e-6; /*plus 1e-6 to avoid division by zero*/
    double maxHSR =  double (max (max(probH, probS), probR ));
    double rawConf = maxHSR/sumHSR;
    return rawConf;
}/*}}}*/
double GetRawS3Confidence(int *probShape8, int n )/*{{{*/
/*Determine raw confidence of the prediction based on the probability on * SRUVKATG * , Nanjiang*/
{
    int probSRUV = probShape8[0] +  probShape8[1] + probShape8[2] + probShape8[3] ; 
    int probKA =   probShape8[4] +  probShape8[5];
    int probGT =   probShape8[6] +  probShape8[7];
    double sumProb3 = double( probSRUV+probKA+probGT) + 1e-6;
    double maxProb3 = double (max (max(probSRUV, probKA), probGT ));
    double rawConf = maxProb3/sumProb3;
    return rawConf;
}/*}}}*/
double GetRawS8Confidence(int *probShape8, int n )/*{{{*/
/*Determine raw confidence of the prediction based on the probability on * SRUVKATG * , Nanjiang*/
{
    double sumProb8 = double(probShape8[0] +  probShape8[1] + probShape8[2] + probShape8[3] + probShape8[4] +  probShape8[5] + probShape8[6] +  probShape8[7])+ 1e-6; /*plus 1e-6 to avoid division by zero*/
    
    double maxProb8 =  double (max_element(probShape8, 0, n-1));
    double rawConf = maxProb8/sumProb8;
    return rawConf;
}/*}}}*/
double GetNormHSRConfidence(int probH, int probS, int probR)/*{{{*/
/*Determine normalized confidence of the prediction based on the probability on H, S
 * and R, parameters were obtained by polynormial regression (order 2) of the
 * raw confidence and two Q3 at different raw confidence bins
 * 2009-08-02 
 * */
{
    double sumHSR =  double(probH+probS+probR) + 1e-6; /*plus 0.00001 to avoid division by zero*/
    double maxHSR =  double (max (max(probH, probS), probR ));
    double rawConf =(maxHSR/sumHSR)*100.0;

    double a1 = 0.8399;
    double a2 = 19.15;
    double normConf = a1 * rawConf + a2;

    normConf /= 100.0;
    if (normConf < 0.0)
    {
        normConf = 0.0;
    }
    else if (normConf> 1.0)
    {
        normConf = 1.0;
    }
    return normConf;

}/*}}}*/
double GetNormS3Confidence(int *probShape8, int n )/*{{{*/
/*Determine normalized confidence of the prediction based on the probability on * SRUVKATG * , Nanjiang*/
{
    int probSRUV = probShape8[0] +  probShape8[1] + probShape8[2] + probShape8[3] ; 
    int probKA =   probShape8[4] +  probShape8[5];
    int probGT =   probShape8[6] +  probShape8[7];
    double sumProb3 = double( probSRUV+probKA+probGT) + 1e-6;
    double maxProb3 = double (max (max(probSRUV, probKA), probGT ));

    double rawConf = maxProb3/sumProb3;
    rawConf *= 100.0;

    double a1 = 0.8602; 
    double a2 = 15.223;
    double normConf;

    if (rawConf >= 30.0)
    {
        normConf = a1 * rawConf + a2;
    }
    else
    {
        normConf = rawConf; /*2009-10-30, the polynomial regression is only approriate for the majority region, not for some weread points*/
    }

    normConf /= 100.0;
    if (normConf <0.0)
    {
        normConf = 0.0;
    }
    else if (normConf >= 1.0)
    {
        normConf = 1.0;
    }
    return normConf;
}/*}}}*/
double GetNormS8Confidence(int *probShape8, int n )/*{{{*/
/*Determine raw confidence of the prediction based on the probability on * SRUVKATG * , Nanjiang*/
{
    double sumProb8 =  Sum(probShape8,0, n-1) + 1e-6; /*plus 1e-6 to avoid division by zero*/
    double maxProb8 =  double (max_element(probShape8, 0, n-1));

    double rawConf = maxProb8/sumProb8;
    rawConf *= 100.0;

    double normConf = 0.0;

    double a1 = 0.9648;
    double a2 = 3.6344;

    if (normConf >= 20.0 )
    {
        normConf = a1 * rawConf + a2;
    }
    else
    {
        normConf = rawConf;
    }

    normConf /= 100.0;
    if (normConf <0.0)
    {
        normConf = 0.0;
    }
    else if (normConf >= 1.0)
    {
        normConf = 1.0;
    }

    return normConf;
}/*}}}*/

int GetHSRState(int probH, int probS, int probR, int method_HSR /*= 1*/)/*{{{*/
/*Determine the state of secondary structure based on the probability on H, S
 * and R
 * the return value is 0 or 1 or 2, representing H, S and R respectively
 * 2009-07-13, Nanjiang*/
{
    int hsrState = 2;
    if (method_HSR == 0)
    {
        if  (    (probS>=probH) && (probS>=probR) ) { hsrState = 1; }
        else if ((probH>=probS) && (probH>=probR) ) { hsrState = 0; }
        else                                        { hsrState = 2; }
    }
    else if (method_HSR == 1)
    {
        if  (    (probR>=probH) && (probR>=probS) ) { hsrState = 2; }
        else if ((probS>=probH) && (probS>=probR) ) { hsrState = 1; }
        else { hsrState = 0; }
    }
    else 
    {
        if  (    (probR>=probH) && (probR>=probS) ) { hsrState = 2; }
        else if ((probH>=probS) && (probH>=probR) ) { hsrState = 0; }
        else { hsrState = 1; }
    }

    return hsrState;
}/*}}}*/
int GetShape8State(int *probShape8, int method_SHAPE /*= 1*/)/*{{{*/
/*Determine the state of shape string based on the probability on 
 * SRUVKATG
 * the return value is 0-7, representing SRUVKATG respectively
 * 2009-07-24, Nanjiang*/
{
    int j = 0;
    int shape8State = 7;
    int sumProb = Sum(probShape8,0,7);
    int adjust_prob = sumProb/15;

    int probSRUV = probShape8[0]+  probShape8[1]+  probShape8[2]+ probShape8[3]; 
    int probKA = probShape8[4]+  probShape8[5];
    int probTG = probShape8[6]+  probShape8[7] + adjust_prob;

    int maxProbSRUV = 0;

    if (method_SHAPE == 0)/*{{{*/
    {
        if  (    (probTG>=probSRUV) && (probTG>=probKA) ) /*GT*/
        {
            if (probShape8[7] >= (probShape8[6]-adjust_prob))
            {
                shape8State = 7; 
            }
            else
            {
                shape8State = 6; 
            }
        }
        else if ((probSRUV>=probKA) && (probSRUV>=probTG) )/*SRUV*/
        {
            shape8State = 0;
            maxProbSRUV = probShape8[0];
            int probSRUV[4] = { 0 };
            for(j = 0; j < 4; j ++)
            {
                probSRUV[j] = probShape8[j];
            }
            probSRUV[2] += adjust_prob; /*U*/
            probSRUV[3] += adjust_prob; /*V*/ 
            for (j=1; j<4; j++)
            {   
                if (  probSRUV[j] >= maxProbSRUV  )
                {
                    maxProbSRUV = probSRUV[j];
                    shape8State = j;
                }
            }
        }
        else /*KA*/
        {
            if (  probShape8[4] >= (probShape8[5]-adjust_prob)  )
            {
                shape8State = 4;
            }
            else
            {
                shape8State = 5;
            }
        }
    }/*}}}*/
    else if (method_SHAPE == 1)/*{{{*/
    {
        if (    (probKA>=probSRUV) && (probKA>=probTG) )/*KA*/
        {
            if (  probShape8[4] >= (probShape8[5]-adjust_prob)  )
            {
                shape8State = 4;
            }
            else
            {
                shape8State = 5;
            }
        }
        else if  ((probSRUV>=probKA) && (probSRUV>=probTG) )/* SRUV*/
        {
            shape8State = 0;
            maxProbSRUV = probShape8[0];
            int probSRUV[4] = { 0 };
            for(j = 0; j < 4; j ++)
            {
                probSRUV[j] = probShape8[j];
            }
            probSRUV[2] += adjust_prob; /*U*/
            probSRUV[3] += adjust_prob; /*V*/ 
            for (j=1; j<4; j++)
            {   
                if (  probSRUV[j] >= maxProbSRUV  )
                {
                    maxProbSRUV = probSRUV[j];
                    shape8State = j;
                }
            }
        }
        else //if  (    (probTG>=probSRUV) && (probTG>=probKA) ) /*GT*/
        {
            if (probShape8[7] >= (probShape8[6]-adjust_prob))
            {
                shape8State = 7; 
            }
            else
            {
                shape8State = 6; 
            }
        }
    }/*}}}*/
    else if (method_SHAPE == 2)/*{{{*/
    {
        if (    (probKA>=probSRUV) && (probKA>=probTG) )/*KA*/
        {
            if (  probShape8[4] >= (probShape8[5]-adjust_prob)  )
            {
                shape8State = 4;
            }
            else
            {
                shape8State = 5;
            }
        }
        else if  (    (probTG>=probSRUV) && (probTG>=probKA) ) /*GT*/
        {
            if (probShape8[7] >= (probShape8[6]-adjust_prob))
            {
                shape8State = 7; 
            }
            else
            {
                shape8State = 6; 
            }
        }
        else /* if ((probSRUV>=probKA) && (probSRUV>=probTG) )// SRUV*/
        {
            shape8State = 0;
            maxProbSRUV = probShape8[0];
            int probSRUV[4] = { 0 };
            for(j = 0; j < 4; j ++)
            {
                probSRUV[j] = probShape8[j];
            }
            probSRUV[2] += adjust_prob; /*U*/
            probSRUV[3] += adjust_prob; /*V*/ 
            for (j=1; j<4; j++)
            {   
                if (  probSRUV[j] >= maxProbSRUV  )
                {
                    maxProbSRUV = probSRUV[j];
                    shape8State = j;
                }
            }
        }
    }/*}}}*/
    else if (method_SHAPE == 3)/*{{{*/
    { /*decrease the wieght on A and increase the weight on R , U and V and T */
        probSRUV = probShape8[0]+  probShape8[1]+  probShape8[2]+ probShape8[3]; 
        probKA = probShape8[4]+  probShape8[5];
        probTG = probShape8[6]+  probShape8[7];
        if (    (probKA>=probSRUV) && (probKA>=probTG) )/*KA*/
        {
            if (  probShape8[4] >= (probShape8[5])  )
            {
                shape8State = 4;
            }
            else
            {
                shape8State = 5;
            }
        }
        else if  (    (probTG>=probSRUV) && (probTG>=probKA) ) /*GT*/
        {
            if (probShape8[7] >= (probShape8[6]))
            {
                shape8State = 7; 
            }
            else
            {
                shape8State = 6; 
            }
        }
        else /* if ((probSRUV>=probKA) && (probSRUV>=probTG) )// SRUV*/
        {
            shape8State = 0;
            maxProbSRUV = probShape8[0];
            int probSRUV[4] = { 0 };
            for(j = 0; j < 4; j ++)
            {
                probSRUV[j] = probShape8[j];
            }
            //probSRUV[2] += adjust_prob; [>U<]
            //probSRUV[3] += adjust_prob; [>V<] 
            for (j=1; j<4; j++)
            {   
                if (  probSRUV[j] >= maxProbSRUV  )
                {
                    maxProbSRUV = probSRUV[j];
                    shape8State = j;
                }
            }
        }
    }/*}}}*/
    else /*{{{*/
    {
        if (    (probKA>=probSRUV) && (probKA>=probTG) )/*KA*/
        {
            if (  probShape8[4] >= (probShape8[5])  )
            {
                shape8State = 4;
            }
            else
            {
                shape8State = 5;
            }
        }
        else if  ((probSRUV>=probKA) && (probSRUV>=probTG) )/* SRUV*/
        {
            shape8State = 0;
            maxProbSRUV = probShape8[0];
            int probSRUV[4] = { 0 };
            for(j = 0; j < 4; j ++)
            {
                probSRUV[j] = probShape8[j];
            }
            probSRUV[2] +=0  ; /*U*/
            probSRUV[3] +=0 ; /*V*/ 
            for (j=1; j<4; j++)
            {   
                if (  probSRUV[j] >= maxProbSRUV  )
                {
                    maxProbSRUV = probSRUV[j];
                    shape8State = j;
                }
            }
        }
        else //if  (    (probTG>=probSRUV) && (probTG>=probKA) ) /*GT*/
        {
            if (probShape8[7] >= (probShape8[6]))
            {
                shape8State = 7; 
            }
            else
            {
                shape8State = 6; 
            }
        }
    }/*}}}*/

    return shape8State;
}/*}}}*/

int ReadInSecPredFile(const char *infile, char *aaSeq, char *obsSec, char *predSec,int fileFormat)/*{{{*/
{
    int linesize;
    int maxline = 400;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    char aa;
    char obs_sec;
    char pred_sec;
    int status_sscanf = 0;

    int seqLength = 0;
    FILE *fp = fopen(infile, "r");
    checkfilestream(fp, infile, "r", true);

    if (fileFormat == 0) /*Res file format*//*{{{*/
    {
        int num ; 
        char shape;
        int probH;
        int probS;
        int probR;
        int predHSRstate = 0; /*predicted HSR state, 0 -- H, 1 -- S, 2 -- R*/
        int cntRes = 0;
        while((linesize=fgetline(fp, line, maxline)) != EOF)
        { 
            status_sscanf = sscanf(line,"%d %c %c %c %d %d %d" , &num,&aa,&shape, &obs_sec,&probH,&probS,&probR);
            if (status_sscanf != 7)
            {
                fprintf(stderr,"sscanf error!");
                fprintf(stderr,"File: %s\n", infile);
                fprintf(stderr,"line: %s\n", line);
                assert (status_sscanf == 7);
            }
            predHSRstate = GetHSRState(probH, probS, probR, method_HSR);

            aaSeq[cntRes] = aa;
            obsSec[cntRes] = obs_sec;
            predSec[cntRes] = HSRalphabet[predHSRstate];
            cntRes ++;
        } 
        seqLength = cntRes;
    }/*}}}*/
    else if (fileFormat == 1) /*AA OSEC PSEC NUM*//*{{{*/
    {
        /*AA OSEC PSEC NUM*/
        linesize=fgetline(fp,line,maxline); /*neglect the header line*/
        int cntRes = 0;
        while((linesize=fgetline(fp, line, maxline)) != EOF)
        { 
            status_sscanf = sscanf(line," %c %c %c" , &aa,&obs_sec,&pred_sec);
            if (status_sscanf != 3)
            {
                fprintf(stderr,"sscanf error!");
                fprintf(stderr,"File: %s\n", infile);
                fprintf(stderr,"line: %s\n", line);
                assert (status_sscanf == 3);
            }
            if (obs_sec == 'E' ) { obs_sec = 'S'; }
            else if (obs_sec == 'C' || obs_sec == 'L') { obs_sec = 'R'; }

            if (pred_sec == 'E' ) { pred_sec = 'S'; }
            else if (pred_sec == 'C' || pred_sec == 'L') { pred_sec = 'R'; }

            aaSeq[cntRes] = aa;
            obsSec[cntRes] = obs_sec;
            predSec[cntRes] = pred_sec;
            cntRes ++;
        } 
        seqLength = cntRes;
    }/*}}}*/
    if (fp != NULL) { fclose(fp); }
    return seqLength;
}/*}}}*/
int ReadInShapePredFile(const char *infile, char *aaSeq, char *obsShape8, char *predShape8, int fileFormat)/*{{{*/
{
    int linesize;
    int maxline = 400;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    char aa;
    char obs_shape8;
    char pred_shape8;
    int status_sscanf = 0;

    int seqLength = 0;
    FILE *fp = fopen(infile, "r");
    checkfilestream(fp, infile, "r", true);

    if (fileFormat == 0) /*Res file format*//*{{{*/
    {
        int num ; 
        char sec;
        int tmpi;
        int probShape8[8] = { 0 }; /*S, R, U, V, K, A, T, G*/
        int predShape8State = 0; /*predicted HSR state, 0 -- H, 1 -- S, 2 -- R*/
        int cntRes = 0;
        while((linesize=fgetline(fp, line, maxline)) != EOF)
        { 
            if (linesize < 1) { continue; }
            status_sscanf = sscanf(line,"%d %c %c %c %d %d %d %d %d %d %d %d %d %d %d %d" , &num,&aa,&obs_shape8, &sec,&tmpi,&tmpi,&tmpi, &tmpi, &probShape8[0], &probShape8[1],&probShape8[2],&probShape8[3],&probShape8[4],&probShape8[5],&probShape8[6],&probShape8[7]);
            if (status_sscanf != 16)
            {
                fprintf(stderr,"sscanf error!");
                fprintf(stderr,"File: %s\n", infile);
                fprintf(stderr,"line: %s\n", line);
                assert (status_sscanf == 16);
            }
            if ((!isHSRrestricted ) || IsInCharSet(sec, restrictHSRList))
            {
                predShape8State = GetShape8State(probShape8, method_SHAPE);

                aaSeq[cntRes] = aa;
                obsShape8[cntRes] = obs_shape8;
                predShape8[cntRes] = S8Alphabet[predShape8State];

                cntRes ++;
            }
        } 
        seqLength = cntRes;
    }/*}}}*/
    else if (fileFormat == 1) /*AA OSHAPE PSHAPE NUM*//*{{{*/
    {
        /*AA OSHAPE PSHAPE NUM*/
        linesize=fgetline(fp,line,maxline); /*neglect the header line*/
        int cntRes = 0;
        while((linesize=fgetline(fp, line, maxline)) != EOF)
        { 
            if (linesize < 1) { continue; }
            status_sscanf = sscanf(line," %c %c %c" , &aa,&obs_shape8,&pred_shape8);
            if (status_sscanf != 3)
            {
                fprintf(stderr,"sscanf error!");
                fprintf(stderr,"File: %s\n", infile);
                fprintf(stderr,"line: %s\n", line);
                assert (status_sscanf == 3);
            }
            if (obs_shape8 == 'E' ) { obs_shape8 = 'S'; }
            else if (obs_shape8 == 'C' || obs_shape8 == 'L') { obs_shape8 = 'R'; }

            if (pred_shape8 == 'E' ) { pred_shape8 = 'S'; }
            else if (pred_shape8 == 'C' || pred_shape8 == 'L') { pred_shape8 = 'R'; }

            aaSeq[cntRes] = aa;
            obsShape8[cntRes] = obs_shape8;
            predShape8[cntRes] = pred_shape8;
            cntRes ++;
        } 
        seqLength = cntRes;
    }/*}}}*/
    if (fp != NULL) { fclose(fp); }
    return seqLength;
}/*}}}*/

int ProofReading(char *predSec, int seqLength, int method_ProofReading = 0)/*{{{*/
{
    //proof reading for predicted result, for 1 long strand (sheet)
    int Part1;
    int Part2;
    int ik;
    int ij;
    if (method_ProofReading == 0)/*{{{*/
    {
        int NShe = 0;
        for (ik=0; ik<seqLength; ik++)
        {
            if (  predSec[ik] == 'S'  )
            {
                NShe++;
            }
            else
            {
                if (  NShe == 1  )
                {
                    //previous
                    Part1 = ik-3;
                    Part2 = ik +1 ;
                    if (  (Part1>=0) && (predSec[Part1]=='S')  )
                    {
                        predSec[Part1+1] = 'S';
                    }
                    else if (  (Part2<seqLength) && (predSec[Part2]=='S')  )
                    {
                        predSec[Part2-1] = 'S';
                    }
                    else
                    {
                        for (ij=ik-NShe; ij<ik; ij++)
                        {
                            predSec[ij] = 'R';
                        }
                    }
                }
                NShe = 0;
            }
        }


        //proof reading for predicted result, for 1 or 2 long helix
        int NHel = 0;
        for (ik=0; ik<seqLength; ik++)
        {
            if (  predSec[ik] == 'H'  )
            {
                NHel++;
            }
            else
            {
                if  (   (NHel>=1) &&  (NHel<=2)   )
                {
                    //previous
                    Part1 = ik-NHel-2;
                    Part2 = ik +1 ;
                    if (  (Part1>=0) && (predSec[Part1]=='H')  )
                    {
                        predSec[Part1+1] = 'H';
                    }
                    else if (  (Part2<seqLength) && (predSec[Part2]=='H')  )
                    {
                        predSec[Part2-1] = 'H';
                    }
                    else
                    {
                        for (ij=ik-NHel; ij<ik; ij++)
                        {
                            predSec[ij] = 'R';
                        }
                    }
                }
                NHel = 0;
            }
        } 
    }/*}}}*/
    else if (method_ProofReading == 1)/*{{{*/
    {  /*do not proof reading single residue sheet, when using the DSSP8to3 scheme BE --> Sheet, 2009-10-29*/
        //proof reading for predicted result, for 1 or 2 long helix
        int NHel = 0;
        for (ik=0; ik<seqLength; ik++)
        {
            if (  predSec[ik] == 'H'  )
            {
                NHel++;
            }
            else
            {
                if  (   (NHel>=1) &&  (NHel<=2)   )
                {
                    //previous
                    Part1 = ik-NHel-2;
                    Part2 = ik +1 ;
                    if (  (Part1>=0) && (predSec[Part1]=='H')  )
                    {
                        predSec[Part1+1] = 'H';
                    }
                    else if (  (Part2<seqLength) && (predSec[Part2]=='H')  )
                    {
                        predSec[Part2-1] = 'H';
                    }
                    else
                    {
                        for (ij=ik-NHel; ij<ik; ij++)
                        {
                            predSec[ij] = 'R';
                        }
                    }
                }
                NHel = 0;
            }
        } 
    }/*}}}*/
    return seqLength;
}/*}}}*/
int ReadConf(const char *infile, double *predConfSec, double *predConfS3, double *predConfS8, int typeConfidence)/*{{{*/
{
    int linesize;
    int maxline = 400;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;
    char aa;
    char obs_sec;
    char obs_shape8;
    int status_sscanf = 0;

    int seqLength = 0;
    FILE *fp = fopen(infile, "r");
    checkfilestream(fp, infile, "r", true);
    int num ; 
    int tmpi;
    int probH;
    int probS;
    int probR;
    int probShape8[8] = { 0 }; /*S, R, U, V, K, A, T, G*/
    int cntRes = 0;
    while((linesize=fgetline(fp, line, maxline)) != EOF)
    { 
        status_sscanf = sscanf(line,"%d %c %c %c %d %d %d %d %d %d %d %d %d %d %d %d" , &num,&aa,&obs_shape8, &obs_sec,&probH, &probS, &probR, &tmpi, &probShape8[0], &probShape8[1],&probShape8[2],&probShape8[3],&probShape8[4],&probShape8[5],&probShape8[6],&probShape8[7]);
        if (status_sscanf != 16)
        {
            fprintf(stderr,"sscanf error!");
            fprintf(stderr,"File: %s\n", infile);
            fprintf(stderr,"line: %s\n", line);
            assert (status_sscanf == 16);
        }

        if (typeConfidence == 0) /*Raw confidence*/
        {
            predConfSec[cntRes] = GetRawHSRConfidence(probH, probS, probR);
            predConfS3[cntRes] = GetRawS3Confidence(probShape8, NUM_SHAPE8_STATE);
            predConfS8[cntRes] = GetRawS8Confidence(probShape8, NUM_SHAPE8_STATE);
        }
        else 
        {
            predConfSec[cntRes] = GetNormHSRConfidence(probH, probS, probR);
            predConfS3[cntRes] = GetNormS3Confidence(probShape8, NUM_SHAPE8_STATE);
            predConfS8[cntRes] = GetNormS8Confidence(probShape8, NUM_SHAPE8_STATE);
        }

        cntRes ++;
    } 
    seqLength = cntRes;
    if (fp != NULL) { fclose(fp); }
    return seqLength;
}/*}}}*/
int Shape8to3(char *shape8, char *shape3, int seqLength)/*{{{*/
{
    for (int i = 0; i< seqLength; i ++) {
        if(shape8[i] == 'S' ||  shape8[i] == 'R' || shape8[i] == 'U' || shape8[i] == 'V' ) {
            shape3[i] = 'S';
        } else if(shape8[i] == 'K' ||  shape8[i] == 'A') {
            shape3[i] = 'H';
        } else if(shape8[i] == 'T' ||  shape8[i] == 'G') {
            shape3[i] = 'T';
        } else {
            shape3[i] = shape8[i];
        }
    }
    return seqLength;
}/*}}}*/
void WriteResultFile(FILE* fpout, char *aaSeq, char *predSec,char *predS3, char *predS8, double *predConfSec, double *predConfS3, double *predConfS8,int seqLength, const char *id, int counter)/*{{{*/
{
    /*header file*/
    fprintf(fpout, "//BEGIN query %d, id = %s\n", counter, id);
    fprintf(fpout,"#%4s %2s %3s %6s %3s %6s %3s %6s\n", "Num", "AA", "Sec", "ConfSec", "S8", "ConfS8", "S3", "ConfS3");
    for (int i=0;i<seqLength; i++) {
        fprintf(fpout,"%5d %2c %3c %6.3lf %3c %6.3lf %3c %6.3lf\n", i+1, aaSeq[i], predSec[i], predConfSec[i], predS8[i], predConfS8[i], predS3[i], predConfS3[i] );
    }
    fprintf(fpout, "//END\n");
}/*}}}*/

void WriteFailMessage(FILE *fpout,  const char *id, int counter)/*{{{*/
{
    fprintf(fpout, "//BEGIN query %d, id = %s\n", counter, id);
    fprintf(fpout,"Frag1D failed for this query");
    fprintf(fpout, "//END\n");
}/*}}}*/
int ReadInIDList(string file, set <string> &idSet)/*{{{*/
{
    ifstream ifp (file.c_str(), ios::binary);
    if (ifp.is_open()){
        ifp.seekg (0, ios::end);
        int length = ifp.tellg();
        ifp.seekg (0, ios::beg);
        char *buffer = new char [ length+1];
        ifp.read(buffer,length);
        buffer[length] = '\0';
        ifp.close();
        char *pch;
        pch = strtok(buffer, "\n");
        while (pch != NULL){
            if ( strlen(pch)> 0 ){
                idSet.insert(pch);
            }
            pch = strtok(NULL, "\n");
        }
        delete [] buffer;
    }else{
        cerr << "Failed to open idlistfile "
            << file 
            << endl;
        return -1;
    }
    return 0;
}/*}}}*/

void GetFrag1DVersion(string rundir)/*{{{*/
{
    string versionFile = rundir + string("/../version");
    ifstream ifp (versionFile.c_str(), ios::binary);
    if (ifp.is_open()){
        getline(ifp, frag1d_version);
        ifp.close();
    }
}/*}}}*/
int main(int argc, char** argv)/*{{{*/
{
    bool isNonOptionArg = false;

    if(argc < 2) {
        PrintHelp();
        return 1;
    }
    int i,j;
    char outname[MAX_PATH+1] = "query";
    char outpath[MAX_PATH+1] = "./";
    char resultPath[MAX_PATH+1] = "";
    char idListFile[MAX_PATH+1] = "";
    const char control_option[] = ""; //options which control the program, and does not take parameters
    int fileFormat = 0;
    
    set <string> ::iterator iss;
    set <string> idList_set;

    i = 1;
    while(i < argc)/*{{{*/ {
        if(argv[i][0] == '-' && !isNonOptionArg) /*options*/ {
            isNonOptionArg = false;
            if(IsInCharSet(argv[i][1], control_option))//if argv[i][1] is in control_option, it might be used as -aqs
            {
                for(j = 1 ; j < int(strlen(argv[i])); j++)
                {
                    switch (argv[i][j])
                    {
                        default : fprintf(stderr,"Invalid option, non-control option '%c' can be used together with contorl-option\n", argv[i][j]); return -1;
                    }
                }
                i ++;
            } else if(strcmp(argv[i],"-h") == 0 ||strcmp(argv[i],"--help")==0 ) {
                PrintHelp(); 
                return 0;
            } else if(strcmp(argv[i],"-H") == 0 ) {
                PrintVerboseHelp();
                return 0;
            } else if( (strcmp(argv[i],"--outname") == 0) || (strcmp(argv[i], "-outname") == 0))  {
                if( ( i = option_parser_filename_old(argc, argv, i, outname)) == -1)
                    return -1;
            } else if( (strcmp(argv[i],"--outpath") == 0) || (strcmp(argv[i], "-outpath") == 0))  {
                if( ( i = option_parser_filename_old(argc, argv, i, outpath)) == -1)
                    return -1;
            } else if( (strcmp(argv[i],"--datapath") == 0) || (strcmp(argv[i],"-datapath") == 0))  {
                if( ( i = option_parser_filename_old(argc, argv, i, resultPath)) == -1)
                    return -1;
            } else if( (strcmp(argv[i],"-l") == 0) || (strcmp(argv[i],"--l") == 0) || (strcmp(argv[i], "--list") == 0))  {
                if( ( i = option_parser_filename_old(argc, argv, i, idListFile)) == -1)
                    return -1;
            } else if (strcmp(argv[i], "--") == 0)//next item is non option argument
            {
                isNonOptionArg = true;
                i ++;
                continue;
            } else {
                fprintf(stderr,"Error! Invalid argument '%s'\n", argv[i]);
                return -1;
            }
        } else /*non-option argument*/ {
            idList_set.insert(argv[i]);
            i ++;
        }
    }/*}}}*/
    struct stat st;
    bool isFolderExist= false;

    if (!(isFolderExist = (stat(resultPath, &st) == 0))){
        cerr << "Datapath not set or does not exist. Exit." << endl;
        return -1;
    }

    if (strcmp(idListFile , "") != 0) {
        ReadInIDList(idListFile, idList_set);
    }
    /*============================= main procedure goes here ================ */
    string rundir;
    Array1D <char> dir(strlen(argv[0]+1));
    getfilepath(argv[0], dir.array1D);
    rundir = dir.array1D;
    GetFrag1DVersion(rundir);

    Array1D <char> aaSeq_1darray(maxSeqLength);
    Array1D <char> predSec_1darray(maxSeqLength);
    Array1D <char> predS3_1darray(maxSeqLength);
    Array1D <char> predS8_1darray(maxSeqLength);
    aaSeq_1darray.Init('-');
    predSec_1darray.Init('-');
    predS3_1darray.Init('-');
    predS8_1darray.Init('-');
    char *aaSeq = aaSeq_1darray.array1D;
    char *predSec = predSec_1darray.array1D;
    char *predS3 = predS3_1darray.array1D;
    char *predS8 = predS8_1darray.array1D;

    Array1D <double> predConfSec_1darray(maxSeqLength);
    Array1D <double> predConfS3_1darray(maxSeqLength);
    Array1D <double> predConfS8_1darray(maxSeqLength);
    double *predConfSec = predConfSec_1darray.array1D;
    double *predConfS3 = predConfS3_1darray.array1D;
    double *predConfS8 = predConfS8_1darray.array1D;

    Array1D <char> obsSec_1darray(maxSeqLength);
    Array1D <char> obsShape8_1darray(maxSeqLength);
    char *obsSec = obsSec_1darray.array1D;
    char *obsShape8 = obsShape8_1darray.array1D;

    int seqLength = 0;
    char resFile[600] = "";
    FILE *fpout = NULL;

    /*output the result*/
    string outfile="";
    outfile += string(outpath) + string("/") + string(outname) + ".predfrag1d";
    fpout = fopen(outfile.c_str(), "w");
    if (fpout == NULL ) {
        cerr << "Failed to output to "
            << outfile 
            << ". Exit."
            << endl;
        return -1;
    }

    WriteHeader(fpout);
    int cnt = 0;
    for(iss = idList_set.begin(); iss != idList_set.end(); iss ++) {
        cnt++;
        /*read in secondary structure prediction file*/
        sprintf(resFile,"%s/%s/Res_%s.txt", resultPath, "chkRstRound2", (*iss).c_str());
        if ((seqLength = ReadInSecPredFile(resFile, aaSeq, obsSec, predSec, fileFormat) < 0)){
            WriteFailMessage(fpout,(*iss).c_str(), cnt);
            continue;
        }
        if (isProofReading) {
            ProofReading(predSec, seqLength, method_ProofReading);
        }
        /*read in shape string prediction file*/
        sprintf(resFile,"%s/%s/Res_%s.txt", resultPath, "chkFstRound2", (*iss).c_str());
        if ((seqLength = ReadInShapePredFile(resFile, aaSeq, obsShape8, predS8, fileFormat)) < 0){
            WriteFailMessage(fpout,(*iss).c_str(), cnt);
            continue;
        }
        Shape8to3(predS8,predS3, seqLength);

        /*get the confidence*/
        sprintf(resFile,"%s/%s/Res_%s.txt", resultPath, "chkFstRound1", (*iss).c_str());
        if ((ReadConf(resFile, predConfSec, predConfS3,predConfS8, typeConfidence)) < 0) {
            WriteFailMessage(fpout,(*iss).c_str(), cnt);
            continue;
        }

        WriteResultFile(fpout, aaSeq, predSec,predS3, predS8,predConfSec,predConfS3, predConfS8,seqLength, (*iss).c_str(), cnt);
    }
    fprintf(stdout,"Frag1D successed. Result has been output to %s\n", outfile.c_str());
    fclose(fpout);
    /*============================= main procedure ends here ================ */
    return 0;
}
/*}}}*/

