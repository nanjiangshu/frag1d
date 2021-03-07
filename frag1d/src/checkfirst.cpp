// checkfirst.cpp : Extracting the primary result from searched fragments files

/*
 * =====================================================================================
 * 
 *        Program:  checkfirst
 *        Description:  predicting the secondary structure, shape string and homology detection.
 *                     1) the source fragment files are from program search_frag;
 *                     2) this routine extracts the primary predicted result which are used for the program checksecond.exe;
          Paraneters:
 *                  -ltst   check_list;
 *                  -qtst   location_qijmatrix;
 *                  -mtst   location_modmatrix;
 *                  -npte   Database_pertemp;
 *                  -sloc   sourcefile;
 *                  -ploc   HRS_profile;
 *                  -rloc   check_Result_location;
 *                  -ltrn   database_list;
 *                  -qtrn   database_qijmatrix;
 *                  -mtrn   database_modmatrix;
 *                  -atrn   database_facc;
  *                 -nptr   NPer_Frag_Database;
 *                  -phsr   HSR_frag;
 *                  --rb    read binary database
 *                  --ratioscheme    set the ratio scheme
 *                  --beg   beg number of the id of test set
 *                  --end   end number of the id of test set
 * =====================================================================================
 */

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include "myfunc.h"
#include "mypro.h"
#include "mytemplate.h"
#include "array.h"
#include "base.h"

using namespace std;

bool isNotPredUnDefined = false; /*whether not predict those without secondary structure definition, default = false*/
bool isTreatZeroNSUM = true; /*whether set the zero NSUM to 1, default = true*/
bool isUseHomologyInfo = true; /*whether use homology info, default = true*/
int threshold_shift = 0;/*ignore those with Number_class  < threshold_shift, default = 0*/ 
int NSpecialBeg = 3;  /*treat specially the first NSpecialBeg residues*/
int NSpecialEnd = 3;  /*treat specially the last NSpecialEnd residues*/
#ifdef DEBUG_NUMCLASS
int debugNumberClass_threshold = 2;/*print out the prediction result for those with Number_class <= debugNumberClass_threshold; 2009-07-11 , nanjiang*/
#endif
int method_HSR = 1; /*using method 1 to determine the predicted secondary structure from three probabilities*/
float min_HSR_confidence = 0.;  /*merge HSR profile only when confidence is between min_HSR_confidence, max_HSR_confidence, 2009-07-14*/
float max_HSR_confidence = 200.; /**/

int dsspMapMethod = 1; /* added 2009-10-14, using the commom dssp8to3 mapping scheme, GHI->Helis, EB->sheet, other->coil*/

int8 typeProfile = 1; /*2009-11-08*/
using namespace std;

char usage[]/*{{{*/="\n\
usage: checkfirst [options]\n\
Description: extract first round result\n\
Options:\n\
  -dbtype INT  (default: 1)\n\
          0: database stored as one id one file\n\
          1: database stored as dumped file\n\
  -ltst   The file path and name for testing set(text file, only one columns:\n\
          chain ID). \n\
  -qtst   Location of qijmatrix for testing set\n\
  -mtst   Location of modmatrix for testing set\n\
  -npte   The percentage of hsrfrag in the testing set (default 30)\n\
  -sloc   Location for source file(searcehed frag files)\n\
  -rloc   Location for predicted result\n\
  -ploc   Location for HSR profiles\n\
  -ltrn   The file path and name for training set(text file, two columns:\n\
          series number and chain ID). \n\
  -qtrn   Location of qijmatrix for training set\n\
  -mtrn   Location of modjmatrix for training set\n\
  -atrn   Location of acc frag matrix for training set\n\
  -nptr   The percentage of accfrag in the training database (default 40)\n\
  -phsr   Whether using the HSR profiles. 1: yes, 0: no, (default: 1)\n\
  --rb    Read binary database\n\
  --nhomo|--not-use-homology   \n\
          Do not use homology information\n\
  --nudef|--not-pred-undefined\n\
          Do not predict residues with secondary structure defined\n\
  --tshift|--threshold-shift INT\n\
          Ignore the first and last N residues, (default: 0)\n\
  --specialhead INT\n\
          First N residues to treat specially, (default: 3)\n\
  --specialtail INT\n\
          Last N residues to treat specially, (default: 3)\n\
  --ratioscheme 0|1|2  \n\
          Set the ratio scheme \n\
  --beg INT  \n\
          Set the starting number of id to be run \n\
  --end INT \n\
          Set the end number of id to be run, beg = 0, end = 1 will run only\n\
          the first one \n\
  --minconf FLOAT\n\
          Set the minimum HSR_confidence above which to merge the HSR_profile \n\
  --maxconf FLOAT\n\
          Set the maximum HSR_confidence below which to merge the HSR_profile\n\
  --dsspmap 0|1     \n\
          Set the dssp mapping method, (default: 1)\n\
          Scheme 0, GHI->Helix, E->Sheet, other->coil \n\
          Scheme 1, GHI->Helix, BE->Sheet, other->coil \n\
\n\
Created 2009-07-01, Updated 2011-10-17\n\
\n\
DEBUG_NUMCLASS\n\
    --thcls INT \n\
           set the threshold for debug number class, so that the program will\n\
           print out result for those with Number_class <= N, default = 3\n\
";/*}}}*/
void PrintHelp() { 
    fprintf(stdout,"%s\n", usage);
} 
 
float GetHSRConfidence(int probH, int probS, int probR)/*{{{*/
/*Determine confidence of the prediction based on the probability on H, S
 * and R
 * 2009-07-21, Nanjiang*/
{
    float sumHSR =  float(probH+probS+probR) + 0.00001; /*plus 0.00001 to avoid division by zero*/
    float maxHSR =  float (max (max(probH, probS), probR ));
    float rawConf = maxHSR/sumHSR;
    float adjustConf = sumHSR / (10 * 100);  /*add roughly 10% confidence to when sumHSR == 100*/
    return min(rawConf + adjustConf, float(100.0)); /*2009-07-21*/
}/*}}}*/
int main(int argc, char* argv[])
{
    if(argc < 2) {
        PrintHelp();
        return 0;
    }
    int dbtype=1; /*0: database stored as one file for each id, 
                    1: database stored as dumped files. 
                    added 2011-10-17*/

    int linesize;
    int maxline = 1000;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    char CSubname[1000] = ""; /*string for seq ID*/

    char first_non_blank_char = ' '; /*2009-07-15, Nanjiang, using first_non_blank_char to determine whether the matrix are valid data lines*/
    string OpenName; /*string for file name*/
    FILE *fpread,*wfp,*wfpper, *fptest, *wppp;
    char  Camino,Cshp, CSAS[10],CCvalue;
    int Sumunit,ik,ij,im,ip, HSR,NAAS,NSeries, temp1;
    NSeries = 0;
    int TPer00;
    char CTempAASseq[LONGEST_SEQ];
    char CTempShapeseq[LONGEST_SEQ];
    char CTempHSRseq[LONGEST_SEQ];
    char CTempDSSPseq[LONGEST_SEQ];

    int  *LengthList,  NPer[25];
    //int iksub,iklen, ScoreSum, RECONTINUE;
    int ScoreSum = 0;
    int Target_Naas[LONGEST_SEQ];
    int GGroup[5], Numberofgroup, NNshape, TotalNumber,Sumnumber;
    int SCore_Sample;
    //int Save_Sample;
    int Consens_Sample;
    float Totalvalue, Svalue, Sumvalue,param1,param2, value, Back_comp[21], Log_per[105],Per_ten_one[105],Log_Back_Comp[21];
    int totalNumValiedPredResidue = 0; /*the number of valied predictions, that is, those predicted residue positions with secondary structure definition, 2009-07-07 */
    float valuelog2,Sumvalue1, Totalvalue1, Zero_deviation;

    Array2D <int> Psi_blosum_2darray(MAX_SEQ_LENGTH, 20);
    Array2D <int> Psi_Frag_2darray(MAX_SEQ_LENGTH, 20);
    int ** Psi_blosum = Psi_blosum_2darray.array2D;
    int ** Psi_Frag = Psi_Frag_2darray.array2D;
    //int Psi_blosum[LONGEST_SEQ][20], Psi_Frag[LONGEST_SEQ][20];


    int Weight_per, NTG_COM, NPer_Frag_Database;
    //candidate
    int  *Candidate_Score_iksub[LONGEST_SEQ],*Candidate_Score_iklen[LONGEST_SEQ], *Candidate_Score_Maxline[LONGEST_SEQ];
    float *Candidate_Score[LONGEST_SEQ], Sum_value;
    int Len_target, NFragment;
    //int Halv_NFragment;
    int HSRsum1,HSRcorrect1,HSRsum,HSRcorrect;
//     int AAS_Code[27], Shape_Code[27];
    int iktest;
    int Pred_Result[LONGEST_SEQ][8],NSUM,SRUV[8], PHSR, Nun_Cand, Num_Candhomo;

    int Pred_Resulthsr[LONGEST_SEQ][4]; /*Pred_Resulthsr[][3] stores the observed secondary structure, 2009-07-03*/

    //char CNameHSR[] = "HSRR" ;
    char CNameHSR[] = "HSR-" ;  /*changed by Nanjiang 2009-07-07*/
    char CShape[] = "SRUVKATG-";
    int *Candidate_Distribution; /*store the number of dots for a candidate chain*/
    int *Candidate_Distribution_per, SFrag, Number_class;

    float weight[LONGEST_SEQ], Xmol[21],   *Pair_data;

    int Number_fam[LONGEST_SEQ], Homology_number, Seg_Wedth, Nkeep;
    int Database_pertemp, Len1_dotline, Len2_dotline,Len_dotline_ID;
    float Homology_Scale;
    //int database_Size, HSR_Cont[3][2], Per_homo ;
    int database_Size = 0;
    int Per_homo = 0;
    Array2D <int> HSR_Cont_2darray(4,2); /*HSR_Cont[4][2], [ij][0] for observed sec, [ij][1] for predicted sec, [ij] for H, S ,R and '-', 2009-07-24*/
    HSR_Cont_2darray.Init(0);
    int ** HSR_Cont = HSR_Cont_2darray.array2D;
    //int Check_file_path = 0;

    /*filenames*/
    char *check_list            = NULL;
    char *check_Result_location = NULL;
    char *sourcefile            = NULL;
    char *location_qijmatrix    = NULL;
    char *location_modmatrix    = NULL;
    char *HRS_profile           = NULL;
    char *database_list         = NULL;
    char *database_modmatrix    = NULL;
    char *database_qijmatrix    = NULL;
    char *database_facc         = NULL;

    int HSR_frag ;/*frag for whether using HSR profile*/

    //check for help

    //default value
    Database_pertemp = 30;//for the HSR fragment files together with qigmatrix profiles
    Homology_Scale = 1;//the weights given to the homologue
    HSR_frag = 1; /*by default, use HSR_frag*/
    SCore_Sample = 100; /*The number of high scoring fragment to use from the target sequence
                        for each fragment of the query sequence*/
    //Save_Sample = 100;
    database_Size = 6100;
    NPer_Frag_Database = 40;
    NTG_COM = 2;
    Num_Candhomo = SCore_Sample;
    Nun_Cand = SCore_Sample;

    bool isReadBinaryFile = false;
    int ratioScheme = 0; /*by default the ratioscheme = 0*/
    int begNum = 0;
    int endNum = 5000000;


    bool isNonOptionArg = false;
    int i,j;
    const char control_option[] = ""; //options which control the program, and does not take parameters

    i = 1;
    while(i < argc)/*{{{*/
    {
        if(argv[i][0] == '-' && !isNonOptionArg) //options
        {
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
            }
            else if(strcmp(argv[i],"-h") == 0 ||strcmp(argv[i],"--help")==0 )
            {
                PrintHelp(); 
                return 0;
            }
            else if( (strcmp(argv[i],"-ltst") == 0) || (strcmp(argv[i], "--idlist-test") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, &check_list)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "-dbtype") == 0) || (strcmp(argv[i], "--dbtype")== 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, dbtype, true, 0, 1)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-qtst") == 0) || (strcmp(argv[i], "--qij-test") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, &location_qijmatrix)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-mtst") == 0) || (strcmp(argv[i], "--modm-test") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, &location_modmatrix)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-npte") == 0) || (strcmp(argv[i], "--nperHSR") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, Database_pertemp, 0, 100)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-sloc") == 0) || (strcmp(argv[i], "--fragpath") == 0)) /*frag files*/ 
            {
                if( ( i = option_parser_filename(argc, argv, i, &sourcefile)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-rloc") == 0) || (strcmp(argv[i], "--outpath") == 0)) /*location for the result files*/ 
            {
                if( ( i = option_parser_filename(argc, argv, i, &check_Result_location)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-ploc") == 0) || (strcmp(argv[i], "--HSRprofile") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, &HRS_profile)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-ltrn") == 0) || (strcmp(argv[i], "--trainlist") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, &database_list)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-qtrn") == 0) || (strcmp(argv[i], "--qij-train") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, &database_qijmatrix)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-mtrn") == 0) || (strcmp(argv[i], "--modm-train") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, &database_modmatrix)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-atrn") == 0) || (strcmp(argv[i], "--fragacc-train") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, &database_facc)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "-nptr") == 0) || (strcmp(argv[i], "--nper")== 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, NPer_Frag_Database, true, 0, 100)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "-phsr") == 0) || (strcmp(argv[i], "--isUseHSRProfile")== 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, HSR_frag, true, 0, 1)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--rb") == 0) || (strcmp(argv[i], "--read-binary")== 0))  
            {
                isReadBinaryFile = true;
                i++;
            }
            else if( (strcmp(argv[i], "--nhomo") == 0) || (strcmp(argv[i], "--not-use-homology")== 0))  
            {
                isUseHomologyInfo = false;
                i++;
            }
            else if( (strcmp(argv[i], "--nudef") == 0) || (strcmp(argv[i], "--not-pred-undefined")== 0))  
            {
                isNotPredUnDefined = true;
                i++;
            }
            else if( (strcmp(argv[i], "--tshift") == 0) || (strcmp(argv[i], "--threshold-shift")== 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, threshold_shift, true, 0, 500)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--specialhead") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, NSpecialBeg, true, 0, 500)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--specialtail") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, NSpecialEnd, true, 0, 500)) == -1)
                    return -1;
            }
#ifdef DEBUG_NUMCLASS
            else if( (strcmp(argv[i], "--thcls") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, debugNumberClass_threshold, true, 0, 500)) == -1)
                    return -1;
            }
#endif
            else if( (strcmp(argv[i], "--ratioscheme") == 0) )  
            {
                if( ( i = option_parser_numeric(argc, argv, i, ratioScheme, true, 0, 10)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--beg") == 0) )  
            {
                if( ( i = option_parser_numeric(argc, argv, i, begNum, true, 0, 5000000)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--end") == 0) )  
            {
                if( ( i = option_parser_numeric(argc, argv, i, endNum, true, 0, 5000000)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--dsspmap") == 0) )  
            {
                if( ( i = option_parser_numeric(argc, argv, i, dsspMapMethod, true, 0, 1)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--minconf") == 0) )  
            {
                if( ( i = option_parser_numeric(argc, argv, i, min_HSR_confidence, true, float(0.0), float(5000.0))) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--maxconf") == 0) )  
            {
                if( ( i = option_parser_numeric(argc, argv, i, max_HSR_confidence, true, float(0.0), float(5000.0))) == -1)
                    return -1;
            }
            else if (strcmp(argv[i], "--") == 0)//next item is non option argument
            {
                isNonOptionArg = true;
                i ++;
                continue;
            }
            else
            {
                fprintf(stderr,"Error! Invalid argument '%s'\n", argv[i]);
                return -1;
            }
        }
        else //non-option argument
        {
            i ++;
        }
    }/*}}}*/
    VerifyFolder(check_Result_location);

    //fprintf(stdout,"finish argument parsing\n");

    //check datasize
//     fp = fopen(database_list,"r");
//     sumunit = 0;
//     while (  fgets(line,300,fp) != null  )
//     {
//         sumunit++;
//     }
//     fclose(fp);
    Sumunit = fgetlinecnt(database_list);
    database_Size = Sumunit + 1;

    char *(*Namelistunit),*(*CAASSequence), *(*CShapesequence),*(*CHSRsequence);
    int *(*Chain_code_AAS), *(*blosum_chain);

    Namelistunit = new  ( char (*[database_Size]) );
    CAASSequence = new  ( char (*[database_Size]) );
    CShapesequence = new  ( char (*[database_Size]) );
    CHSRsequence = new  ( char (*[database_Size]) );
    Chain_code_AAS = new ( int (*[database_Size]) );
    blosum_chain = new ( int (*[database_Size]) );
    LengthList = new int [database_Size];
    Candidate_Distribution = new int [database_Size];
    Candidate_Distribution_per = new int [database_Size];

    NFragment = 9;
    //Halv_NFragment = NFragment/2;
    Nkeep = 0;

    //sprintf(OpenName,"%s/detailed.txt",check_Result_location);
    OpenName = string(check_Result_location) + "/detailed.txt";
    wppp = fopen(OpenName.c_str(),"w");
    fclose(wppp);

    TotalNumber = 0;
    Totalvalue = 0;
    Totalvalue1 = 0;
    HSRsum = 0;
    HSRcorrect = 0;

    Zero_deviation = float ( 0.25 );
    valuelog2 = float ( log(2) );
    for (ik=0; ik<103; ik++)
    {
        Per_ten_one[ik] = float ( ik + Zero_deviation ); 
        Log_per[ik] = float ( log(Per_ten_one[ik])/valuelog2 );
    }
    //background composition--Back_comp
    Back_comp[0] = float (7.88);
    Back_comp[1] = float (6.73);
    Back_comp[2] = float (9.65);
    Back_comp[3] = float (5.90);
    Back_comp[4] = float (4.83);
    Back_comp[5] = float (3.96);
    Back_comp[6] = float (2.38);
    Back_comp[7] = float (5.93);
    Back_comp[8] = float (5.40);
    Back_comp[9] = float (2.29);
    Back_comp[10] = float (6.96);
    Back_comp[11] = float (6.83);
    Back_comp[12] = float (5.40);
    Back_comp[13] = float (1.50);
    Back_comp[14] = float (3.03);
    Back_comp[15] = float (4.13);
    Back_comp[16] = float (6.67);
    Back_comp[17] = float (1.14);
    Back_comp[18] = float (5.34);
    Back_comp[19] = float (3.95);


    for (ik=0; ik<20; ik++)
    { 
        Log_Back_Comp[ik] = float ( log(Back_comp[ik])/valuelog2 );//see Per_ten_one[ik]
    }

    for (ik=0; ik<Sumunit; ik++)
    {
        LengthList[ik] = 0;
    }
    //read the database into  Namelistunit, CAASSequence, CShapesequence, CHSRsequence, Chain_code_AAS, blosum_chain, LengthList
    //Read_databse(database_list, NPer_Frag_Database, database_qijmatrix, database_modmatrix, database_facc,  check_Result_location,  Namelistunit, CAASSequence,CShapesequence, CHSRsequence, blosum_chain, Chain_code_AAS, LengthList );

    fprintf(stdout,"Read training set in %s format with ratioScheme = %d ...\n", isReadBinaryFile == true ? "binary" : "text", ratioScheme);
    Sumunit = Read_databse_SHU(dbtype, database_list, NPer_Frag_Database, ratioScheme,  database_qijmatrix, database_modmatrix, database_facc,  check_Result_location,  Namelistunit, CAASSequence,CShapesequence, CHSRsequence, blosum_chain, Chain_code_AAS, LengthList, dsspMapMethod, isReadBinaryFile , typeProfile);
    fprintf(stdout,"Read training set finished, %d chains read in.\n", Sumunit);

#ifdef DEBUG_NUMCLASS /*2009-07-11*/
        char tmpDebugNumClassFile[1000] = "";
        sprintf(tmpDebugNumClassFile, "%s/resultNumberClassLE%d.dat", check_Result_location, debugNumberClass_threshold);
        FILE *fwpTmpDebugNumberClass = fopen(tmpDebugNumClassFile, "w");
#endif 

    iktest = -1;
    Homology_Scale = float (5000.0/Sumunit); /*scale to normalize the raw homology score, which is actually the number_of_dots / average_length */
    fptest = fopen(check_list,"r");

    char alphabet[50] = "";
    double parameter[8];
    int seqLength = 0;
    Array1D <ProfileSADByte> profileSADByte_1darray(LONGEST_SEQ);
    ProfileSADByte * profileSADByte = profileSADByte_1darray.array1D;

    fprintf(stdout,"Start secondary structure prediction...\n");
    /*while (  fgets(line,100,fptest) != NULL  )*/
    while ((linesize = fgetline(fptest,line,maxline)) != EOF)
    { 
        if (linesize <= 0) {
            continue;
        }

        sscanf(line, "%s",CSubname);
        iktest++;

        if (iktest < begNum || iktest >= endNum) {
            continue;
        }
        fprintf(stdout,"%d predicting the sequence %s\n", iktest, CSubname);

        /*initialize*/
        for (ik=0; ik<LONGEST_SEQ; ik++)
        {
            for (ij=0; ij<20; ij++)
            {
                Psi_blosum[ik][ij] = 0;
            }
            Target_Naas[ik] = 20;
            //Pred_Resulthsr[ik][3] = 2; [>initialize the observed HSR by 'R'<]
            Pred_Resulthsr[ik][3] = 3; /*initialize the observed HSR by '-', changed by Nanjiang 2009-07-07*/
            CTempHSRseq[ik] = '-';
        }


        if (!isReadBinaryFile)
        {
            //sprintf(OpenName,"%s/%s.Qij",location_qijmatrix,CSubname);
            OpenName = string(location_qijmatrix) + "/" + string(CSubname) + ".Qij";
            fpread = fopen(OpenName.c_str(),"r");

            if (  fpread == NULL  ) {
                continue;
            }
            while((linesize = fgetline(fpread, line,maxline)) != EOF)
            {
                sscanf(line, " %c", &first_non_blank_char);
                if (first_non_blank_char <'0' || first_non_blank_char > '9') /*only when the first_non_blank_char is digit */
                {
                    continue;
                }
                sscanf(line," %d ", &NSeries );
                NSeries = NSeries - 1;
                //percent
                sscanf(line,"%d %c %s %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%f%f",
                        &HSR,&Camino,CSAS,&NPer[0],&NPer[1],&NPer[2],&NPer[3],&NPer[4],&NPer[5],&NPer[6],&NPer[7],&NPer[8],&NPer[9],
                        &NPer[10],&NPer[11],&NPer[12],&NPer[13],&NPer[14],&NPer[15],&NPer[16],&NPer[17],&NPer[18],&NPer[19],&param1,&param2);
                Cshp = CSAS[0];
                NAAS = Camino - 'A';
                if (  (NAAS>=0) && (NAAS<=25)  )
                {
                    NAAS = AAS_Code[NAAS];
                }
                else
                {
                    NAAS = 20;
                }
                Target_Naas[NSeries] = NAAS;

                if (  NAAS >= 20  )
                {
                    //amino acids
                    CTempAASseq[NSeries] = 'X';  /*for non-standard amino acids, set the amino acid type to 'X' and set the shape string and secondary structure to '-', 2009-07-15*/
                    //shape string
                    CTempShapeseq[NSeries] = '-';
                    CTempHSRseq[NSeries] = '-'; /*observed HSR string, '-' --> '-' 2009-07-10,  */
                    CTempDSSPseq[NSeries] = '-'; /**/
                    Pred_Resulthsr[NSeries][3] = 3;
                    continue;
                }
                TPer00 = 0;
                for (ij=0; ij<20; ij++)
                {
                    if (  NPer[ij] >= 1  )
                    {
                        TPer00 = 1;
                        break;
                    }
                }
                if (  TPer00 == 0 )
                {
                    NPer[NAAS] = 100;
                }
                for (ij=0; ij<20; ij++)
                {    
                    Psi_blosum[NSeries][ij] = NPer[ij];
                }

                //amino acids
                CTempAASseq[NSeries] = Camino;
                //shape string
                CTempShapeseq[NSeries] = Cshp;
                //HSR
                /*map DSSP  8 state to 3 state with the scheme:
                 * H, G, I --> H
                 * E       --> S
                 * '-'     --> '-'
                 * rest    --> R
                 * 2009-07-03 */
                if (  (CSAS[2]=='H') || (CSAS[2]=='G') || (CSAS[2]=='I')  )
                    //if (  (CSAS[2]=='H')  )
                {
                    temp1 = 0;
                    CCvalue = 'H';
                }
                //else if (  (CSAS[2]=='E') || (CSAS[2]=='B')    )
                else if (  (CSAS[2]=='E') || (dsspMapMethod == 1 && CSAS[2] == 'B')     )
                {
                    temp1 = 1;
                    CCvalue = 'S';
                }
                else if (  (CSAS[2] != '-')     )
                {
                    temp1 = 2;
                    CCvalue = 'R';
                }
                else             /*added by Nanjiang 2009-07-07*/
                {
                    temp1 = 3;
                    CCvalue = '-';
                }
                //CTempHSRseq[NSeries] = CSAS[0]; [>observed HSR string, '-' --> '-' ??? isn't CASA[0] = shape ?  <]
                CTempHSRseq[NSeries] = CCvalue; /*observed HSR string, '-' --> '-' 2009-07-10,  */
                CTempDSSPseq[NSeries] = CSAS[2]; /**/
                Pred_Resulthsr[NSeries][3] = temp1;
            }
            fclose(fpread);
            NSeries++;
            Len_target = NSeries;
        }
        else/*{{{*/
        {   /*Read In bindary Qijfile for the test set*/
            //sprintf(OpenName,"%s/%s.Qijbin",location_qijmatrix,CSubname);
            OpenName = string(location_qijmatrix) + "/" + string(CSubname) + ".Qijbin";
            if (GetBinaryMODM(OpenName.c_str(), alphabet, seqLength, profileSADByte,
                        parameter, typeProfile) == -1) {
                fprintf(stderr, "can not open QijFile %s\n", OpenName.c_str());
                return -1;
            }
            for (ik = 0; ik<seqLength; ik++)
            {
                Camino = profileSADByte[ik].aa;
                NAAS = Camino - 'A';
                if (  (NAAS>=0) && (NAAS<=25)  )
                {
                    NAAS = AAS_Code[NAAS];
                }
                else
                {
                    NAAS = 20;
                }
                Target_Naas[ik] = NAAS;

                if (  NAAS >= 20  )
                {
                    //amino acids
                    CTempAASseq[ik] = 'X';  /*for non-standard amino acids, set the amino acid type to 'X' and set the shape string and secondary structure to '-', 2009-07-15*/
                    //shape string
                    CTempShapeseq[ik] = '-';
                    CTempHSRseq[ik] = '-'; /*observed HSR string, '-' --> '-' 2009-07-10,  */
                    CTempDSSPseq[ik] = '-'; /**/
                    Pred_Resulthsr[ik][3] = 3;
                }
                TPer00 = 0;
                for (ij=0; ij<20; ij++)
                {
                    if (  profileSADByte[ik].p[ij] >= 1  )
                    {
                        TPer00 = 1;
                        break;
                    }
                }
                if (  TPer00 == 0 )
                {
                    profileSADByte[ik].p[NAAS] = 100;
                }
                for (ij=0; ij<20; ij++)
                {    
                    Psi_blosum[ik][ij] = profileSADByte[ik].p[ij];
                }

                //amino acids
                CTempAASseq[ik] = profileSADByte[ik].aa;
                //shape string
                CTempShapeseq[ik] = profileSADByte[ik].shape;
                char dsspSec = profileSADByte[ik].dsspSec;
                //HSR
                /*map DSSP  8 state to 3 state with the scheme:
                 * H, G, I --> H
                 * E       --> S
                 * '-'     --> '-'
                 * rest    --> R
                 * 2009-07-03 */
                if (  (dsspSec=='H') || (dsspSec=='G') || (dsspSec=='I')  )
                    //if (  (dsspSec=='H')  )
                {
                    temp1 = 0;
                    CCvalue = 'H';
                }
                //else if (  (dsspSec=='E') || (dsspSec=='B')    )
                else if (  (dsspSec=='E') || (dsspMapMethod == 1 && dsspSec == 'B')     )
                {
                    temp1 = 1;
                    CCvalue = 'S';
                }
                else if (  (dsspSec != '-')     )
                {
                    temp1 = 2;
                    CCvalue = 'R';
                }
                else             /*added by Nanjiang 2009-07-07*/
                {
                    temp1 = 3;
                    CCvalue = '-';
                }
                //CTempHSRseq[NSeries] = CSAS[0]; [>observed HSR string, '-' --> '-' ??? isn't CASA[0] = shape ?  <]
                CTempHSRseq[ik] = CCvalue; /*observed HSR string, '-' --> '-' 2009-07-10,  */
                CTempDSSPseq[ik] = dsspSec; /**/
                Pred_Resulthsr[ik][3] = temp1;
            }
            Len_target = seqLength;
        }/*}}}*/

#ifdef DEBUG_READQIJ
        /*print the Qij matrix read in 2009-07-24*/
        fprintf(stdout, "Read in %s\n", OpenName.c_str());
        for (ik = 0; ik < Len_target; ik ++)
        {
            fprintf(stdout,"%4d %c ", ik, CTempAASseq[ik]);
            for(ij = 0; ij < 20; ij ++)
            {
                fprintf(stdout, " %3d", Psi_blosum[ik][ij]);
            }
            fprintf(stdout,"\n");
        }
        fprintf(stdout, "\n");
#endif

        //read from aacfrag_HSR
        if (   HSR_frag == 1  )  /*use HSR profile, 2009-07-03*/
        {
            if(!isReadBinaryFile)
            {
                //modmatrix
                //sprintf(OpenName,"%s/%s.modm",location_modmatrix,CSubname);
                OpenName = string(location_modmatrix) + "/" + string(CSubname) + ".modm";
#ifdef DEBUG_READQIJ
                /*print the Qij matrix read in 2009-07-24*/
                fprintf(stdout,"\n");
                fprintf(stdout, "Read in %s\n", OpenName.c_str());
#endif
                fpread = fopen(OpenName.c_str(),"r");
                if (  fpread != NULL  )
                { 
                    /*while (  fgets(line,300,fpread) != NULL  )*/
                    while((linesize = fgetline(fpread, line,maxline)) != EOF)
                    {   
                        sscanf(line, " %c", &first_non_blank_char);
                        if (first_non_blank_char <'0' || first_non_blank_char > '9') /*only when the first_non_blank_char is digit */
                        {
                            continue;
                        }
                        sscanf(line," %d ", &NSeries );
                        NSeries = NSeries - 1;
                        //percent
                        sscanf(line,"%d %c %s %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%f%f",
                                &HSR,&Camino,CSAS,&NPer[0],&NPer[1],&NPer[2],&NPer[3],&NPer[4],&NPer[5],&NPer[6],&NPer[7],&NPer[8],&NPer[9],
                                &NPer[10],&NPer[11],&NPer[12],&NPer[13],&NPer[14],&NPer[15],&NPer[16],&NPer[17],&NPer[18],&NPer[19],&param1,&param2);

                        for (ij=0; ij<20; ij++)
                        {     
                            Psi_Frag[NSeries][ij] = NPer[ij];
                        } 
#ifdef DEBUG_READQIJ
                        /*print the Qij matrix read in 2009-07-24*/
                        fprintf(stdout,"%4d %c ", NSeries, Camino);
                        for(ij = 0; ij < 20; ij ++)
                        {
                            fprintf(stdout, " %3d", Psi_Frag[NSeries][ij]);
                        }
                        fprintf(stdout,"\n");
#endif
                    }   
                    fclose(fpread);
                }
            } 
            else/*{{{*/
            {   /*Read In bindary Qijfile for the test set*/
                //sprintf(OpenName,"%s/%s.modmbin",location_modmatrix,CSubname);
                OpenName = string(location_modmatrix) + "/" + string(CSubname) + ".modmbin";
                if (GetBinaryMODM(OpenName.c_str(), alphabet, seqLength,
                            profileSADByte, parameter, typeProfile) == -1) {
                    fprintf(stderr, "can not open QijFile %s\n", OpenName.c_str());
                    return -1;
                }
                for (ij=0; ij<20; ij++) {
                    Psi_Frag[ik][ij] = profileSADByte[ik].p[ij];
                } 
            }/*}}}*/

            /*Read in HSR profile, note that the name scheme is the same
             * as shape string profile, 2009-07-03 */
            //sprintf(OpenName,"%s/Frag_%s.txt",HRS_profile,CSubname);
            OpenName = string(HRS_profile) + "/Frag_" + string(CSubname) + ".txt";
#ifdef DEBUG_READQIJ
            /*print the Qij matrix read in 2009-07-24*/
            fprintf(stdout,"\n");
            fprintf(stdout, "Read in %s\n", OpenName.c_str());
#endif
            fpread = fopen(OpenName.c_str(),"r");
            if (  fpread != NULL  )
            { 
                //while (  fgets(line,300,fpread) != NULL  )
                while ((linesize = fgetline(fpread, line, maxline)) != EOF)
                {   
                    sscanf(line, " %c", &first_non_blank_char);
                    if (first_non_blank_char <'0' || first_non_blank_char > '9') /*only when the first_non_blank_char is digit */
                    {
                        continue;
                    }
                    sscanf(line," %d ", &NSeries );
                    //for (ij=0; ij<5; ij++)
                    //{ 
                        //str_char1[ij] = line[ij];
                    //} 
                    //str_char1[ij] = '\0';
                    //NSeries = -1; [>2009-07-15 <]
                    //NSeries = atoi(str_char1);
                    //if (   (NSeries<1) || (NSeries>=Max_Length)   )//for the first two lines  
                    //{ 
                        //continue;
                    //} 
                    NSeries = NSeries - 1;
                    //percent
                    /*2009-07-14, The last column of the Frag_$id.txt file is
                     * the confidence of the secondary structure prediction*/
                    sscanf(line,"%d %c %s %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%f%f",
                            &HSR,&Camino,CSAS,&NPer[0],&NPer[1],&NPer[2],&NPer[3],&NPer[4],&NPer[5],&NPer[6],&NPer[7],&NPer[8],&NPer[9],
                            &NPer[10],&NPer[11],&NPer[12],&NPer[13],&NPer[14],&NPer[15],&NPer[16],&NPer[17],&NPer[18],&NPer[19],&param1,&param2);
                    HSR = 0;

#ifdef DEBUG_READQIJ
                    /*print the Qij matrix read in 2009-07-24*/
                    fprintf(stdout,"%4d %c ", NSeries, Camino);
                    for(ij = 0; ij < 20; ij ++)
                    {
                        fprintf(stdout, " %3d", NPer[ij]);
                    }
                    fprintf(stdout,"\n");
#endif
                    for (ij=0; ij<20; ij++)
                    { 
                        if (param2 >= min_HSR_confidence && param2 <=max_HSR_confidence) /*merge the HSRprofile only when the confidence (param2 >= min_HSR_confidence && param2 <= max_HSR_confidence), 2009-07-14 */
                        {
                            NPer[ij] = NPer[ij]*Database_pertemp + Psi_Frag[NSeries][ij]*(100-Database_pertemp);
                        }
                        HSR = HSR + NPer[ij];

                    } 
#ifdef DEBUG_NPTR /*debug the merge ratio of modm and HSRprofile, 2009-07-24*/
                    if(Database_pertemp != 50)
                    {
                        fprintf(stderr,"Database_pertemp changed by illegal write. Database_pertemp = %d ID=%s NSeries=%d\n", Database_pertemp, CSubname, NSeries) ;
                    }
#endif
                    for (ij=0; ij<20; ij++)
                    {     
                        Psi_blosum[NSeries][ij] = int ( NPer[ij]*100.0/HSR + 0.4999 );
                        /*debug for Psi_blosum*/
                        if (Psi_blosum[NSeries][ij] < 0 || Psi_blosum[NSeries][ij] > 100 )
                        {
                            fprintf(stderr,"%s %d Psi_blosum error, Psi_blosum[%d][%d] = %d\n", CSubname, NSeries, NSeries, ij,Psi_blosum[NSeries][ij] );
                        }
                        /*debug for Psi_blosum*/
                    } 
                }   
                fclose(fpread);
            } 
        }

#ifdef DEBUG_READQIJ
        fprintf(stdout,"\n");
        fprintf(stdout, "Psi_blosum %s\n", CSubname);
        for (ik = 0; ik < Len_target; ik ++)
        {
            fprintf(stdout,"%4d %c ", ik, CTempAASseq[ik]);
            for(ij = 0; ij < 20; ij ++)
            {
                fprintf(stdout, " %3d", Psi_blosum[ik][ij]);
            }
            fprintf(stdout,"\n");
        }
        fprintf(stdout, "\n");
#endif
        //sequence weight
        /*the weight for the profile-profile score2, see the manuscript*/
        for (ik=0; ik<Len_target; ik++)
        {
            weight[ik] = 0;
            HSR = 0;
            value = 0;
            for (im=0; im<20; im++)
            { 
                Xmol[im] = Psi_blosum[ik][im]/Back_comp[im];
                value = value + Xmol[im];
                HSR = HSR + Psi_blosum[ik][im];
            }
            if (  HSR <50  )    /*this HSR is the sum of profile, should be around 100*/
            {
                continue;
            }
            for (im=0; im<20; im++)
            { 
                Xmol[im] = Xmol[im]/value;
            }
            value = 0;
            for (im=0; im<20; im++)
            { 
                if (  Xmol[im] >= 0.00001  )
                {
                    value = value + float (Xmol[im]*log(Xmol[im]));
                }
            }
            weight[ik] = (-value+1)*(-value+1);
        }


        /*initialize*/
        for (ik=0; ik<Len_target; ik++)
        {
            Candidate_Score[ik] = new float [SCore_Sample+5];
            Candidate_Score_iksub[ik] = new int [SCore_Sample+5];
            Candidate_Score_iklen[ik] = new int [SCore_Sample+5];
            Candidate_Score_Maxline[ik] = new int [SCore_Sample+5];
        }

        for (ik=0; ik<Len_target; ik++)
        {
            for (ij=0; ij<SCore_Sample+5; ij++)
            {
                Candidate_Score[ik][ij] = 0;
                Candidate_Score_Maxline[ik][ij] = 0;
            }
        }
        for (ik=0; ik<database_Size; ik++)
        {
            Candidate_Distribution[ik] = 0;
            Candidate_Distribution_per[ik] = 0;

        }

        //open the result file
        //sprintf(OpenName,"%s/%s.txtbin",sourcefile,CSubname);
        OpenName = string(sourcefile) + "/" + string(CSubname) + ".txtbin";
        /* Candidate_Score [m][100]
         * Candidate_Score_iksub [m][100]: the index of chain
         * Candidate_Score_iklen [m][100]: the position of the candidate fragment*/
        ReadInPairBinary( OpenName.c_str(),  Candidate_Score, Namelistunit,
                Sumunit, Candidate_Score_iksub, Candidate_Score_iklen);

        if (isUseHomologyInfo)
        {
            for (ik=0; ik<Len_target; ik++)
            {
                for (ij=0; ij<SCore_Sample+1; ij++)
                {
                    if (Candidate_Score[ik][ij] > 1)
                    {
                        Candidate_Distribution[Candidate_Score_iksub[ik][ij]] ++;
                    }
                }
            }
        }


        Homology_number = 0;
        // for (ik=0; ik<database_Size; ik++)
        for (ik=0; ik<Sumunit; ik++)/* change database_size to Sumunit to avoid the error: jump depend on uninsitialized values. 2009-07-12, */
        {
            HSR = int (Candidate_Distribution[ik]*2.0*100.0/(Len_target+LengthList[ik])/Homology_Scale);
            Per_homo = HSR;
            if (  HSR >= 10 )
            {
                //which is longer Len1_dotline >= Len2_dotline, Len_dotline_ID;
                if (  Len_target >= LengthList[ik] )
                {
                    Len1_dotline = Len_target;
                    Len2_dotline = LengthList[ik];
                    Len_dotline_ID = 0;//iktest, target
                }
                else
                {
                    Len1_dotline = LengthList[ik];
                    Len2_dotline = Len_target;
                    Len_dotline_ID = 1;//ik, candidate
                }

                ScoreSum = Len1_dotline*Len2_dotline+1;
                Pair_data = new float [ScoreSum];
                for (im=0; im<ScoreSum; im++)
                {
                    Pair_data[im] = -1;
                }
                for (im=0; im<Len_target; im++)
                {
                    for (ip=0; ip<Num_Candhomo; ip++)
                    {
                        if (  Candidate_Score[im][ip] <= 1   )
                        {
                            continue;
                        }
                        if (  (Candidate_Score_iksub[im][ip]==ik) && (Candidate_Score_iklen[im][ip]>=0) && (Candidate_Score_iklen[im][ip]<LengthList[ik])  )
                        {
                            if (  Len_dotline_ID == 0  )
                            {
                                NAAS = im*Len2_dotline + Candidate_Score_iklen[im][ip];
                            }
                            else
                            {
                                NAAS = Candidate_Score_iklen[im][ip]*Len2_dotline + im;
                            }
                            Pair_data[NAAS] = float (im*Num_Candhomo + ip);  /*this dot is filled and also show the location*/
                        }
                    }
                }

                NAAS = Treat_pair_Percent_homology(Pair_data, Len1_dotline,
                        Len2_dotline, NPer, Candidate_Score_Maxline,
                        Num_Candhomo);
                HSR = int (NAAS*2.0*100.0/(Len_target+LengthList[ik])/Homology_Scale);
                delete [] Pair_data;


                if (  HSR >= 6 )
                {
                    Candidate_Distribution_per[ik] = HSR;
                    //printf("HSR=%d\n",HSR);
                }
                else
                {
                    Candidate_Distribution_per[ik] = 0;
                }

            }
            else
            {
                Candidate_Distribution_per[ik] = 0;
            }


            if (  HSR >= 60  )/*Here HSR = homology score, 2009-07-03*/
            {
                Candidate_Distribution[ik] = 22;  /*Candidate_Distribution is the homology weight*/
            }
            else if (   (HSR<60) && (HSR>=50)  )
            {
                Candidate_Distribution[ik] = 20;
            }
            else if (  (HSR<50) && (HSR>=40)  )
            {
                Candidate_Distribution[ik] = 18;
            }
            else if (  (HSR<40) && (HSR>=30)  )
            {
                Candidate_Distribution[ik] = 16;
            }
            else if (  (HSR<30) && (HSR>=20)  )
            {
                Candidate_Distribution[ik] = 14;
            }
            else if (  (HSR<20) && (HSR>=15)  )
            {
                Candidate_Distribution[ik] = 12;
            }
            else if (  (HSR<15) && (HSR>=10)  )
            {
                Candidate_Distribution[ik] = 11;
            }
            else
            {
                Candidate_Distribution[ik] = 10;
            }

            //Candidate_Distribution[ik] = Candidate_Distribution[ik] + (Candidate_Distribution[ik]-10)/2;

            if ( HSR >= 11 )
            {
                Homology_number = Homology_number + Candidate_Distribution[ik] - 10;
            }

            //Candidate_Distribution[ik] = 10;
            //Candidate_Distribution[ik] = Candidate_Distribution[ik]*Candidate_Distribution[ik]/10;

            if  (   (HSR_frag==0) && ( (HSR>=8)||(Per_homo>=12) )   ) {
                //sprintf(OpenName,"%s/Homologue_list.txt",check_Result_location);
                OpenName = string(check_Result_location) + "/" +  "Homologue_list.txt";
                wppp = fopen(OpenName.c_str(),"a");
                fprintf(wppp,"%5s %5s  %3d %3d\n",CSubname,Namelistunit[ik],Per_homo,HSR);
                fclose(wppp);
            } 
        }

        //*********************************************
        NSUM = Nun_Cand;
        Seg_Wedth = 0;
        for (ik=0; ik<Len_target-NFragment; ik++)
        {

            for (ij=Nkeep; ij<NSUM; ij++)
            {
                if (  Candidate_Score[ik][ij] <= 1  )
                {
                    continue;
                }

                Sum_value = 900000;
                param1 = 0;
                for ( SFrag=-Seg_Wedth; SFrag<NFragment+Seg_Wedth; SFrag++)
                {
                    if (     (  (ik+SFrag) >= Len_target ) ||  (  (ik+SFrag) < 0 ) || (  (Candidate_Score_iklen[ik][ij]+SFrag) >= LengthList[Candidate_Score_iksub[ik][ij]] )  ||  (  (Candidate_Score_iklen[ik][ij]+SFrag) < 0 )  ||  (Target_Naas[ik+SFrag]>=20)    )
                    {
                        continue;
                    }

                    HSR = (Candidate_Score_iklen[ik][ij]+SFrag)*20;


                    value = 0;
                    for (im=0; im<20; im++)
                    { 
                        param1 = Per_ten_one[Psi_blosum[ik+SFrag][im]]*(Log_per[blosum_chain[Candidate_Score_iksub[ik][ij]][HSR+im]]-Log_Back_Comp[im])
                            + Per_ten_one[blosum_chain[Candidate_Score_iksub[ik][ij]][HSR+im]]*(Log_per[Psi_blosum[ik+SFrag][im]]-Log_Back_Comp[im]);
                        //param1 = Per_ten_one[Psi_blosum[ik+SFrag][im]]*(Log_per[blosum_chain[Candidate_Score_iksub[ik][ij]][HSR+im]]-back_ground[Candidate_Score_iksub[ik][ij]][im])
                        //	   + Per_ten_one[blosum_chain[Candidate_Score_iksub[ik][ij]][HSR+im]]*(Log_per[Psi_blosum[ik+SFrag][im]]-back_ground[Candidate_Score_iksub[ik][ij]][im]);
                        value = value + param1;
                    }
                    //temp1 = Target_Naas[ik+SFrag];//AAS
                    //value = float (value*log(2.5 + Psi_blosum[ik+SFrag][temp1]/1000.0));

                     /* value: is the profile-profile score 
                      * weight: is the information content (psudeo); 2009-07-14*/ 
                    Sum_value = Sum_value + value*30*weight[ik+SFrag];  
                }
                //Sum_value = Sum_value*Candidate_Distribution[ik]/10;
                Candidate_Score[ik][ij] = Sum_value;
            }


            //set the order 
            for (ij=Nkeep; ij<NSUM; ij++)
            {
                for (im=ij+1; im<NSUM; im++)
                {
                    if (  Candidate_Score[ik][im] > Candidate_Score[ik][ij] )
                    {
                        Sum_value = Candidate_Score[ik][im];
                        Candidate_Score[ik][im] = Candidate_Score[ik][ij];
                        Candidate_Score[ik][ij] = Sum_value;

                        HSR = Candidate_Score_iksub[ik][im];
                        Candidate_Score_iksub[ik][im] = Candidate_Score_iksub[ik][ij];
                        Candidate_Score_iksub[ik][ij] = HSR;

                        HSR = Candidate_Score_iklen[ik][im];
                        Candidate_Score_iklen[ik][im] = Candidate_Score_iklen[ik][ij];
                        Candidate_Score_iklen[ik][ij] = HSR;

                        HSR = Candidate_Score_Maxline[ik][im];
                        Candidate_Score_Maxline[ik][im] = Candidate_Score_Maxline[ik][ij];
                        Candidate_Score_Maxline[ik][ij] = HSR;
                    }
                }

            }

        }





        //*********************************************    
        if (  Homology_number >= 20  )
        {
            Homology_number = 19;
        }

        if (  Homology_number <= 2  )
        {
            Consens_Sample = 8;
        }
        else
        {
            Consens_Sample = 5;
        }
        //Consens_Sample = 5;




        for (ik=0; ik<LONGEST_SEQ; ik++)
        {
            Number_fam[ik] = 0;
            for (ij=0; ij<8; ij++)
            {
                Pred_Result[ik][ij] = 0;
            }
            for (ij=0; ij<3; ij++)
            {
                Pred_Resulthsr[ik][ij] = 0;
            }
        }

        for (ik=0; ik<Len_target-NFragment; ik++)
        {
            //extract result
            for (im=0; im<NFragment; im++)//for 9 positions of 9-long fragment
            {
                for (ij=0; ij<Consens_Sample; ij++)
                {
                    if (  Candidate_Score[ik][ij] <= 0  )
                    {
                        continue;
                    }
/*debug*/
//if (ik+im == 125)
//{
    //fprintf(stdout,"ik+im=%d\n",ik+im);
//}

/*debug*/

                    //shape of candidate
                    CCvalue = CShapesequence[Candidate_Score_iksub[ik][ij]][Candidate_Score_iklen[ik][ij]+im];
                    HSR = CCvalue - 'A';
                    if  (  (HSR>=0) && (HSR<=25)  )
                    {
                        HSR = Shape_Code[HSR];
                    }
                    else
                    {
                        HSR = 8;
                    }
                    if (  HSR <= 7 )
                    {
                        Pred_Result[ik+im][HSR] = Pred_Result[ik+im][HSR] + Candidate_Distribution[Candidate_Score_iksub[ik][ij]];
                    }
                    //HSR
                    CCvalue = CHSRsequence[Candidate_Score_iksub[ik][ij]][Candidate_Score_iklen[ik][ij]+im];
                    if  (  CCvalue == 'H'   )
                    {
                        HSR = 0;
                    }
                    else if (  CCvalue == 'S'  )
                    {
                        HSR = 1;
                    }
                    else if (  CCvalue == 'R')   /*code added 2009-07-03,*/
                    {
                        HSR = 2;
                    }
                    else /* == '-' or other chars*/
                    {
                        HSR = 3;
                    }
                    if (HSR < 3) /*code added 2009-07-03,*/ 
                    {

                        Pred_Resulthsr[ik+im][HSR] = Pred_Resulthsr[ik+im][HSR] + Candidate_Distribution[Candidate_Score_iksub[ik][ij]];
                        Number_fam[ik+im]++;
                    }
                }
            }

        }
        /*normalization*/
        for (ik=0; ik<Len_target; ik++)
        {
            Number_fam[ik] = ( (Pred_Resulthsr[ik][0] + Pred_Resulthsr[ik][1] + Pred_Resulthsr[ik][2]) - Number_fam[ik]*10 )/4;
            if (  Number_fam[ik] >= 20 ) 
            {
                Number_fam[ik] = 19;
            }

            for (ij=0; ij<8; ij++)
            {
                Pred_Result[ik][ij] = Pred_Result[ik][ij]/10;
            }
            for (ij=0; ij<3; ij++)
            {
                Pred_Resulthsr[ik][ij] = Pred_Resulthsr[ik][ij]/10;

            }
        }


        HSR_Cont_2darray.Init(0);
        //for (ik=0; ik<3; ik++)
        //{
            //for (ij=0; ij<2; ij++)
            //{
                //HSR_Cont[ik][ij] = 0;
            //}
        //}

        Sumvalue = 0;
        Sumvalue1 = 0;
        Sumnumber = 0;
        HSRsum1 = 0;
        HSRcorrect1 = 0;
        

        int numValiedPredResidue = 0; /*the number of valied predictions, that is, those predicted residue positions with secondary structure definition, 2009-07-07 */
        //sprintf(OpenName,"%s/Res_%s.txt",check_Result_location,CSubname);
        OpenName = string(check_Result_location) + "/Res_" + string(CSubname) + ".txt";
        wfp = fopen(OpenName.c_str(),"w");
        //shape and HSR
        for (ij = 0 ; ij <Len_target; ij ++ ){ /*2009-07-09, change 1800 to
                                                 Len_target, nanjiang*/
            NSUM = 0;
            for (im=0; im<8;im++) { 
                NSUM = NSUM + Pred_Result[ij][im];
                NPer[im] = Pred_Result[ij][im];
            }

            if ( isTreatZeroNSUM ) {
                if (NSUM == 0) {
                    NSUM = 1; /*to avoid divide_by_zero error, 2009-07-09*/
                }
            }

            if (  NSUM < 1  ) { /*ignore those without any matched fragment,
                                 2009-07-03, this code can be changed to that
                                 assign this state to 'R' and set the
                                 confidence to -1, yet to be changed
                                 2009-07-08: do not ignore those without any
                                 matched fragments.*/
                continue;
            }
            Number_class = (NSUM-1)/Consens_Sample;
            NTG_COM = NSUM/15;
            GGroup[0] = NPer[0] + NPer[1] + NPer[2] + NPer[3];
            GGroup[1] = NPer[4] + NPer[5];
            GGroup[2] = NPer[6] + NPer[7] + NTG_COM;

            //prediction
            //NNshape--observed shape, Numberofgroup--predicted shape
            //Sumnumber-- the number of AAS for predicted shape
            //Sumvalue-- the score sum for the Sumnumber AAS for predicted 8-shape
            //Sumvalue1-- the score sum for the Sumnumber AAS for predicted 3-shape

            Numberofgroup = 7 ;/*set the initial value, 2009-07-16 */
            //if (  Number_class >= (Halv_NFragment-0) ) [>igore the first 4 residues and the last 4 residues, 2009-07-03 <]
            if ( Number_class >= threshold_shift) /*do not ignore the beg end*/
            {
                if (  (GGroup[2]>=GGroup[0]) && (GGroup[2]>=GGroup[1])   )//GT
                {
                    if (  NPer[7] >= (NPer[6]-NTG_COM)  )
                    {
                        Numberofgroup = 7;
                    }
                    else
                    {
                        Numberofgroup = 6;
                    }
                    Weight_per = (GGroup[2]-NTG_COM)*100/NSUM;
                }
                else if (  (GGroup[0]>=GGroup[1]) && (GGroup[0]>=GGroup[2])   )//SRUV
                {
                    Numberofgroup = 0;
                    HSR = NPer[0];
                    for (im=0; im<4; im++)
                    {
                        SRUV[im] = NPer[im];
                    }
                    SRUV[2] = SRUV[2] + NTG_COM;//U
                    SRUV[3] = SRUV[3] + NTG_COM;//V
                    for (im=1; im<4; im++)
                    {   
                        if (  SRUV[im] >= HSR  )
                        {
                            HSR = SRUV[im];
                            Numberofgroup = im;
                        }
                    }
                    Weight_per = GGroup[0]*100/NSUM;
                }
                else// (  (GGroup[1]>GGroup[0]) && (GGroup[1]>GGroup[2])   ), K A
                {
                    if (  NPer[4] >= (NPer[5]-NTG_COM)  )
                    {
                        Numberofgroup = 4;
                    }
                    else
                    {
                        Numberofgroup = 5;
                    }
                    Weight_per = GGroup[1]*100/NSUM;
                }

                NNshape = CTempShapeseq[ij] - 'A';
                if  (  (NNshape>=0) && (NNshape<=25)  )
                {
                    NNshape = Shape_Code[NNshape];
                }
                else
                {
                    NNshape = 8;
                }

                Svalue = -1;
                Weight_per = 0;
                if (  NNshape < 8 )
                {
                    Svalue = 0;
                    if (  NNshape == Numberofgroup  )
                    {
                        Svalue = 1;
                    }
                    else
                    {
                        if (  (NNshape<=3) && (Numberofgroup<=3)  )
                        {
                            Svalue = 0.5;
                        }
                        else if  (  (NNshape>=6) && (Numberofgroup>=6)  )
                        {
                            Svalue = 1;
                        }
                        else if (  (NNshape>=4) && (NNshape<=5) && (Numberofgroup<=5) && (Numberofgroup>=4)  )
                        {
                            Svalue = 0.5;
                        }
                    }

                    Sumvalue = Sumvalue + Svalue;
                    if (  Svalue > 0.9 )//1
                    {
                        Sumvalue1 = Sumvalue1 + 1;
                    }
                    Sumnumber++;
                }


            }//if (  Number_class >= NFragment )
            else
            {
                Numberofgroup = 0;
                Weight_per = 0;
                Svalue = -1;
            }


            //HSR---
            //NAAS--predicted HSR, PHSR--predicted score
            //HSRsum1-- the number of AAS for predicted HSR
            //HSRcorrect1-- the number of correctly predicted for HSRsum1 AAS of the predicted HSR
            //numValiedPredResidue : [>the number of valied predictions, that is, those predicted residue positions with secondary structure definition, 2009-07-07 <]

            /*======  remove this continue, so that the program
             * will predict even on those residues without the
             * secondary structure definition =============== */

            if (isNotPredUnDefined)    /*2009-07-09 , Nanjiang*/
            {  
                if (  CTempHSRseq[ij] == '-'  )
                {
                    continue;
                }
            }

            NPer[10] = Pred_Resulthsr[ij][0];
            NPer[11] = Pred_Resulthsr[ij][1];
            NPer[12] = Pred_Resulthsr[ij][2];
            NTG_COM = (NPer[10]+NPer[11]+NPer[12])/15;
            NPer[11] = NPer[11] + NTG_COM;
            HSR = NPer[10] + NPer[11] + NPer[12];

            /*special treat of the beginning and end*/
            /*if all NPer[] for secondary structure are zero, predict as 'R'; 2009-07-12*//*{{{*/
            if (NPer[10] == 0 && NPer [11] == 0 && NPer[12] == 0) 
            {
                    NPer[10] += 0;
                    NPer[11] += 1;
                    NPer[12] += 9; /*just set a higher probability to the state "R"*/
            }

            /*when Number_class <= 1, at the beginning top 3 and last 3, increase the weight on 'R', in the middle, increase the weight on 'H' and 'S'*/
            if (Number_class <= 0 )
            {
                if (ij < NSpecialBeg || ij > Len_target-NSpecialEnd-1) /*for the first NSpecialBeg and last NSpecialEnd residues, treat specially*/
                {
                    int sft = min(ij, -ij+Len_target-1);
                    sft = min (sft, 1);
                    NPer[10] += 0;
                    NPer[11] += 1;
                    NPer[12] += 10-sft; /*just set a higher probability to the state "R"*/
                }
                else if (ij > 5 && ij < Len_target - 5)
                {
                    NPer[10] += 5;
                    NPer[11] += 3;
                    NPer[12] += 0; /*just set a higher probability to the state "R"*/
                }
            }
            else if (Number_class <= 1 )
            {
                if (ij < NSpecialBeg || ij > Len_target-NSpecialEnd-1) /*for the first NSpecialBeg and last NSpecialEnd residues, treat specially*/
                {
                    int sft = min(ij, -ij+Len_target-1);
                    sft = min (sft, 0);
                    NPer[10] += 1;
                    NPer[11] += 0;
                    NPer[12] += 8-sft; /*just set a higher probability to the state "R"*/
                }
                else if (ij > 5 && ij < Len_target - 5)
                {
                    NPer[10] += 5;
                    NPer[11] += 3;
                    NPer[12] += 0; /*just set a higher probability to the state "R"*/
                }
            }
            else if (Number_class <= 2 )
            {
                if (ij < NSpecialBeg || ij > Len_target-NSpecialEnd-1) /*for the first NSpecialBeg and last NSpecialEnd residues, treat specially*/
                {
                    int sft = min(ij, -ij+Len_target-1);
                    sft = min (sft, 1);
                    NPer[10] += 0;
                    NPer[11] += 0;
                    NPer[12] += 2-sft; /*just set a higher probability to the state "R"*/
                }
                else if (ij > 5 && ij < Len_target - 5)
                {
                    NPer[10] += 2;
                    NPer[11] += 2;
                    NPer[12] += 0; /*just set a higher probability to the state "R"*/
                }
            }
            else if (Number_class <= 3 )
            {
                if (ij < NSpecialBeg || ij > Len_target-NSpecialEnd-1) /*for the first NSpecialBeg and last NSpecialEnd residues, treat specially*/
                {
                    NPer[10] += 0;
                    NPer[11] += 0;
                    NPer[12] += 0; /*just set a higher probability to the state "R"*/
                }
                else if (ij > 5 && ij < Len_target - 5)
                {
                    NPer[10] += 2;
                    NPer[11] += 1;
                    NPer[12] += 0; /*just set a higher probability to the state "R"*/
                }
            }

            /*force the first 1 and last 1 to 'R', 2009-07-09, Nanjiang*/
            if (ij == 0 || ij == Len_target-1)
            {
                NPer[10] = 1;
                NPer[11] = 2;
                NPer[12] = 13; /*just set a higher probability to the state "R"*/
            }               /*}}}*/
            /* ======== special treat of the beginning and end finished, 2009-07-13 */

            //if (  Number_class >= (Halv_NFragment-0) )
            if (  Number_class >= threshold_shift)/*predict all positions, do not ignore the start and end*/

                //if (  Number_class >= 5 )
            {
                NAAS = GetHSRState(NPer[10], NPer[11], NPer[12], method_HSR); /*2009-07-13, method to determine the state of the secondary structure is changed, old: in the order of S -> H -> R, method(1) R -> S -> H */
                //if  (   (NPer[11]>=NPer[10]) && (NPer[11]>=NPer[12])   )
                //{
                    //NAAS = 1;

                //}
                //else if (   (NPer[10]>=NPer[11]) && (NPer[10]>=NPer[12])   )
                //{
                    //NAAS = 0;
                //}
                //else
                //{
                    //NAAS = 2;
                //}
            }
            else if (  Number_class <= 4 )
            {
                continue;
            }

            else
            {
                continue;
            }
            PHSR = 0; /*frag for the prediction: 0 for  wrong,  1 for correct*/
            if (  NAAS == Pred_Resulthsr[ij][3] )
            {
                HSRcorrect1++;   /*the number of correctly predicted amino acids for this target chain*/
                PHSR = 1;
            }
            if ( Pred_Resulthsr[ij][3] != 3 )
            {
                numValiedPredResidue ++;
            }

            HSRsum1++; /*the total number of predicted amino acid for this target chain*/

            //Pred_Resulthsr[ij][3]--observed HSR from DSSP
            HSR_Cont[Pred_Resulthsr[ij][3]][0]++;  /*observed, Error found! Illegal writting, for the original code Tuping sent to me Pred_Resulthsr[ij][3] can only be 0 1 or 2, and residues with undefined structure is not correctly analyzed. I've changed the bug, but now Pred_Resulthsr[ij][3] can be 3, and HSR_Cont is of the size [3][2], 2009-07-24*/
            HSR_Cont[NAAS][1]++;//predicted

#ifdef DEBUG_NUMCLASS
            if (Number_class <= debugNumberClass_threshold)
            {
                fprintf(fwpTmpDebugNumberClass,"%s %4d %c %c %c %3d %3d %3d %d numclass= %2d NSUM= %3d Len= %4d\n", CSubname, ij, CTempAASseq[ij], CNameHSR[Pred_Resulthsr[ij][3]], CNameHSR[NAAS], NPer[10], NPer[11], NPer[12], PHSR, Number_class , NSUM, Len_target);
            }

#endif

            fprintf(wfp, "%4d %c %c %c ",ij,CTempAASseq[ij],CTempShapeseq[ij],CNameHSR[Pred_Resulthsr[ij][3]]);
            fprintf(wfp, " %3d  %3d  %3d  %1d   ",NPer[10],NPer[11],NPer[12],PHSR);
            for (im=0; im<8; im++)
            {
                fprintf(wfp, "%3d ",NPer[im]);
            }
/*debug CShape[Numberofgroup]; 2009-07-22 */
            if (CShape[Numberofgroup] < 'A' || CShape[Numberofgroup] > 'Z')
            {
                fprintf(stderr,"%s ij=%d, Numberofgroup=%d\n", CSubname, ij, Numberofgroup);
            }

            fprintf(wfp, "  %c %4.1f %3d ",CShape[Numberofgroup],Svalue,Weight_per); /*here Numberofgroup must assign a initial value, otherwise, in case of invalid reading, CShape[Numberofgroup] might be '\n', as in the case of Res_2DXRB.txt, 2009-07-16*/
            fprintf(wfp, "numClass= %3d NSUM= %4d Conf= %.3f", Number_class, NSUM,GetHSRConfidence(NPer[10], NPer[11], NPer[12]) ); /*added 2009-07-22*/
            fprintf(wfp, "\n");
        }
        fclose(wfp);

        if (  HSRsum1 >= 1  )
        {
            Totalvalue = Totalvalue + Sumvalue;
            Totalvalue1 = Totalvalue1 + Sumvalue1;
            TotalNumber = TotalNumber + Sumnumber;
            HSRcorrect = HSRcorrect + HSRcorrect1;
            totalNumValiedPredResidue = totalNumValiedPredResidue + numValiedPredResidue; /*total number of valied predicted residue positions*/
            HSRsum= HSRsum + HSRsum1;
            //sprintf(OpenName,"%s/Test_resultper.txt",check_Result_location);
            OpenName = string(check_Result_location) + "/" +  "Test_resultper.txt";
            wfpper = fopen(OpenName.c_str(),"a");
            //fprintf(wfpper,"%4d %5s %6.1f %4d %6.2f %6.2f    %9.1f %7d %6.2f %6.2f      %4d %4d %6.2f  %8d %8d %6.2f\n",iktest,CSubname,Sumvalue,Sumnumber,Sumvalue*100/Sumnumber,Sumvalue1*100/Sumnumber,Totalvalue,TotalNumber,Totalvalue*100/TotalNumber,Totalvalue1*100/TotalNumber,
                    //HSRcorrect1,HSRsum1,HSRcorrect1*100.0/HSRsum1,HSRcorrect,HSRsum,HSRcorrect*100.0/HSRsum);

            fprintf(wfpper,"%4d %5s %6.1f %4d %6.2f %6.2f    %9.1f %7d %6.2f %6.2f      %4d %4d %6.2f  %8d %8d %6.2f\n",iktest,CSubname,Sumvalue,Sumnumber,Sumvalue*100/Sumnumber,Sumvalue1*100/Sumnumber,Totalvalue,TotalNumber,Totalvalue*100/TotalNumber,Totalvalue1*100/TotalNumber,
                    HSRcorrect1,numValiedPredResidue,HSRcorrect1*100.0/numValiedPredResidue,HSRcorrect,totalNumValiedPredResidue,HSRcorrect*100.0/totalNumValiedPredResidue);
            fclose(wfpper);

            HSR = HSR_Cont[0][0] + HSR_Cont[1][0] + HSR_Cont[2][0];
            NAAS = HSR_Cont[0][1] + HSR_Cont[1][1] + HSR_Cont[2][1];
            //sprintf(OpenName,"%s/detailed.txt",check_Result_location);
            OpenName = string(check_Result_location) + "/detailed.txt";
            wppp = fopen(OpenName.c_str(),"a");
            fprintf(wppp,"%4d %5s  %6.2f  %6.2f  %6.2f      %6.2f  %6.2f  %6.2f %4d  %4d  %6.2f\n",
                    iktest,CSubname,
                    HSR_Cont[0][0]*100.0/HSR, HSR_Cont[1][0]*100.0/HSR,
                    HSR_Cont[2][0]*100.0/HSR , HSR_Cont[0][1]*100.0/NAAS,
                    HSR_Cont[1][1]*100.0/NAAS, HSR_Cont[2][1]*100.0/NAAS,
                    HSRcorrect1,HSR , HSRcorrect1*100.0/HSRsum1);
            fclose(wppp);
        }

        for (ik=0; ik<Len_target; ik++)
        {
            delete [] Candidate_Score[ik];
            delete [] Candidate_Score_iksub[ik];
            delete [] Candidate_Score_iklen[ik];
            delete [] Candidate_Score_Maxline[ik];
        }
    }
    fclose(fptest);
#ifdef DEBUG_NUMCLASS
    fclose(fwpTmpDebugNumberClass);
#endif 

    /*Free memory*/
    /*1. Free memory of filenames*/

    delete [] check_list            ;
    delete [] check_Result_location ;
    delete [] sourcefile            ;
    delete [] location_qijmatrix    ;
    delete [] location_modmatrix    ;
    delete [] HRS_profile           ;
    delete [] database_list         ;
    delete [] database_modmatrix    ;
    delete [] database_qijmatrix    ;
    delete [] database_facc         ;

    delete [] LengthList;
    for (ik=0; ik<Sumunit; ik++)
    {
        delete [] Namelistunit[ik];
        delete [] CAASSequence[ik];
        delete [] CShapesequence[ik];
        delete [] CHSRsequence[ik];
        delete [] blosum_chain[ik];
        delete [] Chain_code_AAS[ik];
    }
    delete [] Namelistunit;
    delete [] CAASSequence;
    delete [] CShapesequence;
    delete [] CHSRsequence;
    delete [] blosum_chain;
    delete [] Chain_code_AAS;
    delete [] Candidate_Distribution;
    delete [] Candidate_Distribution_per;

    return EXIT_SUCCESS ;

}
