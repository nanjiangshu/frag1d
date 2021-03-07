#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include "homolgy.h"
#include "base.h"
#include "myfunc.h"
#include "mypro.h"
#include "array.h"

using namespace std;

//#include <windows.h>

/*valgrind checked 2009-07-13, no error and no memory leak, Nanjiang*/

#if defined(_Windows) || defined(__WINDOWS__) || \
	defined(__WIN32__) || defined(WIN32) || \
defined(__WINNT__) || defined(__NT__)
#   ifndef WINDOWS
#       define WINDOWS
#   endif
#endif

#ifdef DEBUG_NUMCLASS
int debugNumberClass_threshold = 2;/*print out the prediction result for those with Number_class <= debugNumberClass_threshold; 2009-07-11 , nanjiang*/
#endif

bool isNotUseCandidateWithUndefinedStructure = true; /**/
bool isNotPredUnDefined = false; /*whether not predict those without secondary structure definition, default = false, 2009-07-13*/
bool isTreatZeroNSUM = true; /*whether set the zero NSUM to 1, default = true*/
bool isUseHomologyInfo = true; /*whether use homology info, default = true*/
int threshold_shift = 0;/*ignore those with Number_class  < threshold_shift, default = 0*/ 
int NSpecialBeg = 3;  /*treat specially the first NSpecialBeg residues*/
int NSpecialEnd = 3;  /*treat specially the last NSpecialEnd residues*/
int method_HSR = 1; /*using method 1 to determine the predicted secondary structure from three probabilities*/
int dsspMapMethod = 1; /* added 2009-10-14, using the commom dssp8to3 mapping scheme, GHI->Helis, EB->sheet, other->coil*/
int8 typeProfile = 1; /*2009-11-08*/

string usage/*{{{*/="\n\
usage: checkresult [options]\n\
Description: second round 1D structure prediction\n\
Options:\n\
  -dbtype INT  (default: 1)\n\
        0: database stored as one id one file\n\
        1: database stored as dumped file\n\
  --trainlist LISTFILE  \n\
        Set the trainging list, format is #  num id # for each line\n\
  --modm-train   |-mtrn DIR  Set the modmatrix for the training set\n\
  --accfrag-train|-atrn DIR  Set the accfrag for the training set\n\
  --qij-train    |-qtrn DIR  Set the qij path for the training set\n\
  --outpath      |-rloc DIR  \n\
  --fragpath     |-sloc DIR  \n\
  --HSRprofile   |-ploc DIR  Set the path to the HSR profile \n\
  --round1-pred         DIR  Set the result of secondary structure prediction\n\
                             from the previous step\n\
  --idlist-test  |-ltst FILE Set the idlist for the testing set\n\
  --qij-test     |-qtst DIR  Set the qij matrix for the testing set\n\
  --modm-test    |-mtst DIR \n\
  --nper         |-nptr INT  \n\
        The percentage of accfrag in the training database (default 40)\n\
  --nperHSR      |-npte INT  \n\
        The percentage of hsrfrag in the testing set (default 30)\n\
  --ratioscheme 0|1|2|3      Set the ratio scheme \n\
  --beg INT     \n\
        Set the starting number of id to be run \n\
  --end INT  \n\
        Set the end number of id to be run, beg = 0, end = 1 will run only the\n\
        first one \n\
  --rb                       Read binary database\n\
  --not-use-undef            Whether not use candidate fragment with undefined\n\
                             structures\n\
  --tshift|--threshold-shift INT \n\
        Ignore the first and last N residues, default = 0\n\
  --specialhead INT          First N residues to treat specially, (default: 3)\n\
  --specialtail INT          Last N residues to treat specially, (default: 3)\n\
\n\
  --blosum      FILE        Set the blosum file, (default: blosum62_tu.txt)\n\
  --dsspmap 0|1             Sset the dssp mapping method, (default: 1 )\n\
                            scheme 0, GHI->Helix, E->Sheet, other->coil \n\
                            scheme 1, GHI->Helix, BE->Sheet, other->coil \n\
DEBUG_NUMCLASS\n\
  --thcls \n\
        Set the threshold for debug number class, so that the program will\n\
        print out result for those with Number_class <= N, default = 3\n\
\n\
Created 2009-07-01, updated 2011-10-17, Nanjiang Shu\n\
";/*}}}*/

void PrintHelp() {
    fprintf(stdout,"%s\n", usage.c_str());
}

int main(int argc, char** argv)/*{{{*/
{
    bool isNonOptionArg = false;
    if(argc < 2) {
        PrintHelp();
        return 0;
    }

    int dbtype=1; /*0: database stored as one file for each id, 
1: database stored as dumped files. 
added 2011-10-17*/

    char first_non_blank_char ; /*2009-07-15, Nanjiang, using first_non_blank_char to determine whether the matrix are valid data lines*/
    char OpenName[1000];
    FILE *fp, *fpread,*wfp,*wfpper, *fptest;
    int seqIndex = 0; /*the index of the amino acid in sequence*/
    char str_char1[20],Camino,Cshp, CSAS[10],CCvalue;
    int Sumunit,ik,ij,im,ip, NSeries, Max_Length, temp1, temp2;
    int HSR = 0;
    int NAAS = 0;
    int TPer00;
    //char *Namelistunit[6600],*CAASSequence[6600], *CShapesequence[6600],*CHSRsequence[6600]; 
    //int  *LengthList;  
    int NPer[21] = { 0 };
    int iksub,iklen, ScoreSum;
    int NNshape = 25;
    int GGroup[5], Numberofgroup, TotalNumber,Sumnumber, SCore_Sample;
    int Save_Sample, Consens_Sample;
    float Totalvalue, Svalue, Sumvalue,param1,param2, value, Back_comp[21];
    float Log_per[105],Per_ten_one[105],Log_Back_Comp[21];
    int totalNumValiedPredResidue = 0; /*the number of valied predictions, that is, those predicted residue positions with secondary structure definition, 2009-07-07 */
    float valuelog2,Sumvalue1, Totalvalue1, Zero_deviation;
    //int *Chain_code_AAS[6600],*blosum_chain[6600];
    //int  SUM_Component[LONGEST_SEQ];
     
    Array2D <int> Psi_blosum_2darray(LONGEST_SEQ, 20);
    Array2D <int> Psi_blosum_Frag_2darray(LONGEST_SEQ, 20);
    Array2D <int> HSR_confidence_2darray(3, LONGEST_SEQ);
    Array2D <int> Shape_confidence_2darray(8, LONGEST_SEQ);
    int ** Psi_blosum = Psi_blosum_2darray.array2D;
    int ** Psi_blosum_Frag = Psi_blosum_Frag_2darray.array2D;
    int **HSR_confidence = HSR_confidence_2darray.array2D;
    int **Shape_confidence = Shape_confidence_2darray.array2D;

    int *Target_Naas = new int[LONGEST_SEQ];

    //int Nsearch_beg;
    //int NSer_sub;
    int Weight_per, NTG_COM, NPer_Frag_Database;
    //int Number_chain;
    //candidate
    int **Candidate_Score_iksub = new int*[LONGEST_SEQ];
    int **Candidate_Score_iklen = new int*[LONGEST_SEQ];
    float **Candidate_Score = new float*[LONGEST_SEQ];

    float Sum_value;
    int Len_target, NFragment, Halv_NFragment,HSRsum1,HSRcorrect1,HSRsum,HSRcorrect;
    int Database_pertemp = 30; /*ratio to merge HSRprofile and modmProfile NPer[ij] = NPer[ij]*Database_pertemp + Psi_blosum_Frag[NSeries][ij]*(100-Database_pertemp); default = 30, 2009-07-12, Nanjiang */
    //     int AAS_Code[27], Shape_Code[27];
    int iktest;

    Array2D <int> Pred_Result_2darray(LONGEST_SEQ, 8);
    Array2D <int> Pred_Resulthsr_2darray(LONGEST_SEQ, 4);
    int **Pred_Result = Pred_Result_2darray.array2D;
    int **Pred_Resulthsr = Pred_Resulthsr_2darray.array2D;

    int NSUM,SRUV[8], PHSR, Nun_Cand, Num_Candhomo;

    char CNameHSR[] = "HSR-"; 
    char CShape[] = "SRUVKATG-";
    //     int Candidate_Distribution[6600], Candidate_Distribution_per[6600]; 
    int Blosum_table[405], SFrag, Number_class, *Pair_data;
    Array1D <float> weight_1darray(LONGEST_SEQ);
    Array1D <float> statis_1darray(LONGEST_SEQ);
    weight_1darray.Init(INIT_FLOAT);
    statis_1darray.Init(INIT_FLOAT);
    float *weight = weight_1darray.array1D;
    float *statis = statis_1darray.array1D;
    float Xmol[21];
    //float Xmolbak[21]; 
    //float Value_score;
    //int Number_fam_read;
    //int Number_super_read;
    //int Homology_Score;
    int Homology_number, Seg_Wedth,  Nkeep;
    float Homology_Scale;
    //int HSR_DIS[2][3];
    int OLDConfidence, Homology_Total;
    //int (*blast)[21];
    int Realhomology;
    int NonSCOP  = 0;

    //int *homo_chain_number[6600];
    //int Number_fam[6600]; 
    int Dataset_self = 1;   /*test dataset self*/
    //Database_pertemp = 30;
    //Database_pertemp = 40;
    //Dataset_self = 0;//test dataset self


    /*file names*/
    char *blosumFile            = NULL; /* 2009-07-10                             */
    char *CDataBase             = NULL; /* idlist for the training set            */
    char *CBlast_Location_mod   = NULL;
    char *CBlast_Location_Frag  = NULL;
    char *CBlast_Location_Qij   = NULL;
    char *CData_Result_location = NULL; /* output path for the result; 2009-07-12 */
    char *sourcefile            = NULL;
    char *First_round_HSR       = NULL;
    char *first_round_loc       = NULL;
    char *Ctestlist             = NULL;
    char *Ctestqij              = NULL;
    char *Ctestmod              = NULL;


    SCore_Sample = 100;
    Save_Sample = 100;
    Consens_Sample = 5;
    NPer_Frag_Database = 40;
    //Nsearch_beg = 0;

    //Number_chain = 5900;
    Nkeep = 0;

    bool isReadBinaryFile = false;
    int ratioScheme = 0; /*by default the ratioscheme = 0*/
    int begNum = 0;
    int endNum = 5000000;



    const char control_option[] = ""; //options which control the program, and does not take parameters

    int j = 0;
    int i = 1;
    while(i < argc)/*{{{*/
    {
        if(argv[i][0] == '-' && !isNonOptionArg) //options
        {
            isNonOptionArg = false;
            if(IsInCharSet(argv[i][1], control_option))//if argv[i][1] is in control_option, it might be used as -aqs
            {
                for(j = 1 ; j <int( strlen(argv[i])); j++)
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
            else if( (strcmp(argv[i], "-dbtype") == 0) || (strcmp(argv[i], "--dbtype")== 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, dbtype, true, 0, 1)) == -1)
                    return -1;
            } else if( (strcmp(argv[i],"--trainlist") == 0) ||
                    strcmp(argv[i],"-ltrn") == 0 ) {
                if( ( i = option_parser_filename(argc, argv, i, &CDataBase)) == -1){
                    return -1;
                }
            } else if( (strcmp(argv[i],"--modm-train") == 0) ||
                    strcmp(argv[i],"-mtrn")==0) {
                if( ( i = option_parser_filename(argc, argv, i,
                                &CBlast_Location_mod)) == -1){
                    return -1;
                }
            } else if( (strcmp(argv[i],"--accfrag-train") == 0) ||
                    strcmp(argv[i],"-atrn")==0) {
                if( ( i = option_parser_filename(argc, argv, i, &CBlast_Location_Frag)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--qij-train") == 0) || strcmp(argv[i],"-qtrn") == 0) 
            {
                if( ( i = option_parser_filename(argc, argv, i, &CBlast_Location_Qij)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--outpath") == 0) || strcmp(argv[i],"-rloc")==0) /*location of the result file*/
            {
                if( ( i = option_parser_filename(argc, argv, i, &CData_Result_location)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--fragpath") == 0)|| strcmp(argv[i],"-sloc")==0)   /*location of the sourcefile (frag file)*/
            {
                if( ( i = option_parser_filename(argc, argv, i, &sourcefile)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--HSRprofile") == 0) ||strcmp(argv[i], "-ploc") == 0)  
            {
                if( ( i = option_parser_filename(argc, argv, i, &First_round_HSR)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--round1-pred") == 0))  
            {
                if( ( i = option_parser_filename(argc, argv, i, &first_round_loc)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--idlist-test") == 0)|| strcmp(argv[i],"-ltst")==0)  
            {
                if( ( i = option_parser_filename(argc, argv, i, &Ctestlist)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--qij-test") == 0)|| strcmp(argv[i],"-qtst")==0)  
            {
                if( ( i = option_parser_filename(argc, argv, i, &Ctestqij)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--modm-test") == 0)|| strcmp(argv[i],"-mtst")==0)  
            {
                if( ( i = option_parser_filename(argc, argv, i, &Ctestmod)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--blosum") == 0)|| strcmp(argv[i],"-blosum")==0)  /*added 2009-07-28*/
            {
                if( ( i = option_parser_filename(argc, argv, i, &blosumFile)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "-nptr") == 0) || (strcmp(argv[i], "--nper")== 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, NPer_Frag_Database, true, 0, 100)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-npte") == 0) || (strcmp(argv[i], "--nperHSR") == 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, Database_pertemp, 0, 100)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--rb") == 0) || (strcmp(argv[i], "--read-binary")== 0))  
            {
                isReadBinaryFile = true;
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
            else if(  (strcmp(argv[i], "--not-use-undef")== 0))  
            {
                isNotUseCandidateWithUndefinedStructure = false;
                i++;
            }
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
            } else if( (strcmp(argv[i], "--dsspmap") == 0) )  {
                if( ( i = option_parser_numeric(argc, argv, i, dsspMapMethod,
                                true, 0, 1)) == -1){
                    return -1;
                }
            } else if (strcmp(argv[i], "--") == 0){/*next item is non option argument*/
                isNonOptionArg = true;
                i ++;
                continue;
            } else {
                fprintf(stderr,"Error! Invalid argument '%s'\n", argv[i]);
                return -1;
            }
        } else {
            i ++;
        }
    }/*}}}*/

    VerifyFolder(CData_Result_location);
    //modified confidence (y): y = -0.0133x2 + 2.763x - 45.053,  x-- real confidence

    int linesize;
    int maxline = 1000;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    int status_sscanf = 0;

    Nun_Cand = 100;
    Num_Candhomo = 100;

    NFragment = 9;

    Array2D <int> HSR_Cont_2darray(4,2); /*HSR_Cont[4][2], [ij][0] for observed sec, [ij][1] for predicted sec, [ij] for H, S ,R and '-', 2009-07-24*/
    Array2D <int> HSR_DIS_2darray(2,4); /*HSR_Cont[2][4], [][ij], ij for H, S, R and '-' */
    Array3D <int> HSR_HSR_3darray(2,4,4);  /*the second and the thrid iterators are for H, S, R and '-', 2009-07-27*/ 
    HSR_DIS_2darray.Init(0);
    HSR_Cont_2darray.Init(0);
    HSR_DIS_2darray.Init(0);
    int ** HSR_Cont = HSR_Cont_2darray.array2D;
    int ** HSR_DIS = HSR_DIS_2darray.array2D;
    int *** HSR_HSR = HSR_HSR_3darray.array3D;

    float vvmin, vvmax, Evolution;
    //float Target_comp[12][20]; 
    int Dist_conf[10][10]  ;
    //int MExcess;
    int FST_HSR[9][4][2][2]; /*the second iterator is for H, S, R and '-', 2009-07-24*/ 
    //int HSR_Cont[3][2];
    //int HSR_HSR[2][3][3];   
    int Shape_Shape[2][8][8];

    //int Two_position[2][21]; 
    //int Confid_Cand[30]; 
    int  Confi_AAS[103][2],Statis_AAS[253][2], Confi_Statis[253][2],Nconfd;
    //int  Confi_AASbak[103][2];
    //int Statis_AASbak[253][2]; 
    //int Confi_Statisbak[253][2];
    int Nconf_aas, Homology_Per_chain[20][2],Homology_Per_AAS[20][2];
    //int DSSP_Shape[10][10];
    FILE *wppp;
    int Stat_Conf_AAS[101][51][2];

    for (ik=0; ik<100; ik++)
    {
        for (ij=0; ij<50; ij++)
        {
            Stat_Conf_AAS[ik][ij][0] = 0;
            Stat_Conf_AAS[ik][ij][1] = 0;
        }
    }



    vvmin = 10000;
    vvmax = 0;
    //MExcess = 0;
    for (ik=0; ik<9; ik++) {
        for (ij=0; ij<4; ij++) {
            for (im=0; im<2; im++) {
                FST_HSR[ik][ij][im][0] = 0;
                FST_HSR[ik][ij][im][1] = 0;
            }
        }
    }

    for (ik=0; ik<10; ik++) {
        for (ij=0; ij<10; ij++) {
            Dist_conf[ik][ij] = 0;
            //DSSP_Shape[ik][ij] = 0;
        }
    }

    for (ik=0; ik<8; ik++) {
        for (ij=0; ij<8; ij++) {
            Shape_Shape[0][ik][ij] = 0;
            Shape_Shape[1][ik][ij] = 0;
        }
    }

    HSR_HSR_3darray.Init(0);
    //for (ik=0; ik<4; ik++)
    //{
    //for (ij=0; ij<3; ij++)
    //{
    //HSR_HSR[0][ik][ij] = 0;
    //HSR_HSR[1][ik][ij] = 0;
    //}
    //}

    for (ik=0; ik<102; ik++) {
        Confi_AAS[ik][0] = 0;
        Confi_AAS[ik][1] = 0;
    }
    for (ik=0; ik<252; ik++) {
        Statis_AAS[ik][0] = 0;
        Statis_AAS[ik][1] = 0;
        Confi_Statis[ik][0] = 0;
        Confi_Statis[ik][1] = 0;
        Statis_AAS[ik][0] = 0;
        Statis_AAS[ik][1] = 0;
        Confi_Statis[ik][0] = 0;
        Confi_Statis[ik][1] = 0;
    }

    NPer_Frag_Database = 0;
    TotalNumber = 0;
    Totalvalue = 0;
    Totalvalue1 = 0;
    HSRsum = 0;
    HSRcorrect = 0;
    Max_Length = LONGEST_SEQ;
    NTG_COM = 2;

    Halv_NFragment = NFragment/2;
    Zero_deviation = 0.075;

    Zero_deviation = 0.25;
    //Zero_deviation = 1;
    /* ==== read in blosum matrix in tuping's format ==== */
    fp = fopen(blosumFile, "r");
    checkfilestream(fp , blosumFile,"r", true);
    int cntRow = 0;
    while ((linesize = fgetline(fp, line, maxline)) != EOF) /*2009-07-28, fscanf changed, never use fscanf to read text file*/
    {



        status_sscanf = sscanf(line,
                "%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d", &Target_Naas[0],
                &Target_Naas[1], &Target_Naas[2], &Target_Naas[3],
                &Target_Naas[4], &Target_Naas[5], &Target_Naas[6],
                &Target_Naas[7], &Target_Naas[8], &Target_Naas[9],
                &Target_Naas[10], &Target_Naas[11], &Target_Naas[12],
                &Target_Naas[13], &Target_Naas[14], &Target_Naas[15],
                &Target_Naas[16], &Target_Naas[17], &Target_Naas[18],
                &Target_Naas[19] );
        if (status_sscanf != 20)
        {
            fprintf(stderr,"Error! Read Blosum file error, status_sscanf != 20\n");
            assert (status_sscanf == 20);
        }

        for (ik=0; ik<20; ik++)
        {
            Blosum_table[cntRow*20+ik] = Target_Naas[ik];
        }
        cntRow++;
    }
    fclose(fp);

    Blosum_table[400] = -10;
    Blosum_table[401] = -2;
    Blosum_table[402] = -2;
    Blosum_table[403] = -2;

    valuelog2 = float ( log(2) );
    for (ik=0; ik<103; ik++) {
        Per_ten_one[ik] = float ( ik + Zero_deviation ); 
        Log_per[ik] = float ( log(Per_ten_one[ik])/valuelog2 );
    }
    //conser_parameters.txt---search number: sample number:

    //fp = fopen(CParameter_conser_control,"r");
    //while (  fscanf(fp,"%d%d%d%d%d\n", &SCore_Sample, &Save_Sample, &Consens_Sample, &NPer_Frag_Database, &Nsearch_beg ) != EOF  )
    //{
    //}
    //fclose(fp);


    NTG_COM = Consens_Sample/15; /*add an extra 1/15 = 6.7% percent confidence to S state*/

    //Database_pertemp = NPer_Frag_Database;

    //background composition--Back_comp
    //fp = fopen(Back_composition,"r");
    //while (  fscanf(fp,"%d%f\n", &NAAS,&value ) != EOF  )
    //{
    //Back_comp[NAAS] = value;//percentage*10
    //}
    //fclose(fp);

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

    for (ik=0; ik<20; ik++) {
        Log_Back_Comp[ik] = float ( log(Back_comp[ik])/valuelog2 );//see Per_ten_one[ik]
    }

    //fp = fopen(CDataBase,"r");
    //Sumunit = 0;
    //while (  fgets(string,300,fp) != NULL  )
    //{
    //Sumunit++;
    //}
    //rewind(fp);
    //if (  Sumunit > Number_chain  )    
    //{
    //printf("too many files: %d  > %d",Sumunit , Number_chain);
    //return 0;
    //}
    //

    char *CTempAASseq = new char[LONGEST_SEQ+1];
    char *CTempShapeseq = new char[LONGEST_SEQ+1];
    char *CTempHSRseq = new char[LONGEST_SEQ+1];
    char *CTempDSSPseq = new char[LONGEST_SEQ+1];

    int database_Size = fgetlinecnt(CDataBase)+1;
    char **Namelistunit              = new (char(*[database_Size]) );
    char **CAASSequence              = new (char(*[database_Size]) );
    char **CShapesequence            = new (char(*[database_Size]) );
    char **CHSRsequence              = new (char(*[database_Size]) );
    int **Chain_code_AAS             = new (int (*[database_Size]) );
    int **blosum_chain               = new (int (*[database_Size]) );
    int * LengthList                 = new int [database_Size];
    int * Candidate_Distribution     = new int [database_Size];
    int * Candidate_Distribution_per = new int [database_Size];
    int * Homology_List              = new int[database_Size];

    /*read in the training set, 2009-07-12 */
    fprintf(stdout,"Read training set in %s format with ratioScheme = %d ...\n",
            isReadBinaryFile == true ? "binary" : "text", ratioScheme);
    Sumunit = Read_databse_SHU(dbtype, CDataBase, NPer_Frag_Database, ratioScheme,  CBlast_Location_Qij,  CBlast_Location_mod,CBlast_Location_Frag,  CData_Result_location,  Namelistunit, CAASSequence,CShapesequence, CHSRsequence, blosum_chain, Chain_code_AAS, LengthList, dsspMapMethod, isReadBinaryFile , typeProfile);
    fprintf(stdout,"Read training set finished, %d chains read in.\n", Sumunit);

    //Number_fam_read = 0;
    //Number_super_read = 0;

#ifdef DEBUG_NUMCLASS /*2009-07-13*/
    char tmpDebugNumClassFile[MAX_PATH+1] = "";
    sprintf(tmpDebugNumClassFile, "%s/resultNumberClassLE%d.dat", CData_Result_location, debugNumberClass_threshold);
    FILE *fwpTmpDebugNumberClass = fopen(tmpDebugNumClassFile, "w");
#endif 

    sprintf(OpenName,"%s/detailed.txt",CData_Result_location);
    wppp = fopen(OpenName,"w");
    fclose(wppp);

    Homology_Scale = 5000.0/Sumunit;//the weights given to the homologues

    for (ik=0; ik<20; ik++)
    {
        Homology_Per_chain[ik][0] = 0;
        Homology_Per_chain[ik][1] = 0;
        Homology_Per_AAS[ik][0] = 0;
        Homology_Per_AAS[ik][1] = 0;
    }


    char alphabet[50] = "";
    double parameter[8];
    int seqLength = 0;
    Array1D <ProfileSADByte> profileSADByte_1darray(LONGEST_SEQ);
    ProfileSADByte * profileSADByte = profileSADByte_1darray.array1D;

    iktest = -1;
    fptest = fopen(Ctestlist,"r");
    fprintf(stdout,"Start secondary structure prediction (step 2)...\n");
    while((linesize = fgetline(fptest, line ,maxline)) != EOF)// idList for test set
    {

        if (linesize <= 0) {
            continue;
        }
        char *CSubname = new char[linesize+1];
        sscanf(line, "%s",  CSubname);

        iktest++;

        if (iktest < begNum || iktest >= endNum)
        {
            continue;
        }


        for (ik=0; ik<LONGEST_SEQ; ik++)
        {
            for (ij=0; ij<20; ij++)
            {
                Psi_blosum[ik][ij] = 0;
            }
            Target_Naas[ik] = 20;
            Pred_Resulthsr[ik][3] = 2;
            CTempHSRseq[ik] = '-';
        }

        if (!isReadBinaryFile){/*{{{*/ 
            sprintf(OpenName,"%s/%s.Qij",Ctestqij,CSubname);
            fpread = fopen(OpenName,"r");
            if (  fpread == NULL  ) {
                continue;
            }
            while((linesize = fgetline(fpread, line ,maxline)) != EOF)
            {
                sscanf(line, " %c", &first_non_blank_char);
                if (first_non_blank_char <'0' || first_non_blank_char > '9') /*only when the first_non_blank_char is digit 2009-07-15*/
                {
                    continue;
                }
                sscanf(line," %d ", &NSeries );
                NSeries -= 1;
                //percent
                sscanf(line, "%d %c %s %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%f%f",
                        &seqIndex, &Camino, CSAS, &NPer[0], &NPer[1], &NPer[2],
                        &NPer[3], &NPer[4], &NPer[5], &NPer[6], &NPer[7],
                        &NPer[8], &NPer[9], &NPer[10], &NPer[11], &NPer[12],
                        &NPer[13], &NPer[14], &NPer[15], &NPer[16], &NPer[17],
                        &NPer[18], &NPer[19], &param1, &param2);
                statis[NSeries] = param2;
                Cshp = CSAS[0];
                NAAS = int(Camino - 'A');
                if (  (NAAS>=0) && (NAAS<=25)  ) {
                    NAAS = AAS_Code[NAAS];
                } else {
                    NAAS = 20;
                }
                Target_Naas[NSeries] = NAAS;
#ifdef DEBUG_CONDITIONAL_JUMP
                fprintf(stdout,"iktest= %d idtest= %s NSeries= %d Camino= %c NAAS=%d\n", iktest, CSubname, NSeries, Camino, NAAS);
#endif

                if (  NAAS >= 20  ) {
                    CTempAASseq[NSeries] = 'X';/*for non-standard amino acids,
                                                 set the amino acid type to 'X'
                                                 and set the shape string and
                                                 secondary structure to '-',
                                                 2009-07-15
                                                   */
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
                //            if (  (CSAS[2]=='H') || (CSAS[2]=='G')  )
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
                else if (  CSAS[2] == '-'     )
                {
                    temp1 = 3;
                    CCvalue = '-';
                }
                else
                {
                    temp1 = 2;
                    CCvalue = 'R';
                }
                CTempHSRseq[NSeries] = CCvalue;
                CTempDSSPseq[NSeries] = CSAS[2];
                Pred_Resulthsr[NSeries][3] = temp1;
            }
            fclose(fpread);
            NSeries++;
            Len_target = NSeries;
        }/*}}}*/
        else/*{{{*/
        {   /*Read In bindary Qijfile for the test set*/
            sprintf(OpenName,"%s/%s.Qijbin",Ctestqij,CSubname);
            if (GetBinaryMODM(OpenName, alphabet, seqLength, profileSADByte,
                        parameter, typeProfile) == -1) {
                fprintf(stderr, "can not open QijFile %s\n", OpenName);
                return -1;
            }
            for (ik = 0; ik<seqLength; ik++)
            {
                Camino = profileSADByte[ik].aa;
                statis[ik] = profileSADByte[ik].score2;
                NAAS = Camino - 'A';
                if (  (NAAS>=0) && (NAAS<=25)  ) {
                    NAAS = AAS_Code[NAAS];
                } else {
                    NAAS = 20;
                }
                Target_Naas[ik] = NAAS;

                if (  NAAS >= 20  ) {
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

        //read the modmatrix
        if (!isReadBinaryFile)/*{{{*/
        {
            sprintf(OpenName,"%s/%s.modm",Ctestmod,CSubname);
            fpread = fopen(OpenName,"r");
            for (ik=0; ik<LONGEST_SEQ; ik++)
            {
                for (ij=0; ij<20; ij++)
                {
                    Psi_blosum_Frag[ik][ij] = 0;
                }
            }
            if (  fpread != NULL  )
            {
                //while (  fgets(line,300,fpread) != NULL  )
                while((linesize = fgetline(fpread, line ,maxline)) != EOF)
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
                    //NSeries = atoi(str_char1);
                    //if (   (NSeries<1) || (NSeries>=Max_Length)   )//for the first two lines
                    //{
                    //continue;
                    //}
                    NSeries = NSeries - 1;
                    //percent
                    sscanf(line, "%d %c %s %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%f%f", 
                            &seqIndex, &Camino, CSAS, &NPer[0], &NPer[1],
                            &NPer[2], &NPer[3], &NPer[4], &NPer[5], &NPer[6],
                            &NPer[7], &NPer[8], &NPer[9], &NPer[10], &NPer[11],
                            &NPer[12], &NPer[13], &NPer[14], &NPer[15],
                            &NPer[16], &NPer[17], &NPer[18], &NPer[19],
                            &param1, &param2);
                    for (ij=0; ij<20; ij++) {
                        Psi_blosum_Frag[NSeries][ij] =  NPer[ij];
                    }
                }
                fclose(fpread);
            }

        }/*}}}*/
        else/*{{{*/
        {   /*Read In bindary Qijfile for the test set*/
            sprintf(OpenName,"%s/%s.modmbin",Ctestmod,CSubname);
            if (GetBinaryMODM(OpenName, alphabet, seqLength, profileSADByte, parameter, typeProfile) == -1)
            {
                fprintf(stderr, "can not open QijFile %s\n", OpenName);
                return -1;
            }
            for (ij=0; ij<20; ij++)
            {     
                Psi_blosum_Frag[ik][ij] = profileSADByte[ik].p[ij];
            } 
        }/*}}}*/

        //read from aacfrag_HSR
        sprintf(OpenName,"%s/Frag_%s.txt",First_round_HSR,CSubname);
        //sprintf(OpenName,"%s%s.txt",CBlast_Location_Frag,CSubname);
        fpread = fopen(OpenName,"r");
        if (  fpread != NULL  )
        {

            //while (  fgets(line,300,fpread) != NULL  )
            while((linesize = fgetline(fpread, line ,maxline)) != EOF)
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
                //NSeries = atoi(str_char1);
                //if (   (NSeries<1) || (NSeries>=Max_Length)   )//for the first two lines  
                //{
                //continue;
                //}
                NSeries = NSeries - 1;
                //percent
                sscanf(line,"%d %c %s %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%f%f",
                        &seqIndex, &Camino, CSAS, &NPer[0], &NPer[1], &NPer[2],
                        &NPer[3], &NPer[4], &NPer[5], &NPer[6], &NPer[7],
                        &NPer[8], &NPer[9], &NPer[10], &NPer[11], &NPer[12],
                        &NPer[13], &NPer[14], &NPer[15], &NPer[16], &NPer[17],
                        &NPer[18], &NPer[19], &param1, &param2);

                int sumProfile = 0;
                for (ij=0; ij<20; ij++)
                {
                    NPer[ij] = NPer[ij]*Database_pertemp + Psi_blosum_Frag[NSeries][ij]*(100-Database_pertemp);
                    //NPer[ij] = NPer[ij]*Database_pertemp + Psi_blosum[NSeries][ij]*(100-Database_pertemp);
                    sumProfile = sumProfile + NPer[ij];
                }
#ifdef DEBUG_NPTR /*debug the merge ratio of modm and HSRprofile, 2009-07-24*/
                if(Database_pertemp != 50)
                {
                    fprintf(stderr,"Database_pertemp changed by illegal write. Database_pertemp = %d ID=%s NSeries=%d\n", Database_pertemp, CSubname, NSeries) ;
                }
#endif
                for (ij=0; ij<20; ij++)
                {    
                    Psi_blosum[NSeries][ij] = int ( NPer[ij]*100.0/sumProfile + 0.4999 );
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


        //read from the first round result
        sprintf(OpenName,"%s/Res_%s.txt",first_round_loc,CSubname);

        /*2009-07-13, store the step 1 predicted result*/
        int maxSizeLinePred = 100;
        fgetlinecnt(OpenName, maxSizeLinePred);
        Array2D <char> linesStep1Pred_2darray(Len_target, maxSizeLinePred+1);
        char ** linesStep1Pred = linesStep1Pred_2darray.array2D;

        fpread = fopen(OpenName,"r");
        if (  fpread == NULL  )
        {
            continue;
        }
        for (ik=0; ik<LONGEST_SEQ; ik++)
        {
            HSR_confidence[0][ik] = 40;
            HSR_confidence[1][ik] = 20;
            HSR_confidence[2][ik] = 40;
            for (ij=0; ij<8; ij++)
            {
                Shape_confidence[ij][ik] = 0;
            }

        }
        OLDConfidence = 0;
        HSRsum1 = 0;
        //while (  fgets(line,300,fpread) != NULL  )
        while((linesize = fgetline(fpread, line ,maxline)) != EOF)
        {
            sscanf(line,"%d %c %c %c %d %d %d %d %d %d %d %d %d %d %d %d %c %f %d",
                    &NAAS,&str_char1[0],&str_char1[1], &str_char1[2],&NPer[0],&NPer[1],&NPer[2], &NPer[3],
                    &NPer[11],&NPer[12],&NPer[13],&NPer[14],&NPer[15],&NPer[16],&NPer[17],&NPer[18], &str_char1[4],&param2, &NPer[19]);

            int tmpn = 0;
            tmpn = strlen(line);
            strncpy(linesStep1Pred[NAAS], line, tmpn);/*store the step 1 predicted result, 2009-07-13, Nanjiang*/
            linesStep1Pred[NAAS][tmpn] = '\0';

            if  (  (NPer[1]>=NPer[0]) && (NPer[1]>=NPer[2])   )
            {
                temp2 = 1;
            }
            else if (  (NPer[0]>=NPer[1]) && (NPer[0]>=NPer[2])   )
            {
                temp2 = 0;
            }
            else
            {
                temp2 = 2;
            }
            HSR = NPer[0] + NPer[1] + NPer[2];
            for (ij=0; ij<3; ij++)
            { 
                HSR_confidence[ij][NAAS] = int ( NPer[ij]*100.0/HSR + 0.499 );
            } 
            OLDConfidence = OLDConfidence + HSR_confidence[temp2][NAAS];
            HSRsum1++;
            HSR = 0;
            for (ij=0; ij<8; ij++)
            {
                HSR = HSR + NPer[ij+11];
            }

            for (ij=0; ij<8; ij++)
            {
                Shape_confidence[ij][NAAS] = int (NPer[ij+11]*100.0/HSR + 0.499);
            }

        }
        fclose(fpread);
        if (  HSRsum1 >= 1  )
        {
            OLDConfidence = OLDConfidence/HSRsum1;
        }



        Confi_AAS[101][0] = 0;
        Confi_AAS[101][1] = 0;
        Statis_AAS[251][0] = 0;
        Statis_AAS[251][1] = 0;
        Confi_Statis[251][0] = 0;
        Confi_Statis[251][1] = 0;



        //sequence weight
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
            if (  HSR <50  )
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
                    value = value + Xmol[im]*log(Xmol[im]);
                }
            }
            weight[ik] = (-value+1)*(-value+1);
        }




        for (ik=0; ik<Len_target; ik++)
        {
            Candidate_Score[ik] = new float [SCore_Sample+5];
            Candidate_Score_iksub[ik] = new int [SCore_Sample+5];
            Candidate_Score_iklen[ik] = new int [SCore_Sample+5];
        }

        for (ik=0; ik<Len_target; ik++)
        {
            for (ij=0; ij<SCore_Sample+1; ij++)
            {
                Candidate_Score[ik][ij] = 0;
            }
        }
        for (ik=0; ik<database_Size; ik++)  /*2009-07-13, bug fixed, when initializing the array, the index should not exceed the alllocated size. Better to use the Init fuction of the template Array1D to initializing the array. Nanjiang*/
        {
            Candidate_Distribution[ik] = 0;
            Candidate_Distribution_per[ik] = 0;

        }

        //open the result file, changed 2009-07-12, by Nanjiang
        sprintf(OpenName,"%s/%s.txtbin",sourcefile,CSubname);
        /* Candidate_Score [m][100]
         * Candidate_Score_iksub [m][100]: the index of chain
         * Candidate_Score_iklen [m][100]: the position of the candidate fragment*/
        ReadInPairBinary( OpenName,  Candidate_Score, Namelistunit, Sumunit, Candidate_Score_iksub, Candidate_Score_iklen);
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

        // OLD-CODE open the result file
        //sprintf(OpenName,"%s/%s.txt",sourcefile,CSubname);
        //fpread = fopen(OpenName,"r");
        //if (  fpread == NULL  )
        //{
        //continue;
        //}
        //Pair_data,*
        //while (  fscanf(fpread,"%d%d%s%d%d%d\n", &ik,&ij,CSAS,&iksub, &iklen, &ScoreSum ) != EOF  )
        //{
        //Candidate_Score[ik][ij] = ScoreSum;
        //Candidate_Score_iksub[ik][ij] = iksub;
        //Candidate_Score_iklen[ik][ij] = iklen;
        //if (  ij < Num_Candhomo  )
        //{
        //Candidate_Distribution[iksub]++;
        //}
        //if (  iksub >= Number_chain)
        //{
        //printf("%s %d >= %d\n",CSubname,iksub,Number_chain);
        //return 1;
        //}    

        //}
        //fclose(fpread);


        Seg_Wedth = 0;

        //         if (  Dataset_self == 1  )
        //         {
        //             NonSCOP = Number_fam[iktest];
        //         }
        //         else
        //         {
        //             NonSCOP = 1;
        //         }
        Homology_Total = 0;
        Homology_number = 0;
        Realhomology = 0;



        //
        //for (ik=0; ik<Number_chain; ik++)
        for (ik=0; ik<Sumunit; ik++)  /*bug fixed 2009-07-12,Nanjiang, Number_chain will cause invalid read of memory, when ik >= Sumunit*/
        {
            HSR = int ( Candidate_Distribution[ik]*2*100/(Len_target+LengthList[ik])/Homology_Scale );

#ifdef DEBUG_CONDITIONAL_JUMP
            fprintf(stdout,"iktest= %d idtest= %s ik= %d Candidate_Distribution[ik]= %d LengthList[ik]= %d HSR=%d\n", iktest, CSubname, ik,  Candidate_Distribution[ik], LengthList[ik],  HSR);
#endif
            if (  HSR >= 10 )
            {
                ScoreSum = Len_target*LengthList[ik]+1;
                Pair_data = new int [ScoreSum];
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
                        if (  (Candidate_Score_iksub[im][ip]==ik) && (Candidate_Score_iklen[im][ip]>=0) && (Candidate_Score[im][ip]>=1)   )
                        {
                            NAAS = im*LengthList[ik] + Candidate_Score_iklen[im][ip];
                            assert (  NAAS <  (ScoreSum-1)  );
                            Pair_data[NAAS] = 1000;
                        }
                    }
                }

                NAAS = Treat_pair_Percent_homology(Pair_data, Len_target, LengthList[ik]);
                HSR = int ( NAAS*2*100/(Len_target+LengthList[ik])/Homology_Scale );

                delete [] Pair_data;

                if (  HSR >= 10 )
                {
                    Candidate_Distribution_per[ik] = HSR;
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


#ifdef DEBUG_CONDITIONAL_JUMP
            fprintf(stdout,"iktest= %d idtest= %s ik= %d  LengthList[ik]= %d HSR=%d\n", iktest, CSubname, ik, LengthList[ik], HSR);
#endif
            if (  HSR >= 60  ) {
                Candidate_Distribution[ik] = 22;
            } else if (   (HSR<60) && (HSR>=50)  ) {
                Candidate_Distribution[ik] = 20;
            } else if (  (HSR<50) && (HSR>=40)  ) {
                Candidate_Distribution[ik] = 18;
            } else if (  (HSR<40) && (HSR>=30)  ) {
                Candidate_Distribution[ik] = 16;
            } else if (  (HSR<30) && (HSR>=20)  ) {
                Candidate_Distribution[ik] = 14;
            } else if (  (HSR<20) && (HSR>=15)  ) {
                Candidate_Distribution[ik] = 12;
            } else if (  (HSR<15) && (HSR>=11)  ) {
                Candidate_Distribution[ik] = 11;
            } else {
                Candidate_Distribution[ik] = 10;
            }

            //Candidate_Distribution[ik] = Candidate_Distribution[ik] + (Candidate_Distribution[ik]-10)/2;


            if ( HSR >= 12 ) {
                Homology_number = Homology_number + Candidate_Distribution[ik] - 10;
            }

            if ( HSR >= 15 )
            {
                Homology_List[ik] = 1;
            }
            //Candidate_Distribution[ik] = 10;
        }

        Homology_Total = Homology_number;
        Homology_number = Homology_number/4;
        if (  Homology_number >= 20  )
        {
            Homology_number = 19;
        }
        else if (  Homology_number < 0  )
        {
            Homology_number = 0;
        }

        if (  Homology_number < 2 )
        {
            //Seg_Wedth = 1;
            Seg_Wedth = 0;
        }
        else
        {
            Seg_Wedth = 0;
        }

        if (  Homology_number <= 2  )
        {
            Consens_Sample = 8;
        }
        else
        {
            Consens_Sample = 5;
        }













        //*********************************************
        NSUM = Nun_Cand;
        for (ik=0; ik<Len_target-NFragment; ik++)
        {

            NSUM = Nun_Cand;
            for (ij=Nkeep; ij<NSUM; ij++)
            {
                if (  Candidate_Score[ik][ij] <= 1  )
                {
                    continue;
                }
                iksub = Candidate_Score_iksub[ik][ij];
                iklen = Candidate_Score_iklen[ik][ij];


                if (isNotUseCandidateWithUndefinedStructure) /*by default, do not use candidate with undefined structures (shape or HSR)*/
                {
                    for ( SFrag=0; SFrag<NFragment; SFrag++)
                    {
                        if (   (CTempHSRseq[ik+SFrag]=='-')  &&  (CHSRsequence[iksub][iklen+SFrag]=='-')   )
                        {
                            Candidate_Score[ik][ij] = 0;
                            break;
                        }
                    }
                }

                if (  Candidate_Score[ik][ij] <= 1  )
                {
                    continue;
                }


                for ( SFrag=0; SFrag<NFragment; SFrag++)
                {
                    CCvalue = CHSRsequence[iksub][iklen+SFrag];
                    if (  CCvalue == 'H'  )
                    { 
                        temp1 = 0;
                    } 
                    else if (  CCvalue == 'S'  )
                    { 
                        temp1 = 1;
                    } 
                    else
                    { 
                        temp1 = 2;
                    } 
                    NPer[SFrag] = temp1;

                    NPer[10+SFrag] = 8;
                    CCvalue = CShapesequence[Candidate_Score_iksub[ik][ij]][Candidate_Score_iklen[ik][ij]+SFrag];
                    temp1 = CCvalue - 'A';
                    if  ( (temp1<=25) && (temp1>=0)  )
                    { 
                        temp1 = Shape_Code[temp1];
                        NPer[10+SFrag] = temp1;
                    } 
                    else
                    {
                        NPer[10+SFrag] = 8;
                    }
                } 



                Sum_value = 900000;
                param1 = 0;
                param2 = Candidate_Distribution[Candidate_Score_iksub[ik][ij]]/10.0;;
                for ( SFrag=-Seg_Wedth; SFrag<NFragment+Seg_Wedth; SFrag++)
                {
                    if (     (  (ik+SFrag) >= Len_target ) ||  (  (ik+SFrag) < 0 ) || (  (Candidate_Score_iklen[ik][ij]+SFrag) >= LengthList[Candidate_Score_iksub[ik][ij]] )  ||  (  (Candidate_Score_iklen[ik][ij]+SFrag) < 0 )      )
                    {
                        continue;
                    }

                    HSR = (Candidate_Score_iklen[ik][ij]+SFrag)*20;
                    vvmin = 0;
                    if (  NPer[10+SFrag] < 8  )
                    {
                        vvmin = Shape_confidence[NPer[10+SFrag]][ik+SFrag];
                    }

                    value = 0;
                    for (im=0; im<20; im++)
                    { 
                        param1 = Per_ten_one[Psi_blosum[ik+SFrag][im]]*(Log_per[blosum_chain[Candidate_Score_iksub[ik][ij]][HSR+im]]-Log_Back_Comp[im])
                            + Per_ten_one[blosum_chain[Candidate_Score_iksub[ik][ij]][HSR+im]]*(Log_per[Psi_blosum[ik+SFrag][im]]-Log_Back_Comp[im]);
                        value = value + param1;
                    }

                    if (  value > 0  ) {
                        value = value*30*weight[ik+SFrag]*(100+HSR_confidence[NPer[SFrag]][ik+SFrag]+vvmin/2)/10;//secondary
                        //value = value*30*weight[ik+SFrag]*(100+Shape_confidence[NPer[10+SFrag]][ik+SFrag])/10;//shape string
                    } else {
                        value = value*30*weight[ik+SFrag]*(100-HSR_confidence[NPer[SFrag]][ik+SFrag])/100;//secondary
                        //value = value*30*weight[ik+SFrag]*(100-Shape_confidence[NPer[10+SFrag]][ik+SFrag])/100;//shape string
                    }
                    //Sum_value = Sum_value + value*30*weight[ik+SFrag];
                    Sum_value = Sum_value + value;
                }
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
                    }
                }
            }



            //NSUM = 15;
            NSUM = 30;
            for (ij=Nkeep; ij<NSUM; ij++)
            {
                if (  Candidate_Score[ik][ij] <= 1  ) {
                    continue;
                }
                //HSR for candidate
                for ( SFrag=0; SFrag<NFragment; SFrag++)
                {
                    CCvalue = CHSRsequence[Candidate_Score_iksub[ik][ij]][Candidate_Score_iklen[ik][ij]+SFrag];
                    if (  CCvalue == 'H'  ) {
                        temp1 = 0;
                    } else if (  CCvalue == 'S'  ) {
                        temp1 = 1;
                    } else {
                        temp1 = 2;
                    }
                    NPer[SFrag] = temp1;


                    NPer[10+SFrag] = 8;
                    CCvalue = CShapesequence[Candidate_Score_iksub[ik][ij]][Candidate_Score_iklen[ik][ij]+SFrag];
                    temp1 = CCvalue - 'A';
                    if (  (temp1<=25) && (temp1>=0)  )
                    {
                        temp1 = Shape_Code[temp1];
                        if (  (temp1<8) && (temp1>=0)  ) {
                            NPer[10+SFrag] = temp1;
                        } else {
                            NPer[10+SFrag] = 8;
                        }
                    } else {
                        NPer[10+SFrag] = 8;
                    }
                }

                vvmin = 0;
                vvmax = 1;
                for ( SFrag=0; SFrag<NFragment; SFrag++)
                {
                    temp1 = ik+SFrag;
                    if ( (temp1 >= Len_target) || (temp1<0) || 
                         ((Candidate_Score_iklen[ik][ij]+SFrag) >=
                          LengthList[Candidate_Score_iksub[ik][ij]] )  ||
                            (  (Candidate_Score_iklen[ik][ij]+SFrag) < 0 )
                       ) {
                        Candidate_Score[ik][ij] = 0;
                        continue;
                    }
                    if  (   (NPer[10+SFrag]<8) && (NPer[10+SFrag]>=0)   ) {
                        temp2 = NPer[10+SFrag];
                        // Halv point
                        if (  temp2 == 0 ){/*S*/
                            param1 = Shape_confidence[0][ik+SFrag] + 0.5*(Shape_confidence[1][ik+SFrag]+Shape_confidence[2][ik+SFrag]+Shape_confidence[3][ik+SFrag]);
                        } else if (  temp2 == 1 ){//R 
                            param1 = Shape_confidence[1][ik+SFrag] + 0.5*(Shape_confidence[0][ik+SFrag]+Shape_confidence[2][ik+SFrag]+Shape_confidence[3][ik+SFrag]);
                        } else if (  temp2 == 2 ){//U 
                            param1 = Shape_confidence[2][ik+SFrag] + 0.5*(Shape_confidence[1][ik+SFrag]+Shape_confidence[0][ik+SFrag]+Shape_confidence[3][ik+SFrag]);
                        } else if (  temp2 == 3 ){//V 
                            param1 = Shape_confidence[3][ik+SFrag] + 0.5*(Shape_confidence[1][ik+SFrag]+Shape_confidence[2][ik+SFrag]+Shape_confidence[0][ik+SFrag]);
                        } else if (  temp2 == 4 ){//K 
                            param1 = Shape_confidence[4][ik+SFrag] + 0.5*Shape_confidence[5][ik+SFrag];
                        } else if (  temp2 == 5 ){//A 
                            param1 = Shape_confidence[5][ik+SFrag] + 0.5*Shape_confidence[4][ik+SFrag];
                        } else if (  temp2 == 6 ){//T 
                            param1 = Shape_confidence[6][ik+SFrag] + Shape_confidence[7][ik+SFrag];
                        } else if (  temp2 == 7 ){//G 
                            param1 = Shape_confidence[6][ik+SFrag] + Shape_confidence[7][ik+SFrag];
                        } else {
                            param1 = 0;
                        }
                    }

                    NAAS = Target_Naas[temp1];
                    if  (   (NAAS>=20) || (NAAS<0)   )
                    {
                        Candidate_Score[ik][ij] = 0;
                        continue;
                    }
                    //vvmax = vvmax*(100+HSR_confidence[NPer[SFrag]][ik+SFrag]*1.0+param1*0.20)/100;//secondary
                    vvmax = vvmax*(100+HSR_confidence[NPer[SFrag]][ik+SFrag])/100;//secondary----higher
                    //vvmax = vvmax*(100+HSR_confidence[NPer[SFrag]][ik+SFrag]*0.20+param1*1.00)/100;//shape string
                    //vvmax = vvmax*(100+param1)/100;//shape string
                }
                Candidate_Score[ik][ij] = vvmax;


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
                    }
                }
            }


        }



        HSR_DIS_2darray.Init(0);
        //for (ik=0; ik<2; ik++)
        //{
        //for (ij=0; ij<3; ij++)
        //{
        //HSR_DIS[ik][ij] = 0;
        //}
        //}

        for (ik=0; ik<LONGEST_SEQ; ik++)
        {
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
                    else if (  CCvalue == '-'  )
                    {   
                        HSR = 3;
                    }
                    else
                    {
                        HSR = 2;
                    }
                    if (  HSR < 3  )
                    {
                        Pred_Resulthsr[ik+im][HSR] = Pred_Resulthsr[ik+im][HSR] + Candidate_Distribution[Candidate_Score_iksub[ik][ij]];
                    }
                }
            }
        }


        //normalization
        for (ik=0; ik<LONGEST_SEQ; ik++)
        {
            for (ij=0; ij<8; ij++)
            {
                Pred_Result[ik][ij] = Pred_Result[ik][ij]/10;
            }
            for (ij=0; ij<3; ij++)
            {
                Pred_Resulthsr[ik][ij] = Pred_Resulthsr[ik][ij]/10;

            }
        }


        for (ik=0; ik<3; ik++)
        {
            for (ij=0; ij<2; ij++)
            {
                HSR_Cont[ik][ij] = 0;
            }
        }


        Sumvalue = 0;
        Sumvalue1 = 0;
        Sumnumber = 0;
        HSRsum1 = 0;
        HSRcorrect1 = 0;
        Nconfd = 0;
        Evolution = 0;
        int numValiedPredResidue = 0; /*the number of valied predictions, that is, those predicted residue positions with secondary structure definition, 2009-07-13 */

        sprintf(OpenName,"%s/Res_%s.txt",CData_Result_location,CSubname);
        wfp = fopen(OpenName,"w");
        //shape and HSR

        for (ij=0; ij<Len_target;ij++) { /*ij iterating the target sequence*/
            if (isNotPredUnDefined) /*by default, predict even on those
                                      residues without structure definition,
                                      2009-07-13, Nanjiang it seems for
                                      checkresult, predicting on those regions
                                      worse the result*/
            {
                //if((  CTempHSRseq[ij] == '-'  )   ||  (CTempShapeseq[ij]=='-')    )
                if( (  CTempHSRseq[ij] == '-'  )  ) {
                    continue;
                }
            }

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
            if (  NSUM < 1  ) /*when isTreatZeroNSUM = true, no residues will be ignored, this code if NSUM <= 1 may cause the neglection of the last or first residue in the sequence*/
            {
                continue;
            }
            Number_class = (NSUM-1)/Consens_Sample;
            NTG_COM = NSUM/15;
            GGroup[0] = NPer[0] + NPer[1] + NPer[2] + NPer[3];
            GGroup[1] = NPer[4] + NPer[5];
            GGroup[2] = NPer[6] + NPer[7] + NTG_COM;

            //prediction
            //if (  Number_class >= (Halv_NFragment) )
            if ( Number_class >= threshold_shift) /*set the threshold_shift, when threshold_shift == 0, do not ignore those with Number_class <=4, which are mostly located at the beginning and end, and those residues without secondary structure definition, in the sequence, 2009-07-12*/
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
                    Shape_Shape[0][NNshape][Numberofgroup]++;
                    Shape_Shape[1][Numberofgroup][NNshape]++;
                }


            }//if (  Number_class >= NFragment )
            else
            {
                Numberofgroup = 0;
                Weight_per = 0;
                Svalue = -1;
            }

            /* calculate NPer[10], NPer[11], NPer[12] which are probabilities of the residue on H, S and R*/
            if   (     (  CTempHSRseq[ij] == '-'  )  ) /*borrow the result from step 1*/
            {
                int tmpint = 0;
                char tmpchar = ' ';
                sscanf(linesStep1Pred[ij],"%d %c %c %c %d %d %d ", &tmpint,&tmpchar,&tmpchar, &tmpchar,&NPer[10],&NPer[11],&NPer[12]);
            } else {
                //HSR
                Nconf_aas = 0;
                NPer[10] = Pred_Resulthsr[ij][0];
                NPer[11] = Pred_Resulthsr[ij][1];
                NPer[12] = Pred_Resulthsr[ij][2];



                NTG_COM = (NPer[10]+NPer[11]+NPer[12])/15;
                NPer[11] = NPer[11] + NTG_COM;
                HSR = NPer[10] + NPer[11] + NPer[12];
                //if (  Number_class >= (Halv_NFragment) )

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
            }
            /*based on NPer[10], NPer[11], NPer[12] to predict HSR, 2009-07-13*/
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
                Nconf_aas = int ( NPer[10+NAAS]*100.0/HSR + 0.5 );
            } else if (  Number_class <= 4 ) {
                continue;
                NAAS = 2;
                //FST_HSR[Number_class][Pred_Resulthsr[ij][3]][0]++;
                sprintf(OpenName,"%s/detailed.txt",CData_Result_location);
                wppp = fopen(OpenName,"a");
                fprintf(wppp,"%5s %4d %1d %1d\n",CSubname,ij,Pred_Resulthsr[ij][3],Number_class);
                fclose(wppp);
            } else {
                continue;
            }

            PHSR = 0;
            if (  NAAS == Pred_Resulthsr[ij][3] ) {
                HSRcorrect1++;
                PHSR = 1;
            }

            if ( Pred_Resulthsr[ij][3] != 3 ) {
                numValiedPredResidue ++;
            }

            HSRsum1++;

            HSR_DIS[0][Pred_Resulthsr[ij][3]]++;
            HSR_DIS[1][NAAS]++;




            HSR_HSR[0][Pred_Resulthsr[ij][3]][NAAS]++;
            HSR_HSR[1][NAAS][Pred_Resulthsr[ij][3]]++;

            //shape string
            //for eight shape
            if (  (NNshape>=0) && (NNshape<8)   )
            {
                FST_HSR[NNshape][Pred_Resulthsr[ij][3]][0][0]++;//for 8-shape all
                FST_HSR[NNshape][3][0][0]++;//for 8-shape all all
                if (  Svalue  >= 0.6  )
                { 
                    FST_HSR[NNshape][Pred_Resulthsr[ij][3]][0][1]++;//for 8-shape  right
                    FST_HSR[NNshape][3][0][1]++;//for 8-shape  right, all
                } 

                if (NNshape<=3) { 
                    HSR = 0;
                } else if (NNshape<=5) { 
                    HSR = 1;
                } else { 
                    HSR = 2;
                } 

                FST_HSR[HSR][Pred_Resulthsr[ij][3]][1][0]++;//for 3-shape all
                FST_HSR[HSR][3][1][0]++;//for 3-shape all, all
                if (  Svalue  >= 0.1  )
                {
                    FST_HSR[HSR][Pred_Resulthsr[ij][3]][1][1]++;//for 3-shape right
                    FST_HSR[HSR][3][1][1]++;//for 3-shape right,all
                } 
            }
            //for all 




            HSR_Cont[Pred_Resulthsr[ij][3]][0]++;
            HSR_Cont[NAAS][1]++;



            if (  CTempDSSPseq[ij] == 'H' ) {
                HSR = 0;
            } else if (   CTempDSSPseq[ij] == 'G'  ) {
                HSR = 1;
            } else if (   CTempDSSPseq[ij] == 'I'  ) {
                HSR = 2;
            } else if (   CTempDSSPseq[ij] == 'E'  ) {
                HSR = 3;
            } else if (   CTempDSSPseq[ij] == 'B'  ) {
                HSR = 4;
            } else if (   CTempDSSPseq[ij] == 'T'  ) {
                HSR = 5;
            } else if (   CTempDSSPseq[ij] == 'S'  ) {
                HSR = 6;
            } else /* (   CTempDSSPseq[ij] == 'R'  )*/ {
                HSR = 7;
            }
            Dist_conf[HSR][NAAS]++;

            Nconfd = Nconfd + Nconf_aas;
            Evolution = Evolution + statis[ij];

            if (  Nconf_aas <= 0 )
            {
                ik=ik;
            }

            //confidence and right percentage
            if ( Nconf_aas > 100 )
            {
                Nconf_aas = 100;
            }
            else if ( Nconf_aas < 0 )
            {
                Nconf_aas = 0;
            }

            Confi_AAS[Nconf_aas][0]++;
            Confi_AAS[Nconf_aas][1] = Confi_AAS[Nconf_aas][1] + PHSR;

            //statistic and right percentage
            temp1 = Integer (statis[ij]*100);
            if (statis[ij] == INIT_FLOAT){
                fprintf(stderr,"Uninitialized statis for ij=%d\n", ij);
                temp1 = -1;
            }

            if (  temp1 >= 250 ) {
                temp1 = 250;
            } else if (  temp1 < 0 ) {
                temp1 = 0;
            }

            Statis_AAS[temp1][0]++;
            Statis_AAS[temp1][1] = Statis_AAS[temp1][1] + PHSR;
            //confidence and statistic
            Confi_Statis[temp1][0]++;
            Confi_Statis[temp1][1] = Confi_Statis[temp1][1] + Nconf_aas;
            //combination of statistical and confidence
            temp1 = temp1/2;
            if (  temp1 >= 100 ) {
                temp1 = 99;
            } else if (  temp1 < 0 ) {
                temp1 = 0;
            }
            temp2 = Nconf_aas - 51;
            if (  temp2 >= 50 ) {
                temp2 = 49;
            } else if (  temp2 < 0 ) {
                temp2 = 0;
            }

            Stat_Conf_AAS[temp1][temp2][0]++;
            Stat_Conf_AAS[temp1][temp2][1] = Stat_Conf_AAS[temp1][temp2][1] + PHSR;



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
            /*debug CShape[Numberofgroup]; 2009-07-22 */
            fprintf(wfp, "  %c %4.1f %3d\n",CShape[Numberofgroup],Svalue,Weight_per);

        }
        fclose(wfp);

        if (  HSRsum1 >= 1  )
        {
            Homology_Per_chain[Homology_number][0] = Homology_Per_chain[Homology_number][0] + HSRsum1;
            Homology_Per_chain[Homology_number][1] = Homology_Per_chain[Homology_number][1] + HSRcorrect1;



            Totalvalue = Totalvalue + Sumvalue;
            Totalvalue1 = Totalvalue1 + Sumvalue1;
            TotalNumber = TotalNumber + Sumnumber;
            HSRcorrect = HSRcorrect + HSRcorrect1;
            totalNumValiedPredResidue = totalNumValiedPredResidue + numValiedPredResidue; /*total number of valied predicted residue positions*/
            HSRsum= HSRsum + HSRsum1;
            sprintf(OpenName,"%s/Test_resultper.txt",CData_Result_location);
            wfpper = fopen(OpenName,"a");

            //fprintf(wfpper,"%4d %5s %4d %6.1f %4d %6.2f %6.2f    %9.1f %7d %6.2f %6.2f      %4d %4d %6.2f  %8d %8d %6.2f\n",iktest,CSubname,Len_target,Sumvalue,Sumnumber,Sumvalue*100/Sumnumber,Sumvalue1*100/Sumnumber,Totalvalue,TotalNumber,Totalvalue*100/TotalNumber,Totalvalue1*100/TotalNumber,
            //HSRcorrect1,HSRsum1,HSRcorrect1*100.0/HSRsum1,HSRcorrect,HSRsum,HSRcorrect*100.0/HSRsum);

            fprintf(wfpper,"%4d %5s %6.1f %4d %6.2f %6.2f    %9.1f %7d %6.2f %6.2f      %4d %4d %6.2f  %8d %8d %6.2f\n",iktest,CSubname,Sumvalue,Sumnumber,Sumvalue*100/Sumnumber,Sumvalue1*100/Sumnumber,Totalvalue,TotalNumber,Totalvalue*100/TotalNumber,Totalvalue1*100/TotalNumber,
                    HSRcorrect1,numValiedPredResidue,HSRcorrect1*100.0/numValiedPredResidue,HSRcorrect,totalNumValiedPredResidue,HSRcorrect*100.0/totalNumValiedPredResidue);

            //fprintf(wfpper,"%5s %4d %6.1f %4d %6.2f %6.2f    %5.1f %5.1f %5.1f      %4d %4d %6.2f  %5.1f %5.1f %5.1f   %3d  %5.2f\n",CSubname,Len_target,Sumvalue,Sumnumber,Sumvalue*100/Sumnumber,Sumvalue1*100/Sumnumber,HSR_DIS[0][0]*100.0/HSRsum1,HSR_DIS[0][1]*100.0/HSRsum1,HSR_DIS[0][2]*100.0/HSRsum1,
            //HSRcorrect1,HSRsum1,HSRcorrect1*100.0/HSRsum1,HSR_DIS[1][0]*100.0/HSRsum1,HSR_DIS[1][1]*100.0/HSRsum1,HSR_DIS[1][2]*100.0/HSRsum1,OLDConfidence, Evolution/HSRsum1);
            fclose(wfpper);

            HSR = HSR_Cont[0][0] + HSR_Cont[1][0] + HSR_Cont[2][0];
            NAAS = HSR_Cont[0][1] + HSR_Cont[1][1] + HSR_Cont[2][1];
            sprintf(OpenName,"%s/detailed.txt",CData_Result_location);
            wppp = fopen(OpenName,"a");
            fprintf(wppp,"%4d %5s  %6.2f  %6.2f  %6.2f      %6.2f  %6.2f  %6.2f  %4d  %4d  %6.2f  %3d  %3d %3d %3d %1d\n",iktest, CSubname,
                    HSR_Cont[0][0]*100.0/HSR, HSR_Cont[1][0]*100.0/HSR, HSR_Cont[2][0]*100.0/HSR ,
                    HSR_Cont[0][1]*100.0/NAAS, HSR_Cont[1][1]*100.0/NAAS, HSR_Cont[2][1]*100.0/NAAS, HSRcorrect1,HSRsum1 ,
                    HSRcorrect1*100.0/HSRsum1, Nconfd/HSRsum1, int (Evolution*100/HSRsum1),Homology_Total, Realhomology,NonSCOP);
            fclose(wppp);
        }

        //whole protein
        /*
           if (  HSRsum1 >= 10  )
           {
           Confi_AAS[101][0] = int ( Confi_AAS[101][0]*1.0/HSRsum1+0.499);//conf
           if (  Confi_AAS[101][0] > 100  )
           {
           Confi_AAS[101][0] = 100;
           }
           Statis_AAS[251][0] = int (Statis_AAS[251][0]*1.0/HSRsum1+0.499);//statist
           if (  Statis_AAS[251][0] > 240  )
           {
           Statis_AAS[251][0] = 240;
           }
           temp2 = int ( HSRcorrect1*100.0/HSRsum1 + 0.499);
           Confi_AAS[Confi_AAS[101][0]][0]++;
           Confi_AASbak[Confi_AAS[101][0]][0] = Confi_AASbak[Confi_AAS[101][0]][0] + HSRsum1;
           Confi_AASbak[Confi_AAS[101][0]][1] = Confi_AASbak[Confi_AAS[101][0]][1] + HSRcorrect1;

           Statis_AAS[Statis_AAS[251][0]][0]++;
           Statis_AASbak[Statis_AAS[251][0]][0] = Statis_AASbak[Statis_AAS[251][0]][0] + HSRsum1;
           Statis_AASbak[Statis_AAS[251][0]][1] = Statis_AASbak[Statis_AAS[251][0]][1] + HSRcorrect1;

           }
           */


        fprintf(stdout,"%4d %s finished\n",iktest, CSubname);

        /*Free memory*/
        for (ik=0; ik<Len_target; ik++)
        {
            delete [] Candidate_Score[ik];
            delete [] Candidate_Score_iksub[ik];
            delete [] Candidate_Score_iklen[ik];
        }
        delete [] CSubname;
    }

    fclose(fptest);
#ifdef DEBUG_NUMCLASS
    fclose(fwpTmpDebugNumberClass);
#endif 


    sprintf(OpenName,"%s/detailed.txt",CData_Result_location);
    wppp = fopen(OpenName,"a");

    fprintf(wppp,"Sum for HSR %8d  %8d %6.2f\n",HSRcorrect,HSRsum,HSRcorrect*100.0/HSRsum);
    fprintf(wppp,"\n\n");


    fprintf(wppp,"shape string\n");
    fprintf(wppp,"eihgt shape string\n");
    for (ik=0; ik<4; ik++)//7 shape
    {
        FST_HSR[6][ik][0][0] = FST_HSR[6][ik][0][0] + FST_HSR[7][ik][0][0];
        FST_HSR[6][ik][0][1] = FST_HSR[6][ik][0][1] + FST_HSR[7][ik][0][1];
    }
    fprintf(wppp,"1--helix, 2--sheet, 3--random, 4-all\n");
    for (ik=0; ik<4; ik++)
    {
        fprintf(wppp,"%1d-observed ", ik);
        HSR = 0;
        for (ij=0; ij<7; ij++)
        {
            HSR = HSR + FST_HSR[ij][ik][0][0];
        }
        for (ij=0; ij<7; ij++)
        {
            fprintf(wppp,"%6.2f ", FST_HSR[ij][ik][0][0]*100.0/HSR);
        }
        fprintf(wppp,"%8d\n",HSR);

        fprintf(wppp,"%1d-predictd ", ik);
        NAAS = 0;
        for (ij=0; ij<7; ij++)
        {
            fprintf(wppp,"%6.2f ", FST_HSR[ij][ik][0][1]*100.0/FST_HSR[ij][ik][0][0]);
            NAAS = NAAS + FST_HSR[ij][ik][0][1];
        }
        fprintf(wppp,"%8.2f\n",NAAS*100.0/HSR);
    }



    fprintf(wppp,"\n");
    fprintf(wppp,"observed: 3-shape  predicted\n");

    for (ik=0; ik<4; ik++)
    {
        fprintf(wppp,"%1d-observed ", ik);
        HSR = 0;
        for (ij=0; ij<3; ij++)
        {
            HSR = HSR + FST_HSR[ij][ik][1][0];
        }
        for (ij=0; ij<3; ij++)
        {
            fprintf(wppp,"%6.2f ", FST_HSR[ij][ik][1][0]*100.0/HSR);
        }
        fprintf(wppp,"%8d\n",HSR);

        fprintf(wppp,"%1d-predictd ", ik);
        NAAS = 0;
        for (ij=0; ij<3; ij++)
        {
            fprintf(wppp,"%6.2f ", FST_HSR[ij][ik][1][1]*100.0/FST_HSR[ij][ik][1][0]);
            NAAS = NAAS + FST_HSR[ij][ik][1][1];
        }
        fprintf(wppp,"%8.2f\n",NAAS*100.0/HSR);
    }





    fprintf(wppp,"DSSP to HSR\n");
    for (ik=0; ik<8; ik++)
    {
        HSR = 0;
        for (ij=0; ij<3; ij++)
        {
            HSR = HSR + Dist_conf[ik][ij];
        }

        if (  HSR >= 1  )
        {
            fprintf(wppp,"%2d %8d   ",ik,HSR);
            for (ij=0; ij<3; ij++)
            {
                fprintf(wppp,"%6.2f ",Dist_conf[ik][ij]*100.0/HSR);
            }
            fprintf(wppp,"\n");
        }

    }


    fprintf(wppp,"HSR to DSSP\n");
    for (ik=0; ik<3; ik++)
    {
        HSR = 0;
        for (ij=0; ij<8; ij++)
        {
            HSR = HSR + Dist_conf[ij][ik];
        }

        if (  HSR >= 1  )
        {
            fprintf(wppp,"%2d %8d   ",ik,HSR);
            for (ij=0; ij<8; ij++)
            {
                fprintf(wppp,"%6.2f ",Dist_conf[ij][ik]*100.0/HSR);
            }
            fprintf(wppp,"\n");
        }

    }
    fprintf(wppp,"\nHSR to HSR, Observed to predicted \n");
    for (ik=0; ik<3; ik++)
    {
        HSR = 0;
        for (ij=0; ij<3; ij++)
        {
            HSR = HSR + HSR_HSR[0][ik][ij];
        }
        fprintf(wppp,"%8d    ",HSR);
        for (ij=0; ij<3; ij++)
        {
            fprintf(wppp,"%8d ",HSR_HSR[0][ik][ij]);
        }
        fprintf(wppp,"%6.2f\n",HSR_HSR[0][ik][ik]*100.0/HSR);
    }

    fprintf(wppp,"\nHSR to HSR, predicted to Observed \n");
    for (ik=0; ik<3; ik++)
    {
        HSR = 0;
        for (ij=0; ij<3; ij++)
        {
            HSR = HSR + HSR_HSR[1][ik][ij];
        }
        fprintf(wppp,"%8d    ",HSR);
        for (ij=0; ij<3; ij++)
        {
            fprintf(wppp,"%8d ",HSR_HSR[1][ik][ij]);
        }
        fprintf(wppp,"%6.2f\n",HSR_HSR[1][ik][ik]*100.0/HSR);
    }

    fprintf(wppp,"\nShape to Shape, Observed to predicted \n");
    for (ik=0; ik<8; ik++)
    {
        HSR = 0;
        for (ij=0; ij<8; ij++)
        {
            HSR = HSR + Shape_Shape[0][ik][ij];
        }
        fprintf(wppp,"%8d    ",HSR);
        for (ij=0; ij<8; ij++)
        {
            fprintf(wppp,"%8d ",Shape_Shape[0][ik][ij]);
        }
        fprintf(wppp,"%6.2f\n",Shape_Shape[0][ik][ik]*100.0/HSR);
    }

    fprintf(wppp,"\nShape to Shape, predicted to Observed \n");
    for (ik=0; ik<7; ik++)
    {
        HSR = 0;
        for (ij=0; ij<8; ij++)
        {
            HSR = HSR + Shape_Shape[1][ik][ij];
        }
        fprintf(wppp,"%8d    ",HSR);
        for (ij=0; ij<8; ij++)
        {
            fprintf(wppp,"%8d ",Shape_Shape[1][ik][ij]);
        }
        fprintf(wppp,"%6.2f\n",Shape_Shape[1][ik][ik]*100.0/HSR);
    }

    fprintf(wppp,"\n confidence \n");
    for (ik=0; ik<101; ik++)
    {
        if (  Confi_AAS[ik][0] >= 1   )
        {
            fprintf(wppp,"%3d  %9d %9d  %6.2f \n",ik, Confi_AAS[ik][1], Confi_AAS[ik][0],Confi_AAS[ik][1]*100.0/Confi_AAS[ik][0] );
        }
    }

    fprintf(wppp,"\n Statistical qaulity VS confidence and right percentage\n");
    for (ik=0; ik<251; ik++)
    {
        if  (   (Statis_AAS[ik][0]>=1) && (Confi_Statis[ik][0]>=1)   )
        {
            fprintf(wppp,"%3d  %9d  %9d ",ik,Statis_AAS[ik][1],Statis_AAS[ik][0] );
            fprintf(wppp," %6.2f ", Confi_Statis[ik][1]*1.0/Confi_Statis[ik][0] );
            fprintf(wppp," %6.2f\n", Statis_AAS[ik][1]*100.0/Statis_AAS[ik][0] );
        }
    }


    fprintf(wppp,"\n Statistical qaulity ,confidence,right percentage for single amino acids\n");
    for (ik=0; ik<100; ik++)
    {
        for (ij=0; ij<50; ij++)
        {
            if (  Stat_Conf_AAS[ik][ij][0] >= 1 )
            {
                fprintf(wppp,"%3d %3d  %7d %7d  %6.2f\n",ik,ij,Stat_Conf_AAS[ik][ij][1],Stat_Conf_AAS[ik][ij][0], Stat_Conf_AAS[ik][ij][1]*100.0/Stat_Conf_AAS[ik][ij][0]);
            }
        }
    }


    fprintf(wppp,"\n result based on the chain homology \n");
    for (ik=0; ik<20; ik++)
    {
        fprintf(wppp,"%8d %8d  %6.2f\n",Homology_Per_chain[ik][1], Homology_Per_chain[ik][0], Homology_Per_chain[ik][1]*100.0/Homology_Per_chain[ik][0]);
    }

    fprintf(wppp,"\n result based on the segment(AAS) homology \n");
    HSR = 0;
    NAAS = 0;
    for (ik=0; ik<20; ik++)
    {
        fprintf(wppp,"%8d %8d  %6.2f\n",Homology_Per_AAS[ik][1], Homology_Per_AAS[ik][0], Homology_Per_AAS[ik][1]*100.0/Homology_Per_AAS[ik][0]);
        HSR = HSR + Homology_Per_AAS[ik][1];
        NAAS = NAAS + Homology_Per_AAS[ik][0];
    }
    fprintf(wppp,"result:  %8d %8d  %6.2f\n",HSR,NAAS,HSR*100.0/NAAS);


    fclose(wppp);
    //fclose(fp4416);

    /* Free memory*/
    delete [] blosumFile            ;
    delete [] CDataBase             ;
    delete [] CBlast_Location_mod   ;
    delete [] CBlast_Location_Frag  ;
    delete [] CBlast_Location_Qij   ;
    delete [] CData_Result_location ;
    delete [] sourcefile            ;
    delete [] First_round_HSR       ;
    delete [] first_round_loc       ;
    delete [] Ctestlist             ;
    delete [] Ctestqij              ;
    delete [] Ctestmod              ;

    delete []  CTempAASseq;
    delete []  CTempShapeseq;
    delete []  CTempHSRseq;
    delete []  CTempDSSPseq;

    delete []  Target_Naas;


    delete [] LengthList;
    delete [] Candidate_Distribution;
    delete [] Candidate_Distribution_per;
    delete [] Homology_List;
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

    delete [] Candidate_Score_iksub ;
    delete [] Candidate_Score_iklen ;
    delete [] Candidate_Score  ;
    //     for (ik=0; ik<=NonSCOP; ik++)
    //     {
    //         if (  Number_fam[ik] >= 1   )
    //         {
    //             delete homo_chain_number[ik];
    //         }
    //     }

    return EXIT_SUCCESS;
}/*}}}*/
