// build_hsrfrag.cpp : building profiles based on predicted secondary structures

/*
 * =====================================================================================
 * 
 *        Program:  build_hsrfrag
 *        Description:  building profiles based on predicted secondary structures
 *                     1) the trainging set is the same as that used for search_fragment program;
 *                     2) the frgament blocks are built on the base of the amino acid similarity and the similarity between predicted and observed secondary structure;
Paraneters:
 *                  -ltst   CHSR_list;
 *                  -qtst   CHSR_list_qij;
 *                  -sloc   First_Round_Prediction_Location;
 *                  -ploc   CResult_Location;
 *                  -bloc   Clocation_blosum;
 *                  -ltrn   CDatabase;
 *                  -qtrn   database_qijmatrix;
 *                  -sbeg   Searchbeg;
 *                  -send   Searchend;
 * =====================================================================================
 */


#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <sys/types.h>
#include <sys/stat.h>
#include <algorithm>
#include "base.h"
#include "array.h"
#include "mypro.h"
#include "myfunc.h"
#include "mytemplate.h"

using namespace std;
/*valgrind checked 2009-07-13, no error and no memory leak, Nanjiang*/
/*2009-07-14, updated, the last column of the Frag_ file added with the
 * confidence of the predicted residue position*/

int method_HSR = 1; /*using method 1 to determine the predicted secondary structure from three probabilities*/

int isReadBinaryFile = true;

string usage="\n\
usage: build_hsrfrag [options]\n\
Description: building profiles from predicted secondary structure\n\
Options:\n\
  -dbtype INT  (default: 1)\n\
          0: database stored as one id one file\n\
          1: database stored as dumped file\n\
  -ltst   The file path and name for testing set(text file, only one columns:\n\
          chain ID). \n\
  -qtst   Location of qijmatrix for testing set\n\
  -rloc   Location for predicted secondary structure\n\
  -ploc   Location for coming HSR profiles\n\
  -bloc   The location and file name for blosum table(20 X 20 table)\n\
  -ltrn   The file path and name for training set(text file, two columns:\n\
          series number and chain ID). \n\
  -qtrn   Location of qijmatrix for training set\n\
  -sbeg   The start series number of chains to be made for profiles\n\
  -send   The end series number of chains to be made for profiles\n\
Created 2009-07-24, updated 2011-10-17, Nanjiang Shu\n\
";
void PrintHelp() { 
    fprintf(stdout,"%s\n",usage.c_str());
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
    if(argc < 2) 
    {
        PrintHelp();
        return 0;
    }

    int dbtype=1; /*0: database stored as one file for each id, 
                    1: database stored as dumped files. 
                    added 2011-10-17*/

    //Build_profile_HSR_Predition(char *CHSR_list, char *First_Round_Prediction_Location, char *programpath, char *CHSR_list_qij, char *CDatabase, char *database_qijmatrix)
    //make the profiles which are based on the first round secondary structure (HSR) predicted result
    char first_non_blank_char = ' '; /*2009-07-15, Nanjiang, using first_non_blank_char to determine whether the matrix are valid data lines*/
    FILE *fp, *wfp,*fpread, *fplist;
    char  str_char1[20], str_char2[20];
    char CSubunit[200];
    int NSeries,ik,ij, Max_Length, NNPer[21], HSR, NAAS; 

    Array1D <char> CTempShape_1darray(LONGEST_SEQ);
    Array1D <char> CTempAAS_1darray(LONGEST_SEQ);
    Array1D <char> CTempHSR_1darray(LONGEST_SEQ);
    char *CTempShape = CTempShape_1darray.array1D;
    char *CTempAAS = CTempAAS_1darray.array1D;
    char *CTempHSR = CTempHSR_1darray.array1D;

    float param1, param2;
    bool TCount;
    int Sumunit, *LengthList,iij,iikdata,iijdata;
    int  Blosum62_table[21][21];
    int NScore, High_Score[300], Pos_Chain[300], Pos_length[300];
    int NSamplehigh,low_score, NSore_Fragment;
    int INVALID_ACC;
    int HSR_HSR[4][4], Part1, Part2, NScore_AAS ,NScore_HSR;
    int Target_Length, Ntar, Searchbeg, Searchend;
    int *(*CCODEAAS), *(*CCODEHSR);

    Array2D <int> Fragment_profile_2darray(LONGEST_SEQ, 20);
    Array1D <int> CTarHSR_1darray(LONGEST_SEQ);
    Array1D <int> CTarAAS_1darray(LONGEST_SEQ);
    int **Fragment_profile = Fragment_profile_2darray.array2D;
    int *CTarHSR = CTarHSR_1darray.array1D;
    int *CTarAAS = CTarAAS_1darray.array1D;

    int Nfragment;
    //int AAS_Code[26];

    /*File names*/
    /*
    char *CResult_Location = NULL;
    char *CHSR_list = NULL;
    char *First_Round_Prediction_Location = NULL;
    char *CHSR_list_qij = NULL ;
    char *CDatabase = NULL ;
    char *database_qijmatrix = NULL;
    char *Clocation_blosum = NULL ;
    */

    char CResult_Location[MAX_PATH+1];
    char openfile[MAX_PATH+1];
    char CHSR_list[MAX_PATH+1] = ""; 
    char First_Round_Prediction_Location[MAX_PATH+1] = "";  
    char CHSR_list_qij[MAX_PATH+1] =""; 
    char CDatabase[MAX_PATH+1] =""; 
    char database_qijmatrix[MAX_PATH+1] = "";
    char Clocation_blosum[MAX_PATH+1] = "";

    Searchbeg = 0;
    Searchend = 60000;

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
                if( ( i = option_parser_filename_old(argc, argv, i, CHSR_list)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "-dbtype") == 0) || (strcmp(argv[i], "--dbtype")== 0))  
            {
                if( ( i = option_parser_numeric(argc, argv, i, dbtype, true, 0, 1)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-qtst") == 0) || (strcmp(argv[i], "--qij-test") == 0))  
            {
                if( ( i = option_parser_filename_old(argc, argv, i, CHSR_list_qij)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-rloc") == 0) )  
            {
                if( ( i = option_parser_filename_old(argc, argv, i, First_Round_Prediction_Location)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-ploc") == 0) || (strcmp(argv[i], "--outpath") == 0))  
            {
                if( ( i = option_parser_filename_old(argc, argv, i, CResult_Location)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"-bloc") == 0) || (strcmp(argv[i], "--blosum") == 0))  
            {
                if( ( i = option_parser_filename_old(argc, argv, i, Clocation_blosum)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--trainlist") == 0) || strcmp(argv[i],"-ltrn") == 0 )
            {
                if( ( i = option_parser_filename_old(argc, argv, i, CDatabase)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i],"--qij-train") == 0) || strcmp(argv[i],"-qtrn") == 0) 
            {
                if( ( i = option_parser_filename_old(argc, argv, i, database_qijmatrix)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "--rb") == 0) || (strcmp(argv[i], "--read-binary")== 0))  
            {
                isReadBinaryFile = true;
                i++;
            }
            else if( (strcmp(argv[i], "-sbeg") == 0) || (strcmp(argv[i], "--beg") == 0) )  
            {
                if( ( i = option_parser_numeric(argc, argv, i, Searchbeg, true, 0, 5000000)) == -1)
                    return -1;
            }
            else if( (strcmp(argv[i], "-send") == 0) || (strcmp(argv[i], "--end") == 0) )  
            {
                if( ( i = option_parser_numeric(argc, argv, i, Searchend, true, 0, 5000000)) == -1)
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

    VerifyFolder(CHSR_list );
    VerifyFolder(First_Round_Prediction_Location );
    VerifyFolder(CHSR_list_qij );
    VerifyFolder(CDatabase);
    VerifyFolder(database_qijmatrix );
    VerifyFolder(CResult_Location);

    int linesize;
    int maxline = 1000;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D; 
    int status_sscanf = 0;


    INVALID_ACC = 10;
    NSore_Fragment = 15;
    Nfragment = 9;
    //HSR table --- for the comparison between secondary structures
    HSR_HSR[0][0] = 4;
    HSR_HSR[0][1] = -3;
    HSR_HSR[0][2] = -1;
    HSR_HSR[1][0] = -3;
    HSR_HSR[1][1] = 6;
    HSR_HSR[1][2] = 3;
    HSR_HSR[2][0] = -1;
    HSR_HSR[2][1] = -1;
    HSR_HSR[2][2] = 3;
    NSamplehigh = 100;
    HSR = 0;

    fp = fopen(Clocation_blosum,"r");
    while ((linesize = fgetline(fp,line,maxline)) != EOF)
    {
        status_sscanf = sscanf(line,"%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d",
                &NNPer[0],&NNPer[1],&NNPer[2],&NNPer[3],&NNPer[4],&NNPer[5],&NNPer[6],&NNPer[7],&NNPer[8],&NNPer[9],
                &NNPer[10],&NNPer[11],&NNPer[12],&NNPer[13],&NNPer[14],&NNPer[15],&NNPer[16],&NNPer[17],&NNPer[18],&NNPer[19] );
        if (status_sscanf != 20)
        {
            fprintf(stderr, "Format error for blosum file: %s\n", Clocation_blosum);
            fprintf(stderr, "line: %s\n", line);
            assert(status_sscanf == 20);
        }
        for (ik=0; ik<20; ik++)
        {
            Blosum62_table[HSR][ik] = NNPer[ik];
        }
        HSR++;
    }
    fclose(fp);

    Sumunit = 0;
    fp = fopen(CDatabase,"r");
    while ((linesize = fgetline(fp,line,maxline)) != EOF)
    {
        Sumunit++;
    }
    rewind(fp);
    fprintf(stdout, "Reading training set with %d sequences...\n", Sumunit);
    CCODEAAS = new ( int (*[Sumunit]) );
    CCODEHSR = new ( int (*[Sumunit]) );
    LengthList = new int [Sumunit];
    Max_Length = LONGEST_SEQ;
    Sumunit = 0;

    char alphabet[50] = "";
    double parameter[8];
    int8 typeProfile = 1;
    int seqLength = 0;
    Array1D <ProfileSADByte> profileSADByte_1darray(LONGEST_SEQ);
    ProfileSADByte * profileSADByte = profileSADByte_1darray.array1D;

    map <string, dbindex> dbindexqij;
    int maxDBIndexNumber_qij = 0;
    vector <FILE*> fpList_qij;
    if (dbtype==1){ /*read index file of the dumped database, added 2011-10-13*/
        if (ReadDatabaseIndex(string(database_qijmatrix), dbindexqij, maxDBIndexNumber_qij) == -1 ){
            fprintf(stderr,"Read db qij failed for %s\n", database_qijmatrix);
            exit(1);
        }
        GetDBFPList( fpList_qij, string(database_qijmatrix), maxDBIndexNumber_qij);
    }

    FILE *fp_qij=NULL;
    int status_fseek = 0;

    char id[200];

    while ((linesize = fgetline(fp,line,maxline)) != EOF)
    {
        status_sscanf = sscanf(line,"%d %s",&HSR,CSubunit);
        if (status_sscanf != 2)     /*2009-07-24*/
        {
            fprintf(stderr,"Format error for training list file: %s\n", CDatabase);
            fprintf(stderr,"line(%d): %s\n", Sumunit, line);
            assert(status_sscanf == 2);
        }
        strcpy(id,CSubunit);

        if (dbtype==1){
            if (dbindexqij.find(id) == dbindexqij.end()){
                fprintf(stderr,"Can not find id %s in db %s. Ignore this id.\n",id, database_qijmatrix);
                continue;
            }

            fp_qij = fpList_qij[dbindexqij[id].dbfileindex];
            if ((status_fseek = fseek(fp_qij, dbindexqij[id].offset, SEEK_SET)) != 0){
                fprintf(stderr,"Fatal! fseek failed for id %s in db qij. Exit.\n",id);
                exit(1);
            }
        }

        if (!isReadBinaryFile)
        {
            if (dbtype==0){
                sprintf(openfile,"%s/%s.Qij",database_qijmatrix,CSubunit);
                fpread = fopen(openfile,"r");
                if (  fpread == NULL  ) { continue; }
            }else{
                fpread=fp_qij;
            }
            NSeries = 0;
            while ((linesize = fgetline(fpread,line,maxline)) != EOF)
                //while (  fgets(string,300,fpread) != NULL  )
            {
                sscanf(line, " %c", &first_non_blank_char);
                if (first_non_blank_char <'0' || first_non_blank_char > '9') /*only when the first_non_blank_char is digit */
                {
                    continue;
                }
                sscanf(line," %d ", &NSeries );
                NSeries = NSeries - 1;
                //percent
                status_sscanf = sscanf(line,"%d %c  %s%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%f%f",
                        &HSR,&str_char1[0],str_char2,&NNPer[0],&NNPer[1],&NNPer[2],&NNPer[3],&NNPer[4],&NNPer[5],&NNPer[6],&NNPer[7],&NNPer[8],&NNPer[9],
                        &NNPer[10],&NNPer[11],&NNPer[12],&NNPer[13],&NNPer[14],&NNPer[15],&NNPer[16],&NNPer[17],&NNPer[18],&NNPer[19],&param1,&param2);
                if (status_sscanf != 25)     /*2009-07-24*/
                {
                    fprintf(stderr,"Format error for the file: %s\n", openfile);
                    fprintf(stderr,"line(%d): %s\n", NSeries, line);
                    assert(status_sscanf == 25);
                }
                //amino acids
                CTempAAS[NSeries] = str_char1[0];
                //HSR
                CTempHSR[NSeries] = str_char2[2];

            }
            if (dbtype==0 && fpread != NULL){
                fclose(fpread);
            }
            NSeries++;
        } else {   /*Read In bindary Qijfile*/
            if (dbtype==0){
                sprintf(openfile,"%s/%s.Qijbin",database_qijmatrix,CSubunit);
                if (GetBinaryMODM(openfile, alphabet, seqLength, profileSADByte, parameter, typeProfile) == -1) {
                    fprintf(stderr, "can not open QijFile %s\n", openfile);
                    return -1;
                }
            }else if (dbtype==1){
                if (GetBinaryMODM_fp(fp_qij, dbindexqij[id].size, alphabet, seqLength, profileSADByte, parameter, typeProfile) == -1) {
                    fprintf(stderr, "Reading qij matrix failed for %s. offset=%ld, size=%ld. __LINE__=%d\n", id, dbindexqij[id].offset, dbindexqij[id].size, __LINE__);
                    return -1;
                }
            } else {
                fprintf(stderr,"dbtype = %d has not been implemented yet. Exit.\n", dbtype);
                exit(1);
            }
            for (ik = 0; ik<seqLength; ik++)
            {
                CTempAAS[ik] = profileSADByte[ik].aa;
                CTempHSR[ik] = profileSADByte[ik].dsspSec;
            }

            NSeries = seqLength;
        }

        CCODEAAS[Sumunit] = new int[NSeries+1];
        CCODEHSR[Sumunit] = new int [NSeries+1];
        LengthList[Sumunit] = NSeries;

        for (ik=0; ik<NSeries; ik++)
        {
            HSR = CTempAAS[ik] - 'A';
            if  (  (HSR>=0) && (HSR<=25)  )
            {
                HSR = AAS_Code[HSR];
            }
            else
            {
                HSR = 20;
            }
            CCODEAAS[Sumunit][ik] = HSR;
            //HSR
            if  (  (CTempHSR[ik]=='H') || (CTempHSR[ik]=='G') || (CTempHSR[ik]=='I')   )
            {
                HSR = 0;
            }
            else if (  CTempHSR[ik] == 'E'  )
            {
                HSR = 1;
            }
            else
            {
                HSR = 2;
            }
            CCODEHSR[Sumunit][ik] = HSR;
        }
        Sumunit++;
    }
    if (dbtype==1){
        int i;
        for (i=0;i<=maxDBIndexNumber_qij;i++){ fclose(fpList_qij[i]); }
    }
    fclose(fp);
    fprintf(stdout, "Read training set finished, %d sequences are read in\n", Sumunit);



    fprintf(stdout, "Start HSR profile building...\n");
    NSamplehigh = 100;
    Ntar = -1;
    fplist = fopen(CHSR_list,"r");
    while ((linesize = fgetline(fplist,line,maxline)) != EOF)
	{
		if (linesize <=0)
		{
			continue;
		}

		sscanf(line, "%s",CSubunit);
		Ntar++;

        if  (   (Ntar<Searchbeg) || (Ntar>=Searchend)   )
        {
            continue;
        }

        if (!isReadBinaryFile)
        {
            sprintf(openfile, "%s/%s.Qij",CHSR_list_qij,CSubunit); /*test Qij file*/
            fp = fopen(openfile,"r");
            if ( fp == NULL  ) { continue; }
            //while (  fgets(string,300,fp) != NULL  )
            while ((linesize = fgetline(fp,line,maxline)) != EOF)
            {
                sscanf(line, " %c", &first_non_blank_char);
                if (first_non_blank_char <'0' || first_non_blank_char > '9') /*only when the first_non_blank_char is digit */
                {
                    continue;
                }
                sscanf(line," %d ", &NSeries );
                NSeries = NSeries - 1;
                //percent
                status_sscanf = sscanf(line,"%d %c  %s %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%f%f",
                        &HSR,&str_char1[0],str_char2,&NNPer[0],&NNPer[1],&NNPer[2],&NNPer[3],&NNPer[4],&NNPer[5],&NNPer[6],&NNPer[7],&NNPer[8],&NNPer[9],
                        &NNPer[10],&NNPer[11],&NNPer[12],&NNPer[13],&NNPer[14],&NNPer[15],&NNPer[16],&NNPer[17],&NNPer[18],&NNPer[19],&param1,&param2);
                if (status_sscanf != 25)     /*2009-07-24*/
                {
                    fprintf(stderr,"Format error for file: %s\n",openfile );
                    fprintf(stderr,"line(%d): %s\n", NSeries, line);
                    assert(status_sscanf == 25);
                }
                //amino acids
                CTempAAS[NSeries] = str_char1[0];
                //shape string
                CTempShape[NSeries] = str_char2[0];
                //HSR
                CTempHSR[NSeries] = str_char2[2];
                //access surface
            }
            fclose(fp);
            Target_Length = NSeries + 1;
        }
        else/*{{{*/
        {   /*Read In bindary Qijfile for the test set*/
            sprintf(openfile, "%s/%s.Qijbin",CHSR_list_qij,CSubunit); /*test Qij file*/
            if (GetBinaryMODM(openfile, alphabet, seqLength, profileSADByte, parameter, typeProfile) == -1)
            {
                fprintf(stderr, "can not open QijFile %s\n", openfile);
                return -1;
            }
            for (ik = 0; ik<seqLength; ik++)
            {
                CTempAAS[ik] = profileSADByte[ik].aa;
                CTempHSR[ik] = profileSADByte[ik].dsspSec;
                CTempShape[ik] = profileSADByte[ik].shape;
            }
            Target_Length = seqLength;
        }/*}}}*/

        for (ik=0; ik<Target_Length; ik++)
        {
            HSR = CTempAAS[ik] - 'A';
            if  (  (HSR>=0) && (HSR<=25)  )
            {
                HSR = AAS_Code[HSR];
            }
            else
            {
                HSR = 20;
            }
            CTarAAS[ik] = HSR;
            CTarHSR[ik] = 3;
        }

        Array1D <float> HSR_confidence_1darray(Target_Length+1);
        float *HSR_confidence = HSR_confidence_1darray.array1D;
        HSR_confidence_1darray.Init(float(0.0));

        //read the HSR from first-round prediction
        sprintf(openfile, "%s/Res_%s.txt",First_Round_Prediction_Location,CSubunit);
        fpread =  fopen(openfile,"r");
        if (  fpread == NULL  ) { continue; }
        //while (  fgets(string,300,fpread) != NULL  )
        while ((linesize = fgetline(fpread,line,maxline)) != EOF)
        { 
            status_sscanf = sscanf(line,"%d %c %c %c %d %d %d %d %d %d %d %d %d %d %d %d %c %f %d",
                    &NAAS,&str_char1[0],&str_char1[1], &str_char1[2],&NNPer[0],&NNPer[1],&NNPer[2], &NNPer[3],
                    &NNPer[11],&NNPer[12],&NNPer[13],&NNPer[14],&NNPer[15],&NNPer[16],&NNPer[17],&NNPer[18], &str_char1[4],&param2, &NNPer[19]);


            if (status_sscanf != 19)     /*2009-07-24*/
            {
                fprintf(stderr,"Format error for the file: %s\n", openfile);
                fprintf(stderr,"line(%d): %s\n", NAAS,line );
                assert(status_sscanf == 19);
            }
            HSR = GetHSRState (NNPer[0], NNPer[1], NNPer[2], method_HSR);
            /*calculate the confidence 2009-07-14, Nanjiang*/
            HSR_confidence[NAAS] = GetHSRConfidence(NNPer[0], NNPer[1], NNPer[2]) ;
            
            //if  (   (NNPer[1]>=NNPer[0]) && (NNPer[1]>=NNPer[2])   )
            //{
                //HSR = 1;
            //}
            //else if (   (NNPer[0]>=NNPer[2]) && (NNPer[0]>=NNPer[1])   )
            //{
                //HSR = 0;
            //}
            //else
            //{
                //HSR = 2;
            //}
            CTarHSR[NAAS] = HSR;
        } 
        fclose(fpread);


        //profiles
        for (iikdata=0; iikdata<Target_Length; iikdata++)
        {
            for (iijdata=0; iijdata<20; iijdata++)
            {
                Fragment_profile[iikdata][iijdata] = 0;
            }
        }
        for (iij=0; iij<Target_Length-9; iij++)
        {
            TCount = 1;
            for (ik=iij; ik<iij+9; ik++)
            {
                if (  (CTarAAS[ik]>=20) || (CTarHSR[ik]>=3)  )
                {
                    TCount = 0;
                    break;
                }
            } 
            if (  TCount == 0  )
            {
                continue;
            }

            for (ik=0; ik<NSamplehigh; ik++)
            {   
                High_Score[ik] = -5;
            }   
            for (iikdata=0; iikdata<Sumunit; iikdata++)
            { 
                for (iijdata=0; iijdata<LengthList[iikdata]-9; iijdata++)
                {
                    TCount = 1;
                    for (ij=iijdata; ij<iijdata+9; ij++)
                    {
                        if (  CCODEAAS[iikdata][ij] >= 20  )
                        {
                            TCount = 0;
                            break;
                        }
                    }
                    if (  TCount == 0  )
                    {
                        continue;
                    }

                    //compare aas of fragment to check the similarity
                    NScore_AAS = 0;
                    for (ik=0; ik<Nfragment; ik++)
                    {
                        NScore_AAS = NScore_AAS + Blosum62_table[CTarAAS[iij+ik]][CCODEAAS[iikdata][iijdata+ik]];
                    }
                    if (   NScore_AAS < NSore_Fragment  )
                    {
                        continue;
                    }
                    NScore_HSR = 0;
                    //compare the HSR between predicted from target and observed from database
                    for (ik=0; ik<Nfragment; ik++)
                    {
                        Part1 = CTarHSR[iij+ik];
                        Part2 = CCODEHSR[iikdata][iijdata+ik];
                        NScore_HSR = NScore_HSR + HSR_HSR[Part1][Part2];
                    }

                    NScore = NScore_HSR + NScore_AAS;
                    if (   NScore < High_Score[0]  )
                    {
                        continue;
                    }
                    High_Score[0] = NScore;
                    Pos_Chain[0] = iikdata;
                    Pos_length[0] = iijdata;
                    low_score = NScore;
                    HSR = 0;
                    for (ik=1; ik<NSamplehigh; ik++)
                    {
                        if ( High_Score[ik] < low_score )
                        {
                            low_score = High_Score[ik];
                            HSR = ik;
                        }
                    }

                    if (  HSR > 0 )
                    {
                        //score
                        Part2 = High_Score[0];
                        High_Score[0] = High_Score[HSR];
                        High_Score[HSR] = Part2;
                        //pos-chaim
                        Part2 = Pos_Chain[0];
                        Pos_Chain[0] = Pos_Chain[HSR];
                        Pos_Chain[HSR] = Part2;
                        //Pos-length
                        Part2 = Pos_length[0];
                        Pos_length[0] = Pos_length[HSR];
                        Pos_length[HSR] = Part2;
                    }
                } 
            } 
            //building profiles
            NScore = NSore_Fragment + 20;
            for (ik=0; ik<NSamplehigh; ik++)
            {
                if ( High_Score[ik] < NScore  )
                {
                    continue;
                }
                for (ij=0; ij<Nfragment; ij++)
                {
                    Part1 = CCODEAAS[Pos_Chain[ik]][Pos_length[ik]+ij];
                    Fragment_profile[iij+ij][Part1]++;
                }
            }
        }
        sprintf(openfile,"%s/Frag_%s.txt",CResult_Location,CSubunit);
        wfp = fopen(openfile,"w");
        for (iikdata=0; iikdata<Target_Length; iikdata++)
        {
            HSR = 0;
            for (iijdata=0; iijdata<20; iijdata++)
            {
                HSR = HSR + Fragment_profile[iikdata][iijdata];
            }
            if (  HSR < 100  )
            {
                continue;
            }
            fprintf(wfp,"%4d %c A9R  ",iikdata+1, CTempAAS[iikdata]);
            for (iijdata=0; iijdata<20; iijdata++)
            {
                fprintf(wfp,"%3d ",int(Fragment_profile[iikdata][iijdata]*100.0/HSR+0.499) );
            }
            fprintf(wfp," 1.00 %.2f\n", HSR_confidence[iikdata]*100);
        }
        fclose(wfp);
        fprintf(stdout,"%d %s HSR profile built, output to %s\n", Ntar, CSubunit, openfile);
    }
    fclose(fplist);


    delete [] LengthList;
    for (ik=0; ik<Sumunit; ik++)
    {
        delete [] CCODEAAS[ik];
        delete [] CCODEHSR[ik];
    }
    delete [] CCODEAAS;
    delete [] CCODEHSR;

    /*
    delete [] CResult_Location ;
    delete [] openfile ;
    delete [] CHSR_list ;
    delete [] First_Round_Prediction_Location ;
    delete [] CHSR_list_qij ;
    delete [] CDatabase ;
    delete [] database_qijmatrix ;
    delete [] Clocation_blosum ;
    */

    return EXIT_SUCCESS;
}
