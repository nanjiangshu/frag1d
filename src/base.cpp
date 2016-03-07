#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include "base.h"
#include "mytemplate.h"
#include "myfunc.h"
#include "mypro.h"
#include "array.h"

using namespace std;

/*****************************************************************************
 * ChangeLog 2009-10-14
 *      The DSSP8to3 mapping scheme updated, see the annotation of this
 *      function
 ****************************************************************************/

using namespace std;

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

INLINE int GetChainCodeAAS(char aa)/*{{{*/
{
    int tint = 0;
    tint = aa - 'A';
    if (  (tint>=0) && (tint<=25)  )
    {
        return AAS_Code[tint];
    }
    else
    {
        return 20;
    }
}/*}}}*/
bool IsInCharSet(const char ch, const char *charSet, int n /*= 0 */)/*{{{*/
    /*****************************************************************************
     * check if the character "ch" is in charSet,
     ****************************************************************************/
{
    if(n == 0)
        n = strlen(charSet);
    int i;
    for(i = 0 ;i < n ; i ++)
    {
        if(ch == charSet[i])
            return true;
    }
    return false;
}/*}}}*/
int Read_databse(char *database_list, int NPer_Frag_Database, char *database_qijmatrix,char *database_modmatrix, char *database_facc,  char *CResult_location,  char *(*Namelistunit), char *(*CAASSequence),char *(*CShapesequence),char *(*CHSRsequence),int *(*blosum_chain),int *(*Chain_code_AAS), int *LengthList )/*{{{*/
{

    FILE *fp, *wfp, *fpread;
    char Camino, Cshp;
    char CSAS[10] = "";
    char OpenName[1000]; /*general string for filenames*/
    char string[1000];
    char CSubname[100]; /*for ID name*/
    char CTempAASseq[LONGEST_SEQ], CTempShapeseq[LONGEST_SEQ];
    char CTempDSSPseq[LONGEST_SEQ], CTempHSRseq[LONGEST_SEQ];
    int ik, ij, im, Sumunit, NSer_sub, NPer[21], HSR, NAAS, temp1;
    int Max_Length;
    int NSeries = 0;
    int Psi_blosum_Frag[LONGEST_SEQ][20], Psi_blosum[LONGEST_SEQ][20];
    int SUM_Component[LONGEST_SEQ];
    int Target_Naas[30], AAS_Code[30], TPer00;
    float param1,param2;
    char first_non_blank_char ; /*2009-07-15, Nanjiang, using first_non_blank_char to determine whether the matrix are valid data lines*/

    AAS_Code[0] = 0;//A--ALA
    AAS_Code[1] = 20;//B--Non
    AAS_Code[2] = 13;//C--CYS
    AAS_Code[3] = 18;//D--ASP
    AAS_Code[4] = 16;//E--GLU
    AAS_Code[5] = 5;//F--PHE
    AAS_Code[6] = 10;//G--GLY
    AAS_Code[7] = 9;//H--HIS
    AAS_Code[8] = 3;//I--ILE
    AAS_Code[9] = 20;//J--Non
    AAS_Code[10] = 7;//K--LYS
    AAS_Code[11] = 2;//L--LEU
    AAS_Code[12] = 6;//M--MET
    AAS_Code[13] = 15;//N--ASN
    AAS_Code[14] = 20;//O--Non
    AAS_Code[15] = 4;//P--PRO
    AAS_Code[16] = 19;//Q--GLN
    AAS_Code[17] = 8;//R--ARG
    AAS_Code[18] = 11;//S--SER
    AAS_Code[19] = 12;//T--THR
    AAS_Code[20] = 20;//U--Non
    AAS_Code[21] = 1;//V--VAL
    AAS_Code[22] = 17;//W--TRP
    AAS_Code[23] = 20;//X--Non
    AAS_Code[24] = 14;//Y--TYR
    AAS_Code[25] = 20;//Z--Non



    Max_Length = LONGEST_SEQ;
    Sumunit = 0;
    fp = fopen(database_list,"r");
    sprintf(OpenName,"%s/subname_list.txt",CResult_location);
    wfp = fopen(OpenName,"w");
    while (  fscanf(fp,"%d%s\n",&NSer_sub,CSubname) != EOF  )
    {
        //open frag database for profiles
        sprintf(OpenName,"%s/frag_%s.txt",database_facc,CSubname);
        fpread = fopen(OpenName,"r");
        for (ik=0; ik<LONGEST_SEQ; ik++)
        {
            for (ij=0; ij<20; ij++)
            {
                Psi_blosum_Frag[ik][ij] = 0;
            }
            SUM_Component[ik] = 0;
        }
        if (  fpread != NULL  )
        {
            while (  fgets(string,300,fpread) != NULL  )
            {
                sscanf(string, " %c", &first_non_blank_char);
                if (first_non_blank_char <'0' || first_non_blank_char > '9') /*only when the first_non_blank_char is digit */
                {
                    continue;
                }
                sscanf(string," %d ", &NSeries );
                //for (ij=0; ij<5; ij++)
                //{
                //str_char1[ij] = string[ij];
                //}
                //str_char1[ij] = '\0';
                //NSeries = atoi(str_char1);
                //if (   (NSeries<1) || (NSeries>=Max_Length)   )//for the first two lines
                //{
                //continue;
                //}
                NSeries = NSeries - 1;
                //percent
                sscanf(string,"%d %c %s %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%f%f",
                        &HSR, &Camino, CSAS,
                        &NPer[0], &NPer[1], &NPer[2], &NPer[3], &NPer[4],
                        &NPer[5], &NPer[6], &NPer[7], &NPer[8], &NPer[9],
                        &NPer[10], &NPer[11], &NPer[12], &NPer[13], &NPer[14],
                        &NPer[15], &NPer[16], &NPer[17], &NPer[18], &NPer[19],
                        &param1, &param2);
                for (ij=0; ij<20; ij++) {
                    Psi_blosum_Frag[NSeries][ij] =  NPer[ij];
                    SUM_Component[NSeries] = SUM_Component[NSeries] +
                        Psi_blosum_Frag[NSeries][ij];
                }
            }
            fclose(fpread);
        }
        //open Qij profiles for reading rows whose sum of components are zero
        sprintf(OpenName,"%s/Qijmatrix_%s.txt",database_qijmatrix,CSubname);
        fpread = fopen(OpenName,"r");
        if (  fpread == NULL  )
        {
            continue;
        }
        while (  fgets(string,300,fpread) != NULL  )
        {
            sscanf(string, " %c", &first_non_blank_char);
            if (first_non_blank_char <'0' || first_non_blank_char > '9') /*only when the first_non_blank_char is digit */
            {
                continue;
            }
            sscanf(string," %d ", &NSeries );
            //for (ij=0; ij<5; ij++)
            //{
            //str_char1[ij] = string[ij];
            //}
            //str_char1[ij] = '\0';
            //NSeries = atoi(str_char1);
            //if (   (NSeries<1) || (NSeries>=Max_Length)   )//for the first two lines
            //{
            //continue;
            //}
            NSeries = NSeries - 1;
            if (  SUM_Component[NSeries] >= 10  )
            {
                continue;
            }
            //percent
            sscanf(string,"%d %c %s %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%f%f",
                    &HSR, &Camino, CSAS, &NPer[0], &NPer[1], &NPer[2],
                    &NPer[3], &NPer[4], &NPer[5], &NPer[6], &NPer[7], &NPer[8],
                    &NPer[9], &NPer[10], &NPer[11], &NPer[12], &NPer[13],
                    &NPer[14], &NPer[15], &NPer[16], &NPer[17], &NPer[18],
                    &NPer[19], &param1, &param2);
            for (ij=0; ij<20; ij++) {
                Psi_blosum_Frag[NSeries][ij] =  NPer[ij];
            }
        }
        fclose(fpread);

        //open normal psi-blast profiles--modmatirx
        sprintf(OpenName,"%s/modmatrix_%s.txt",database_modmatrix,CSubname); 
        fpread = fopen(OpenName,"r");
        if (  fpread == NULL  )
        {
            continue;
        }
        for (ij=0; ij<LONGEST_SEQ; ij++)
        {
            CTempShapeseq[ij] = '-';
            CTempHSRseq[ij] = '-';
            CTempAASseq[ij] = '-';
            CTempDSSPseq[ij] = 'R';
        }

        while (  fgets(string,300,fpread) != NULL  )
        {
            sscanf(string, " %c", &first_non_blank_char);
            if (first_non_blank_char <'0' || first_non_blank_char > '9') /*only when the first_non_blank_char is digit */
            {
                continue;
            }
            sscanf(string," %d ", &NSeries );
            //for (ij=0; ij<5; ij++)
            //{
            //str_char1[ij] = string[ij];
            //}
            //str_char1[ij] = '\0';
            //NSeries = atoi(str_char1);
            //if (   (NSeries<1) || (NSeries>=Max_Length)   )//for the first two lines
            //{
            //continue;
            //}
            NSeries = NSeries - 1;
            //percent
            TPer00 = 0;
            sscanf(string,"%d %c %s %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%f%f",
                    &HSR,&Camino,CSAS,&NPer[0],&NPer[1],&NPer[2],&NPer[3],&NPer[4],&NPer[5],&NPer[6],&NPer[7],&NPer[8],&NPer[9],
                    &NPer[10],&NPer[11],&NPer[12],&NPer[13],&NPer[14],&NPer[15],&NPer[16],&NPer[17],&NPer[18],&NPer[19],&param1,&param2);

            Cshp= CSAS[0];
            for (ij=0; ij<20; ij++)
            {    
                Psi_blosum[NSeries][ij] =  NPer[ij]; /*Psi_blosum is storing mod matrix, 2009-07-06*/
                if (  NPer[ij] >= 1  )
                {
                    TPer00 = 1;
                }
            }
            //AAS
            if (  TPer00 == 0  )
            {
                NAAS = Camino - 'A';
                if (  (NAAS>=0) && (NAAS<25)  )
                {
                    NAAS = AAS_Code[NAAS];
                }
                else
                {
                    NAAS = 20;
                }
                if (  NAAS < 20  )
                {
                    Psi_blosum[NSeries][NAAS] = 100;
                }
            }
            //profile combination of Frag- and psi-blast
            if (  SUM_Component[NSeries] >= 10  )
            {
                temp1 = 0;
                for (ij=0; ij<20; ij++)
                {
                    Target_Naas[ij] = (100-NPer_Frag_Database)*Psi_blosum[NSeries][ij] + NPer_Frag_Database*Psi_blosum_Frag[NSeries][ij];
                    temp1 = temp1 + Target_Naas[ij];
                }
                for (ij=0; ij<20; ij++)
                {
                    Psi_blosum[NSeries][ij] = int ( Target_Naas[ij]*100.0/temp1 + 0.5  );
                }
            }
            else //if all components are zero in frag profiles, use the Qij profiles
            {
                for (ij=0; ij<20; ij++)
                {
                    Psi_blosum[NSeries][ij] = Psi_blosum_Frag[NSeries][ij];
                }
            }

            //amino acids
            CTempAASseq[NSeries] = Camino;
            //shape string
            CTempShapeseq[NSeries] = Cshp;
            //DSSP
            CTempDSSPseq[NSeries] = CSAS[2];
            //HSR
            if (  (CSAS[2]=='H') || (CSAS[2]=='G') || (CSAS[2]=='I')   )
                //if (  (CSAS[2]=='H')    )
            {
                CTempHSRseq[NSeries] = 'H';
            }
            //else if (   (CSAS[2]=='E') || (CSAS[2]=='B')  )
            else if (   (CSAS[2]=='E')  )
            {
                CTempHSRseq[NSeries] = 'S';
            }
            else
            {
                CTempHSRseq[NSeries] = 'R';
            }
            //CTempHSRseq[NSeries] = CSAS[2];
        }
        fclose(fpread);
        NSeries++;

        Namelistunit[Sumunit] = new char[6];
        CAASSequence[Sumunit] = new char[NSeries+3];
        CShapesequence[Sumunit] = new char[NSeries+3];
        CHSRsequence[Sumunit] = new char[NSeries+1];
        blosum_chain[Sumunit] = new int [NSeries*20+3];
        Chain_code_AAS[Sumunit] = new int [NSeries+3];

        for (ij=0; ij<NSeries; ij++)
        {
            CAASSequence[Sumunit][ij] = CTempAASseq[ij];
            NAAS = CTempAASseq[ij] - 'A';
            if (  (NAAS>=0) && (NAAS<=25)  )
            {
                NAAS = AAS_Code[NAAS];
            }
            else
            {
                NAAS = 20;
            }
            Chain_code_AAS[Sumunit][ij] = NAAS;
            CShapesequence[Sumunit][ij] = CTempShapeseq[ij];
            CHSRsequence[Sumunit][ij] = CTempHSRseq[ij];
            HSR = ij*20;
            for (im=0; im<20; im++)
            {
                blosum_chain[Sumunit][HSR+im] = Psi_blosum[ij][im];
            }
        }
        CAASSequence[Sumunit][ij] = '\0';
        CShapesequence[Sumunit][ij] = '\0';
        CHSRsequence[Sumunit][ij] = '\0';

        LengthList[Sumunit] = NSeries;
        HSR = strlen(CSubname);
        for (ij=0; ij<HSR; ij++)
        {
            Namelistunit[Sumunit][ij] = CSubname[ij];
        }
        Namelistunit[Sumunit][ij] = '\0';
        fprintf(wfp,"%4d %5s\n",Sumunit,Namelistunit[Sumunit]);
        Sumunit++;
    }
    fclose(fp);
    fclose(wfp);

    return 1;
}/*}}}*/
int Read_databse_SHU(int dbtype, char *database_list, int NPer_Frag_Database,
        int ratioScheme, char *database_qijmatrix,char *database_modmatrix,
        char *database_facc,  char *CResult_location,  char *(*Namelistunit),
        char *(*CAASSequence),char *(*CShapesequence),char *(*CHSRsequence),int
        *(*blosum_chain),int *(*Chain_code_AAS), int *LengthList ,int
        dsspMapMethod /*= 0*/, bool isReadBinaryFile /*= true*/, int8
        typeProfile  /*=1*/)/*{{{*/
    /*Modified by Nanjiang, using the function ReadInDatabase*/
{
    int linesize;
    int maxline = 300;
    Array1D <char> line_1darray(maxline+1);
    char *line = line_1darray.array1D;

    Array2D <int> matQij_2darray(LONGEST_SEQ, NUM_20_AA);
    Array2D <int> matFrag_2darray(LONGEST_SEQ, NUM_20_AA);
    Array2D <int> matMODM_2darray(LONGEST_SEQ, NUM_20_AA);
    Array2D <int> matMerged_2darray(LONGEST_SEQ, NUM_20_AA);
    int ** matQij = matQij_2darray.array2D; /*Qij profile*/
    int ** matFrag= matFrag_2darray.array2D; /*fragacc profile*/
    int ** matMODM = matMODM_2darray.array2D; /*modm profile, modm and Qij might be the same, 2008-01-18, Nanjiang*/
    int ** matMerged = matMerged_2darray.array2D;

    matQij_2darray.Init(0);
    matMODM_2darray.Init(0);
    matMerged_2darray.Init(0);

    Array1D <float> matFragScore1_1darray(LONGEST_SEQ);
    Array1D <float> matFragScore2_1darray(LONGEST_SEQ);
    Array1D <float> matQijScore1_1darray(LONGEST_SEQ);
    Array1D <float> matQijScore2_1darray(LONGEST_SEQ);
    Array1D <float> matMODMScore1_1darray(LONGEST_SEQ);
    Array1D <float> matMODMScore2_1darray(LONGEST_SEQ);
    float * matFragScore1 = matFragScore1_1darray.array1D;
    float * matFragScore2 = matFragScore2_1darray.array1D;
    float * matQijScore1 = matQijScore1_1darray.array1D;
    float * matQijScore2 = matQijScore2_1darray.array1D;
    float * matMODMScore1 = matMODMScore1_1darray.array1D;
    float * matMODMScore2 = matMODMScore2_1darray.array1D;

    Array1D <int> digitAASeq_1darray(LONGEST_SEQ);
    Array1D <char> aaSeq_1darray(LONGEST_SEQ+1);
    Array1D <char> shapeSeq_1darray(LONGEST_SEQ+1);
    Array1D <char> dsspSecSeq_1darray(LONGEST_SEQ+1);
    int *digitAASeq = digitAASeq_1darray.array1D;
    char *aaSeq = aaSeq_1darray.array1D;
    char *shapeSeq = shapeSeq_1darray.array1D;
    char *dsspSecSeq = dsspSecSeq_1darray.array1D;

    int qijformat =  QIJ_FORMAT_NANJIANG;
    int modmformat = MODM_FORMAT_NANJIANG;
    int fragaccformat = FRAGACC_FORMAT_NANJIANG;
    int type_dataset = TRAINING_SET;
    int mergeSide = TRAINING_SET;

    map <string, dbindex> dbindexfragacc;
    map <string, dbindex> dbindexqij;
    map <string, dbindex> dbindexmodm;
    int maxDBIndexNumber_fragacc = 0;
    int maxDBIndexNumber_qij = 0;
    int maxDBIndexNumber_modm = 0;
    vector <FILE*> fpList_fragacc;
    vector <FILE*> fpList_qij;
    vector <FILE*> fpList_modm;
    string filename;

    if (dbtype==1){ /*read index file of the dumped database, added 2011-10-13*/
        if (ReadDatabaseIndex(string(database_facc), dbindexfragacc, maxDBIndexNumber_fragacc) == -1 ){
            fprintf(stderr,"Read db fragacc failed for %s\n", database_facc);
            exit(1);
        }
        if (ReadDatabaseIndex(string(database_qijmatrix), dbindexqij, maxDBIndexNumber_qij) == -1 ){
            fprintf(stderr,"Read db qij failed for %s\n", database_qijmatrix);
            exit(1);
        }
        if (ReadDatabaseIndex(string(database_modmatrix), dbindexmodm, maxDBIndexNumber_modm) == -1 ){
            fprintf(stderr,"Read db modm failed for %s\n", database_modmatrix);
            exit(1);
        }
        GetDBFPList( fpList_fragacc, string(database_facc), maxDBIndexNumber_fragacc);
        GetDBFPList( fpList_qij, string(database_qijmatrix), maxDBIndexNumber_qij);
        GetDBFPList( fpList_modm, string(database_modmatrix), maxDBIndexNumber_modm);
    }


    FILE *fpTrainList = fopen(database_list,"r");
    string OpenName;
    //sprintf(OpenName,"%s/subname_list.txt",CResult_location);
    OpenName = string(CResult_location) + "/subname_list.txt";
    FILE *wfp = fopen(OpenName.c_str(),"w");
    int num = 0;
    int lengthSeq = 0;

    int cntID = 0;
    while((linesize = fgetline(fpTrainList, line ,maxline)) != EOF) //idList, main interation
    {
        if(linesize <= 0) { continue; }
        char *id = new char[linesize+1];
        sscanf(line, "%d %s", &num, id);
        if (dbtype==0){
            lengthSeq = ReadInDatabase(id, database_facc, database_qijmatrix,database_modmatrix, qijformat, modmformat, fragaccformat, matFrag, matQij, matMODM, matMerged, matFragScore1,matFragScore2, matQijScore1,matQijScore2,matMODMScore1,matMODMScore2, aaSeq, digitAASeq, shapeSeq, dsspSecSeq, type_dataset, isReadBinaryFile, NPer_Frag_Database, ratioScheme, mergeSide, dsspMapMethod, typeProfile);
        } else if (dbtype == 1){ /*added 2011-10-17 */
            lengthSeq = ReadInDatabase_dumpedfile(id, dbindexfragacc, dbindexqij, dbindexmodm, fpList_fragacc, fpList_qij, fpList_modm, matFrag, matQij, matMODM, matMerged, matFragScore1,matFragScore2, matQijScore1,matQijScore2,matMODMScore1,matMODMScore2, aaSeq, digitAASeq, shapeSeq, dsspSecSeq, type_dataset, isReadBinaryFile, NPer_Frag_Database, ratioScheme, mergeSide, dsspMapMethod, typeProfile);
        } else {
            fprintf(stderr,"dbtype = %d has not been implemented yet. Exit.\n", dbtype);
            exit(1);
        }
        if (lengthSeq <= 0) {
            fprintf(stderr,"%s: ReadInDatabase error, lengthSeq = %d\n", id, lengthSeq);
            continue;
        }

        int sizeid = strlen(id);
        Namelistunit[cntID] = new char[sizeid+1];
        CAASSequence[cntID] = new char[lengthSeq+3];
        CShapesequence[cntID] = new char[lengthSeq+3];
        CHSRsequence[cntID] = new char[lengthSeq+3];
        blosum_chain[cntID] = new int [lengthSeq*20+3];
        Chain_code_AAS[cntID] = new int [lengthSeq+3];
        /*assign the values to array in the argument list*/
        strcpy(Namelistunit[cntID], id);
        strcpy(CAASSequence[cntID], aaSeq);
        strcpy(CShapesequence[cntID], shapeSeq);
        strcpy(CHSRsequence[cntID], dsspSecSeq);

        int i,j;
        int begpos;
        for (i = 0; i < lengthSeq; i ++)
        {
            begpos = i*NUM_20_AA;
            for (j = 0 ; j < NUM_20_AA; j ++) {
                blosum_chain[cntID][begpos+j] = matMerged[i][j];
            }
        }

        for (i = 0; i < lengthSeq; i ++) {
            Chain_code_AAS[cntID][i] = GetChainCodeAAS(aaSeq[i]);
        }
        LengthList[cntID] = lengthSeq;
        cntID ++;
        /*free memory*/
        delete [] id;
    }
    if (dbtype==1){
        int i;
        for (i=0;i<=maxDBIndexNumber_fragacc;i++){ fclose(fpList_fragacc[i]); }
        for (i=0;i<=maxDBIndexNumber_qij;i++){ fclose(fpList_qij[i]); }
        for (i=0;i<=maxDBIndexNumber_modm;i++){ fclose(fpList_modm[i]); }
    }
    fclose(fpTrainList);
    fclose(wfp);

    return cntID; /*return the number of chains actually read in*/
}/*}}}*/

int Treat_pair_Percent_homology(float *Datapair, int Lenlong, int Lenshort, int *Max_sum, int *(*Candidate_Maxline), int Num_Candhomo)/*{{{*/
{
    /*  Input:
     *      Datapair: including all the hits for two sequences. It is a
     *                matrix(Lenshort, Lenlong), but in a one-dimension variable
     *      Lenlong : length of the query sequence
     *      Lenshort: length of the target sequence
     *      Max_sum: The set of returned data;
     Max_sum[0]: the maximum Consecutive hits;
     Max_sum[1]: the total effective hits (by taking away the sporatic points);
     Max_sum[2]: Maximum points for finally optimizide, the returned number;
     Max_sum[3]: The number of matched fragments;
     Max_sum[4]: Max_point_single_line;
     Max_sum[5]: Max_point_single_line_id;
     *  Candidate_Maxline: the data from source file, matrix(target_length, Number_candidate)
     *  Num_Candhomo: The number of candidates for homology analysis
     */
    int Max_point_line, *Max_Line_X, *Max_Line_Y, Sum_seq2;

    int Total_Len, *Total_lines,  ik,ij, im,  Max_Numbere_EffPair, Max_Numbere_Line;
    int beg, end, Lbeg, Lend, Lsegment, NFragment, Npos, Npos_cur;
    int Turn_long_short;
    int Max_point_single_line = 0; 
    int Max_point_single_line_id= 0;
    int Nhomopos1 = 0;
    int Nhomopos2 = 0;
    int *Max_Line_X_point,*Max_Line_Y_point, LenX, LenY, datareturn[5], Max_pos, ntemp, Nper ;
    int Max_Consecutive, Len_consecutive, PrevX, PrevY, All_data;
    float *copy_Datapair;


    Sum_seq2 = Lenlong + Lenshort +2;
    Max_Line_X = new int[Sum_seq2];
    Max_Line_Y = new int[Sum_seq2];

    for (ik=0; ik<Sum_seq2; ik++)
    {
        Max_Line_X[ik] = -1;
        Max_Line_Y[ik] = -1;
    }

    Max_Line_X_point = new int [Lenlong + 2];
    Max_Line_Y_point = new int [Lenshort +2];
    for (ik=0; ik<Lenlong+2; ik++)
    {
        Max_Line_X_point[ik] = 0;
    }
    for (ij=0; ij<Lenshort+2; ij++)
    {
        Max_Line_Y_point[ij] = 0;
    }

    All_data = Lenlong*Lenshort + 1;
    copy_Datapair = new float[All_data];
    for (ij=0; ij<All_data; ij++)
    {
        copy_Datapair[ij] = 0;
    }

    //Max_point_line = Treat_pair_to_optimisized_Line(Datapair, Lenlong, Lenshort, Max_Line_X, Max_Line_Y);
    //Max_sum[1] = Max_point_line;

    Max_Numbere_EffPair = 0;
    Max_point_single_line = 0;
    NFragment = 0;
    Turn_long_short = Lenlong - Lenshort;

    //set the numbr of points on one doagonal
    Total_Len = Lenlong + Lenshort + 2;
    Total_lines = new int[Total_Len];
    for (ik=0; ik<Total_Len ; ik++)
    {
        Total_lines[ik] = 0;
    }


    /*
    //copy data from datapair to copy_datapair/
    //suppose: X-long sequence, Y--short sequence
    //upp-left
    beg = 0;
    end = Lenshort;
    for (ik=beg; ik<end; ik++)
    {
    Lend = Lenshort - ik;
    Lsegment = 0;
    for (ij=0; ij<Lend; ij++)
    {
    Npos = ij*Lenshort + ik + ij;
    if (  Datapair[Npos] >= 0  )
    {
    Lsegment++;
    }
    else
    {
    if (  Lsegment >= 3  )
    {
    for (im=ij-Lsegment; im<ij; im++)
    {
    Npos = im*Lenshort + ik + im;
    copy_Datapair[Npos] = Datapair[Npos];
    }
    }
    Lsegment = 0;
    }
    }
    }

    //from the begining of the longer one to the end of the longer one, skip last 10;
    Lbeg = 1;
    Lend = Lenlong;
    for (ik=Lbeg; ik<Lend; ik++)
    {
    if (  ik <= Turn_long_short  )
    {
    beg = 0;
    end = Lenshort;
    }
    else
    {
    beg = 0;
    end = Lenshort - (ik - Turn_long_short);
    }
    Lsegment = 0;

    for (ij=beg; ij<end; ij++)
    {
    Npos = (ik+ij)*Lenshort + ij;
    if (  Datapair[Npos] >= 0  )
    {
    Lsegment++;
    }
    else
    {
    if (  Lsegment >= 3  )
    {
    for (im=ij-Lsegment; im<ij; im++)
    {
    Npos = (ik+im)*Lenshort + im;
    copy_Datapair[Npos] = Datapair[Npos] ;
    }
    }
    Lsegment = 0;
    }
    }
    }
    */



    //*******************
    int Nsegment_line, Ntemp1, Ntemp2,pos1,pos2;
    float xxx;
    //*******************

    //suppose: X-long sequence, Y--short sequence
    //upp-left
    beg = 0;
    end = Lenshort;
    for (ik=beg; ik<end; ik++)
    {
        Lend = Lenshort - ik;
        Max_Numbere_Line = 0;
        Lsegment = 0;
        //*******************
        for (ij=0; ij<Lend; ij++)
        {
            copy_Datapair[ij] = 0;
        }
        Nsegment_line = 0;
        //*******************
        for (ij=0; ij<Lend; ij++)
        {
            Npos = ij*Lenshort + ik + ij;
            if (  Datapair[Npos] >= 0  )
            {
                Lsegment++;
            }
            else
            {
                if (  Lsegment >= 3  )
                {
                    for (im=ij-Lsegment; im<ij; im++)
                    {
                        Npos = im*Lenshort + ik + im;
                        ntemp = int ( Datapair[Npos] + 0.1 );
                        Nhomopos1 = ntemp/Num_Candhomo;
                        Nhomopos2 = max(0, ntemp - Nhomopos1*Num_Candhomo);
                        Candidate_Maxline[Nhomopos1][Nhomopos2] = Lsegment;
                        //*******************
                        copy_Datapair[im] = float (Nsegment_line*1000) + float (Lsegment);
                        //*******************
                    }
                    //*******************
                    Nsegment_line++;
                    //*******************
                    Max_Numbere_Line = Max_Numbere_Line + Lsegment;
                    NFragment++;
                }
                else if (  Lsegment >= 1  )
                {
                    for (im=ij-Lsegment; im<ij; im++)
                    {
                        Npos = im*Lenshort + ik + im;
                        Datapair[Npos] = -2;
                    }
                }
                Lsegment = 0;
            }
        }
        Max_Numbere_EffPair = Max_Numbere_EffPair + Max_Numbere_Line;
        Total_lines[Lenshort-ik-1] = Max_Numbere_Line;
        if (  Max_Numbere_Line > Max_point_single_line  )
        {
            Max_point_single_line = Max_Numbere_Line;
            Max_point_single_line_id = Lenshort-ik-1;
        }
        //******************* set on same line
        if (  Nsegment_line >= 2  )
        {
            for (ij=0; ij<Lend; ij++)
            { 
                if (  copy_Datapair[ij] < 2  )
                { 
                    continue;
                }
                Ntemp1 = int (copy_Datapair[ij]/1000);
                Npos = ij*Lenshort + ik + ij;
                ntemp = int ( Datapair[Npos] + 0.1 );
                pos1 = ntemp/Num_Candhomo;
                pos2 = ntemp - pos1*Num_Candhomo;
                if (Datapair[Npos]<0)
                {
                    printf("1-1 %4d %4d  %4d %4d  %f  %4d %4d\n",ik,ij, Lenshort, Lenlong,Datapair[Npos], pos1,pos2);
                }

                xxx = 0;
                for (im=0; im<Lend; im++)//
                { 
                    if (  copy_Datapair[im] < 2  )
                    { 
                        continue;
                    } 
                    Ntemp2 = int (copy_Datapair[im]/1000);
                    if (  Ntemp1 == Ntemp2  )
                    {
                        continue;
                    }
                    xxx = xxx + float (copy_Datapair[im]*0.1);
                    if (Datapair[Npos]<0)
                    {
                        printf("1-2 %4d %4d  %4d %4d  %f  %4d %4d\n",ik,ij, Lenshort, Lenlong,Datapair[Npos], Nhomopos1,Nhomopos2);
                    }
                } 
                Candidate_Maxline[pos1][pos2] = Candidate_Maxline[pos1][pos2] + int (xxx);

            } 
        }
        //*******************
    }


    //from the begining of the longer one to the end of the longer one, skip last 10;
    Lbeg = 1;
    Lend = Lenlong;
    for (ik=Lbeg; ik<Lend; ik++)
    {
        if (  ik <= Turn_long_short  )
        {
            beg = 0;
            end = Lenshort;
        }
        else
        {
            beg = 0;
            end = Lenshort - (ik - Turn_long_short);
        }
        Lsegment = 0;
        Max_Numbere_Line = 0;
        //*******************
        for (ij=beg; ij<end; ij++)
        {
            copy_Datapair[ij] = 0;
        }
        Nsegment_line = 0;
        //*******************

        for (ij=beg; ij<end; ij++)
        {
            Npos = (ik+ij)*Lenshort + ij;
            if (  Datapair[Npos] >= 0  )
            {
                Lsegment++;
            }
            else
            {
                if (  Lsegment >= 3  )
                {
                    for (im=ij-Lsegment; im<ij; im++)
                    {
                        Npos = (ik+im)*Lenshort + im;
                        ntemp = int ( Datapair[Npos] + 0.1 );
                        Nhomopos1 = ntemp/Num_Candhomo;
                        Nhomopos2 = max(0, ntemp - Nhomopos1*Num_Candhomo);
                        Candidate_Maxline[Nhomopos1][Nhomopos2] = Lsegment;
                        //*******************
                        copy_Datapair[im] = float (Nsegment_line*1000) + float (Lsegment);
                        //*******************

                    }

                    //*******************
                    Nsegment_line++;
                    //*******************
                    Max_Numbere_Line = Max_Numbere_Line + Lsegment;
                    NFragment++;
                }
                else if (  Lsegment >= 1  )
                {
                    for (im=ij-Lsegment; im<ij; im++)
                    {
                        Npos = (ik+im)*Lenshort + im;
                        Datapair[Npos] = -2;
                    }
                }
                Lsegment = 0;
            }
        }
        Max_Numbere_EffPair = Max_Numbere_EffPair + Max_Numbere_Line;
        Total_lines[Lenshort+ik-1] = Max_Numbere_Line;
        if (  Max_Numbere_Line > Max_point_single_line  )
        {
            Max_point_single_line = Max_Numbere_Line;
            Max_point_single_line_id = Lenshort+ik-1;
        }
        //******************* set on same line
        if (  Nsegment_line >= 2  )
        {
            for (ij=beg; ij<end; ij++)
            { 
                if (  copy_Datapair[ij] < 2  )
                { 
                    continue;
                }
                Ntemp1 = int (copy_Datapair[ij]/1000);
                Npos = (ik+ij)*Lenshort + ij;
                ntemp = int ( Datapair[Npos] + 0.1 );
                pos1 = ntemp/Num_Candhomo;
                pos2 = ntemp - pos1*Num_Candhomo;
                if (Datapair[Npos]<0)
                {
                    printf("2-1 %4d %4d  %4d %4d %f  %4d %4d\n",ik,ij, Lenshort, Lenlong,Datapair[Npos],pos1,pos2);
                }

                xxx = 0;
                for (im=beg; im<end; im++)
                { 
                    if (  copy_Datapair[im] < 2  )
                    { 
                        continue;
                    } 
                    Ntemp2 = int (copy_Datapair[im]/1000);
                    if (  Ntemp1 == Ntemp2  )
                    {
                        continue;
                    }
                    xxx = xxx + float( copy_Datapair[im]*0.1 );
                    if (Datapair[Npos]<0)
                    {
                        printf("2-2 %4d %4d  %4d %4d %f  %4d %4d\n",ik,ij, Lenshort, Lenlong,Datapair[Npos],Nhomopos1,Nhomopos2);
                    }
                } 
                Candidate_Maxline[pos1][pos2] = Candidate_Maxline[pos1][pos2] + int (xxx);

            } 
        }
        //*******************

    }
    /*Get the index of the diagonal line with maximum number of points on the
     * dot plot*/
    Max_point_line = Treat_pair_to_optimisized_Line(Datapair, Lenlong, Lenshort, Max_Line_X, Max_Line_Y, datareturn, Candidate_Maxline, Num_Candhomo);

    //max consecutive
    Max_Consecutive = 1;
    Len_consecutive = 1;
    PrevX = Max_Line_X[0];
    PrevY = Max_Line_Y[0];

    Max_pos = datareturn[1];
    for (ik=0; ik<Max_pos; ik++)
    {
        if  (  (Max_Line_X[ik]<0) || (Max_Line_Y[ik]<0)  )
        {
            printf("(Max_Line_X[ik]<0) || (Max_Line_Y[ik]<0)\n");
        }
        Npos = Max_Line_X[ik]*Lenshort + Max_Line_Y[ik] ;
        if (  Datapair[Npos] >= 0 )
        {

            Max_Line_X_point[Max_Line_X[ik]] = 1;
            Max_Line_Y_point[Max_Line_Y[ik]] = 1;
            if  (  (abs(Max_Line_X[ik]-PrevX)==1) &&  (abs(Max_Line_Y[ik]-PrevY)==1)   )
            {
                Len_consecutive++;
            }
            else
            {
                if (  Len_consecutive > Max_Consecutive  )
                {
                    Max_Consecutive = Len_consecutive;
                }
                Len_consecutive = 1;
            }
        }
        else
        {
            if (  Len_consecutive > Max_Consecutive  )
            {
                Max_Consecutive = Len_consecutive;
            }
            Len_consecutive = 1;
        }
        PrevX = Max_Line_X[ik];
        PrevY = Max_Line_Y[ik];

    }

    LenX = 0;
    for (ik=0; ik<Lenlong; ik++)
    {
        if (  Max_Line_X_point[ik] == 1  )
        {
            LenX++;
        }
    }

    LenY = 0;
    for (ij=0; ij<Lenshort; ij++)
    {
        if (  Max_Line_Y_point[ij] == 1  )
        {
            LenY++;
        }
    }
    Max_point_line = min(LenX,LenY);





    Nper = Max_point_line*100*2/(Lenshort+Lenlong);
    if (  Nper >= 8  )
    {
        //set the max line : 2 away
        for (ik=0; ik<Max_pos; ik++)
        { 
            Npos = Max_Line_X[ik]*Lenshort + Max_Line_Y[ik] ;
            if (  Datapair[Npos] >= 0 )//current point
            {
                ntemp = int ( Datapair[Npos] + 0.1 );
                Nhomopos1 = ntemp/Num_Candhomo;
                Nhomopos2 = max(0, ntemp - Nhomopos1*Num_Candhomo);
                if (  Candidate_Maxline[Nhomopos1][Nhomopos2] >= 1 )
                {
                    Candidate_Maxline[Nhomopos1][Nhomopos2] = Candidate_Maxline[Nhomopos1][Nhomopos2]*1 + Nper;
                }
            }
            // 2 more down
            beg = max(0, Max_Line_Y[ik]-2);
            for (ij=beg; ij<Max_Line_Y[ik]; ij++ )
            {  
                Npos_cur = ik + Max_Line_Y[ik] - ij;
                if  (   (Npos_cur>=Max_pos) ||  (Max_Line_X[Npos_cur]==Max_Line_X[ik])  )
                {  
                    continue;
                }  
                Npos = Max_Line_X[ik]*Lenshort + ij;
                if (  Datapair[Npos] >= 0 )
                {  
                    ntemp = int ( Datapair[Npos] + 0.1 );
                    Nhomopos1 = ntemp/Num_Candhomo;
                    Nhomopos2 = max(0, ntemp - Nhomopos1*Num_Candhomo);
                    if (  Candidate_Maxline[Nhomopos1][Nhomopos2] >= 1 )
                    {   
                        Candidate_Maxline[Nhomopos1][Nhomopos2] = Candidate_Maxline[Nhomopos1][Nhomopos2] + Nper;
                    }   
                }  
            }  
            // 2 more up
            beg = Max_Line_Y[ik]+1;
            end = min(Max_Line_Y[ik]+3, Lenshort);
            for (ij=beg; ij<end; ij++ )
            {  
                Npos_cur = ik + ij - Max_Line_Y[ik];
                if  (   (Npos_cur>=Max_pos) ||  (Max_Line_X[Npos_cur]==Max_Line_X[ik])  )
                {  
                    continue;
                }  
                Npos = Max_Line_X[ik]*Lenshort + ij;
                if (  Datapair[Npos] >= 0 )
                {  
                    ntemp = int ( Datapair[Npos] + 0.1 );
                    Nhomopos1 = ntemp/Num_Candhomo;
                    Nhomopos2 = max(0, ntemp - Nhomopos1*Num_Candhomo);
                    if (  Candidate_Maxline[Nhomopos1][Nhomopos2] >= 1 ) {   
                        Candidate_Maxline[Nhomopos1][Nhomopos2] = Candidate_Maxline[Nhomopos1][Nhomopos2] +  Nper;
                    }
                }
            }
        }
    }







    Max_sum[0] = Max_Consecutive;
    Max_sum[1] = Max_Numbere_EffPair;
    //only fragments left 

    Max_sum[2] = Max_point_line;
    Max_sum[3] = NFragment;
    Max_sum[4] = Max_point_single_line;
    Max_sum[5] = Max_point_single_line_id;


    delete [] Max_Line_X;
    delete [] Max_Line_Y;
    delete [] Max_Line_X_point;
    delete [] Max_Line_Y_point;
    delete [] Total_lines;
    delete [] copy_Datapair;

    return Max_point_line;
}/*}}}*/


int Treat_pair_to_optimisized_Line(float *Datapair, int Lenlong, int Lenshort, int *Max_Line_X, int *Max_Line_Y, int *datareturn, int *(*Candidate_Maxline), int Num_Candhomo )/*{{{*/
    /*
     * Description: 
     *      Scan the dot plot (of high scoring fragment pairs between the query
     *      and target sequence) to find the index of the diagonal line with maximum
     *      number of points
     * Input:
     *      Datapair: array of scores for all point on the dot plot
     *      lenlong : length of the query sequence
     *      lenshort: length of the target sequence
     *      Num_Candhomo: number of high scoring fragment to be used for each
     *                    query fragment (default: 100)
     * */
{
    float xgap, Xpoint,Xpoint_Halv, xline_score, xscore, *Score_table, score_diagonal, score_down, score_left, val_big;
    int Ntotal_Point,Max_Point_single_line, ik, ij, Npos, Npos_cur, Max_pos;
    float Max_Score;
    int beg, Max_Score_pos_long, Max_Score_pos_short, Point_long, Point_short, SUM_point;
    int Nhomopos1, Nhomopos2, Ntemp;

    Npos = 0;

    xgap = float (-0.1);
    Xpoint = 1;
    Xpoint_Halv = Xpoint/2;
    xline_score = float (0.1);
    Ntotal_Point = Lenlong*Lenshort + 1;
    Score_table = new float [Ntotal_Point];


    //fill the score table Score_table[]
    //fill the first column
    for (ik=0; ik<Lenshort; ik++) {/*{{{*/
        Npos = ik;
        if (  Datapair[Npos] >= 0 ) {
            Ntemp = int ( Datapair[Npos] + 0.1 );
            Nhomopos1 = Ntemp/Num_Candhomo;
            Nhomopos2 = max(0, Ntemp - Nhomopos1*Num_Candhomo);
            xscore = Candidate_Maxline[Nhomopos1][Nhomopos2]*xline_score;
            xscore = xscore + Xpoint_Halv;
        } else {
            xscore = 0;
        }
        Score_table[Npos] = xscore;
    }/*}}}*/

    //fill the first row
    for (ij=0; ij<Lenlong; ij++) {/*{{{*/
        Npos = ij*Lenshort;
        if (  Datapair[Npos] >= 0 ) {
            Ntemp = int ( Datapair[Npos] + 0.1 );
            Nhomopos1 = Ntemp/Num_Candhomo;
            Nhomopos2 = max(0, Ntemp - Nhomopos1*Num_Candhomo);
            xscore = Candidate_Maxline[Nhomopos1][Nhomopos2]*xline_score;
            xscore = xscore + Xpoint_Halv;
        } else {
            xscore = 0;
        }
        Score_table[Npos] = xscore;
    }/*}}}*/

    //fill the other places in table
    for (ik=1; ik<Lenshort; ik++) {
        for (ij=1; ij<Lenlong; ij++)
        {
            // diagonal score
            Npos_cur = ij*Lenshort + ik ;
            if (  Datapair[Npos_cur] >= 0 ) {
                Ntemp = int ( Datapair[Npos] + 0.1 );
                Nhomopos1 = Ntemp/Num_Candhomo;
                Nhomopos2 = max(0, Ntemp - Nhomopos1*Num_Candhomo);
                xscore = Candidate_Maxline[Nhomopos1][Nhomopos2]*xline_score;
                xscore = xscore + Xpoint;
            } else { 
                xscore = 0;
            }
            Npos = (ij-1)*Lenshort + ik-1 ;
            score_diagonal = Score_table[Npos] + xscore;
            // left score
            Npos = (ij-1)*Lenshort + ik ;
            if  (   (Datapair[Npos]<0) && (Datapair[Npos_cur]>=0)   )
            {
                Ntemp = int ( Datapair[Npos] + 0.1 );
                Nhomopos1 = Ntemp/Num_Candhomo;
                Nhomopos2 = max(0, Ntemp - Nhomopos1*Num_Candhomo);
                xscore = Candidate_Maxline[Nhomopos1][Nhomopos2]*xline_score;
                xscore = xscore + Xpoint_Halv;
            } else {
                xscore = xgap;
            }
            score_left = Score_table[Npos] + xscore;
            // down score
            Npos = ij*Lenshort + ik - 1 ;
            if  (   (Datapair[Npos]<0) && (Datapair[Npos_cur]>=0)   ) {
                Ntemp = int ( Datapair[Npos] + 0.1 );
                Nhomopos1 = Ntemp/Num_Candhomo;
                Nhomopos2 = max(0, Ntemp - Nhomopos1*Num_Candhomo);
                xscore = Candidate_Maxline[Nhomopos1][Nhomopos2]*xline_score;
                xscore = xscore + Xpoint_Halv;
            } else {
                xscore = xgap;
            }
            score_down = Score_table[Npos] + xscore;
            //decide which one is biggest
            val_big = max(score_down, score_left);
            Score_table[Npos_cur] = max(val_big, score_diagonal);
        } 
    }
    //search the maxium pass, the top line (axis) can be divided into two parts. For top line (long sequence), first part-- upper-left triangle,
    //second part--top line of square;  For right line (short sequence), there is only one, that is right-down triangle
    //search the squares (  Lenshort X Lenshort  ), on second part of top line, long sequence

    //for second part of long sequence
    Max_Score = 0;
    Max_Score_pos_long = Lenshort -1;//from left to right for long sequence, then from top to low for short sequence
    Max_Score_pos_short = Lenshort -1;//from left to right for long sequence, then from top to low for short sequence
    beg = Lenshort -1;
    for (ik=beg; ik<Lenlong; ik++)
    {  
        Npos = ik*Lenshort + Lenshort - 1;
        if  (  Score_table[Npos] > Max_Score  )
        {
            Max_Score = Score_table[Npos];
            Max_Score_pos_long = ik;
            Max_Score_pos_short = Lenshort -1;
        }
    }
    //for the first part of long sequence
    beg = Lenshort -2;
    for (ik=beg; ik>=10; ik--)
    {  
        Npos = ik*Lenshort + Lenshort - 1;
        if  (  Score_table[Npos] > Max_Score  )
        {
            Max_Score = Score_table[Npos];
            Max_Score_pos_long = ik;
            Max_Score_pos_short = Lenshort -1;
        }
    }
    //for short sequence
    beg = Lenlong - 1;
    for (ik=Lenshort -1; ik>=10; ik--)
    {  
        Npos = beg*Lenshort + ik;
        if  (  Score_table[Npos] > Max_Score  )
        {
            Max_Score = Score_table[Npos];
            Max_Score_pos_long = Lenlong - 1;
            Max_Score_pos_short = ik;
        }
    }

    //trace back from the biggest score point, the "Max_Score_pos" with "Max_Score"
    SUM_point = int (Max_Score + 0.4999);
    Max_Point_single_line = 0;
    Max_pos = 0;
    if (  Max_Score >= 5  )
    {
        Point_long = Max_Score_pos_long;
        Point_short = Max_Score_pos_short;
        Max_pos = 0;
        Max_Line_X[Max_pos] = Point_long;
        Max_Line_Y[Max_pos] = Point_short;

        Max_Point_single_line = 0;
        while (  (Point_long>=1) && (Point_short>=1)  )
        {
            //count point
            Npos_cur = Point_long*Lenshort + Point_short ;
            //extract the path with maximium points
            // diagonal score
            Npos = (Point_long-1)*Lenshort + Point_short-1 ;
            score_diagonal = Score_table[Npos];
            // left score
            Npos = (Point_long-1)*Lenshort + Point_short ;
            score_left = Score_table[Npos];
            // down score
            Npos = Point_long*Lenshort + Point_short - 1 ;
            score_down = Score_table[Npos];
            //decide which one is biggest
            if (  (score_diagonal>=score_left) && (score_diagonal>=score_down)  )
            {
                Point_long--;
                Point_short--;
                if (  Datapair[Npos_cur] >= 0  )
                {   
                    Max_Point_single_line++;
                }   
            }
            else if (  (score_left>=score_diagonal) && (score_left>=score_down)  )
            {
                Npos = (Point_long-1)*Lenshort + Point_short ;
                if  (  (Datapair[Npos_cur]>=0) && (Datapair[Npos]<0)  )
                {   
                    Max_Point_single_line++;
                }   
                Point_long--;
            }
            else//if (  (score_down>=score_diagonal) && (score_down>=score_left)  )
            {
                Npos = Point_long*Lenshort + Point_short - 1 ;
                if  (  (Datapair[Npos_cur]>=0) && (Datapair[Npos]<0)  )
                {   
                    Max_Point_single_line++;
                }   
                Point_short--;
            }
            Max_pos++;
            Max_Line_X[Max_pos] = Point_long;
            Max_Line_Y[Max_pos] = Point_short;
        }
    }
    else//too few points
    {
        Max_Point_single_line = SUM_point;
    }
    if (  Max_Point_single_line != SUM_point  )
    {
        //printf("Max_Point_single_line(%d) != SUM_point(%d)",Max_Point_single_line,SUM_point);
    }
    datareturn[0] = Max_Point_single_line;
    datareturn[1] = Max_pos;

    delete [] Score_table;
    return Max_Point_single_line;
}/*}}}*/

int ReadInPairBinary( const char *binfragPairFile,  float *(*fragScore), char *(*proIDCan), int datasize, int *(*subikk), int *(*sublen))/*{{{*/
    /*****************************************************************************
     * Read in pairs from binary frag files
     * 2009-06-25
     ****************************************************************************/
{
    /*check whether the input file is a vaid fragbinary file*/
    int lenFileName = strlen(binfragPairFile);
    if(strcmp (&(binfragPairFile[lenFileName-3]), "bin" ) != 0)
    {
        fprintf(stderr,"Error, the binfragfile \"%s\" is not end with \"bin\"\n", binfragPairFile);
        return -1;
    }

    int numID = 0;
    int maxSizeID = 0;
    int length = 0;
    int totalFragCan = 0;
    int fragformat, ncan, index;

    ReadInFragParameter(binfragPairFile, fragformat, numID, maxSizeID, length, totalFragCan);
    Array2D <char> idList_2darray(numID, maxSizeID+1);
    idList_2darray.Init('\0');
    char **idList = idList_2darray.array2D;

    Array1D <short> tmpposTar_1darray(length);
    Array1D <short> tmpnumCan_1darray(length);
    tmpposTar_1darray.Init(0);
    tmpnumCan_1darray.Init(0);
    short *tmpposTar = tmpposTar_1darray.array1D;
    short *tmpnumCan = tmpnumCan_1darray.array1D;
    FragCanShort5 *tmpfragCan5 = NULL;
    FragCanShort6 *tmpfragCan6 = NULL;
    if (fragformat == 5)
    { tmpfragCan5 = new FragCanShort5[totalFragCan+1]; }
    else 
    { tmpfragCan6 = new FragCanShort6[totalFragCan+1]; }

    if (fragformat == 5)
    {
        GetBinaryFragShort5(binfragPairFile, fragformat, idList, numID, maxSizeID, length, tmpposTar, tmpnumCan, totalFragCan, tmpfragCan5);
    }
    else
    {
        GetBinaryFragShort6(binfragPairFile, fragformat, idList, numID, maxSizeID, length, tmpposTar, tmpnumCan, totalFragCan, tmpfragCan6);
    }
    int cntPair = 0;
    int i,j;
    for (i = 0; i < length; i ++)
    {
        ncan = min(tmpnumCan[i], short(100));//in case more than 100
        for (j = 0; j <ncan ; j ++)
        {
            if (fragformat == 5)
            {
                index = Quick_Find_str_whole(proIDCan, datasize, idList[tmpfragCan5[cntPair].idxInner]);
                if (  index >= 0  )
                {
                    subikk[tmpposTar[i]][j] = index;
                    sublen[tmpposTar[i]][j] = tmpfragCan5[cntPair].posCan;
                    fragScore[tmpposTar[i]][j] = tmpfragCan5[cntPair].score;
                }
            }
            else
            {
                subikk[tmpposTar[i]][j] = tmpfragCan6[cntPair].idxOuter;
                sublen[tmpposTar[i]][j] = tmpfragCan6[cntPair].posCan;
                fragScore[tmpposTar[i]][j] = tmpfragCan6[cntPair].score;
            }
            cntPair ++;
        }
    }

    /*free memory*/
    if (fragformat == 5)
    { delete [] tmpfragCan5; }
    else 
    { delete [] tmpfragCan6; }

    return cntPair;
}/*}}}*/

int ReadInFragParameter(const char *file, int &fragformat, int &numID, int &maxSizeID, int &length, int &totalFragCan)/*{{{*/
    /* read in the fragformat
     * totalFragCan: the total number of candidate fragments for all target positions
     * */
{
    FILE *fpin = fopen(file,"rb");
    if (  fpin == NULL  )
    {
        return -1;
    }
    size_t nread = 0;/*return value for fread*/
    /*read the fragFileType*/
    nread=fread(&fragformat, sizeof(int), 1, fpin); /*read the fragFileType*/ 
    if (nread != 1)
    { fprintf(stderr,"fread error! File:%s, function:%s, var:%s\n",file,"ReadInFragParameter","fragformat"); }
    nread=fread(&numID, sizeof(int), 1, fpin); /*read the number of unique ids*/
    if (nread != 1)
    { fprintf(stderr,"fread error! File:%s, function:%s, var:%s\n",file,"ReadInFragParameter","numID"); }
    nread=fread(&maxSizeID, sizeof(int),1,fpin);/*read the maximal size of ids*/
    if (nread != 1)
    { fprintf(stderr,"fread error! File:%s, function:%s, var:%s\n",file,"ReadInFragParameter","maxSizeID"); }
    nread=fread(&length, sizeof(int), 1, fpin);  /*read length of the target sequence*/ 
    if (nread != 1)
    { fprintf(stderr,"fread error! File:%s, function:%s, var:%s\n",file,"ReadInFragParameter","length"); }
    nread=fread(&totalFragCan, sizeof(int), 1, fpin);   /*read the total number of candidate fragments*/ 
    if (nread != 1)
    { fprintf(stderr,"fread error! File:%s, function:%s, var:%s\n",file,"ReadInFragParameter","totalFragCan"); }
    fclose(fpin);
    return fragformat;
}/*}}}*/

int Quick_Find_str_whole(char *(*database), int nsize, char *Cquery)/*{{{*/
{
    int bound_low, bound_high,bound_between, nindex,pre_bound_between;
    bool booldo;
    nindex = -1;

    pre_bound_between = 0;
    if (  nsize <= 0 )
    {
        return nindex;
    }
    bound_low = 0;
    bound_high = nsize-1;
    booldo = 1;
    if (   strcmp(database[bound_low], Cquery) == 0   )
    {
        nindex = bound_low;
        booldo = 0;
    }
    else if (   strcmp(database[bound_high], Cquery) == 0   )
    {
        nindex = bound_high;
        booldo = 0;
    }

    while (   (bound_high>bound_low) && (booldo)   )
    {
        bound_between = (bound_low+bound_high)/2;

        if (   strcmp(database[bound_between], Cquery) == 0   )
        {
            nindex = bound_between;
            booldo = 0;
            continue;
        }
        else if (   strcmp(database[bound_between], Cquery) > 0   )
        {
            bound_high = bound_between;
        }
        else
        {
            bound_low = bound_between;
        }
        if (  bound_between == pre_bound_between  )
        {
            booldo = 0;
        }
        pre_bound_between = bound_between;
    }
    return nindex;

}/*}}}*/
