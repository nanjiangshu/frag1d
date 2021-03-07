#ifndef  _HAS_BASE_H
#define _HAS_BASE_H

#include <map>
#include <vector>
#include <string>
using namespace std;
const int Confid_Cand[] =/*{{{*/
{
    67,
    68,
    71,
    71,
    71,
    73,
    74,
    77,
    78,
    74,
    78,
    81,
    84,
    85,
    85,
    87,
    86,
    88,
    88,
    88,
    88,
    89,
    90,
    95 
};/*}}}*/
const int Shape_Code[] =/*{{{*/
{
		5,//A--A
		8,//B-- -
		8,//C-- -
		8,//D-- -
		8,//E-- -
		8,//F-- -
		7,//G-- G
		8,//H-- -
		8,//I-- -
		8,//J-- -
		4,//K-- K
		8,//L-- -
		8,//M-- -
		8,//N-- -
		8,//O-- -
		8,//P-- -
		8,//Q-- -
		1,//R-- R
		0,//S-- S
		6,//T-- T
		2,//U-- U
		3,//V-- V
		8,//W-- -
		8,//X-- -
		8,//Y-- -
		8 //Z-- -
};/*}}}*/
typedef int                 BOOL;
typedef unsigned char       BYTE;

typedef signed char         int8;
typedef unsigned char       unit8;
typedef signed short        int16;
typedef unsigned short      unit16;
typedef signed int          int32;
typedef unsigned int        unit32;
typedef signed long long    int64;
typedef unsigned long long  unit64;

typedef int8       __int8;
typedef int8        _int8;
typedef unit8      __unit8;
typedef unit8       _unit8;
typedef int16      __int16;  
typedef int16       _int16;
typedef unit16     __unit16;
typedef unit16      _unit16;
typedef int32      __int32;
typedef int32       _int32;
typedef unit32     __unit32;
typedef unit32      _unit32;
typedef int64      __int64;
typedef int64       _int64;
typedef unit64     __unit64;
typedef unit64      _unit64;


#ifndef NULL
#define NULL 0
#endif


#ifndef MAX_PATH
#define MAX_PATH 513
#endif

#ifndef MAX_COMMAND_LINE
#define MAX_COMMAND_LINE 500
#endif

#ifndef QIJ_NAME_FORMAT
#define QIJ_NAME_FORMAT
#define QIJ_FORMAT_TUPING 0
#define QIJ_FORMAT_NANJIANG 1
#endif

#ifndef MODM_NAME_FORMAT
#define MODM_NAME_FORMAT
#define MODM_FORMAT_TUPING 0
#define MODM_FORMAT_NANJIANG 1
#endif

#ifndef FRAGACC_NAME_FORMAT
#define FRAGACC_NAME_FORMAT
#define FRAGACC_FORMAT_TUPING 0
#define FRAGACC_FORMAT_NANJIANG 1
#endif

#ifndef FRAGFORMAT
#define FRAGFORMAT
#define FRAGFORMAT_TUPING 0
#define FRAGFORMAT_NANJIANG 1
#endif

#ifndef PROFILESCORETYPE
#define PROFILESCORETYPE
#define PROFILESCORETYPE_TUPING 0
#define PROFILESCORETYPE_NANJIANG 1
#endif

#define TRAINING_SET 0
#define TEST_SET 1

#undef SIZE_ID
#define SIZE_ID 100


#ifndef INLINE
#  if __GNUC__
#    define INLINE extern inline
#  else
#    define INLINE inline
#  endif
#endif

int GetHSRState(int probH, int probS, int probR, int method_HSR = 1);

int Read_databse_SHU(int dbtype, char *database_list, int NPer_Frag_Database, int ratioScheme, char *database_qijmatrix,char *database_modmatrix, char *database_facc,  char *CResult_location,  char *(*Namelistunit), char *(*CAASSequence),char *(*CShapesequence),char *(*CHSRsequence),int *(*blosum_chain),int *(*Chain_code_AAS), int *LengthList , int dsspMapMethod = 0, bool isReadBinaryFile = true, int8 typeProfile=1);
int Read_databse(char *database_list, int NPer_Frag_Database, char *database_qijmatrix,char *database_modmatrix, char *database_facc,  char *CResult_location,  char *(*Namelistunit), char *(*CAASSequence),char *(*CShapesequence),char *(*CHSRsequence),int *(*blosum_chain),int *(*Chain_code_AAS), int *LengthList );
int Treat_pair_Percent_homology(float *Datapair, int Lenlong, int Lenshort, int *Max_sum, int *(*Candidate_Maxline), int Num_Candhomo);
int Treat_pair_to_optimisized_Line(float *Datapair, int Lenlong, int Lenshort, int *Max_Line_X, int *Max_Line_Y, int *datareturn, int *(*Candidate_Maxline), int Num_Candhomo );
/*int GetBinaryFragShort5( char *file, int &fragFileType, char **idList, int &numID, int &maxSizeID, int &length, short *posTar, short *numCan, int &totalFragCan, FragCanShort5 *fragCan);*/
int ReadInPairBinary(const char *binfragPairFile,  float *(*fragScore) , char *(*proIDCan), int datasize, int *(*subikk), int *(*sublen));
/*int GetBinaryFragShort6( char *file, int &fragFileType, char **idList, int &numID, int &maxSizeID, int &length, short *posTar, short *numCan, int &totalFragCan, FragCanShort6 *fragCan);*/
int ReadInFragParameter(const char *file, int &fragformat, int &numID, int &maxSizeID, int &length, int &totalFragCan);

int Quick_Find_str_whole(char *(*database), int nsize, char *Cquery);



#endif
