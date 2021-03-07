int asm_enum ( int n );
void asm_triangle ( int n, int a[] );
void bell ( int n, int b[] );
void bell_values ( int *n_data, int *n, int *c );
bool bvec_add ( int n, int ivec[], int jvec[], int kvec[] );
bool bvec_check ( int n, int ivec[] );
void bvec_complement2 ( int n, int ivec[], int jvec[] );
void bvec_not ( int n, int ivec[], int jvec[] );
void bvec_print ( int n, int a[], char *title );
void bvec_print_unsigned ( int n, unsigned int a[], char *title );
void bvec_reverse ( int n, int ivec[], int jvec[] );
void bvec_sub ( int n, int ivec[], int jvec[], int kvec[] );
int bvec_to_i_signed ( int n, int ivec[] );
unsigned int bvec_to_i_unsigned ( int n, unsigned int ivec[] );
void bvec_xor ( int n, int ivec[], int jvec[], int kvec[] );
void bvec_xor_unsigned ( int n, unsigned int ivec[], unsigned int jvec[], 
  unsigned int kvec[] );
void catalan ( int n, int c[] );
void catalan_row_next ( bool next, int n, int irow[] );
void catalan_values ( int *n_data, int *n, int *c );
void cfrac_to_rat ( int n, int a[], int p[], int q[] );
void cfrac_to_rfrac ( int m, double g[], double h[], double p[], double q[] );
char ch_cap ( char c );
void change_greedy ( int total, int coin_num, int coin_value[], int *change_num, 
  int change[] );
void change_next ( int total, int coin_num, int coin_value[], int *change_num, 
  int change[], bool *done  );
bool chinese_check ( int n, int m[] );
int chinese_to_i ( int n, int m[], int r[] );
void comb_next ( int n, int k, int a[], bool *done );
void comb_row ( bool next, int n, int row[] );
void comb_unrank ( int m, int n, int rank, int a[] );
double combin ( int n, int k );
int combin2 ( int n, int k );
void comp_next ( int n, int k, int a[], bool *more );
void comp_random ( int n, int k, int *seed, int a[] );
void compnz_next ( int n, int k, int a[], bool *more );
void compnz_random ( int n, int k, int *seed, int a[] );
int congruence ( int a, int b, int c, bool *error );
void count_pose_random ( int *seed, int blocks[], int*goal );
double d_agm ( double a, double b );
double d_epsilon ( void );
double d_huge ( void );
int d_nint ( double x );
double d_pi ( void );
void d_swap ( double *x, double *y );
void d_to_cfrac ( double r, int n, int a[], int p[], int q[] );
void d_to_dec ( double rval, int dec_digit, int *mantissa, int *exponent );
void d_to_rat ( double a, int ndig, int *iatop, int *iabot );
double d_uniform ( double rlo, double rhi, int *seed );
double d_uniform_01 ( int *seed );
void debruijn ( int m, int n, int string[] );
void dec_add ( int mantissa1, int exponent1, int mantissa2, int exponent2, 
  int dec_digit, int *mantissa, int *exponent );
void dec_div ( int mantissa1, int exponent1, int mantissa2, int exponent2, 
  int dec_digit, int *mantissa, int *exponent, bool *error );
void dec_mul ( int mantissa1, int exponent1, int mantissa2, int exponent2, 
  int dec_digit, int *mantissa, int *exponent );
void dec_round ( int mantissa1, int exponent1, int dec_digit, 
  int *mantissa2, int *exponent2 );
double dec_to_d ( int mantissa, int exponent );
void dec_to_rat ( int mantissa, int exponent, int *rat_top, int *rat_bot );
char *dec_to_s ( int mantissa, int exponent );
int dec_width ( int mantissa, int exponent );
void decmat_det ( int n, int atop[], int abot[], int dec_digit, 
  int *dtop, int *dbot );
void decmat_print ( int m, int n, int a[], int b[], char *title );
void derange_back_candidate ( int n, int a[], int k, int *nstack, int stack[], 
  int ncan[] );
void derange_back_next ( int n, int a[], bool *more );
bool derange_check ( int n, int a[] );
int derange_enum ( int n );
void derange_enum2 ( int n, int d[] );
int derange_enum3 ( int n );
void derange_weed_next ( int n, int a[], bool *more );
char digit_to_ch ( int digit );
void digraph_arc_euler ( int nnode, int nedge, int inode[], int jnode[], 
  bool *success, int trail[] );
void digraph_arc_print ( int nedge, int inode[], int jnode[], char *title );
void diophantine ( int a, int b, int c, bool *error, int *x, int *y );
void diophantine_solution_minimize ( int a, int b, int *x, int *y );
void equiv_next ( int n, int *npart, int jarray[], int iarray[], bool *more );
double dmat_det ( int n, double a[] );
void dmat_perm ( int n, double a[], int p[] );
void dmat_perm2 ( int m, int n, double a[], int p[], int q[] );
double dmat_permanent ( int n, double a[] );
void dmat_print ( int m, int n, double a[], char *title );
void dmat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi, 
  int jhi, char *title );
void dmat_transpose_print ( int m, int n, double a[], char *title );
void dpoly ( int n, double a[], double x0, int iopt, double *val );
int dpoly_degree ( int na, double a[] );
void dpoly_div ( int na, double a[], int nb, double b[], int *nq, double q[], 
  int *nr, double r[] );
void dpoly_f2p ( int n, double a[] );
double dpoly_fval ( int n, double a[], double x );
void dpoly_mul ( int na, double a[], int nb, double b[], double c[] );
void dpoly_n2p ( int n, double a[], double xarray[] );
double dpoly_nval ( int n, double a[], double xarray[], double x );
void dpoly_nx ( int n, double a[], double xarray[], double x );
void dpoly_p2f ( int n, double a[] );
void dpoly_p2n ( int n, double a[], double xarray[] );
void dpoly_p2t ( int n, double a[], double x );
void dpoly_power ( int na, double a[], int p, double b[] );
void dpoly_print ( int n, double a[], char *title );
double dpoly_pval ( int n, double a[], double x );
void dpoly_t2p ( int n, double a[], double x );
void dvec_backtrack ( int n, int maxstack, int stack[], double x[], int *indx, 
  int *k, int *nstack, int ncan[] );
double dvec_frac ( int n, double a[], int k );
void dvec_indicator ( int n, double a[] );
bool dvec_mirror_next ( int n, double a[] );
void dvec_print ( int n, double a[], char *title );
void dvec_uniform ( int n, double alo, double ahi, int *seed, double a[] );
void equiv_next2 ( bool *done, int iarray[], int n );
void equiv_print ( int n, int iarray[], char *title );
void equiv_random ( int n, int *seed, int *npart, int a[], double b[] );
void euler ( int n, int ieuler[] );
double fall ( double x, int n );
double gamma_log ( double x );
void gamma_log_values ( int *n_data, double *x, double *fx );
unsigned long get_seed ( void );
void gray_next ( int n, int *change );
int gray_rank ( int gray );
int gray_rank2 ( int gray );
int gray_unrank ( int rank );
int gray_unrank2 ( int rank );
int i_bset ( int i, int bit );
bool i_btest ( int i, int bit );
void i_factor ( int n, int maxfactor, int *nfactor, int factor[], 
  int power[], int *nleft );
int i_factorial ( int n );
int i_gcd ( int i, int j );
int i_huge ( void );
int i_log_10 ( int i );
int i_max ( int i1, int i2 );
int i_min ( int i1, int i2 );
int i_modp ( int i, int j );
int i_moebius ( int n );
void i_partition_conj ( int n, int iarray1[], int mult1[], int npart1, 
  int iarray2[], int mult2[], int *npart2 );
void i_partition_count ( int n, int p[] );
int *i_partition_count2 ( int n );
void i_partition_count_values ( int *n_data, int *n, int *c );
void i_partition_next ( bool *done, int a[], int mult[], int n, int *npart );
void i_partition_next2 ( int n, int a[], int mult[], int *npart, bool *more );
void i_partition_print ( int n, int npart, int a[], int mult[] );
void i_partition_random ( int n, int table[], int *seed, int a[], int mult[], 
  int *npart );
int i_uniform ( int ilo, int ihi, int *seed );
int i_sign ( int i );
void i_sqrt ( int n, int *q, int *r );
void i_sqrt_cf ( int n, int max_term, int *n_term, int b[] );
void i_swap ( int *i, int *j );
void i_to_bvec_signed ( int i, int n, int bvec[] );
void i_to_bvec_unsigned ( unsigned int i, int n, unsigned int bvec[] );
void i_to_chinese ( int j, int n, int m[], int r[] );
void i_to_ipoly ( int intval, int base, int degree_max, int *degree, int a[] );
double i_to_van_der_corput ( int seed, int base );
void imat_01_rowcolsum ( int m, int n, int r[], int c[], int a[], bool *error );
void imat_01_rowcolsum2 ( int m, int n, int r[], int c[], int a[], bool *error );
void imat_u1_inverse ( int n, int a[], int b[] );
void imat_perm ( int n, int a[], int p[] );
void imat_perm2 ( int m, int n, int a[], int p[], int q[] );
void imat_print ( int m, int n, int a[], char *title );
void imat_print_some ( int m, int n, int a[], int ilo, int jlo, int ihi, 
  int jhi, char *title );
void index_box2_next_2d ( int n1, int n2, int ic, int jc, int *i, int *j, 
  bool *more );
void index_box2_next_3d ( int n1, int n2, int n3, int ic, int jc, int kc, 
  int *i, int *j, int *k, bool *more );
void index_box_next_2d ( int n1, int n2, int *i, int *j, bool *more );
void index_box_next_3d ( int n1, int n2, int n3, int *i, int *j, int *k, 
  bool *more );
void index_next0 ( int n, int hi, int a[], bool *more );
void index_next1 ( int n, int hi[], int a[], bool *more );
void index_next2 ( int n, int lo[], int hi[], int a[], bool *more );
int index_rank0 ( int n, int hi, int a[] );
int index_rank1 ( int n, int hi[], int a[] );
int index_rank2 ( int n, int lo[], int hi[], int a[] );
void index_unrank0 ( int n, int hi, int rank, int a[] );
void index_unrank1 ( int n, int hi[], int rank, int a[] );
void index_unrank2 ( int n, int lo[], int hi[], int rank, int a[] );
void ins_perm ( int n, int ins[], int p[] );
void involute_enum ( int n, int s[] );
void ipoly ( int n, int a[], int x0, int iopt, int *val );
void ipoly_cyclo ( int n, int phi[] );
int ipoly_degree ( int na, int a[] );
void ipoly_div ( int na, int a[], int nb, int b[], int *nq, int q[], 
  int *nr, int r[] );
void ipoly_mul ( int na, int a[], int nb, int b[], int c[] );
void ipoly_print ( int n, int a[], char *title );
int ipoly_to_i ( int n, int a[], int x );
bool ivec_ascends ( int n, int x[] );
void ivec_backtrack ( int n, int maxstack, int stack[], int x[], int *indx, 
  int *k, int *nstack, int ncan[] );
bool ivec_descends ( int n, int x[] );
int ivec_frac ( int n, int a[], int k );
void ivec_indicator ( int n, int a[] );
int ivec_index ( int n, int a[], int aval );
int ivec_max ( int n, int a[] );
int ivec_maxloc_last ( int n, int x[] );
int ivec_min ( int n, int a[] );
bool ivec_pairwise_prime ( int n, int a[] );
int ivec_product ( int n, int a[] );
void ivec_reverse ( int n, int a[] );
void ivec_sort_bubble_a ( int n, int a[] );
int *ivec_sort_heap_index_a ( int n, int a[] );
int *ivec_sort_heap_index_d ( int n, int a[] );
int ivec_sum ( int n, int a[] );
void ivec_uniform ( int n, int alo, int ahi, int *seed, int a[] );
void ivec0_print ( int n, int a[], char *title );
void ivec1_print ( int n, int a[], char *title );
void jfrac_to_rfrac ( int m, double r[], double s[], double p[], double q[] );
int josephus ( int n, int m, int k );
void ksub_next ( int n, int k, int a[], bool *more );
void ksub_next2 ( int n, int k, int a[], int *in, int *iout );
void ksub_next3 ( int n, int k, int a[], bool *more, int *in, int *iout );
void ksub_next4 ( int n, int k, int a[], bool *done );
void ksub_random ( int n, int k, int *seed, int a[] );
void ksub_random2 ( int n, int k, int *seed, int a[] );
void ksub_random3 ( int n, int k, int *seed, int a[] );
void ksub_random4 ( int n, int k, int *seed, int a[] );
void ksub_rank ( int k, int a[], int *rank );
void ksub_unrank ( int k, int rank, int a[] );
void matrix_product_opt ( int n, int rank[], int *cost, int order[] );
void moebius_matrix ( int n, int a[], int mu[] );
int morse_thue ( int i );
int multinomial_coef1 ( int nfactor, int factor[] );
int multinomial_coef2 ( int nfactor, int factor[] );
void network_flow_max ( int nnode, int nedge, int iendpt[], int icpflo[], 
  int source, int sink, int cut[], int node_flow[] );
unsigned int nim_sum ( int i, int j );
void padovan ( int n, int p[] );
void pell_basic ( int d, int *x0, int *y0 );
void pell_next ( int d, int x0, int y0, int xn, int yn, int *xnp1, int *ynp1 );
int pent_enum ( int n );
void perm_ascend ( int n, int a[], int *length, int sub[] );
int perm_break_count ( int n, int p[] );
void perm_canon_to_cycle ( int n, int p1[], int p2[] );
bool perm_check ( int n, int p[] );
void perm_cycle ( int n, int p[], int *isgn, int *ncycle, int iopt );
void perm_cycle_to_canon ( int n, int p1[], int p2[] );
void perm_cycle_to_index ( int n, int p1[], int p2[] );
int perm_distance ( int n, int a[], int b[] );
int perm_fixed_enum ( int n, int m );
void perm_free ( int ipart[], int npart, int nfree, int ifree[] );
void perm_index_to_cycle ( int n, int p1[], int p2[] );
void perm_ins ( int n, int p[], int ins[] );
void perm_inv ( int n, int p[] );
void perm_inv2 ( int n, int p[] );
void perm_lex_next ( int n, int p[], bool *more );
void perm_mul ( int n, int p1[], int p2[], int p3[] );
void perm_next ( int n, int p[], bool *more, bool *even );
void perm_next2 ( int n, int p[], bool *done );
void perm_next3 ( int n, int p[], bool *more );
void perm_print ( int n, int p[], char *title );
void perm_random ( int n, int *seed, int p[] );
void perm_random2 ( int n, int *seed, int p[] );
void perm_random3 ( int n, int *seed, int p[] );
int perm_rank ( int n, int p[], int invers[] );
int perm_sign ( int n, int p[] );
void perm_to_equiv ( int n, int p[], int *npart, int jarray[], int iarray[] );
void perm_to_ytb ( int n, int p[], int lambda[], int a[] );
void perm_unrank ( int n, int rank, int p[] );
void perrin ( int n, int p[] );
bool pord_check ( int n, int a[] );
int power_mod ( int a, int n, int m );
void power_series1 ( int n, double alpha, double a[], double b[] );
void power_series2 ( int n, double a[], double b[] );
void power_series3 ( int n, double a[], double b[], double c[] );
void power_series4 ( int n, double a[], double b[], double c[] );
int prime ( int n );
void pythag_triple_next ( int *i, int *j, int *a, int *b, int *c );
unsigned long rand_initialize ( unsigned long seed );
unsigned long random_initialize ( unsigned long seed );
void rat_add ( int itop1, int ibot1, int itop2, int ibot2, int *itop, int *ibot, 
  bool *error );
void rat_div ( int itop1, int ibot1, int itop2, int ibot2, int *itop, 
  int *ibot, bool *error );
void rat_farey ( int n, int max_frac, int *num_frac, int a[], int b[] );
void rat_farey2 ( int n, int a[], int b[] );
void rat_mul ( int itop1, int ibot1, int itop2, int ibot2, int *itop, 
  int *ibot, bool *error );
void rat_normalize ( int *a, int *b );
void rat_sum_formula ( int n, int a[], int b[] );
void rat_to_cfrac ( int ip, int iq, int m, int *n, int a[], bool *error );
void rat_to_dec ( int rat_top, int rat_bot, int *mantissa, int *exponent );
double rat_to_d ( int top, int bot );
int rat_width ( int a, int b );
void ratmat_det ( int n, int iatop[], int iabot[], int *idtop, int *idbot, 
  bool *error );
void ratmat_print ( int m, int n, int a[], int b[], char *title );
void regro_next ( bool *done, int n, int v[], int vmax[] );
void rfrac_to_cfrac ( int m, double p[], double q[], double t[], bool *error );
void rfrac_to_jfrac ( int m, double p[], double q[], double r[], double s[], 
  bool *error );
double rise ( double x, int n );
void s_blank_delete ( char *s );
void s_blanks_delete ( char *s );
bool s_eqi ( char *s1, char *s2 );
int s_len_trim ( char* s );
void schroeder ( int n, int s[] );
void sort_heap_external ( int n, int *indx, int *i, int *j, int isgn );
void sub_by_size_next ( int n, int a[], int *size, bool *more );
void sub_gray_next ( int n, int a[], bool *more, int *ncard, int *iadd );
int sub_gray_rank ( int n, unsigned int a[] );
void sub_gray_unrank ( int rank, int n, unsigned int a[] );
void sub_lex_next ( int n, bool jmp, int ndim, int *k, int a[] );
void sub_random ( int n, int *seed, int a[] );
void subcomp_next ( int n, int k, int a[], bool *more );
void subcompnz_next ( int n, int k, int a[], bool *more );
void subcompnz2_next ( int n_lo, int n_hi, int k, int a[], bool *more );
void thue_binary_next ( int *n, int thue[] );
void thue_ternary_next ( int *n, int thue[] );
void timestamp ( void );
void triang ( int n, int zeta[], int p[] );
void tuple_next ( int m1, int m2, int n, int *rank, int x[] );
void tuple_next_fast ( int m, int n, int rank, int x[] );
void tuple_next_ge ( int m, int n, int *k, int x[] );
void tuple_next2 ( int n, int xmin[], int xmax[], int x[], int *rank );
void vec_next ( int n, int base, int a[], bool *more );
void vec_next2 ( int n, int base[], int a[], bool *done );
void vec_random ( int n, int base, int *seed, int a[] );
int vec_rank ( int n, int base[], int a[] );
void vec_unrank ( int n, int base[], int a[], int rank );
int ytb_enum ( int n );
void ytb_next ( int n, int lambda[], int a[], bool *more );
void ytb_print ( int n, int a[], char *title );
void ytb_random ( int n, int lambda[], int *seed, int a[] );