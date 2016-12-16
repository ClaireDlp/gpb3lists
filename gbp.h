#include<math.h>
#include<stdio.h>
#include<inttypes.h>
#include<m4ri/m4ri.h>

/*
 * Implementation of algorithm for the Generalized Birthday
 * Problem with 3 lists, using Wagner algorithm and
 * various optimizations including ours.
 * Algorithms to be implemented here:
 *
 * -> Basic Wagner Algorithm (Wagner).
 * -> Nicolic & Sasaki (NS) (+ a Parallel version with more targets).
 * -> Joux Algorithm (Joux)
 * -> Our new algorithms:
 * ---> 4-Collision search in L3 (2 methods : Naive one (col4) or Using Mitzenmacher results (Mitz))
 * ---> Using ISD
 *
 * Implementated using the m4ri library.
 */


//New type definition.
typedef mzd_t F2list; //<-- list of n-bit elements seen as vectors over F_2

typedef struct{
  size_t n; //size of the vector.
  size_t n_word; //number of word (uint64_t) in the vector.
  word *v; //vector (v_0....v_n).
} F2vector;

typedef struct{
  size_t l; // Partial ordering on l LSB
  size_t *index; //index[low] : first index such that low_l(x) = low.
  F2list *L; //partialy ordered list.
} po_F2list;


//MACROS
#define BIT_AT_END_OF_STRING(x, b) ((x << 1) ^ b)

//Functions

/* Initialize and free memory */
F2list* F2list_rand_init(size_t n, size_t l); //Get memory and initialize a list with random elements.
po_F2list* po_F2list_init(size_t nentry, size_t nbits, size_t l); //Initilise a partialy ordered list.
void po_F2list_free(po_F2list *T);
F2vector* F2vector_init(size_t n);
void F2vector_free(F2vector *x);
void F2list_copy_row(F2list *B, size_t ib, size_t jb, F2list *A, size_t ia, size_t j0, size_t j1);

/*Partial collision*/
uint32_t F2list_entry_low(F2list *L, size_t i, size_t l); //get low_l(L[i]).
po_F2list* F2list_partial_ordering(F2list *L, size_t l, size_t start, size_t end); //partialy sort list L

/*Parameters estimation*/
size_t FindlValue(size_t n, int prec); //find value of l such that |L| = 2^l

