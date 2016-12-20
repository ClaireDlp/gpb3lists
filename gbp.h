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
 * -> Joux Algorithm (Joux) (Parallel version ??)
 * -> Our new algorithms:
 * ---> 4-Collision search in L3 (2 methods : Naive one (col4) or Using Mitzenmacher results (Mitz))
 * ---> Using ISD
 *
 * Implementated using the m4ri library.
 */


//New type definition.
typedef mzd_t F2list; //<-- list of n-bit elements seen as vectors over F_2

typedef word F2vector; //<-- F2vector is uint64_t (assuming that n <=64).

typedef struct{ //<-- Partial ordering on a F2list.
  size_t l; // Partial ordering on l LSB
  size_t *index; //index[low] : first index such that low_l(L[index[low]]) = low.
  F2list *L; //partialy ordered list.
  size_t *pinv; //inverse of permutation (optional).
    
} po_F2list;



//MACROS
#define BIT_AT_END_OF_VECT(x, b) ((x << 1) ^ b) //<-- add bit b at the end of x.
#define GET_LAST_BIT(x) (x & 0x1)

//Functions

/* Initialize and free memory */
F2list* F2list_rand_init(size_t n, size_t l); //Get memory and initialize a list with random elements.
po_F2list* po_F2list_init(size_t nentry, size_t nbits, size_t l, int keep_pinv); //Initilise a partialy ordered list.
void po_F2list_free(po_F2list *T);
F2vector* F2vector_init(size_t n);
void F2vector_free(F2vector *x);
void F2list_copy_row(F2list *B, size_t ib, size_t jb, F2list *A, size_t ia, size_t j0, size_t j1);

/*Partial collision*/
F2vector F2list_entry_low(F2list *L, size_t i, size_t l); //get low_l(L[i]).
po_F2list* F2list_partial_ordering(F2list *L, size_t l, size_t start, size_t end, int keep_pinv); //partialy sort list L

/*Parameters estimation*/
size_t FindlValue(size_t n, int prec); //find value of l such that |L| = 2^l

