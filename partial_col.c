#include"gbp.h"

//TODO : Move functions into other files.

/*Etant donné n, donne une estimation de l ( en théorie l = 1/2 (n - log(4/3 (n-l))) ) */
size_t FindlValue(size_t n, int prec){
  size_t tmp, l, n_half;
  float p, logp;
  if(prec < 1){
    prec = 5; // choix arbitraire.
  }

  /*initialisation*/
  n_half = n/2; // Division exacte car n puissance de 2.
  l = n/2; 
  tmp = 2*l; // tmp  = 4 * (n-l) avec l = n/2
  p = tmp*1.00/3; 

  /*calcul de la valeur de l*/
  while(prec > 0){
    logp = log2f(p);
    tmp = (size_t)floor(logp/2); // <-- valeur arrondie de 1/2 log(p).
    l = n_half - tmp; // l <-- n/2 - round(1/2 log(p)).
    tmp = 4*(n-l);
    p = tmp*1.00/3;
    prec--;
  }

  return l;
} //TODO change this function to include other cases.

/*Paramètres choisis grâce à la fonction FindlValue.
 *Pour n = 16, l = 7, N = 64.
 *Pour n = 32, l = 14, N = 16384.
 */

/*
 * Function malloc with test.
 */
void * F2malloc(size_t size){
  void *x = malloc(size);
  if (x == NULL) {
    perror("malloc failed");
    exit(1);
  }
  return x;
}

/*Alloue une liste de taille 2^l d'élément de F2^n et l'initialise avec du random*/
F2list* F2list_init(size_t n, size_t l){
  size_t N = (0x1 << l); 
  F2list *L = mzd_init(N,n); 

  mzd_randomize(L);

  return L;
}

/*initialize a po_F2list, where ordering is on the l LSB.*/
po_F2list* po_F2list_init(size_t nentry, size_t nbits, size_t l, int keep_pinv){
  po_F2list *T;
  T = F2malloc(sizeof(po_F2list));
  T->l = l;
  T->index = F2malloc(((1 << l)+1) * sizeof(size_t));
  T->L = mzd_init(nentry, nbits);
  T->pinv = (keep_pinv) ? F2malloc(nentry * sizeof(size_t)) : NULL;

  return T;
}


//copy row ia from matrix A between cols j0 and j1 into row ib matrix B starting at jb.
void F2list_copy_row(F2list *B, size_t ib, size_t jb, F2list *A, size_t ia, size_t j0, size_t j1){
  size_t px;

  assert(j1 >= j0); //check inputs.
  for(px = j0; px < j1; px++){
    mzd_write_bit(B, ib, jb, mzd_read_bit(A, ia, px));
    jb++;
  }
  
}


void po_F2list_free(po_F2list *T){
  if(T== NULL){
    return;
  }
  if(T->index != NULL){
    free(T->index);
  }
  if(T->L != NULL){
    mzd_free(T->L);
  }
  if(T->pinv != NULL){
    free(T->pinv);
  }

  free(T);
}

/*Fonction qui récupère les l bits de poids faibles de la ligne i de la liste*/
F2vector F2list_entry_low(F2list *L, size_t i, size_t l){
  F2vector low = 0;
  size_t j;

  assert(l <= (size_t)L->ncols); //check inputs.

  /*On considère les elements de la liste comme des vecteurs sur F2
    Les bits de poids faibles sont à gauche. */

  for(j = 0; j < l; j++){
    low = BIT_AT_END_OF_VECT(low, mzd_read_bit(L, i, j));
  }

  return low;
}


/*
 * Function that get bits j0 to j1 of row i dans store them in a F2vector.
 */
F2vector F2list_get_entry(F2list *L, size_t i, size_t j0, size_t j1){
  F2vector x = 0;
  size_t j;

  assert(j1-j0 < 65); //check size of output.

  for(j = j0; j< j1; j++){
    x = BIT_AT_END_OF_VECT(x, mzd_read_bit(L, i, j));
  }

  return x;
}

/*Fonction qui parcourt L entre start et end et qui pour chacune des 2^l 
 *valeurs possibles de x , compte le nombre d'éléments de L dont les l bits 
 *de poids faibles valent x
 *Le resultat est stoké dans le tableau w*/
void F2list_count_low(size_t *w, F2list *L, size_t l, size_t start, size_t end){
  F2vector low;
  size_t i;
  size_t N;

  assert(l <= (size_t)L->ncols);
  N = (0x1 << l);
  //assert(N <= (size_t)L->nrows);

  /*initialisation de w*/
  for(i = 0; i < N; i++){
    w[i] = 0;
  }

  /* on parcoure les entrées de L*/
  for(i = start; i < end; i++){
    low = F2list_entry_low(L, i, l);
    w[low]++;
    fprintf(stderr,"\r[count low = %llu] : %zu ...", low, w[low]);
    fflush(stderr);
  }
  return;
}

/*Fonction qui parcourt la liste L entre start et end, et permute les éléments de
 *sorte qu'ils soient ordonnés sur les l-bits de poids faibles.
 *le resultat est stocké dans la po_F2list T.*/
po_F2list* F2list_partial_ordering(F2list *L, size_t l, size_t start, size_t end, int keep_pinv){
  F2vector low;
  size_t i, N, NL = end - start;
  size_t *ind, *w;
  F2list *TL;
  po_F2list *T = po_F2list_init(NL, L->ncols, l, keep_pinv);

  ind = T->index;
  TL = T->L;
  N = (0x1 << l);

  //get workspace
  w = F2malloc(N * sizeof(size_t));

  //count low for each value
  fprintf(stderr,"[Partial ordering]: count low start...\n");
  F2list_count_low(w, L, l, start, end);
  

  //index table.
  ind[0] = 0;
  for(i = 0; i < N; i++){
    ind[i+1] = ind[i] + w[i];
    w[i] = 0; // reset w.
  }
  fprintf(stderr,"\n");
  fprintf(stderr,"[Partial ordering]: re-order list...\n");

  //Partialy sort list L in TL.
  /*we look at each entry in L between*/
  for(i = start; i < end; i++){
    low = F2list_entry_low(L, i, l);
    mzd_copy_row(TL, ind[low]+w[low], L, i); //copy entry i of L in good position in TL
    if(keep_pinv){
      T->pinv[ind[low]+w[low]] = i; // row ind[low]+w[low] is row i of L.
    }
    w[low]++; // incrementation of the counter.
  }

  //free workspace.
  free(w);
  return T;
}

/*
 * Function that iterate over a given list L, and for each entry L[i] search in the 
 * Partialy sorted list T if there is an entry T->L[j] such that
 * low_l(T->L[j]) = low_l(L[i]). 
 * The return value is the number of partial collisions ths found
 */
size_t F2listx2_partial_col_count(F2list *L, po_F2list *T, size_t start, size_t end){
  size_t count = 0, i, l;
  size_t *ind;
  F2vector low;

  l = T->l;
  ind = T->index;


  // Look at each entry of L
  for(i = start; i < end; i++){
    low = F2list_entry_low(L, i, l); // get the l LSB.
    
    if(ind[low+1] == ind[low]){
      continue; // no partial collision for entry L[i].
    }
    
    assert(ind[low+1] > ind[low]); // at least one collision.
    count += ind[low+1] - ind[low]; // number of entries in T->L such that
    fprintf(stderr,"\r[Wagner step 1: Partial collisons] count = %zu...", count);
    fflush(stderr);
    //low_l(T->L[j]) = low_l(L[i]).    
  }

  return count;
}

/*
 * Iterate over the list L and search partial collisions in T. For each entry such that
 * low_l(L[i]) = low_l(TL[j]), store the vector (L[i], TL[j]) in the list L12.
 * return L12.
 */
F2list * F2list_partial_col_list(F2list *L, po_F2list *T, size_t start, size_t end){
  F2list *L12, *TL; // Wagner's Notation.
  size_t count, i, l, n, j;
  size_t *ind;
  F2vector low;

  //initialize values.
  l = T->l;
  ind = T->index;
  TL = T->L;
  n = L->ncols;

  assert(n == (size_t)TL->ncols); // check inputs.

  //count the number of element in L12.
  fprintf(stderr,"[Wagner step 1: Partial collisions]: start count...\n");
  count = F2listx2_partial_col_count(L, T, start, end);

  //allocate memory
  L12 = mzd_init(count, 2*n);
  count = 0; //reset variable count.

  fprintf(stderr,"\n");
  fprintf(stderr,"[Wagner step 1: Partial collisions]: fill list...\n");
  //For each entry in L, search partial collision in TL
  for(i = start; i < end; i++){
   low = F2list_entry_low(L, i, l); // get the l LSB.
    
    if(ind[low+1] == ind[low]){
      continue; // no partial collision for entry L[i].
    }
    
    assert(ind[low+1] > ind[low]); // at least one collision.
    for(j = ind[low]; j < ind[low+1]; j++){
      F2list_copy_row(L12, count, 0, L, i, 0, n); // copy row i of L in the first n cols of L12[count].
      F2list_copy_row(L12, count, n, TL, j, 0, n); //copy row j of TL in the last n cols of L12[count].
      count++;
    }
  }
  return L12; 
}

/*
 * given a list L12 s.t. entries are 2n-bit long, and a row i, return a vector x s.t.
 * x = low_n(L12[i]) ^ high_n(L12[i]).
 */
F2vector F2vector_xor_lr(F2list *L12, size_t i){
  F2vector x = 0;
  int n;

  n = L12->ncols;

  //check inputs.
  assert(((n >> 1) << 1) == n); // n is even.

  x = F2list_get_entry(L12, i, 0, (n>>1)); // store the first n/2 bits of L12[i] in x.
  x ^= F2list_get_entry(L12, i, (n>>1), n); //xore with the last n/2 bits of L12[i].

  return x;
}

/*
 * given a list L12 of 2n-bit long entries, create a F2list L1xor2, s.t.
 * L1xor2[i] = low_n(L12[i]) ^ high_n(L12[i]), for all start <= i < end.
 * As long as low_l(L1xor2[i]) = 0 in our case, we don't store the coefficients.
 */
F2list *F2list_xor_lr_entries(F2list *L12, size_t start, size_t end, int l){
  F2list *L1xor2;
  size_t i, N;
  int j, n;
  F2vector x;

  assert(end <= (size_t)L12->nrows); // check inputs.
  N = end - start;
  n = (L12->ncols)>>1;

  assert(n > l); //check inputs.
  
  L1xor2 = mzd_init(N, (size_t)(n-l)); // allocation of memory.

  fprintf(stderr,"[Wagner prepare step 2]: Xor left to right part of list...\n");

  for(i = start; i < end; i++){
    x  = F2vector_xor_lr(L12, i); //get the xor in x.
    for(j = n-l-1; j>=0; j--){
      mzd_write_bit(L1xor2, i-start, j, GET_LAST_BIT(x)); //write from right to left.
      x = (x >> 1); // update x so that new last bit would be bit j-1.
      }
  }
  return L1xor2;
}

/*
 * Given a vector x, check whether x is in the partialy ordered list TL.
 * xn is the size of x.
 * return 1 if x is in TL and 0 otherwise.
 */
int F2list_search_vector(po_F2list *T, F2vector x, size_t nx){
  F2vector low, *vec;
  size_t i, l;
  size_t *ind;
  F2list *TL;

  l = T->l;
  ind = T->index;
  TL = T->L;

  assert(nx >= l); //check inputs.

  //low <-- low_l(x).
  low = (x>>(nx-l));

  if(ind[low+1]-ind[low] == 0){ // no entry in TL that collide with x on the l LBS.
    fprintf(stderr,"\r[Check vector]: No collision found...");
    fflush(stderr);
    return 0;
  }

  for(i = ind[low]; i< ind[low+1]; i++){ // For all entry that collide with x on the l LSB, check if x = TL[i].
    vec = mzd_row(TL,i); //get pointer on the first row of TL (we assume TL has only one m4ry word).
    if(x == vec[1]){
      fprintf(stderr,"\r[Check vector]: Collision found. End of research.\n");
      fflush(stderr);
      return 1;
    }
  }

  fprintf(stderr,"\r[Check vector]: No collision found...\n");
  fflush(stderr);
  return 0;
  
}

int main(){
  /*Calcul de la valeur de l théorique pour n in {16, 32, 64} */
  size_t l16, l32, l64;
  //uint32_t i;
  l16 = FindlValue(16,5);
  l32 = FindlValue(32,5);
  l64 = FindlValue(64,5);

  //Test: printf("%zu : %zu : %zu\n", l16, l32, l64);

  /*initialisation d'une liste L*/
  F2list *L = F2list_init(16, l16);

  //mzd_print(L);
  //F2vector x;
  //x = F2list_entry_low(L, 0, 4);

  //printf("\n");
  po_F2list *T = F2list_partial_ordering(L, l16, 0, (0x1 << l16), 0);
  
  //mzd_print(T->L);
  
  /* for(i = 0; i < (0x1 << (l16-1)); i++){ */
  /*   printf("%zu\n", w[i]); */
  /* } */

  F2list *L12 = F2list_partial_col_list(L, T, 0, 5);
  po_F2list_free(T); //free memory.
  
  F2list *L1xor2 = F2list_xor_lr_entries(L12, 0, L12->nrows, l16);
  T = F2list_partial_ordering(L1xor2, l16, 0, L12->nrows, 1); //sort list L1xor2.
  F2list_search_vector(T, 0, L1xor2->ncols); 
  
  //mzd_print(L12);
  //printf("\n");
  //mzd_print(L1xor2);

  //free memory
  mzd_free(L12);
  mzd_free(L1xor2);
  mzd_free(L);
  po_F2list_free(T);
  
  return 0;
}
