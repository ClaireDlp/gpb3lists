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
}

/*Paramètres choisis grâce à la fonction FindlValue.
 *Pour n = 16, l = 7, N = 64.
 *Pour n = 32, l = 14, N = 16384.
 */


/*Alloue une liste de taille 2^l d'élément de F2^n et l'initialise avec du random*/
F2list* F2list_init(size_t n, size_t l){
  size_t N = (0x1 << l); 
  F2list *L = mzd_init(N,n); 

  mzd_randomize(L);

  return L;
}

/*initialize a po_F2list, where ordering is on the l LSB.*/
po_F2list* po_F2list_init(size_t nentry, size_t nbits, size_t l){
  po_F2list *T;
  T = malloc(sizeof(po_F2list));
  T->l = l;
  T->index = malloc((nentry+1) * sizeof(size_t));
  T->L = mzd_init(nentry, nbits);

  return T;
}

//Vector allocation.
F2vector* F2vector_init(size_t n){
  F2vector *x;
  x = malloc(sizeof(F2vector));
  x->n = n;
  x->n_word = (n >> 6) + 1; //ceil(n/64).
  x->v = malloc(x->n_word*sizeof(word));
  return x;
}

//Free vector.
void F2vector_free(F2vector *x){
  if(x == NULL){
    return;
  }
  if(x->v != NULL){
    free(x->v);
  }
  free(x);
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

  free(T);
}

/*Fonction qui récupère les l bits de poids faibles de la ligne i de la liste*/
uint32_t F2list_entry_low(F2list *L, size_t i, size_t l){
  uint32_t low = 0;
  size_t j;

  assert(l <= (size_t)L->ncols); //check inputs.

  /*On considère les elements de la liste comme des vecteurs sur F2
    Les bits de poids faibles sont à gauche. */

  for(j = 0; j < l; j++){
    low = BIT_AT_END_OF_STRING(low, mzd_read_bit(L, i, j));
  }

  return low;
}

/*Fonction qui parcourt L entre start et end et qui pour chacune des 2^l 
 *valeurs possibles de x , compte le nombre d'éléments de L dont les l bits 
 *de poids faibles valent x
 *Le resultat est stoké dans le tableau w*/
void F2list_count_low(size_t *w, F2list *L, size_t l, size_t start, size_t end){
  uint32_t low;
  size_t i;
  size_t N;

  assert(l <= (size_t)L->ncols);
  N = (0x1 << l);
  assert(N <= (size_t)L->nrows); //

  /*initialisation de w*/
  for(i = 0; i < N; i++){
    w[i] = 0;
  }

  /* on parcoure les entrées de L*/
  for(i = start; i < end; i++){
    low = F2list_entry_low(L, i, l);
    w[low]++;
  }
  return;
}

/*Fonction qui parcourt la liste L entre start et end, et permute les éléments de
 *sorte qu'ils soient ordonnés sur les l-bits de poids faibles.
 *le resultat est stocké dans la po_F2list T.*/
po_F2list* F2list_partial_ordering(F2list *L, size_t l, size_t start, size_t end){
  uint32_t low;
  size_t i, N, NL = end - start;
  size_t *ind, *w;
  F2list *TL;
  po_F2list *T = po_F2list_init(NL, L->ncols, l);

  ind = T->index;
  TL = T->L;
  N = (0x1 << l);

  //get workspace
  w = malloc(N * sizeof(size_t));

  //count low for each value
  F2list_count_low(w, L, l, start, end);
  

  //index table.
  ind[0] = 0;
  for(i = 0; i < N; i++){
    ind[i+1] = ind[i] + w[i];
    w[i] = 0; // reset w.
  }

  //Partialy sort list L in TL.
  /*we look at each entry in L between*/
  for(i = start; i < end; i++){
    low = F2list_entry_low(L, i, l);
    mzd_copy_row(TL, ind[low]+w[low], L, i); //copy entry i of L in good position in TL
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
  uint32_t low;

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
  uint32_t low;

  //initialize values.
  l = T->l;
  ind = T->index;
  TL = T->L;
  n = L->ncols;

  assert(n == (size_t)TL->ncols); // check inputs.

  //count the number of element in L12.
  count = F2listx2_partial_col_count(L, T, start, end);

  //allocate memory
  L12 = mzd_init(count, 2*n);
  count = 0; //reset variable count.

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
  uint32_t x;
  x = F2list_entry_low(L, 0, 4);

  printf("\n");
  po_F2list *T = F2list_partial_ordering(L, l16, 0, (0x1 << l16));
  
  //mzd_print(T->L);
  
  /* for(i = 0; i < (0x1 << (l16-1)); i++){ */
  /*   printf("%zu\n", w[i]); */
  /* } */

  F2list *L12 = F2list_partial_col_list(L, T, 0, 5);

  mzd_print(L12);

  //free memory
  mzd_free(L);
  po_F2list_free(T);
  
  return 0;
}
