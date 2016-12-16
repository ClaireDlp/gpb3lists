#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <m4ri/m4ri.h>

static inline void swap(rci_t *A, size_t i, size_t j) {
	rci_t x = A[i];
	A[i] = A[j];
	A[j] = x;
}

/* randomizes the permutation P */
void mzp_randomize(mzp_t *P) {
	// mzp_set_ui(P, 1);
	size_t n = P->length;
	for (size_t i = 0; i < n - 1; i++) {
		size_t j = random() % (n - i);
		swap(P->values, i, i + j);
	}
} 

/* binary entropy function */
double H(double x) {
	if (x == 0)
		return 0;
	if (x == 1)
		return 0;
	return -(x * log(x) + (1 - x) * log(1 - x)) / M_LN2;
}


/* reciprocal binary entropy -- bissection method */
double H_inv(double y) {
	double a = 0;
	double b = 0.5;

	while (b - a > 1e-9) {
		const double x = (a + b) / 2;
		const double Hx = H(x);
		if (Hx < y)
			a = x;
		else
			b = x;
	}
	return (a + b) / 2;
}

/* Gilbert-Varshamov bound on the relative distance */
double GV_S(double R) {
	return H_inv(1 - R);
}

/* find a low-weight code word -- single-core, but easly parallelizable */
mzd_t * LeeBrickel_ISD(size_t n, size_t k, size_t target_d, const mzd_t *H) {
	/* setup */
	mzp_t *P = mzp_init(n);
	mzd_t *H_copy = mzd_init(n-k, n);
	mzd_t *A = mzd_init(n-k, n-k);
	mzd_t *I = mzd_init(n-k, n-k);
	mzd_t *Ht = mzd_init(n, n-k);
	mzd_t *Solution = mzd_init(n, 1);

	mzd_set_ui(I, 1);
	const size_t words_per_row = ceil(((double) (n-k)) / m4ri_radix);

	/* one iteration */
	int found = 0;
	size_t iterations = 0;
	size_t best_d = n;
	while (!found) {
		fprintf(stderr, "\rIteration %zu, best distance %zu", iterations++, best_d);
		/* permute phase */
		mzp_randomize(P);
		mzd_copy(H_copy, H);
		mzd_apply_p_right(H_copy, P);
		mzd_echelonize(H_copy, 1);
		
		/* check for leading I submatrix */
		mzd_submatrix(A, H_copy, 0, 0, n-k, n-k);
		if (!mzd_equal(A, I))
			continue;

		/* search phase */
		mzd_transpose(Ht, H_copy);
		for (size_t i = n-k; i < n; i++) {
			/* compute weight of column i */
			size_t w = 0;
			word *row = mzd_row(Ht, i);
			for (size_t j = 0; j < words_per_row; j++)
				w += __builtin_popcountll(row[j]);
			
			/* improvement ? */
			if (w + 1 < best_d) {
				best_d = w + 1;
				mzd_set_ui(Solution, 0);
				for (size_t j = 0; j < n-k; j++)
					mzd_write_bit(Solution, j, 0, mzd_read_bit(Ht, i, j));	
				mzd_write_bit(Solution, i, 0, 1);

 				if (best_d <= target_d)
					found = 1;
				break;
			}
		}
	}

	/* cleanup */
	mzd_apply_p_left(Solution, P);
	mzp_free(P);
	mzd_free(H_copy);
	mzd_free(Ht);
	mzd_free(A);
	mzd_free(I);
	return Solution;
}


int main() {
	srandom(42);

	size_t n = 512;
	size_t k = 32;

	double R = ((double) k) / ((double) n);
	double S = GV_S(R);
	int d = S*n + 8;

	printf("R = %.2f\n", R);
	printf("S = %.2f\n", S);
	printf("d = %d (Gilbert-Varshamov bound)\n", d);
 
        double e = H(S) - (1 - R) * H((S - 1.0/n) / (1-R)); //- log(k) / M_LN2;
        printf("decoding exponent = %.4f\n", e);
        printf("expected #iterations = 2^%.1f\n", n * e);

	/* random code setup */
	mzd_t *H = mzd_init(n-k, n);
	mzd_randomize(H);

	/* action ! */
	mzd_t *w = LeeBrickel_ISD(n, k, d, H);
	printf("\n");

	/* tests */
	mzd_t *check = mzd_mul(NULL, H, w, 0);
	assert (mzd_is_zero(check));
	mzd_free(check);

	int h = 0;
	for (size_t i = 0; i < n; i++)
		h += mzd_read_bit(w, i, 0);
 				
 	assert (h <= d);
 	printf("Correct solution found, weight %d : \n", h);
 	mzd_t *wT = mzd_transpose(NULL, w);
 	mzd_print(wT);
	
	return 0;
}


