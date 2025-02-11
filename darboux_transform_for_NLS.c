#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include "gmp.h"
#include "mpfr.h"

#define PRECISION 500

typedef struct complex_nbr
{
	mpfr_t re, im;
} complex_nbr;

void init_complex(complex_nbr *z)
{
	mpfr_init_set_ui(z->re, 0, MPFR_RNDN);
	mpfr_init_set_ui(z->im, 0, MPFR_RNDN);
}

void clear_complex(complex_nbr *z)
{
	mpfr_clear(z->re);
	mpfr_clear(z->im);
}

void complex_mult(complex_nbr *out, complex_nbr *z1, complex_nbr *z2, mpfr_t d1)
{

	//(a+ib)*(c+id)=ac-bd+i(ad+bc)
	// store ac and ad, then replace out.re	 with bc

	// ac ad bc bd
	mpfr_fmms(d1, z1->re, z2->re, z1->im, z2->im, MPFR_RNDN);
	mpfr_fmma(out->im, z1->re, z2->im, z1->im, z2->re, MPFR_RNDN);

	mpfr_set(out->re, d1, MPFR_RNDN);
}

void complex_mult_d(complex_nbr *out, complex_nbr *z1, complex double z2, complex_nbr *dum)
{
	//(a+ib)*(c+id)=ac-bd+i(ad+bc)
	// store ac and ad, then replace out.re	 with bc

	// ac ad bc bd
	mpfr_mul_d(dum[2].re, z1->re, creal(z2), MPFR_RNDN);
	mpfr_mul_d(dum[2].im, z1->re, cimag(z2), MPFR_RNDN);
	mpfr_mul_d(out->re, z1->im, creal(z2), MPFR_RNDN);
	mpfr_mul_d(dum[1].re, z1->im, cimag(z2), MPFR_RNDN);

	mpfr_add(out->im, dum[2].im, out->re, MPFR_RNDN);

	mpfr_sub(out->re, dum[2].re, dum[1].re, MPFR_RNDN);
}

void complex_add(complex_nbr *out, complex_nbr *z1, complex_nbr *z2)
{
	mpfr_add(out->re, z1->re, z2->re, MPFR_RNDN);
	mpfr_add(out->im, z1->im, z2->im, MPFR_RNDN);
}

void complex_conj(complex_nbr *out, complex_nbr *z)
{
	if (!mpfr_zero_p(z->im))
	{
		mpfr_neg(out->im, z->im, MPFR_RNDN);
	}
	else
	{
		mpfr_set_zero(out->im, 1);
	}
	mpfr_set(out->re, z->re, MPFR_RNDN);
}

void complex_mtrx_add(complex_nbr *out, complex_nbr *M1, complex_nbr *M2)
{
	for (int i = 0; i < 4; i++)
	{
		complex_add(out + i, M1 + i, M2 + i);
	}
}

void complex_mtrx_mul_2x2(complex_nbr *out, complex_nbr *M1, complex_nbr *M2, complex_nbr *dum)
{
	for (int i = 0; i < 2; i++)
	{
		// dum[0]=m11a11 dum[1]=m12a21
		complex_mult(dum + 0, &(M1[2 * i]), &(M2[0]), dum[1].re);
		complex_mult(&(out[2 * i]), &(M1[2 * i + 1]), &(M2[2]), dum[1].re);
		complex_add(&(out[2 * i]), dum + 0, &(out[2 * i]));

		// dum[0]=m11a12 dum[1]=m12a22
		complex_mult(dum + 0, &(M1[2 * i]), &(M2[1]), dum[1].re);
		complex_mult(&(out[2 * i + 1]), &(M1[2 * i + 1]), &(M2[3]), dum[1].re);
		complex_add(&(out[2 * i + 1]), dum + 0, &(out[2 * i + 1]));
	}
}

void complex_d_to_mpfr(complex_nbr *out, complex double z)
{
	mpfr_set_d(out->re, creal(z), MPFR_RNDN);
	mpfr_set_d(out->im, cimag(z), MPFR_RNDN);
}

void complex_exp(complex_nbr *out, complex_nbr *z)
{
	mpfr_t exp_re;
	mpfr_init2(exp_re, PRECISION);

	mpfr_exp(exp_re, z->re, MPFR_RNDN);

	mpfr_cos(out->re, z->im, MPFR_RNDN);
	mpfr_sin(out->im, z->im, MPFR_RNDN);
	mpfr_mul(out->re, out->re, exp_re, MPFR_RNDN);
	mpfr_mul(out->im, out->im, exp_re, MPFR_RNDN);

	mpfr_clear(exp_re);
}

void complex_abs_square(mpfr_t out, complex_nbr *z1, mpfr_t dum)
{
	mpfr_sqr(dum, z1->re, MPFR_RNDN);
	mpfr_sqr(out, z1->im, MPFR_RNDN);
	mpfr_add(out, out, dum, MPFR_RNDN);
}

void complex_cpy(complex_nbr *out, complex_nbr *in)
{
	mpfr_set(out->re, in->re, MPFR_RNDN);
	mpfr_set(out->im, in->im, MPFR_RNDN);
}

void complex_mult_by_i(complex_nbr *out, complex_nbr *z, mpfr_t d)
{
	mpfr_set(d, z->re, MPFR_RNDN);
	if (!mpfr_zero_p(z->im))
	{
		mpfr_neg(out->re, z->im, MPFR_RNDN);
	}
	else
	{
		mpfr_set_zero(out->re, 1);
	}
	mpfr_set(out->im, d, MPFR_RNDN);
}

void complex_mult_by_re(complex_nbr *out, complex_nbr *z, mpfr_t r)
{
	mpfr_mul(out->re, z->re, r, MPFR_RNDN);
	mpfr_mul(out->im, z->im, r, MPFR_RNDN);
}

void complex_inverse(complex_nbr *out, complex_nbr *z, complex_nbr *dum)
{
	if (!mpfr_zero_p(z->im) || !mpfr_zero_p(z->re))
	{
		complex_abs_square(dum->re, z, dum->im);
		mpfr_div(out->im, z->im, dum->re, MPFR_RNDN);
		mpfr_div(out->re, z->re, dum->re, MPFR_RNDN);
		complex_conj(out, out);
	}
}

void complex_printf(complex_nbr *z)
{
	mpfr_printf("%.50Rf%+.50Rf", z->re, z->im);
	printf("i");
}

void complex_set_zero(complex_nbr *z)
{
	mpfr_set_zero(z->re, 1);
	mpfr_set_zero(z->im, 1);
}

void get_init_matrix(complex_nbr **out, double x, complex double *eigVs, int n, complex_nbr *dum)
{
	for (int i = 0; i < n; i++)
	{
		complex double zi = eigVs[i];

		// Get M^0(z*)
		// M[0]=exp(i*zi*x), M[3]=e^(-i*zi*x)
		complex_d_to_mpfr(dum, -I * zi * x);
		complex_exp(&(out[i][0]), dum);

		if (!mpfr_zero_p(dum->re))
		{
			mpfr_neg(dum->re, dum->re, MPFR_RNDN);
		}
		if (!mpfr_zero_p(dum->re))
		{
			mpfr_neg(dum->im, dum->im, MPFR_RNDN);
		}

		complex_exp(&(out[i][3]), dum);

		complex_set_zero(&(out[i][1]));
		complex_set_zero(&(out[i][2]));
	}
}

void darboux_transform(complex_nbr *sol, complex double *eigVs, double *coeffs, complex_nbr **init_matrix, double t, int n, complex_nbr *dum, FILE *log_ptr)
{
	// complex_nbr **M;
	complex_nbr q[2];
	complex_nbr ci;
	complex_nbr Q[4];
	complex_nbr **M;

	init_complex(&Q[0]);
	init_complex(&Q[1]);
	init_complex(&Q[2]);
	init_complex(&Q[3]);
	init_complex(&q[0]);
	init_complex(&q[1]);
	init_complex(&ci);

	complex_set_zero(sol);

	// // Store entries M11 M12 M21 M22
	M = calloc(n, sizeof(complex_nbr *));

	for (int i = 0; i < n; i++)
	{
		M[i] = calloc(4, sizeof(complex_nbr));

		init_complex(&(M[i][0]));
		init_complex(&(M[i][1]));
		init_complex(&(M[i][2]));
		init_complex(&(M[i][3]));

		complex_cpy(&(M[i][0]), init_matrix[i] + 0);
		complex_cpy(&(M[i][1]), init_matrix[i] + 1);
		complex_cpy(&(M[i][2]), init_matrix[i] + 2);
		complex_cpy(&(M[i][3]), init_matrix[i] + 3);
	}

	for (int i = 0; i < n; i++)
	{
		// get c_i(t) = c_i e^-2itz^2
		complex double zi = eigVs[i];
		complex_d_to_mpfr(&ci, -2 * I * zi * zi * t);

		complex_exp(&ci, &ci);
		complex_mult_d(&ci, &ci, coeffs[i], dum);

		// q[0]=M*(zi*)[0]+ci*M*(zi*)[1]
		complex_conj(&(q[0]), &(M[i][1]));
		complex_mult(&(q[0]), &(q[0]), &ci, dum[2].re);
		complex_conj(dum, &(M[i][0]));
		complex_add(&(q[0]), dum, &(q[0]));

		// q[1]=M(zi*)[2]+ci*M(zi*)[3]
		complex_conj(&(q[1]), &(M[i][3]));
		complex_mult(&(q[1]), &(q[1]), &ci, dum[2].re);
		complex_conj(dum, &(M[i][2]));
		complex_add(&(q[1]), &(q[1]), dum);

		// set dum.re=q_1^2, dum.im=q_2^2
		complex_abs_square(dum->re, &(q[0]), dum[2].re);
		complex_abs_square(dum->im, &(q[1]), dum[2].re);

		// set c1.re=2Im(z)/(q1^2+q2^2)
		mpfr_add(ci.re, dum->re, dum->im, MPFR_RNDN);
		mpfr_d_div(ci.re, 2 * cimag(zi), ci.re, MPFR_RNDN);

		// Q[0]=i*ci.re*|q_1|^2, Q[3]=|q_2|^2
		mpfr_set_zero(Q[0].re, 1);
		mpfr_mul(Q[0].im, dum->re, ci.re, MPFR_RNDN);
		mpfr_set_zero(Q[3].re, 1);
		mpfr_mul(Q[3].im, dum->im, ci.re, MPFR_RNDN);

		// set dum=c1.re* q_1^* q_2. Get Q_{12}
		complex_conj(&(q[0]), &(q[0]));
		complex_mult(dum, &(q[0]), &(q[1]), dum[2].re);

		complex_mult_by_re(dum, dum, ci.re);
		complex_mult_by_i(&(Q[1]), dum, dum[2].re);

		// Get Q_{21}
		complex_conj(dum, dum);
		complex_mult_by_i(&(Q[2]), dum, dum[2].re);

		complex_add(sol, sol, &(Q[1]));

		for (int k = i + 1; k < n; k++)
		{
			complex double zk_zi = conj(eigVs[k]) - zi;
			// get 1/(zk^*-zi)
			complex_d_to_mpfr(dum, zk_zi);

			complex_inverse(dum, dum, dum + 1);

			// M[0]=Qn/(zk-zn)

			complex_mtrx_mul_2x2(M[0], Q, M[k], dum + 1);

			complex_mult(&(M[0][0]), &(M[0][0]), dum, dum[2].re);
			complex_mult(&(M[0][1]), &(M[0][1]), dum, dum[2].re);
			complex_mult(&(M[0][2]), &(M[0][2]), dum, dum[2].re);
			complex_mult(&(M[0][3]), &(M[0][3]), dum, dum[2].re);

			complex_mtrx_add(M[k], M[k], M[0]);
		}
	}

	// sol*2i
	mpfr_set(dum->re, sol->re, MPFR_RNDN);
	mpfr_mul_si(sol->re, sol->im, -2, MPFR_RNDN);
	mpfr_mul_si(sol->im, dum->re, 2, MPFR_RNDN);

	for (int i = 0; i < n; i++)
	{
		clear_complex(&(M[i][0]));
		clear_complex(&(M[i][1]));
		clear_complex(&(M[i][2]));
		clear_complex(&(M[i][3]));

		free(M[i]);
	}
	free(M);

	clear_complex(&Q[0]);
	clear_complex(&Q[1]);
	clear_complex(&Q[2]);
	clear_complex(&Q[3]);
	clear_complex(&q[0]);
	clear_complex(&q[1]);
	clear_complex(&ci);
}

// Input format:
// n_eigVs
// eig1
// eig2
//...
// coeff1
// coeff2
//...
// x_start x_step x_end
// t_start t_step t_end

int main(void)
{
	mpfr_set_default_prec(PRECISION);

	char output_path[] = "/Users/ivanpedrocr/Documents/soliton_computation.txt";
	complex double *eigVs;
	double *eigRe, *eigIm, *coeffs;
	int n_eigVs, t_length, x_length;
	double t_step, x_step, t_start, t_end, x_start, x_end;
	complex_nbr dum[3], sol;
	complex_nbr **init_matrix;

	FILE *in_ptr;
	FILE *out_ptr;
	FILE *log_ptr;

	clock_t init_time = clock();

	log_ptr = fopen("/Users/ivanpedrocr/Documents/soliton_log.txt", "w");

	init_complex(&sol);
	init_complex(dum);
	init_complex(dum + 1);
	init_complex(dum + 2);

	in_ptr = fopen("input.txt", "r");
	out_ptr = fopen(output_path, "w");

	fscanf(in_ptr, "%d", &n_eigVs);

	eigVs = calloc(n_eigVs, sizeof(complex double));
	eigRe = calloc(n_eigVs, sizeof(double));
	eigIm = calloc(n_eigVs, sizeof(double));
	coeffs = calloc(n_eigVs, sizeof(double));

	for (int i = 0; i < n_eigVs; i++)
	{
		fscanf(in_ptr, "%lf + %lfi", &(eigRe[i]), &(eigIm[i]));
		// printf("eig %d:%lf + %lfi\n", i + 1, eigRe[i], eigIm[i]);
		eigVs[i] = CMPLX(eigRe[i], eigIm[i]);
	}

	for (int i = 0; i < n_eigVs; i++)
	{
		fscanf(in_ptr, " %lf", coeffs + i);
		// printf("coeff %d:%lf\n", i + 1, coeffs[i]);
	}

	fscanf(in_ptr, "%lf %lf %lf", &x_start, &x_step, &x_end);
	fscanf(in_ptr, "%lf %lf %lf", &t_start, &t_step, &t_end);

	init_matrix = calloc(n_eigVs, sizeof(complex_nbr *));

	for (int i = 0; i < n_eigVs; i++)
	{
		init_matrix[i] = calloc(4, sizeof(complex_nbr));

		init_complex(&(init_matrix[i][0]));
		init_complex(&(init_matrix[i][1]));
		init_complex(&(init_matrix[i][2]));
		init_complex(&(init_matrix[i][3]));
	}

	// printf("x_start: %lf,x_step: %lf,x_end: %lf\n", x_start, x_step, x_end);
	// printf("t_start: %lf,t_step: %lf,t_end: %lf\n", t_start, t_step, t_end);

	x_length = floor((x_end - x_start) / x_step);
	t_length = floor((t_end - t_start) / t_step);

	for (int i = 0; i <= t_length; i++)
	{
		double t = t_start + i * t_step;

		for (int j = 0; j <= x_length; j++)
		{
			double x = x_start + j * x_step;

			get_init_matrix(init_matrix, x, eigVs, n_eigVs, dum);

			// ERROR WHEN t=10,x=30

			// printf("x,t in loop: %lf, %lf\n", x, t);

			darboux_transform(&sol, eigVs, coeffs, init_matrix, t, n_eigVs, dum, log_ptr);
			mpfr_fprintf(out_ptr, "%.16Rf%+.16Rf", sol.re, sol.im);
			fprintf(out_ptr, "i ");
		}
		fprintf(out_ptr, "\n");
	}

	clock_t end_time = clock();

	printf("TIME ELAPSED: %lf\n", (double)((end_time - init_time)) / CLOCKS_PER_SEC);

	for (int i = 0; i < n_eigVs; i++)
	{
		clear_complex(&(init_matrix[i][0]));
		clear_complex(&(init_matrix[i][1]));
		clear_complex(&(init_matrix[i][2]));
		clear_complex(&(init_matrix[i][3]));

		free(init_matrix[i]);
	}

	free(init_matrix);
	fclose(in_ptr);
	fclose(log_ptr);
	fclose(out_ptr);
	free(eigVs);
	free(coeffs);
	free(eigRe);
	free(eigIm);

	clear_complex(&sol);
	clear_complex(dum);
	clear_complex(dum + 1);
	clear_complex(dum + 2);

	return 0;
}
