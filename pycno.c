/******************* Smooth Pycnophylactic Interpolation *******************/

/* factor to use for underrelaxation */
#define UNDERRELAX (0.25)
/* number of iterations */
#define PYC_ITER (100)
/* new densities can be at minimum this proportion of their original value */
#define MIN_DENS_RATIO (0.50)

/* calculate the adjustment */
double smooth_adj(double center, double right, double left, double down,
		double up) {
	double delta;
	delta = -center + 0.25 * (right + left + down + up);
	/* underrelax the adjustment made */
	delta *= UNDERRELAX;
	return (delta);
}

/* process the adjustments for each point */
void pycno_adj(double rho_adj[], double rho_init[], int lx, int ly, int k) {
	int i, j;
	for (i = 1; i < lx - 1; i++)
		for (j = 1; j < ly - 1; j++) {
			rho_adj[i * ly + j] = smooth_adj(rho_init[i * ly + j],
					rho_init[i * ly + j + 1], rho_init[i * ly + j - 1],
					rho_init[(i + 1) * ly + j], rho_init[(i - 1) * ly + j]);
		}
	return;
}

/* cumulate the adjustments by region */
void pycno_sum(double rho_adj[], int n_reg, double s_k[]) {
	int i, j, k;
	int area;
	for (k = -1; k < n_reg; k++) {
		area = 0;
		s_k[k + 1] = 0;
		for (i = 1; i < lx - 1; i++) {
			for (j = 1; j < ly - 1; j++) {
				if (inside[i][j] == k) {
					/* cumulate the region's adjustment values in s_sum */
					s_k[k + 1] += rho_adj[i * ly + j];
					area++;
				}
			}
		}
		/* compute the decrementing factor */
		if (area == 0)
			s_k[k + 1] = 0;
		else
			s_k[k + 1] = -(s_k[k + 1] / area);
	}
	return;
}

/* apply the adjustments with regard to non-negativity constraint */
void pycno_apply(double rho_adj[], double rho_init[], double rho_new[],
		int lx, int ly, int n_reg, double s_k[], double rho_master[]) {
	int i, j;
	for (i = 0; i < lx; i++)
		for (j = 0; j < ly; j++) {
			if (i == 0 || i == lx - 1 || j == 0 || j == ly - 1)
				rho_new[i * ly + j] = rho_init[i * ly + j];
			else if (rho_init[i * ly + j] + rho_adj[i * ly + j]
					+ s_k[inside[i][j] + 1]
					>= MIN_DENS_RATIO * rho_master[i * ly + j])
				rho_new[i * ly + j] = rho_init[i * ly + j]
						+ rho_adj[i * ly + j] + s_k[inside[i][j] + 1];
			else
				rho_new[i * ly + j] = rho_init[i * ly + j];
		}
	return;
}

/* compute difference between initial and new region populations */
void pycno_diff(double rho_init[], double rho_new[], int n_reg, double l_k[],
		int lx, int ly) {
	int i, j, k;
	int area;
	double cum_pop_init, cum_pop_new;
	for (k = -1; k < n_reg; k++) {
		area = 0;
		cum_pop_init = 0;
		cum_pop_new = 0;
		l_k[k + 1] = 0;
		for (i = 0; i < lx; i++) {
			for (j = 0; j < ly; j++) {
				if (inside[i][j] == k) {
					cum_pop_init += rho_init[i * ly + j];
					cum_pop_new += rho_new[i * ly + j];
					area++;
				}
			}
		}
		/* average difference between initial density and new density*/
		if (area == 0)
			l_k[k + 1] = 0;
		else {
			l_k[k + 1] = (cum_pop_init - cum_pop_new) / area;
		}
	}
}

/* apply the corrections to maintain pycnophylactic property */
void pycno_apply_final(double rho_new[], int lx, int ly, double l_k[],
		int examined[], double residual[], double rho_master[]) {
	int i, j, k;
	for (k = -1; k < n_reg; k++) {
		residual[k + 1] = 0;
		examined[k + 1] = 0;
	}
	for (i = 1; i < lx - 1; i++) {
		for (j = 1; j < ly - 1; j++) {
			if (rho_new[i * ly + j] + l_k[inside[i][j] + 1]
					>= MIN_DENS_RATIO * rho_master[i * ly + j]) {
				rho_new[i * ly + j] += l_k[inside[i][j] + 1];
				/* count of points that are non-zero */
				examined[inside[i][j] + 1]++;
			} else {
				/* add to the residual correction */
				residual[inside[i][j] + 1] += l_k[inside[i][j] + 1]
						+ rho_new[i * ly + j]
						- MIN_DENS_RATIO * rho_master[i * ly + j];
				rho_new[i * ly + j] = MIN_DENS_RATIO
						* rho_master[i * ly + j];
			}
		}
	}
	return;
}

/* distribute the residual correction across the region */
void pycno_residual(double rho_new[], int lx, int ly, double l_k[],
		int examined[], double residual[], double rho_master[]) {
	int i, j, k;
	/* initialise l_k */
	for (k = -1; k < n_reg; k++) {
		/* calculate the new corrections */
		if (examined[k + 1] == 0)
			l_k[k + 1] = 0;
		else
			l_k[k + 1] = residual[k + 1] / examined[k + 1];
		/* reset counters for residual */
		residual[k + 1] = 0;
		examined[k + 1] = 0;
	}
	for (i = 1; i < lx - 1; i++) {
		for (j = 1; j < ly - 1; j++) {
			/* residual subtraction is distributed across non-zero points */
			if (rho_new[i * ly + j] + l_k[inside[i][j] + 1]
					>= MIN_DENS_RATIO * rho_master[i * ly + j]) {
				rho_new[i * ly + j] += l_k[inside[i][j] + 1];
				/* count of points that are non-zero */
				examined[inside[i][j] + 1]++;
			} else {
				/* add to the residual correction */
				residual[inside[i][j] + 1] += l_k[inside[i][j] + 1]
						+ rho_new[i * ly + j]
						- MIN_DENS_RATIO * rho_master[i * ly + j];
				rho_new[i * ly + j] = MIN_DENS_RATIO
						* rho_master[i * ly + j];
			}
		}
	}
	return;
}

/* copy the modified density back into original array */
void pycno_copy(double rho_new[], double rho_init[], int lx, int ly) {
	int i, j;
	for (i = 0; i < lx; i++)
		for (j = 0; j < ly; j++)
			rho_init[i * ly + j] = rho_new[i * ly + j];
	return;
}

/* conduct the pycnophylactic interpolation */
void pycno(double *rho_master, double *rho_adj, double *s_k, double *rho_new,
		double *l_k, double *residual, int *examined) {
	/* arrays should already be initialized according to usage purpose */
	int i, j, k;
	/* master copy of rho_init */
	for (i = 0; i < lx; i++)
		for (j = 0; j < ly; j++)
			rho_master[i * ly + j] = rho_init[i * ly + j];
	for (k = 0; k < PYC_ITER; k++) {
		printf("##### Pyc. Int. iter. %d #####\n", k);
		/* STEP 1 */
		pycno_adj(rho_adj, rho_init, lx, ly, k);
		/* STEP 2 */
		pycno_sum(rho_adj, n_reg, s_k);
		/* STEP 3 */
		pycno_apply(rho_adj, rho_init, rho_new, lx, ly, n_reg, s_k,
				rho_master);
		/* STEP 4 */
		pycno_diff(rho_init, rho_new, n_reg, l_k, lx, ly);
		/* STEP 5 */
		pycno_apply_final(rho_new, lx, ly, l_k, examined, residual,
				rho_master);
		pycno_residual(rho_new, lx, ly, l_k, examined, residual, rho_master);
		pycno_copy(rho_new, rho_init, lx, ly);
		/* repeat until stopping rules are met */
	}
}
