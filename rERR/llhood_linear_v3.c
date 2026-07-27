#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>

/*  CALL FROM R
res <- .Call("llhood_linear_v3",beta_2,length(beta_2),length(rsets),
             rsets,nrseti,data,rows_cases,n_lin_vars,n_loglin_vars,constr_ind)
*/

SEXP llhood_linear_v3(SEXP PAR, SEXP NPAR, SEXP NCASES,
                      SEXP RSETS, SEXP NRSETI, SEXP DATA, SEXP ROWS_CASES,
                      SEXP NLIN, SEXP NLOGLIN, SEXP CONSTR_IND)
{
  double *res, sum, sumC, sum_lin, sum_loglin, sum_lin_j, *nrseti, *par, *nrow_case, cons = -9999999.00, doseC, doseRS;
  int *ncases = INTEGER(NCASES), *npar = INTEGER(NPAR), *nlin = INTEGER(NLIN), *nloglin = INTEGER(NLOGLIN), *cons_ind = INTEGER(CONSTR_IND);
  int i, j, ind_in_riskset, n_row;

  SEXP RES, parl;
  PROTECT(RES  = allocVector(REALSXP, *ncases + 1));
  PROTECT(parl = allocVector(REALSXP, *npar));

  res = REAL(RES);
  par = REAL(parl);

  /* save parameter vector into a c vector */
  for (j = 0; j < *npar; j++) {
    par[j] = REAL(PAR)[j];
  }

  /* number of cases loop */
  for (i = 0; i < *ncases; i++) {
    sum       = 0.0;
    sumC      = 0.0;
    sum_lin   = 0.0;
    sum_loglin= 0.0;
    doseC     = 0.0;
    doseRS    = 0.0;
    nrseti    = REAL(VECTOR_ELT(NRSETI, i));
    nrow_case = REAL(VECTOR_ELT(ROWS_CASES, i));

    /* Subject in the relative case risk set loop */
    for (ind_in_riskset = 0; ind_in_riskset < *nrseti; ind_in_riskset++) {
      n_row = REAL(VECTOR_ELT(RSETS, i))[ind_in_riskset] - 1;

      sum_lin    = 0.0;
      sum_loglin = 0.0;
      sum_lin_j  = 0.0;

      /* linear variables */
      for (j = 0; j < *nlin; j++) {
        sum_lin += par[j] * REAL(VECTOR_ELT(DATA, j + 7))[n_row];
        if (j == *cons_ind)
          sum_lin_j = par[j] * REAL(VECTOR_ELT(DATA, j + 7))[n_row];
      }

      /* loglinear variables */
      for (j = 0; j < *nloglin; j++) {
        sum_loglin += par[*nlin + j] * REAL(VECTOR_ELT(DATA, 7 + *nlin + j))[n_row];
      }

      /* linERR model */
      sum += exp(sum_loglin) * (1 + sum_lin);

      doseRS += sum_lin_j / par[*cons_ind];

      /* when it is the case */
      if (n_row == *nrow_case - 1) {
        sumC  = exp(sum_loglin) * (1 + sum_lin);
        doseC = sum_lin_j / par[*cons_ind];
      }
    }

    /* Each risk set contribution to the log likelihood function */
    res[i] = sumC / sum;

    /* constraint of the model with respect the CONSTR_IND-th linear parameter */
    if (cons <= -1 / doseC)
      cons = -1 / doseC;
    if (cons <= -(*nrseti) / doseRS)
      cons = -(*nrseti) / doseRS;
  }

  res[*ncases] = cons;
  UNPROTECT(2);
  return RES;
}