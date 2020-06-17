#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "utils.h"

SEXP r2d2(SEXP refdata, SEXP bc1d, SEXP iconddim, SEXP lag_search, SEXP lag_keep){
  int nrow_sbc, nrow_bc, ncol_bc, nrow_ref, ncol_ref, nconddim;
  int lsearch, lkeep;
  int *icdim;
  double scaling_factor;

  SEXP dim_bc = getAttrib(bc1d, R_DimSymbol);
  nrow_sbc = nrow_bc = INTEGER(dim_bc)[0];
  ncol_bc = INTEGER(dim_bc)[1];
  bc1d = coerceVector(bc1d, REALSXP);

  SEXP dim_ref = getAttrib(refdata, R_DimSymbol);
  nrow_ref = INTEGER(dim_ref)[0];
  ncol_ref = INTEGER(dim_ref)[1];
  refdata = coerceVector(refdata, REALSXP);

  lag_search = coerceVector(lag_search, INTSXP);
  lag_keep = coerceVector(lag_keep, INTSXP);
  lsearch = INTEGER(lag_search)[0];
  lkeep = INTEGER(lag_keep)[0];

  iconddim = coerceVector(iconddim, INTSXP);
  icdim = INTEGER(iconddim);
  nconddim = length(iconddim);
  
  int* ranks_ref = malloc(sizeof *ranks_ref * nrow_ref * ncol_ref);
  int* ranks_bc = malloc(sizeof *ranks_bc * nrow_bc * ncol_bc);
  double* sorted_bc = malloc(sizeof *sorted_bc *nrow_bc * ncol_bc);
  
  rank_bycol(nrow_bc, ncol_bc, REAL(bc1d), ranks_bc);
  rank_bycol(nrow_ref, ncol_ref, REAL(refdata), ranks_ref);
  sort_bycol(nrow_bc, ncol_bc, REAL(bc1d), sorted_bc);

  
  if (nrow_bc < nrow_ref) {
    scaling_factor =  ((double) nrow_bc - 1) / ((double) nrow_ref - 1);
    scale_rank(nrow_ref, ncol_ref, ranks_ref, scaling_factor);
  } else if (nrow_bc > nrow_ref) {
    scaling_factor =  ((double) nrow_ref - 1) / ((double) nrow_bc - 1);
    scale_rank(nrow_bc, ncol_bc, ranks_bc, scaling_factor);
    shrink_bcsorted(nrow_ref, nrow_bc, ncol_bc, &sorted_bc);
    nrow_sbc = nrow_ref;
  }
 
  double* r2d2_bc = malloc(sizeof *r2d2_bc * nrow_bc * ncol_bc);
  int *time_bestanalogue = NULL;
  double *dist_bestanalogue = NULL;
  int* visited_time = malloc(sizeof *visited_time * nrow_ref);
  int nsearch;
  search_byblock(
    &nsearch,
    r2d2_bc,
    &time_bestanalogue,
    &dist_bestanalogue,
    visited_time,
    nconddim,
    icdim,
    nrow_bc,
    ranks_bc,
    sorted_bc,
    nrow_ref,
    ncol_ref,
    ranks_ref,
    lsearch,
    lkeep
    );

  SEXP dim_sbc = PROTECT(allocVector(INTSXP, 2));
  INTEGER(dim_sbc)[0] = nrow_sbc;
  INTEGER(dim_sbc)[1] = ncol_bc;

  SEXP ranks_ref_R = PROTECT(allocVector(INTSXP, nrow_ref * ncol_ref));
  SEXP ranks_bc_R = PROTECT(allocVector(INTSXP, nrow_bc * ncol_bc));
  SEXP sorted_bc_R = PROTECT(allocVector(REALSXP, nrow_sbc * ncol_bc));
  SEXP r2d2_bc_R = PROTECT(allocVector(REALSXP, nrow_bc * ncol_bc));
  setAttrib(ranks_ref_R, R_DimSymbol, dim_ref);
  setAttrib(ranks_bc_R, R_DimSymbol, dim_bc);
  setAttrib(sorted_bc_R, R_DimSymbol, dim_sbc);
  setAttrib(r2d2_bc_R, R_DimSymbol, dim_bc);

  SEXP time_bestanalogue_R = PROTECT(allocVector(INTSXP, nsearch));
  SEXP dist_bestanalogue_R = PROTECT(allocVector(REALSXP, nsearch));
  SEXP visited_time_R = PROTECT(allocVector(INTSXP, nrow_ref));

  for(int i = 0; i < (nrow_ref * ncol_ref); i ++){
    INTEGER(ranks_ref_R)[i] = ranks_ref[i];
  }
  for(int i = 0; i < (nrow_bc * ncol_bc); i ++){
    INTEGER(ranks_bc_R)[i] = ranks_bc[i];
    REAL(r2d2_bc_R)[i] = r2d2_bc[i];
  }
  for(int i = 0; i < (nrow_sbc * ncol_bc); i ++){
    REAL(sorted_bc_R)[i] = sorted_bc[i];
  }
  
  for(int i = 0; i < nsearch; i ++){
    INTEGER(time_bestanalogue_R)[i] = time_bestanalogue[i];
    REAL(dist_bestanalogue_R)[i] = dist_bestanalogue[i];
  }
  for(int i = 0; i < nrow_ref; i ++){
    INTEGER(visited_time_R)[i] = visited_time[i];
  }
  
  SEXP output = PROTECT(allocVector(VECSXP, 7));
  SET_VECTOR_ELT(output, 0, ranks_ref_R);
  SET_VECTOR_ELT(output, 1, ranks_bc_R);
  SET_VECTOR_ELT(output, 2, sorted_bc_R);
  SET_VECTOR_ELT(output, 3, r2d2_bc_R);
  
  SET_VECTOR_ELT(output, 4, time_bestanalogue_R);
  SET_VECTOR_ELT(output, 5, dist_bestanalogue_R);
  SET_VECTOR_ELT(output, 6, visited_time_R);
  
  UNPROTECT(9);
  free(ranks_ref);
  free(ranks_bc);
  free(sorted_bc);
  free(r2d2_bc);
  free(time_bestanalogue);
  free(dist_bestanalogue);
  free(visited_time);
  return output;

}


