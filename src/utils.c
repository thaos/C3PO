#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


void rank_bycol(int nrow, int ncol, double *datamat, int*rankmat){
  double* onecolumn = malloc(sizeof *onecolumn *nrow);
  int* index = malloc(sizeof *index * nrow);


  for (int j = 0; j < ncol; j++) {
    for (int i = 0; i < nrow; i++) {
      onecolumn[i] = datamat[i + nrow * j];
      index[i] = i;
    }
    rsort_with_index (onecolumn, index, nrow);
    for(int i = 0; i < nrow; i++){
      int iend = i;
      while ((iend < nrow - 1) && onecolumn[index[iend]] == onecolumn[index[iend + 1]]) { 
        iend++;
      }
      for (int k = i; k <= iend; k++){
        rankmat[index[k] + nrow * j] = i+1;
      }
    }
  }
}

void sort_bycol(int nrow, int ncol, double* datamat, double* sortedmat){
  double* onecolumn = malloc(sizeof *onecolumn *nrow);

  for (int j = 0; j < ncol; j++) {
    for (int i = 0; i < nrow; i++) {
      onecolumn[i] = datamat[i + nrow * j];
    }
    R_rsort (onecolumn, nrow);
    for(int i = 0; i < nrow; i++){
      sortedmat[i + nrow * j] = onecolumn[i];
    }
  }
}

void scale_rank(int nrow, int ncol, int* mat, double scale_factor){
  for (int i = 0; i < (nrow * ncol); i++) {
      mat[i] = (int) (fround((mat[i] - 1) * scale_factor, 0) + 1);
  }
}

void shrink_bcsorted(int nrow_ref, int nrow_bc, int ncol_bc, double *bc_sorted[]){

  int oldrank = 0, newrank;
  double scale_factor = ((double) (nrow_ref - 1)) / ((double) (nrow_bc - 1));
  double* bc_shrunk = malloc(sizeof *bc_shrunk * nrow_ref * ncol_bc);
  for (int i = 1; i < nrow_bc; i++) {
    newrank = (int) fround(i * scale_factor, 0);
    if(newrank > oldrank){
      for (int j = 0; j < ncol_bc; j ++) {
        bc_shrunk[oldrank + nrow_ref * j] = (*bc_sorted)[i - 1 + nrow_bc * j];
      }
    }
    oldrank = newrank;
  }
  for (int j = 0; j < ncol_bc; j ++) {
    bc_shrunk[oldrank + nrow_ref * j] = (*bc_sorted)[nrow_bc - 1 + nrow_bc * j];
  }
  free(*bc_sorted);
  *bc_sorted = bc_shrunk;

}


void find_closestrank(
    int* bestmatch,
    double* bestdistance,
    int isearchstart,
    int ntimesteps,
    int nconddim,
    int iconddim[],
    int nrow_bc,
    int ranks_bc[],
    int nrow_ref,
    int ncol_ref,
    int ranks_ref[]
) {
  int ans = -1;
  double bestdist = -1;
  double currdist;
  for (int i = 0; i < (nrow_ref - ntimesteps + 1); i++) {
    currdist = 0;
    for (int j = 0; j < nconddim; j++) {
        for (int k = 0; k < ntimesteps; k++) { 
          currdist += pow(ranks_ref[k + i + iconddim[j] * nrow_ref] - ranks_bc[k + isearchstart + iconddim[j] * nrow_bc], 2);
        }
    }
    if(i == 0){
      bestdist = currdist;
      ans = 0;
    } else {
      if (currdist < bestdist){
        ans = i;
        bestdist = currdist;
      }
    }
  }
  *bestmatch=ans;
  *bestdistance=bestdist;
}

void search_byblock(
    int* nsearch,
    double r2d2_bc[],
    int *time_bestanalogue[],
    double *dist_bestanalogue[],
    int visited_time[],
    int nconddim,
    int iconddim[],
    int nrow_bc,
    int ranks_bc[],
    double sorted_bc[],
    int nrow_ref,
    int ncol_ref,
    int ranks_ref[],
    int lag_search,
    int lag_keep
){

  *nsearch =
    ((nrow_bc - lag_search - 1) / (lag_keep + 1)) +
    (((nrow_bc - lag_search - 1) % (lag_keep + 1)) > 0) +
    1;
  int *time_bestana = malloc(sizeof *time_bestana * *nsearch);
  double *dist_bestana = malloc(sizeof *dist_bestana * *nsearch);
  for(int i = 0; i < nrow_ref; i++){
    visited_time[i] = 0;
  }
  int blocksearch_start, blockkeep_start, blockkeep_end = 0;
  int rank_tmp, index_tmp, offset;
  double value_tmp;
  for (int isearch = 0; isearch < *nsearch; isearch++) {
    blockkeep_start = blockkeep_end;
    if (isearch == 0) {
      blockkeep_end +=  lag_search + 1;
    } else if (isearch == (*nsearch - 1)) {
      blockkeep_end = nrow_bc;
    } else {
      blockkeep_end += lag_keep + 1;
    }
    blocksearch_start = blockkeep_end - lag_search - 1;
    
    find_closestrank(
        time_bestana + isearch,
        dist_bestana + isearch,
        blocksearch_start,
        lag_search + 1,
        nconddim,
        iconddim,
        nrow_bc,
        ranks_bc,
        nrow_ref,
        ncol_ref,
        ranks_ref
    );

    offset = blockkeep_start - blocksearch_start;
    for(int i = 0; i < (blockkeep_end - blockkeep_start); i ++){
      for (int j = 0; j < ncol_ref; j ++) {
        rank_tmp = ranks_ref[time_bestana[isearch] + offset + i + j * nrow_ref] - 1;
        value_tmp = sorted_bc[rank_tmp + j * (nrow_ref < nrow_bc ? nrow_ref : nrow_bc) ];
        index_tmp = blockkeep_start + i + j * nrow_bc;
        r2d2_bc[index_tmp] = value_tmp;
      }
      visited_time[time_bestana[isearch] + offset + i]++;
    }
    time_bestana[isearch] += lag_search;
  }
  free(*time_bestanalogue);
  free(*dist_bestanalogue);
  *time_bestanalogue = time_bestana;
  *dist_bestanalogue = dist_bestana;
}
