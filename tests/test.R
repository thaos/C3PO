library(C3PO)

mat <- matrix(
  rep(1:10, 2),
  ncol  = 2 
)
print(mat)

mat2 <- matrix(
  rep(1:15, 2),
  ncol  = 2 
)
print(mat2)

r2d2_bc = r2d2(mat, mat, 1, 0, 0)
stopifnot(
  all.equal(mat, r2d2_bc$r2d2_bc)
)
stopifnot(
  all.equal(r2d2_bc$visited_time, rep(1, nrow(mat)))
)

r2d2_bc = r2d2(mat, mat, 2, 0, 0)
stopifnot(
  all.equal(mat, r2d2_bc$r2d2_bc)
)
stopifnot(
  all.equal(r2d2_bc$visited_time, rep(1, nrow(mat)))
)

r2d2_bc = r2d2(mat, mat2, 1, 0, 0)
stopifnot(
  all.equal(
    matrix(rep(c(1, 3, 3, 4, 6, 6, 8, 8, 9, 11, 11, 12, 14, 14, 15), 2), ncol = 2),
    r2d2_bc$r2d2_bc
  )
)

r2d2_bc = r2d2(mat2, mat, 1, 0, 0)
stopifnot(
  all.equal(mat, r2d2_bc$r2d2_bc)
)

