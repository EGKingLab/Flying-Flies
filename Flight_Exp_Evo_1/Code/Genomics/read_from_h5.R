
#' REad RILhap from h5 file
#'
#' @param fname String path to h5 file
#' @param test.positions Numeric vector of test positions (rows)
#' @param ww Numeric vector of RILs to select
#'
#' @return Array (dim = 3)
#' 
read_from_h5 <- function(fname,
                         test.positions,
                         ww) {
  require(rhdf5)
  require(abind)
  
  h5disableFileLocking()
  
  # Empty list to hold matrices
  m <- list()
  
  # Open for reading
  h5f <- H5Fopen(fname)
  
  # Loop through RILs, loading the matrix for each
  for (ii in 1:length(ww)) {
    iic <- paste0("RIL_", ww[ii])
    p <- h5f&iic
    
    # Subset p matrix by test.position rows
    m[[ii]] <- p[test.positions, ]
  }
  
  # Extract row, column (A1-A8), and RIL ID information
  dim1 <- as.character(h5f$dim1)[test.positions]
  dim2 <- as.character(h5f$dim2)
  dim3 <- as.character(h5f$dim3)[ww]
  
  # Bind matrices to 3d array
  a <- abind(m, along = 3)
  
  # Assign dim names
  dimnames(a) <- list(dim1, dim2, dim3)
  
  # Cloase h5 file
  h5closeAll()
  
  return(a)
}
