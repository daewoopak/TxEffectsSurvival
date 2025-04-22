fn_ds <- function(hcen, dcen, z) {
  
  uni_z <- unique(z)
  h_only <- c(sum(hcen[z == uni_z[2]]), sum(hcen[z == uni_z[1]]), sum(hcen))
  h_and_d <- c(sum((hcen * dcen)[z == uni_z[2]]), sum((hcen * dcen)[z == uni_z[1]]), sum(hcen * dcen))
  d_only <- c(sum(dcen[z == uni_z[2]]), sum(dcen[z == uni_z[1]]), sum(dcen))

  ds_table <- rbind(h_only, h_and_d, d_only)
  colnames(ds_table) <- c(paste("z=", uni_z[2], sep = ""), paste("z=", uni_z[1], sep = ""), "all")  
  ds_table
  
}