#' @title epiView function for scatter plot
#' @description TO ADD
#' @description This is a new line ...
#' @details What's this?
#' @return TO ADD
#' @export


get_comb <- function(comb_ix, y, u_ct, combs) {
  comp = paste(u_ct[combs[1,comb_ix]], "Vs.", u_ct[combs[2,comb_ix]])
  lhs = filter(y, Var1==u_ct[combs[1,comb_ix]])
  rhs = filter(y, Var1==u_ct[combs[2,comb_ix]])
  return(data.frame(comp, merge(lhs, rhs, by="Var2")[,c(1,3,5)]))
}

