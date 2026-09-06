ZY.LCPA <- function(..., method.3step = c("BCH", "ML")) {
  method.3step <- match.arg(method.3step)
  if(method.3step == "BCH"){
    BCH.ZY.LCPA(...)
  }else{
    ML.ZY.LCPA(...)
  }
}
