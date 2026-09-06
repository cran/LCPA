XZ.LCPA <- function(..., method.3step = c("ML", "BCH")) {
  method.3step <- match.arg(method.3step)
  if(method.3step == "ML"){
    ML.XZ.LCPA(...)
  }else{
    BCH.XZ.LCPA(...)
  }
}
