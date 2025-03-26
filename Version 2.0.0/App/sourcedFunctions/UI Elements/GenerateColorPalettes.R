gen_color_pallette <- function(var, seed){
  
  base <- if((seed %% 2) == 0)
  {c("cornsilk3", "orchid")} else {c("cadetblue", "indianred")}
  
  u <- unique(var[,1])
  
  if(length(u) > 2){
    return( setNames( colorRampPalette(c(base[1], "white", base[2]))(length(u)), u ))
  }else{
    return( setNames( colorRampPalette(base)(length(u)), u ))
  }
  
}
