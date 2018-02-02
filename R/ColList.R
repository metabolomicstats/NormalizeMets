# Colour list function
#
# @n An integer
# 

ColList<-function (n) 
{
  n <- round(n)
  if (n <= 15) {
    col_list <- c("#ee3333", "#3366aa", "#009872", "#982187", 
                  "#faa200", "#267aa4", "#910000", "#b56cfe", "#00b7ec", 
                  "#f36a18", "#534731", "#fdb5da", "#064650", "#b5dafe", 
                  "#000000")
  }
  else col_list <- NULL
  return(col_list)
}
