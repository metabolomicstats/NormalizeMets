RainbowCols <- function(n)
{
    
    chlist <- as.hexmode(seq(0,n/5,1)*255)
    redylw<-paste('#ff', chlist, '00', sep='')
    ylwgrn<-paste('#', rev(chlist), 'ff', '00', sep='')
    grncyn<-paste('#00', 'ff', chlist, sep='')
    cynblu<-paste('#00', rev(chlist), 'ff', sep='')
    bluprp<-paste('#', chlist, '00', 'ff', sep='')

return(c(redylw,ylwgrn,grncyn,cynblu,bluprp))

}