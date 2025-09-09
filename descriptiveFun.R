# Determine RSV activity by AAP
getAAPEpi <- function(site, RSVn, threshold) {
  site$AAP <- site[[RSVn]]/sum(site[[RSVn]])*100
  site <- site[order(site$AAP, decreasing = TRUE),]
  site$cumAAP <- cumsum(site$AAP)
  site <- site[order(site$Month, decreasing = FALSE),]
  site$Epi <- site$cumAAP <= threshold
  return(site)
}


