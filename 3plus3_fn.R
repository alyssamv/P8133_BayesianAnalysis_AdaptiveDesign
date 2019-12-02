
performance_3plus3 = function(coh.size, dlt) {
  
  dose = 1
  esc = 0
  mtd = 0
  
  while (esc == 0 && dose <= length(dlt)) {
    dlt.dose_3 = rbinom(1, coh.size, dlt[[dose]])
    
    if (dlt.dose_3 == 0) {
      dose = dose + 1
    } else if (dlt.dose_3 > 1) {
      mtd = dose - 1
      esc = 1
    } else if (dlt.dose_3 == 1) {
      dlt.dose_6 = rbinom(1, coh.size, dlt[[dose]])
      dlt_total = dlt.dose_3 + dlt.dose_6
      
      if (dlt_total == 1 ) {
        dose = dose + 1
      } else {
        mtd = dose - 1
        esc = 1
      }
    }
    
    if (dose >= length(dlt) && mtd == 0){
      mtd = length(dlt)
    }
    
  }
  return(mtd)

}

performance_3plus3(coh.size = 3, dlt = c(0.02, 0.04, 0.10, 0.25, 0.5))

N = 1e3
#set.seed(1)
sims = sapply(1:N, function(i){
  performance_3plus3(coh.size = 3, dlt = c(0.017, 0.043, 0.10, 0.22, 0.41))
})

table(sims)/N
