nacho_fastCrw <- function(k, w, ns, maxx) {
  data_frame(vtm = rexp(ns,w), 
             phi = rvonmises(ns, 0, k),
             d   = cumsum(phi) ) %>%
    mutate(
      t = cumsum(vtm), 
      aux.x = cos(d)*vtm, 
      aux.y = sin(d)*vtm,
      x = cumsum(aux.x), 
      y = cumsum(aux.y)
    ) %>% 
    filter( t <= maxx ) %>%
    select(-aux.x, -aux.y)
}
  