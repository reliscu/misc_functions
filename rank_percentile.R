rank_percentile <- function(num_vec){
 (1-(rank(-num_vec, na.last="keep")/max(rank(-num_vec, na.last="keep"))))*100
}
