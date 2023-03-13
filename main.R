library(psych)
library(tidyverse)
library(vegan)
library(twinspan)

data(varespec)
varspe <- varespec
cl <- factor(cut(twinspan(varespec), 2))

for( i in unique(cl)){
  A <- colSums((varspe == 0)*(cl != i)) # ABSENT outside
  B <- colSums((varspe == 0)*(cl == i)) # ABSENT IN CL
  C <- colSums((varspe != 0)*(cl != i)) # PRESENT outside
  D <- colSums((varspe != 0)*(cl == i)) # PRESENT IN CL
  phi_iter <- (((A*D)-(B*C))/sqrt((A+B)*(C+D)*(A+C)*(B+D)))

  if(i == unique(cl)[1]){
    phi_fin <- phi_iter
  }
  else{
    phi_fin <- cbind(phi_fin, phi_iter)
  }
}

colnames(phi_fin) <- unique(cl)
twintable <- mutate_if(as_tibble(rownames_to_column(data.frame(phi_fin), 'species')),
                       is.numeric,
                       function(x){
                         x = replace_na(x, 0)
                         x[x < 0] <- 0
                         x
                       }
)
th <- .2

no_n_gr <- bind_cols(species = names(varspe), bind_cols(map(split(data.frame(varspe != 0), cl), colSums))) %>%
  pivot_longer(-1, values_to = 'no') %>%
  mutate(name = paste0('X', name))
#order <-
  twintable %>%
  pivot_longer(cols = -c(species)) %>%
  arrange(value < th, name, -value) %>%
  left_join(no_n_gr) %>%
#    mutate(value = paste0(no, ' (', round(value, 2), ')'),
#           value = ifelse(grepl('0 \\(', value), '.', .)) %>%
    select(species, name, value) %>%
  pivot_wider() %>%
    print(n = 100)



sapply(varspe != 0, colSums)


syn_sort <- function(syntable, threshold = .2){
output_step <- syntable %>% pivot_longer(-1) %>% filter(abs(value) > threshold) %>% arrange(name, -abs(value)) %>%
  distinct(species)
left_join(output_step, syntable) %>% print(n = nrow(output_step))
}


syn_sort(twintable)