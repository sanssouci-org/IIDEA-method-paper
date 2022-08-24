format_power <- function(x) {
    as.data.frame(x) %>% 
        tibble::rownames_to_column(var = "selection") %>%
        tidyr::pivot_longer(cols = !`selection`, values_to = "power") %>%
        dplyr::select(!name) %>% 
        dplyr::filter(is.finite(power))  ## only keep those entries with non 0 oracle TP (ie condition on |S \cap H1|>0)
}
