format_power <- function(x, values_to = "power") {
    as.data.frame(x) %>% 
        tibble::rownames_to_column(var = "selection") %>%
        tidyr::pivot_longer(cols = !`selection`, values_to = values_to) %>%
        dplyr::select(!name) %>% 
        dplyr::filter(is.finite(.data[[values_to]]))  ## only keep those entries with non 0 oracle TP (ie condition on |S \cap H1|>0)
}
