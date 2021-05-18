#' @import tidyverse
#' @import reshape2
#' @export
medianIqr_scaling = function(X, bx, refbatch = NULL){
  rn = rownames(X)
  cn = colnames(X)

  colnames(X) = 1:ncol(X)


  if(length(bx) != nrow(X)) stop("dimension mismatch between batch variable and data matrix")

  df = X %>%
    as.data.frame() %>%
    mutate(.bx.indicator = factor(bx), .r = 1:nrow(X)) %>%
    pivot_longer(cols = !all_of(c(".bx.indicator", ".r")), names_to = ".c")

  if(!is.null(refbatch)){
    df$.bx.indicator = relevel(df$.bx.indicator, ref = refbatch)
  }

  df.sum = df %>%
    group_by(.bx.indicator, .c) %>%
    summarise(medval = quantile(value, .5),
              iqrval = diff(quantile(value, c(.25,.75))) ) %>%
    ungroup()

  ref.sum = df.sum %>%
    filter(.bx.indicator == levels(.bx.indicator)[1]) %>%
    rename(refmed = medval, refiqr = iqrval) %>%
    select(-.bx.indicator)

  df.sum = df.sum %>%
    left_join(ref.sum, by = ".c")

  Z = df %>%
    left_join(df.sum, by = c(".bx.indicator", ".c")) %>%
    mutate(scaled.val = refiqr*((value - medval)/iqrval) + refmed) %>%
    acast(.r ~ .c, value.var = "scaled.val")

  rownames(Z) = rn
  colnames(Z) = cn

  return(Z)
}
