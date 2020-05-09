library(ggplot2)
library(cowplot)
library(matrixStats)
library(ggpointdensity)

hypervar <- function(data, span = 0.5, showplot = TRUE, font_size = 14){

  gene_mean_all <- rowMeans(data)
  gene_var_all <- rowVars(data)

  data_filter <- data[gene_mean_all > 0 & gene_var_all > 0,]

  gene_mean <- rowMeans(data_filter)
  gene_var <- rowVars(data_filter)

  data_fit <- data.frame(X=gene_mean,Y=gene_var)
  fit_model <- loess(formula = log2(x=Y) ~ log2(x=X),
                     data = data_fit,
                     span = span)

  gene_var_expect <- fit_model$fitted
  gene_hyper_var <- log2(gene_var) - gene_var_expect

  result <- data.frame(feature=row.names(data_filter), mean=gene_mean, var=gene_var,
                       var_expect_log2=gene_var_expect,hypervar_log2=gene_hyper_var)

  p1 <- ggplot(result, aes(log2(mean), log2(var))) +
    geom_pointdensity() +
    scale_color_viridis_c() +
    geom_point(data=result,aes(log2(mean),var_expect_log2),color="red") +
    theme_bw() +
    theme(axis.text = element_text(size=font_size), axis.title = element_text(size=font_size))

  p2 <- ggplot(result, aes(log2(mean), hypervar_log2)) +
    geom_pointdensity() +
    scale_color_viridis_c() +
    theme_bw() +
    theme(axis.text = element_text(size=font_size), axis.title = element_text(size=font_size))

  combined_plot <- plot_grid(p1, p2, labels = c('A', 'B'))

  if(showplot){
    print(combined_plot)
  }

  return(list(data=result,plot=combined_plot))
}

