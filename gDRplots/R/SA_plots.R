
plot_SA_byCLs = function(SE, CL_names = NULL, Row_names = NULL, sel_normalization = 'GR', 
    colormap = NULL, plot_averaged = T, plot_fit = T) {
  
  if (is.null(Row_names)) {
    Row_names = rownames(SE)
  } else {
    stopifnot(length(Row_names) == nrow(SE))
    rownames(SE) = Row_names
  }  

  plt_list = list()
  for (iR in rownames(SE)) {

    if (is.null(CL_names)) {
      CL_names = colnames(SE)
    } else {
      stopifnot(length(CL_names) == ncol(SE))
      colnames(SE) = CL_names      
    }

    plt_list[[iR]] = grob_SA(
        convert_se_assay_to_dt(SE[iR,], 'Metrics'), 
        grouping = cId,
        df_average = convert_se_assay_to_dt(SE[iR,], 'Averaged'), 
        group_names = CL_names,
        sel_normalization = sel_normalization,
        colormap = colormap,
        plot_averaged = plot_averaged,
        plot_fit = plot_fit
      ) +
      ggtitle(iR) 
    
  }

  return(plt_list)

}


plot_SA_1CL = function(SE, sel_normalization = 'GR', 
    colormap = NULL, plot_averaged = T, plot_fit = T) {
  
  stopifnot(ncol(SE) == 1)

  
  grob_SA(
        convert_se_assay_to_dt(SE, 'Metrics'), 
        grouping = 'DrugName',
        df_average = convert_se_assay_to_dt(SE, 'Averaged'), 
        sel_normalization = sel_normalization,
        colormap = colormap,
        plot_averaged = plot_averaged,
        plot_fit = plot_fit
      ) 
  
}


grob_SA = function(df_metrics, grouping, df_average = NULL, group_names = NULL, sel_normalization = 'GR', 
    colormap = NULL, plot_averaged = !is.null(df_average), plot_fit = T) {

  df_metrics = df_metrics %>% filter(normalization_type == sel_normalization)
  df_average$NormData = if(sel_normalization == 'GR') {df_average$GRvalue} else {df_average$RelativeViability}

  data_range = c( min(min(df_average$NormData),0) -.05, max(max(df_average$NormData), 1) +.05 )
  min_conc = min(df_average$Concentration[df_average$Concentration>0])
  conc_range = .5*c(floor(2*log10(min_conc)-.5), ceiling(2*log10(max(df_average$Concentration))+.3) )    
  sel_conc = 10**(seq(conc_range[1], conc_range[2], .05))
  df_average$Concentration[df_average$Concentration==0] = min_conc/10

  if (is.null(group_names)) group_names = unique(df_metrics[[grouping]])

  df_fit = data.frame()
  if (plot_fit) {
    for (icol in group_names) {
      sel_metrics = as.data.frame(df_metrics[ df_metrics[[grouping]] == icol, ])
      df_fit = rbind(df_fit, cbind( sel_metrics[, grouping, drop = F],
        data.frame(
          Concentration = sel_conc,
          NormData = predict_efficacy_from_conc(sel_conc, sel_metrics$x_inf, sel_metrics$x_0, sel_metrics$ec50, sel_metrics$h)
        )))
    }
  }

  if (is.null(colormap)) {
    colormap = colorRampPalette(c("#dd0000", "#bbdd33", "#0000dd"))(length(group_names))
    names(colormap) = group_names
  }

  df_avg = data.frame()
  if (plot_averaged) {
    stopifnot(!is.null(df_average))
    df_avg = df_average    
  }

  df_avg$grouping = factor(df_avg[[grouping]], levels = group_names)
  df_fit$grouping = factor(df_fit[[grouping]], levels = group_names)

  ggplot(mapping = aes(x = log10(Concentration), y = NormData, color = grouping, group = grouping)) +
    geom_hline(yintercept = 0, color = '#555555') +
    geom_hline(yintercept = 1, color = '#555555') +
    geom_vline(xintercept = 0, color = '#555555') +
    geom_point(data = df_avg) +
    geom_line(data = df_fit) +
    scale_color_manual(values = colormap[names(colormap) %in% unique(df_fit$grouping)]) +
    coord_cartesian(xlim = conc_range, ylim = data_range) +
    scale_x_continuous(breaks = -5:2, labels = c('1e-5', '1e-4', 10**(-3:2))) +
    xlab('Concentration (µM)') +
    ylab(paste(sel_normalization, 'values')) +
    theme_bw()

}