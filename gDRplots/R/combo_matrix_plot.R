

CL = 'MCF-7'

DrugPair = 'GDC−9545 x Everolimus'

sel_row = 1
sel_CL = 1
normalization_type = 'GR'

se1 = mae[['matrix']][sel_row, sel_CL]

y_drug = rowData(se1)$DrugName 
x_drug = rowData(se1)$DrugName_2

mx_names = array(c("SmoothMatrix", 'HSAExcess', 'BlissExcess'), dimnames = list(c("SmoothMatrix", 'HSAExcess', 'BlissExcess')))

df_mx = lapply(mx_names, function(mx_name) {
  df_ = convert_se_assay_to_dt(se1, mx_name) %>% select(-rId, -cId, -parental_identifier)
  if (mx_name == 'SmoothMatrix') {
    df_$value = df_[[normalization_type]]
  } else {
    df_$value = df_$excess
    df_ = df_ %>% filter(normalization_type == normalization_type)
  }  
  df_$value = pmin(1.1,df_$value)

})
df_iso = convert_se_assay_to_dt(se1, "isobolograms")

plts = lapply(mx_names, function(mx_name) {
  df_ = df_mx[[mx_name]]

  x_conc = sort(unique(df_$Concentration_2[df_$Concentration_2>0]))
  x_gap = diff(log10(x_conc))[2]
  df_$x_pos = pmax( log10(df_$Concentration_2), log10(min(df_$Concentration_2[df_$Concentration_2>0])) - x_gap -.1)
  
  y_conc = sort(unique(df_$Concentration[df_$Concentration>0]))
  y_gap = diff(log10(y_conc))[2]
  df_$y_pos = pmax( log10(df_$Concentration), log10(min(df_$Concentration[df_$Concentration>0])) - y_gap -.1)
  
  ggplot(df_, aes(x = x_pos, y = y_pos)) +
    geom_tile(aes(fill=value),
          height = y_gap, width = x_gap) +
    labs(x = paste(x_drug, '[µM]'),
          y = paste(y_drug, '[µM]'),
          title = paste0(colData(se1)$CellLineName, ' (', colData(se1)$clid, '):',
                    mx_name, ' for ', normalization_type)) +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
                      axis.text.y=element_text(size=9),
                      plot.title=element_text(size=11)) +
    scale_x_continuous(breaks = drug2_axis$pos_x, labels = drug2_axis$marks_x,
        expand = c(0,0)) +
    scale_y_continuous(breaks = drug1_axis$pos_y, labels = drug1_axis$marks_y,
        expand = c(0,0)) +
    scale_shape_discrete(name = paste0(
      ifelse(norm_method == 'GRvalue', 'GR', 'IC'), '50'))
}

})
  plt = ggplot(df_, aes(x = pos_x, y = pos_y)) +
    geom_tile(aes(fill=value),
          height = diff(drug1_axis$pos_y[3:4]), width = diff(drug2_axis$pos_x[3:4])) +
    geom_abline(slope = 1, intercept = log10(ref_x50['conc_1']/ref_x50['conc_2'])) +
    labs(x = paste(condition['DrugName_2'], '[µM]'),
          y = paste(condition['DrugName'], '[µM]'),
          title = paste0(condition['CellLineName'], ' (', condition['CLID'], '):',
                    gsub('hsa', 'HSA', gsub('_', ' ', i)), ' for ', 
                    ifelse(norm_method == 'GRvalue', 'GR', 'RV'),
                    ', T=', condition['Duration'], 'h')) +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
                      axis.text.y=element_text(size=9),
                      plot.title=element_text(size=11)) +
    scale_x_continuous(breaks = drug2_axis$pos_x, labels = drug2_axis$marks_x,
        expand = c(0,0)) +
    scale_y_continuous(breaks = drug1_axis$pos_y, labels = drug1_axis$marks_y,
        expand = c(0,0)) +
    scale_shape_discrete(name = paste0(
      ifelse(norm_method == 'GRvalue', 'GR', 'IC'), '50'))
}


    

      all_mx = all_combo_variables[[idc, iCL]]$all_mx
      all_iso = all_combo_variables[[idc, iCL]]$all_iso
      df_CI_100x = all_combo_variables[[idc, iCL]]$df_CI_100x
      drug1_axis = all_combo_variables[[idc, iCL]]$drug1_axis
      drug2_axis = all_combo_variables[[idc, iCL]]$drug2_axis
      condition = all_combo_variables[[idc, iCL]]$condition
      norm_method = all_combo_variables[[idc, iCL]]$norm_method
      ref_x50 = all_combo_variables[[idc, iCL]]$ref_x50
      
      iso_cutoff = as.numeric(names(all_iso))
      if (norm_method == 'GRvalue') {
          iso_colors = sapply(iso_cutoff, function(x)
            sprintf('#%s',paste(as.hexmode(c(70, round((.85-x*.7)*170), round((1.1-x*.7)*200))),
                collapse = '')))
        } else {
            iso_colors = sapply(iso_cutoff, function(x)
              sprintf('#%s',paste(as.hexmode(c(70, round((1-x*.85)*170), round((1.1-x*.85)*240))),
                  collapse = '')))
        }
      names(iso_colors) = iso_cutoff

      # plots the different matrices of interest (measured data,  mean data, HSA excess and Bliss excess)
      plts = lapply(c("mx_response", "mx_full", 'hsa_excess', 'bliss_excess'), function(i) {
        df_ = melt(all_mx[[i]])
        df_$value = pmin(1.1,df_$value)
        colnames(df_)[1:2] = c('conc_1', 'conc_2')
        df_[,1] = round(as.numeric(df_[,1]),5)
        df_[,2] = round(as.numeric(df_[,2]),5)
        df_ = merge(df_, drug1_axis, by = 'conc_1', all.x = F)
        df_ = merge(df_, drug2_axis, by = 'conc_2', all.x = F)
        
        plt = ggplot(df_, aes(x = pos_x, y = pos_y)) +
          geom_tile(aes(fill=value),
                height = diff(drug1_axis$pos_y[3:4]), width = diff(drug2_axis$pos_x[3:4])) +
          geom_abline(slope = 1, intercept = log10(ref_x50['conc_1']/ref_x50['conc_2'])) +
          labs(x = paste(condition['DrugName_2'], '[µM]'),
                y = paste(condition['DrugName'], '[µM]'),
                title = paste0(condition['CellLineName'], ' (', condition['CLID'], '):',
                          gsub('hsa', 'HSA', gsub('_', ' ', i)), ' for ', 
                          ifelse(norm_method == 'GRvalue', 'GR', 'RV'),
                          ', T=', condition['Duration'], 'h')) +
          theme_bw() + theme(axis.text.x=element_text(size=9, angle=45, vjust=1, hjust = 1),
                            axis.text.y=element_text(size=9),
                            plot.title=element_text(size=11)) +
          scale_x_continuous(breaks = drug2_axis$pos_x, labels = drug2_axis$marks_x,
              expand = c(0,0)) +
          scale_y_continuous(breaks = drug1_axis$pos_y, labels = drug1_axis$marks_y,
              expand = c(0,0)) +
          scale_shape_discrete(name = paste0(
            ifelse(norm_method == 'GRvalue', 'GR', 'IC'), '50'))

        if ( !(i %in% c('hsa_excess', 'bliss_excess'))) { # heatmaps with readout values
          if(norm_method == 'GRvalue') {
            plt = plt +
              scale_fill_gradientn(colors = c('black', '#b06000', '#c07700', 'white'),
                  values = c(0, .59/1.7, .61/1.7, 1), limits = c(-.6,1.1), name = 'GR val',
                oob = scales::squish)
          } else {
            plt = plt +
              scale_fill_gradientn(colors=c("#440000", '#ff5500', "white"), values = c(0, .4, 1),
                limits = c(0,1.1), name = 'RV')
            }
        } else { # bliss/hsa excess matrix
            plt = plt +
              scale_fill_gradientn(colors = c('black', '#ffffaa', 'white', 'white', '#aaffff', 'blue'),
                  values = c(0, .35, .48, .51, .65, 1), limits = c(-.6,.6), 
                  name = gsub('bliss', 'Bliss', gsub('hsa', 'HSA', gsub('_', ' ', i))),
                  oob = scales::squish) +
              annotate('text', x = min(df_$pos_x), y = min(df_$pos_y), hjust = 0,
                  label = sprintf('%s score = %.2f',
                  gsub('bliss', 'Bliss', gsub('hsa', 'HSA', gsub('_', ' ', i))), 
                  ifelse(i == 'hsa_excess',
                    agg_results$hsa_q10[idc, iCL], agg_results$bliss_q10[idc, iCL] )))
        }
        if ('0.5' %in% names(all_iso)) { # points of the isobologram at GR/IC50
          plt = plt + geom_point(data = as.data.frame(
              all_iso[[max(1,which(names(all_iso) == '0.5'))]]$df_iso),
            aes(shape = fit_type), show.legend = FALSE)
        }
        if (length(all_iso)>0) { # three isobolograms as lines
          select_iso = names(all_iso)[sort(unique(c(max(1,which(names(all_iso) == '0.5')),
                          which(names(all_iso) %in% ifelse(norm_method[c(1,1,1)] == 'GRvalue',
                              c(0.5, 0.25, 0),  c(0.75, 0.5, 0.25) )))), T)]
          for (iso in select_iso) {
            plt = plt +
              geom_path(mapping = aes(color = color),
                data = cbind(all_iso[[iso]]$df_iso_curve, color = iso), size = 1) +
              geom_path(mapping = aes(x = pos_x_ref, y = pos_y_ref, color = color),
                data = cbind(all_iso[[iso]]$df_iso_curve, color = iso), linetype=2, size=.5)
          }

          plt = plt +
            scale_color_manual(values = iso_colors[select_iso],
              breaks = names(iso_colors[select_iso]),
              labels = paste0(ifelse(norm_method == 'GRvalue', 'GR', 'IC'),
                  100 -100*as.numeric(select_iso)),
              name = 'Isobol')
          }
        return(plt)
      })

      # isobolograms across range of concentration ratios
      plt_iso_compare = ggplot(mapping = aes(x = log10_ratio_conc, y = log2_CI)) +
        geom_line(data = data.frame(log10_ratio_conc = c(-2,2), log2_CI = c(0,0))) +
        geom_hline(yintercept = 0)
      for (iso in names(all_iso)) {
        plt_iso_compare = plt_iso_compare +
            geom_path(mapping = aes(color = color),
              data = cbind(all_iso[[iso]]$df_iso_curve, color = as.numeric(iso)))
      }
      plt_iso_compare = plt_iso_compare  +
        scale_y_continuous(breaks = -5:4, labels = c(paste0('1/', 2**(5:1)), 2**(0:4))) +
        scale_x_continuous(breaks = -3:3, labels = c(paste0('1/', 10**(3:1)), 10**(0:3))) +
        coord_cartesian(ylim = c(-5, 4)) +
        ylab('CI') +
        xlab(paste(condition['DrugName_2'], '/', condition['DrugName'], 'ratio')) +
        theme_bw()
      if (length(all_iso)>1) {
        plt_iso_compare = plt_iso_compare +
          scale_color_gradientn(colors = iso_colors,
            values = scales::rescale(as.numeric(names(iso_colors))),
            breaks = quantile(as.numeric(names(iso_colors)),seq(0,1,.25)),
            labels = 100-round(100*(quantile(as.numeric(names(iso_colors)),seq(0,1,.25)))),
            name = ifelse(norm_method == 'GRvalue', 'GR', 'IC'))
      }

      plt_iso_max = ggplot(df_CI_100x, aes(x=level, y = log2_CI)) +
        geom_hline(yintercept = 0) +
        geom_vline(xintercept = .5) +
        scale_y_continuous(breaks = -5:4, labels = c(paste0('1/', 2**(5:1)), 2**(0:4))) +
        scale_x_continuous(breaks = seq(-1,1,.2), labels = 100-seq(-100,100,20),
            expand = expand_scale(mult = c(.03, 0.13))) +
        coord_cartesian(ylim = c(-5, 5)) +
        ylab('CI') +
        xlab(paste(ifelse(norm_method == 'GRvalue', 'GR', 'IC'), 'value')) +
        theme_bw()

      if (length(all_iso)>0) {
        plt_iso_max = plt_iso_max +
          geom_line(size = 1) +
          annotate('label', x=.5,
              y=df_CI_100x$log2_CI[df_CI_100x$level == .5]+.3,
              hjust = .5, vjust = 0, fill = 'white',
              label = sprintf('CI @%s50 = %.2g', ifelse(norm_method == 'GRvalue', 'GR', 'IC'),
                      2**df_CI_100x$log2_CI[df_CI_100x$level == .5]))
        }

      grid.arrange(grobs = c(plts, list(plt_iso_compare, plt_iso_max)), ncol = 3)


    }
    print('done')
    dev.off()
  }

  # plot the summary of the SE 
  pdf(paste0(folder, norm_method, '_summary_results.pdf'), 7+ncol(agg_results[[1]])/8, 8)

  for (i in 1:4) {
    if (i == 1){
      mx = agg_results$hsa_q10
      mx_title = 'HSA excess'
    } else if (i == 2) {
      mx = agg_results$bliss_q10
      mx_title = 'Bliss excess'
    } else if (i == 3) {
      mx = agg_results$CI_100x_50
      mx_title = 'CI @50'
    } else if (i == 4) {
      mx = agg_results$CI_100x_80
      mx_title = 'CI @80'
    }
    # temp to get proper cell line name
    # mx_CLs = gCLs[match (colData(SE)$clid, gCLs$clid),]
    # stopifnot(all(mx_CLs$canonicalname == colnames(mx)))
    # colnames(mx) = mx_CLs$celllinename
    # # end temp

    if (mx_title  == 'HSA excess') {
      col_fun = colorRamp2(  (c(0, .35, .48, .51, .65, 1)-.5) ,
      c('black', '#ffffaa', 'white', 'white', '#aaffff', 'blue'))
      color_legend = list(
        title = mx_title
      )
    } else if (mx_title  == 'Bliss excess') {
      col_fun = colorRamp2(  .6*(c(0, .35, .48, .51, .65, 1)-.5) ,
      c('black', '#ffffaa', 'white', 'white', '#aaffff', 'blue'))
      color_legend = list(
        title = mx_title
      )
    } else {
      mx = log2(mx)
      col_fun = colorRamp2(  8*(c(0, .35, .48, .51, .65, 1)-.5),
        c('black', '#ffffaa', 'white', 'white', '#aaffff', 'blue'))
      color_legend = list(
          title = mx_title,
          at = -4:4, labels = c(paste0('1/', 2**(4:1)), 2**(0:4))
      )
    }

    mx = mx[order(rowSums(!is.na(mx)) + rowSums(mx,na.rm=T)/1e3),
                order(colSums(!is.na(mx)) + colSums(mx,na.rm=T)/1e3)]
    

    print(Heatmap(mx,
      cluster_rows = ( (!is.null(dim(mx)) && ncol(mx)>1) || length(mx)>1) && all(rowSums(is.na(cor(t(mx), use='pair')))==0),
      clustering_distance_rows = 'spearman',
      cluster_columns =  (!is.null(dim(mx)) && nrow(mx)>1) && all(rowSums(is.na(cor(mx, use='pair')))==0),
      clustering_distance_columns = 'spearman',
      col = col_fun,
      heatmap_legend_param = color_legend,
      column_title = mx_title,
      column_title_gp = gpar(fontsize = 16, fontface = "bold")
    ))

  }
  # Scatter plot of metrics
  metrics = list( c('hsa_q10', 'HSA excess'), c('bliss_q10', 'Bliss excess'),
        c('CI_100x_50', 'CI @50'), c('CI_100x_80', 'CI @80'))
  plts = list()
  for (i_pair in list( c(1,2), c(1,3), c(2,3), c(1,4))) {
    mx_x = agg_results[[ metrics[[i_pair[1]]][1] ]]
    mx_y = agg_results[[ metrics[[i_pair[2]]][1] ]]
    if (grepl('CI', metrics[[i_pair[2]]][1])) { mx_y = log2(mx_y) }
    df_mx = merge(melt(mx_x), melt(mx_y), by = 1:2)
    colnames(df_mx) = c('DrugPair', 'CellLineName', metrics[[i_pair[1]]][1], metrics[[i_pair[2]]][1])

    plt = ggplot(df_mx, 
      aes_string(x = metrics[[i_pair[1]]][1], y = metrics[[i_pair[2]]][1], color = 'DrugPair')) +
      geom_point(size = 2, alpha = .7) +
      xlab(metrics[[i_pair[1]]][2]) + ylab(metrics[[i_pair[2]]][2])

    if (grepl('CI', metrics[[i_pair[2]]][1])) {
      plt = plt + scale_y_continuous(breaks = -5:4, labels = c(paste0('1/', 2**(5:1)), 2**(0:4)))
    }
    plts = c(plts, list( plt ))
  }
  grid.arrange(grobs = plts, ncol = 2)
