
gDR_matrix_plot = function(SE, normalization = 'GR') {
  
sel_SE = SE[1,1]

drug1_axis =  assay(sel_SE, 'isobolograms')

plts = lapply(c("SmoothMatrix", "BlissExcess", 'HSAExcess'), function(i_assay) {
  df_ = assay(sel_SE, i_assay)
  df_$value = pmin(1.1,df_$value)

  #gDRcore::define_matrix_grid_positions(df_$Concentration, df_$Concentration_2) 
  drug1_axis = data.frame(Concentration = df_$Concentration,
                  pos_y = log10(df_$Concentration)
  )
  drug2_axis = data.frame(Concentration_2 = df_$Concentration_2,
                  pos_x = log10(df_$Concentration_2)
  )
        df_ = merge(df_, drug1_axis, by = 'Concentration', all.x = F)
        df_ = merge(df_, drug2_axis, by = 'Concentration_2', all.x = F)
        
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
