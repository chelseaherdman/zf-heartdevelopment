# Functions to create ggplot style principal components analysis plots.

library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)

#-------------------------------------------------------------------------------

screePlot = function(svec, label_prefix="PC",
                     plot_title="Scree plot. Variance explained by each principal component")
{
  # svec can be values from svd()$d, or prcomp()$sdev
  var_explained = (svec^2) / sum(svec^2)
  component_labels = paste(label_prefix, 1:length(svec), sep="")
  tab = data.table(principal_component=factor(component_labels,
                                              levels=component_labels),
                   var_explained=var_explained)
  p = ggplot(tab, aes(x=principal_component, y=var_explained, group=1)) +
    geom_line(alpha=0.3, size=1.3) +
    geom_point(size=3) +
    labs(title=plot_title) +
    ylim(0, 1.0) +
    ylab("fraction variance explained") +
    theme(axis.title.x=element_blank()) +
    theme(axis.text.x=element_text(angle=90))
  
  return(p)
}


#-------------------------------------------------------------------------------

svdMatPlot = function(wide_tab,
                      use_columns=NULL, sort_by_cols=NA_character_,
                      sample_id_variable="NULL",
                      aes_fill="NULL",
                      point_size=5,
                      draw_the_plot=TRUE,
                      fill_palette="Set3",
                      outline_colour="grey30",
                      plot_title="Pairwise scatter plots of principal components.",
                      facet=TRUE)
{
  
  stopifnot(is.data.table(wide_tab))
  
  # Sorting reorders samples along the x-axis.
  if (is.character(sort_by_cols))
  {
    stopifnot(all(sort_by_cols %in% names(wide_tab)))
    setorderv(wide_tab, cols=sort_by_cols, order=-1L)
    new_id_variable_order = as.character(wide_tab[[sample_id_variable]])
    set(wide_tab,
        j=sample_id_variable,
        value=factor(as.character(wide_tab[[sample_id_variable]]),
                     levels=new_id_variable_order))
  }
  
  y_values_name = "value"
  y_variable_name = "principal_component"
  
  mtab = melt(wide_tab,
              measure.vars=use_columns,
              value.name=y_values_name,
              variable.name=y_variable_name)
  
  p = ggplot(mtab, aes_string(x=sample_id_variable, y=y_values_name)) +
    geom_line(aes_string(group=y_variable_name),
              alpha=0.3, size=1.8) +
    geom_point(aes_string(fill=aes_fill), shape=21,
               size=point_size, colour=outline_colour) +
    scale_fill_brewer(palette=fill_palette) +
    theme(axis.title=element_text(size=8)) +
    theme(axis.text=element_text(size=8)) +
    theme(axis.text.x=element_text(angle=90)) +
    labs(title=plot_title)
  
  if(facet)
  {
    facet_formula = as.formula(paste(y_variable_name, "~ ."))
    p = p + facet_grid(facet_formula, scales="free_y")
  }
  
  return(p)
}

#-------------------------------------------------------------------------------


svdDotPlot = function(wide_tab,
                      use_columns=NULL,
                      aes_fill="NULL", aes_colour="NULL",
                      aes_shape="NULL", aes_label="''",
                      point_size=5, label_size=2, label_alpha=0.5,
                      fill_palette="Set3",
                      outline_colour="grey30",
                      plot_title="Dotplots of principal component values.",
                      facet=TRUE)
{
  
  stopifnot(is.data.table(wide_tab))
  
  y_values_name = "value"
  y_variable_name = "principal_component"
  
  mtab = melt(wide_tab,
              measure.vars=use_columns,
              value.name=y_values_name,
              variable.name=y_variable_name)
  
  p = ggplot(mtab, aes_string(x=y_values_name)) +
    geom_point(y=0, aes_string(fill=aes_fill), shape=21,
               size=point_size, colour=outline_colour) +
    geom_text(y=0, aes_string(label=aes_label),
              size=label_size, alpha=label_alpha) +
    scale_fill_brewer(palette=fill_palette) +
    theme(axis.title=element_text(size=8)) +
    theme(axis.text=element_text(size=8)) +
    theme(axis.text.x=element_text(angle=90)) +
    labs(title=plot_title)
  
  if(facet)
  {
    facet_formula = as.formula(paste("~", y_variable_name))
    p = p +
      facet_wrap(facet_formula, scales="free_x", ncol=1, switch="x") +
      theme(strip.background=element_blank()) +
      theme(axis.title.x=element_blank())
  }
  
  return(p)
}

#-------------------------------------------------------------------------------


svdDensityPlot = function(wide_tab,
                          use_columns=NULL,
                          aes_fill="NULL", aes_colour="NULL",
                          aes_shape="NULL", aes_label="''",
                          point_size=5, label_size=2, label_alpha=0.5,
                          fill_palette="Set3",
                          outline_colour="grey30",
                          plot_title="Dotplots of principal component values.",
                          facet=TRUE)
{
  
  stopifnot(is.data.table(wide_tab))
  
  y_values_name = "value"
  y_variable_name = "principal_component"
  
  mtab = melt(wide_tab,
              measure.vars=use_columns,
              value.name=y_values_name,
              variable.name=y_variable_name)
  
  p = ggplot(mtab, aes_string(x=y_values_name)) +
    geom_density(aes_string(fill=aes_fill),
                 colour=outline_colour) +
    scale_fill_brewer(palette=fill_palette) +
    theme(axis.title=element_text(size=8)) +
    theme(axis.text=element_text(size=8)) +
    theme(axis.text.x=element_text(angle=90)) +
    labs(title=plot_title)
  
  if(facet)
  {
    facet_formula = as.formula(paste("~", y_variable_name))
    p = p +
      facet_wrap(facet_formula, scales="free_x", ncol=1, switch="x") +
      theme(strip.background=element_blank()) +
      theme(axis.title.x=element_blank())
  }
  
  return(p)
}

#-------------------------------------------------------------------------------

makeLayoutMatrix = function(n) {
  t(matrix(seq(n^2), nrow=n, byrow=TRUE))
}

#-------------------------------------------------------------------------------

svdPairsPlot = function(wide_tab, use_columns=NULL,
                        aes_fill="NULL", aes_colour="NULL",
                        aes_shape="NULL", aes_label="''",
                        point_size=5, label_size=2,
                        draw_the_plot=TRUE,
                        fill_palette="Set3",
                        outline_colour="grey30",
                        label_alpha=0.5,
                        point_alpha=1,
                        fill_values="NULL",
                        plot_title="Pairwise scatter plots of principal components.")
{
  stopifnot(is.data.table(wide_tab))
  
  pair_mat = combn(use_columns, 2)
  plist = vector("list", length=ncol(pair_mat))
  
  
  for (i in seq(ncol(pair_mat))) {
    xvar = pair_mat[1, i]
    yvar = pair_mat[2, i]
    
    
    if ((aes_shape != "NULL") & (aes_fill != "NULL")) {
      tmp = ggplot(wide_tab, aes_string(x=xvar, y=yvar)) +
        theme_bw() +
        geom_point(aes_string(fill=aes_fill, shape=aes_shape),
                   size=point_size, colour=outline_colour, alpha=point_alpha) +
        geom_text(aes_string(label=aes_label),
                  size=label_size, alpha=label_alpha) +
        scale_shape_manual(values=c(21, 24, 23, 22, 25)) +
        scale_fill_manual(values=fill_values) +
        scale_x_continuous(expand=c(0.1, 0)) +
        scale_y_continuous(expand=c(0.1, 0)) +
        theme(axis.title=element_text(size=8)) +
        theme(axis.text=element_text(size=8)) +
        theme(legend.position="none")
    } else {
      tmp = ggplot(wide_tab, aes_string(x=xvar, y=yvar)) +
        theme_bw() +
        geom_point(aes_string(fill=aes_fill), alpha=point_alpha,
                   shape=21, size=point_size, colour=outline_colour) +
        geom_text(aes_string(label=aes_label),
                  size=label_size, alpha=label_alpha) +
        scale_fill_manual(values=fill_values) +
        scale_x_continuous(expand=c(0.1, 0)) +
        scale_y_continuous(expand=c(0.1, 0)) +
        theme(axis.title=element_text(size=8)) +
        theme(axis.text=element_text(size=8)) +
        theme(legend.position="none")
    }
    
    plist[[i]] = tmp
  }
  
  n_col = length(use_columns) - 1
  newlist = vector("list", length=n_col^2)
  blank_panel = rectGrob(gp=gpar(col=NA))
  panel_type = c(lower.tri(matrix(nrow=n_col, ncol=n_col), diag=TRUE))
  plot_layout = makeLayoutMatrix(n_col)
  ## Need to fix panel order!!
  j = 1
  for (i in seq(n_col^2)) {
    if(panel_type[i]) {
      newlist[[i]] = plist[[j]]
      j = j + 1
    } else {
      newlist[[i]] = blank_panel
    }
  }
  
  newlist = c(newlist, list(ncol=n_col, top=plot_title, layout_matrix=plot_layout))
  p = do.call(arrangeGrob, newlist)
  #if (draw_the_plot) {grid.newpage(); grid.draw(p)}
  return(p)
}


#===============================================================================
svdPairsGradient = function(wide_tab, use_columns=NULL,
                            aes_fill="NULL", aes_colour="NULL",
                            aes_shape="NULL", aes_label="''",
                            point_size=5, label_size=2,
                            draw_the_plot=TRUE,
                            fill_palette="Set3",
                            outline_colour="grey30",
                            label_alpha=0.5,
                            point_alpha=1,
                            fill_values="NULL",
                            plot_title="Pairwise scatter plots of principal components.")
{
  stopifnot(is.data.table(wide_tab))
  
  pair_mat = combn(use_columns, 2)
  plist = vector("list", length=ncol(pair_mat))
  
  
  for (i in seq(ncol(pair_mat))) {
    xvar = pair_mat[1, i]
    yvar = pair_mat[2, i]
    
    
    if ((aes_shape != "NULL") & (aes_fill != "NULL")) {
      tmp = ggplot(wide_tab, aes_string(x=xvar, y=yvar)) +
        theme_bw() +
        geom_point(aes_string(fill=aes_fill, shape=aes_shape),
                   size=point_size, colour=outline_colour, alpha=0.5) +
        geom_text(aes_string(label=aes_label),
                  size=label_size, alpha=label_alpha) +
        scale_shape_manual(values=c(21, 24, 23, 22, 25)) +
        scale_fill_distiller(palette="YlGnBu") +
        scale_x_continuous(expand=c(0.1, 0)) +
        scale_y_continuous(expand=c(0.1, 0)) +
        theme(axis.title=element_text(size=8)) +
        theme(axis.text=element_text(size=8)) +
        theme(legend.position="none")
    } else {
      tmp = ggplot(wide_tab, aes_string(x=xvar, y=yvar)) +
        theme_bw() +
        geom_point(aes_string(fill=aes_fill), alpha=0.5,
                   shape=21, size=point_size, colour=outline_colour) +
        geom_text(aes_string(label=aes_label),
                  size=label_size, alpha=label_alpha) +
        scale_fill_distiller(palette="YlGnBu") +
        scale_x_continuous(expand=c(0.1, 0)) +
        scale_y_continuous(expand=c(0.1, 0)) +
        theme(axis.title=element_text(size=8)) +
        theme(axis.text=element_text(size=8)) +
        theme(legend.position="none")
    }
    
    plist[[i]] = tmp
  }
  
  n_col = length(use_columns) - 1
  newlist = vector("list", length=n_col^2)
  blank_panel = rectGrob(gp=gpar(col=NA))
  panel_type = c(lower.tri(matrix(nrow=n_col, ncol=n_col), diag=TRUE))
  plot_layout = makeLayoutMatrix(n_col)
  ## Need to fix panel order!!
  j = 1
  for (i in seq(n_col^2)) {
    if(panel_type[i]) {
      newlist[[i]] = plist[[j]]
      j = j + 1
    } else {
      newlist[[i]] = blank_panel
    }
  }
  
  newlist = c(newlist, list(ncol=n_col, top=plot_title, layout_matrix=plot_layout))
  p = do.call(arrangeGrob, newlist)
  #if (draw_the_plot) {grid.newpage(); grid.draw(p)}
  return(p)
}




