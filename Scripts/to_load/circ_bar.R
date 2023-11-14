# Create a circular bar

circ_bar <- function(
    df, 
    x, 
    y,
    id,
    taxa,
    labels = NULL, 
    clr_fill = NA, 
    breaks = 1,
    lgd_rows = 5, 
    lgd_pos = "bottom",
    path = here::here(),
    title = "circ_bar",
    dpi = 1200,
    w = 40,
    h = 40, 
    save = FALSE
) {
  
  x_sym <- dplyr::sym(x)
  y_sym <- dplyr::sym(y)
  id_sym <- dplyr::sym(id)
  taxa_sym <- dplyr::sym(taxa)
  
  
  ## Set gaps between groups ----
  empty_bar <- breaks
  nbreaks <- nlevels(df[[x]]) %>% 
    as.numeric()
  
  nObs <- df %>% 
    dplyr::group_by(Sample) %>% 
    dplyr::tally() %>% 
    dplyr::slice(1) %>% 
    .$n %>% 
    as.numeric()
  
  
  to_add <- data.frame(
    matrix(NA, empty_bar * nbreaks * nObs, ncol(df))) %>% 
    purrr::set_names(colnames(df)) %>% 
    dplyr::mutate(!!x_sym := rep(levels(df[[x]]),
                                 each = empty_bar * nObs))
  
  
  ## Define order of bars ----
  p2 <- df %>%
    dplyr::bind_rows(to_add) %>%
    dplyr::arrange(
      desc(!!x_sym),
      !!id_sym
    ) %>%
    dplyr::mutate(labels = rep(seq(1, nlevels(!!id_sym) + empty_bar * nbreaks),
                               each = nObs
    )) %>%
    dplyr::mutate(!!x_sym := as.character(!!x_sym)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(!!id_sym := as.character(!!id_sym)) %>%
    dplyr::mutate(labels = factor(labels,
                                  levels = unique(labels)
    )) %>%
    dplyr::mutate(!!y_sym := ifelse(is.na(!!y_sym),
                                    0,
                                    !!y_sym
    ))
  
  
  ## Define label positions around polar plot ----
  label_data <- p2 %>%
    dplyr::group_by(labels, !!id_sym) %>%
    dplyr::summarise(tot = sum(!!y_sym)) %>%
    dplyr::mutate(labels = as.numeric(as.character(labels)))
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$labels - 0.5) / number_of_bar
  label_data$hjust <- ifelse(angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle + 180, angle)
  
  ## prepare a data frame for base lines ----
  base_data <- p2 %>%
    dplyr::group_by(!!x_sym) %>%
    dplyr::mutate(labels = as.numeric(as.character(labels))) %>%
    dplyr::summarise(
      start = min(labels),
      end = max(labels) - empty_bar
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(title = mean(c(start, end)))
  
  ## prepare a data frame for grid (scales) ----
  grid_data <- base_data
  grid_data$end <- grid_data$end[c(
    nrow(grid_data),
    1:nrow(grid_data) - 1
  )] + 1
  grid_data$start <- grid_data$start - 1
  grid_data <- grid_data[-1, ]
  
  
  ## plot parameters ----
  plot_theme(text_size = 9)
  
  if (is.null(labels)) {
    lables <- names(clrs)
  }
  
  ## polar plot ----
  polar <- ggplot2::ggplot(p2) +
    
    ## control variables for polar plot
    ggplot2::coord_polar() +
    ggplot2::ylim(-0.8, 1.1) +
    
    
    ## add stacked bar for taxa
    ggplot2::geom_bar(aes(
      x = labels,
      y = !!y_sym,
      fill = !!taxa_sym
    ),
    stat = "identity",
    width = 1,
    position = position_stack(reverse = T)
    ) +
    ggplot2::scale_fill_manual(
      values = clr_fill,
      name = "Species",
      label = labels,
      guide = guide_legend(
        nrow = lgd_rows,
      )
    ) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_tile(
      mapping = aes(
        x = labels,
        y = -0.25,
        fill = ReadCounts
      ),
      stat = "identity", 
      height = 0.1, 
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradientn(
      colours = grDevices::hcl.colors(palette = "TealRose", 100),
      na.value = "transparent",
      name = "Read\ncounts"
    ) +
    
    ## Add text for incremental y-axis tick labels
    ggplot2::annotate("text",
                      x = max(as.numeric(as.character(p2$labels))),
                      y = c(0, .25, .50, .75, 1.0),
                      label = c("0", "25", "50", "75", "100"),
                      color = "grey60",
                      size = 5,
                      angle = 0,
                      fontface = "bold",
                      hjust = 1
    ) +
    
    ## Add samples names on top of bars
    ggplot2::geom_text(
      data = label_data,
      aes(
        x = labels,
        y = 1.05,
        label = !!id_sym,
        hjust = hjust
      ),
      color = "grey60",
      size = 4,
      angle = label_data$angle,
      inherit.aes = FALSE
    ) +
    
    ## start new colour scale
    ggnewscale::new_scale_colour() +
    
    ## Add segment at base of the plot
    ggplot2::geom_segment(
      data = base_data,
      aes(
        x = start,
        y = -0.03,
        xend = end,
        yend = -0.03,
        colour = !!x_sym
      ),
      alpha = 0.8,
      linewidth = 0.8,
      inherit.aes = FALSE,
      show.legend = F
    ) +
    
    ## label segments
    geomtextpath::geom_textpath(
      data = base_data,
      aes(
        x = title,
        y = -0.2,
        label = !!x_sym,
        colour = !!x_sym
      ),
      position = position_stack(vjust = 0.5),
      alpha = 0.8,
      size = 4.8,
      fontface = "bold",
      inherit.aes = FALSE,
      show.legend = F
    ) +
    
    ## Add colour
    ggplot2::scale_colour_manual(values = pal2) +
    
    ## adjust labels and theme
    ggplot2::theme(
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(c(-1, -1, 0, -1), "cm"),
      legend.position = lgd_pos,
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.5, "cm")
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab("")
  
  
  # extract read counts legend
  rd_lgd <- cowplot::get_legend(
    ggplot2::ggplot(p2) +
      
      ggplot2::geom_tile(
        mapping = aes(
          x = labels,
          y = -0.25,
          fill = ReadCounts
        ),
        stat = "identity", 
        height = 0.1
      ) +
      ggplot2::scale_fill_gradientn(
        colours = grDevices::hcl.colors(palette = "TealRose", n = 100),
        na.value = "transparent",
        name = "Read\ncounts"
      ) 
  )
  
  polar <- polar + patchwork::inset_element(
    rd_lgd, 
    left = 0.5, 
    right = 0.5, 
    bottom = 0.5, 
    top = 0.5)
  
  
  if (save) {
    ## Save as pdf
    ggplot2::ggsave(
      polar,
      path = path,
      file = paste0(title, ".pdf"),
      device = "pdf",
      dpi = dpi,
      width = w,
      height = h,
      units = "cm"
    )
    
    ggplot2::ggsave(
      polar,
      path = path,
      file = paste0(title, ".png"),
      device = "png",
      dpi = dpi,
      width = w,
      height = h,
      units = "cm"
    )
  }
  
  polar
}
