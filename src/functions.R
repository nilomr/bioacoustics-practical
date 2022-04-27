# ------------------------------------------------------------------------------
# Source code for bioacoustics practical
# (c) Nilo Merino Recalde 2021-present
# University of Oxford
# ------------------------------------------------------------------------------

# DEPENDENCIES

# Load or install dependencies
dependencies <- c(
  "rprojroot",
  "warbleR",
  "dplyr",
  "tibble",
  "magrittr",
  "stringr",
  "ggrepel",
  "patchwork",
  "filesstrings",
  "pbapply",
  "randomForest",
  "yardstick",
  "glue"
)
packages <- lapply(dependencies, function(y) {
  if (!y %in% installed.packages()[, "Package"]) {
    install.packages(y)
  }
  try(require(y, character.only = T), silent = T)
})


# FUNCTION DEFINITIONS

# ---- Data ingest -------------------------------------------------------------

#' Import frequency annotations from .csv files output in Sonic Visualiser
#'
#' @param dir Path to a directory containing the csv files,
#' which must be named `1.csv`, `2.csv`, ... `15.csv`.
#' @param tutor Who is calling the function? Defaults to FALSE.
#'
#' @return a tibble with columns `mean_frequency` <dbl> and
#'  `file` <int>.
#' @export
import_freqdata <-
  function(dir, tutor = FALSE) {
    freq_data_files <- if (tutor) {
      list.files(file.path(dir, "tutor"),
                 pattern = "[0-9].csv$",
                 full.names = TRUE)
    } else {
      list.files(dir, pattern = "[0-9].csv$", full.names = TRUE)
    }
    
    if (length(freq_data_files) < 15) {
      print(
        glue(
          '{15 - length(freq_data_files)}/{15} ',
          '"*.csv" files are missing in {freqdir}'
        )
      )
    } else {
      freq_list <- list()
      for (file in freq_data_files) {
        # Get the min/max frequency data
        df <- read.csv(file.path(file),
                       header = FALSE)[3:4]
        # Get mean for each element and then global mean
        global_mean <- mean(rowMeans(df, na.rm = TRUE))
        freq_list[gsub("\\.csv$", "", file)] <- global_mean
      }
      
      # Build a dataframe with the frequency and ID info
      i <-
        length(stringr::str_split(freq_data_files[1], "/", simplify = TRUE))
      freq_df <- stack(freq_list) %>%
        tibble::as_tibble() %>%
        dplyr::rename(file = ind, mean_frequency = values) %>%
        dplyr::mutate(file = as.integer(str_split(file, "/", simplify = TRUE)[, i]))
      return(freq_df)
    }
  }

#' Merges body mass, species ID and frequency dataframes.
#'
#' @param dir Path to data folder in the project. Requires that
#' `dir/frequency-data/tutor` exist if tutor=TRUE
#' @param body_mass_df A dataframe with a `body_mass` <numeric> column.
#' @param ids_df A dataframe with `file` <chr> and `name` <chr> columns.
#' @param freq_df A dataframe with columns `mean_frequency` <dbl> and
#'  `file` <int>.
#' @param tutor Who is calling the function? Defaults to FALSE.
#'
#' @returns A merged dataframe
#' @export
merge_data <-
  function(dir, body_mass_df, ids_df, freq_df, tutor = FALSE) {
    if (tutor) {
      full_df <-
        read.csv(file.path(dir,
                           "frequency-data", "tutor", "bird-ids-private.csv"),
                 header = TRUE) %>%
        mutate(file = as.integer(file)) %>%
        left_join(., freq_df) %>%
        mutate(across(where(is.numeric), ~ round(., 2)))
    } else{
      full_df <- cbind(ids_df, body_mass_df) %>%
        left_join(., freq_df) %>%
        mutate(across(where(is.numeric), ~ round(., 2)))
    }
    return(full_df)
  }


#' Downlaod data from xeno-canto
#'
#' @param id A vector of integers matching xeno-canto song codes.
#' @param dir A path-like object pointing to the desired output directory.
#' @return None
download_xc <- function(id, dir) {
  name <- rjson::fromJSON(
    file = paste0(
      "https://www.xeno-canto.org/",
      "api/2/recordings?query=",
      paste0("nr:", id),
      "&page=1"
    )
  )$recordings[[1]]$en
  file <- file.path(dir, paste0(id, "-", name, ".mp3"))
  if (!file.exists(file)) {
    download.file(
      url = paste("https://xeno-canto.org/", id, "/download", sep = ""),
      destfile = file,
      quiet = TRUE,
      mode = "wb",
      cacheOK = TRUE,
      extra = getOption("download.file.extra")
    )
  } else {
    print("File exists in destination")
  }
}


# ---- Plotting ----------------------------------------------------------------
# TODO: Document these

#' Practical theme
#'
bioac_theme <- function() {
  theme(
    panel.background = element_rect(fill = "#f5f5f5",
                                    colour = "#f5f5f5"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    text = element_text(size = 15),
    plot.title = element_text(vjust = 0.1),
  )
}

#' Each variable separately + KDE
#'
plot_vars <- function(df) {
  bmass_plot <- df %>% ggplot(aes(x = body_mass)) +
    geom_density(fill = "#bfbfbf", colour = "#bfbfbf") +
    geom_rug(length = unit(0.045, "npc")) +
    theme(aspect.ratio = 1) +
    bioac_theme() +
    theme(plot.margin = unit(c(0, 30, 0, 0), "pt")) +
    labs(x = "\nBody mass (g)",
         y = "Density\n")
  
  freq_plot <- df %>% ggplot(aes(x = mean_frequency)) +
    geom_density(fill = "#bfbfbf", colour = "#bfbfbf") +
    geom_rug(length = unit(0.045, "npc")) +
    theme(aspect.ratio = 1) +
    bioac_theme() +
    labs(x = "\nMean frequency (Hz)",
         y = element_blank())
  
  endplot <- bmass_plot + freq_plot + plot_annotation(
    title = "Smoothed histograms (KDEs) of our two variables\n",
    caption = paste0(
      "\nKDE = Kernel Density Estimation, not the right method ",
      "to use here but useful for vis. purposes"
    ),
    theme = theme(
      plot.title = element_text(size = 18, hjust = 0.5),
      plot.caption = element_text(size = 7, hjust = 0.5)
    )
  )
  return(endplot)
}


#' Body mass vs frequency linear fit
#'
plot_mass_vs_freq <- function(df) {
  # Calculate the coefficient of determination (R-squared)
  rsq <- function(x, y) {
    cor(x, y) ^ 2
  }
  rsq_ <- format(round(rsq(df$body_mass, df$mean_frequency), 2),
                 nsmall = 2)
  
  # Plot the data
  outplot <- df %>%
    # Remove underscore and capitalise species names
    mutate(name = str_replace(str_to_sentence(name), "_", " ")) %>%
    ggplot(aes(x = body_mass, y = mean_frequency)) +
    geom_smooth(
      method = "lm",
      colour = "grey",
      fill = "grey",
      alpha = 0.21
    ) +
    geom_point(size = 3) +
    scale_x_log10() +
    geom_text_repel(
      aes(label = name),
      seed = 666,
      min.segment.length = 0,
      box.padding = 0.5,
      color = "black",
      size = 3.3,
      force = 10,
      alpha = 0.6,
      max.time = 1,
      segment.curvature   = 0.2,
      segment.alpha       = 0.6,
      fontface = "italic"
    ) +
    labs(title = "Mean frequency vs body mass\n",
         x = "\nLog body mass (g)",
         y = "Frequency (Hz)\n") +
    annotate(
      "text",
      x = 400,
      y = 7000,
      label = stringr::str_interp("italic(R) ^ 2 == ${rsq_}"),
      parse = TRUE
    ) +
    theme(aspect.ratio = 0.6, plot.title = element_text(hjust = 0.5)) +
    bioac_theme()
  return(outplot)
}

#' Plot PCA
#'
plot_pca <- function(pca, labs) {
  as_tibble(pca[[5]]) %>%
    tibble::add_column(label = labs) %>%
    ggplot(aes(x = PC1, y = PC2, colour = label)) +
    geom_point(alpha = 0.8, shape = 16) +
    scale_colour_manual(
      values = c("#585955", "#f0b618"),
      labels = c("Coal tit", "Great tit")
    ) +
    bioac_theme() +
    theme(legend.position = "right", aspect.ratio = 0.6, ) +
    labs(title = "Principal component analysis",
         subtitle = "based on spectral measurements\n")
}



#' Plot measures of modulation
#'
plot_modulation <- function(spectral_features) {
  linecolour <- adjustcolor("#c4692d", alpha.f = .9)
  dfrange_plot <-
    spectral_features %>% ggplot(aes(x = label, y = dfrange)) +
    geom_jitter(width = 0.2,
                alpha = 0.2,
                shape = 16) +
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.43,
      size = 1,
      colour = linecolour
    ) +
    bioac_theme() +
    labs(x = NULL, y = "Dominant frequency range (kHz)\n") +
    scale_x_discrete(labels = c("Coal tit", "Great tit")) +
    theme(aspect.ratio = 1.3, plot.margin = unit(c(0, 30, 0, 0), "pt"))
  
  kurt_plot <-
    spectral_features %>% ggplot(aes(x = label, y = kurt)) +
    geom_jitter(width = 0.2,
                alpha = 0.2,
                shape = 16) +
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.43,
      size = 1,
      colour = linecolour
    ) +
    bioac_theme() +
    labs(x = NULL, y = "Kurtosis\n") +
    scale_x_discrete(labels = c("Coal tit", "Great tit")) +
    theme(aspect.ratio = 1.3)
  
  dfrange_plot + kurt_plot +
    plot_annotation(
      title = "Measures of frequency modulation\n",
      caption = paste0(
        "\nHigher kurtosis = more leptokurtic, that is, ",
        "having a distribution more centered around the mean ",
        "than a normal distribution.
    Note that these variables correspond to two ways of ",
        "measuring the same thing and are therefore correlated"
      ),
      theme = theme(
        plot.title = element_text(size = 18, hjust = 0.5),
        plot.caption = element_text(size = 9, hjust = 0.5)
      )
    )
}


#' Plot some temporal characteristics
#'
plot_durations <- function(spectral_features) {
  linecolour <- adjustcolor("#c4692d", alpha.f = .9)
  element_list <- element_list %>% mutate(
    label = spectral_features$label,
    duration = end - start,
    silence = c(element_list$start[2:nrow(element_list)], NA) - end
  )
  silences_plot <- element_list %>%
    filter(silence > 0 & silence < 0.2) %>%
    ggplot(aes(x = label, y = silence)) +
    geom_jitter(width = 0.2,
                alpha = 0.2,
                shape = 16) +
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.43,
      size = 1,
      colour = linecolour
    ) +
    bioac_theme() +
    labs(x = NULL, y = "Silence duration (s)\n") +
    scale_x_discrete(labels = c("Coal tit", "Great tit")) +
    theme(aspect.ratio = 1.3, plot.margin = unit(c(0, 30, 0, 0), "pt"))
  
  durations_plot <- element_list %>%
    ggplot(aes(x = label, y = duration)) +
    geom_jitter(width = 0.2,
                alpha = 0.2,
                shape = 16) +
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.43,
      size = 1,
      colour = linecolour
    ) +
    bioac_theme() +
    labs(x = NULL, y = "Note duration (s)\n") +
    scale_x_discrete(labels = c("Coal tit", "Great tit")) +
    theme(aspect.ratio = 1.3)
  
  silences_plot + durations_plot + plot_annotation(
    title = "Silence and note duration\n",
    theme = theme(
      plot.title = element_text(size = 18, hjust = 0.5),
      plot.caption = element_text(size = 9, hjust = 0.5)
    )
  )
}
