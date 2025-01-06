library(tidyverse)
library(ggplot2)
library(reshape2)
library(stringr)
library(ggforce)
library(patchwork)


def_theme <- theme_classic() +
    theme(
        axis.text.x = element_text(size = 18, angle = 0, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.title = element_text(size = 26),
    )

rowidlabeller <- function(variable) {
    splist <- strsplit(variable, split = "#")[[1]]
    rswtch <- splist[1]

    magname <- strsplit(splist[2], split = "_")[[1]]
    magname2 <- paste(magname[1:length(magname) - 1], collapse = "_")
    if (magname2 != "") {
        tt <- sprintf("%s|%s|%s%s", rswtch, magname2, splist[3], splist[length(splist)])
    } else {
        tt <- sprintf("%s|%s%s", rswtch, splist[3], splist[length(splist)])
    }
    return(tt)
}

peakplotmaker <- function(plog_path, cstats_path, out_path, name) {
    df <- read.csv(plog_path)
    statdf <- read.csv(cstats_path)

    df <- df[df$decision != "", ]

    df$rowid <- sapply(df$rowid, rowidlabeller)
    statdf$rowid <- sapply(statdf$rowid, rowidlabeller)

    statdf$ann <- statdf$delta_mean < 0 &
        statdf$delta_mean_pval < 0.05 &
        statdf$delta_variance_pval < 0.05 &
        statdf$delta_variance > statdf$noiseset_delta_variance
    statdf <- statdf[statdf$ann, ]
    if (nrow(statdf) > 0) statdf$ann_text <- "*"



    peakplot <- ggplot(df, aes(x = from_riboswitch_end_relative, y = coverage_delta_stable_relative)) +
        geom_vline(xintercept = 0, colour = "#959595", linetype = "dashed") +
        geom_line(aes(color = transcriptome)) +
        geom_point() +
        facet_wrap(~rowid, scales = "fixed") +
        geom_mark_ellipse(aes(fill = as.factor(cluster)), expand = unit(0.5, "mm")) +
        labs(
            title = str_glue("{name} riboswitches across all available conditions"),
            y = "Fractional coverage change",
            x = "Position relative to riboswitch 3′ end as a fraction of riboswitch length"
        ) +
        def_theme +
        theme(
            plot.title = element_text(size = 24),
            legend.position = "none",
            strip.text = element_text(size = 20),
            strip.background = element_rect(fill = "#cdcdcd", linewidth = 0),
        ) +
        coord_cartesian(ylim = c(-1.2, 0.5)) +
        guides(alpha = "none", color = "none")

    if (nrow(statdf) > 0) {
        peakplot <- peakplot + geom_text(data = statdf, mapping = aes(x = pos_mean, y = 0.3, label = ann_text), size = 12, fontface = "bold")
    }

    # peakplot

    print(str_glue("Saving plot: {name}"))
    alph <- 0.7
    ggsave(out_path, plot = peakplot, dpi = 320 * alph, units = "px", width = 7000 * alph, height = 4000 * alph)
}
