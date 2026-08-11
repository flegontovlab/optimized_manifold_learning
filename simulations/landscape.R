library(tidyverse)
library(magrittr)
library(pals)
source("helper_functions_hexagons.R")

######## LANDSCAPE ########
## create hexagonal world (approx circle)
make_landscape <- function(radius = 10) {
    span = -radius:radius

    landscape <- expand.grid(span, span, span) %>%
        as_tibble %>%
        magrittr::set_colnames(c("q","r","s"))

    landscape <- landscape[which(rowSums(landscape) == 0),] # only valid hex coordinates

    landscape <- landscape %>% mutate(
        FID = paste0("deme", sprintf("%03d", 1:nrow(landscape))),
        x = q,
        y = (r-s) * sin(deg2rad(60)) * (2/3),
        colour = color_points_in_circle(x, y, num_sectors, rotation)
    ) %>%
    relocate(FID) %>%
    rowwise %>%
    mutate(
        size = radius - max(abs(c(q,r,s))) + 1
    ) %>%
    ungroup
}

num_sectors = 18
rotation = 25

landscape <- make_landscape(10)
write_csv(landscape, "landscape.csv")

## plot
p <- landscape %>%
    ggplot(aes(x,y)) + geom_point(aes(colour = colour, size = size)) +
    scale_color_identity() + theme_bw() + theme(legend.position = "none")

ggsave("landscape.pdf", p, width = 8, height = 8)

######## NEIGHBOURS ########
find_nbs <- function(world, deme) {
    ## merge coords as chars for matching
    world <- world %>% mutate(coords = paste(q,r,s))
    ## filter single deme by FID
    x <- filter(world, FID == deme)
    ## calculate neighbours
    nbs <- c(
        paste(x$q+1, x$r, x$s-1), paste(x$q+1, x$r-1, x$s),
        paste(x$q-1, x$r, x$s+1), paste(x$q-1, x$r+1, x$s),
        paste(x$q, x$r-1, x$s+1), paste(x$q, x$r+1, x$s-1)
    )
    ## match neighbours
    demes <- world %>% filter(coords %in% nbs) %>% pull(FID)
    return(demes)
}

## make a list of neighbour demes
nbs <- list()
for (i in landscape$FID) {
    nbs[[i]] <- paste(i, find_nbs(landscape, i))
}

## convert to a data frame
neighbour_pairs <- nbs %>%
    unlist(use.names = F) %>%
    as_tibble %>%
    separate_wider_delim(cols = value, delim = " ", names = c("pop1", "pop2"))

write_tsv(neighbour_pairs, "neighbour_pairs.tsv")
