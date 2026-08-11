## helper functions for hexagonal landscapes

## angle conversions
## https://stackoverflow.com/q/32370485/5184574
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}

## coloured coordinates
#num_sectors = 18
#rotation = 25
# Function to generate random points in a circle and color them by sector
color_points_in_circle <- function(x, y, num_sectors = 18, rotation = 25) {
  # Calculate angle of each point relative to the center of the circle
  angles <- ((atan2(y, x) + pi) * (180 / pi)) + rotation  # Convert radians to degrees, and adjust range to [0, 360]
  angles <- ifelse(angles > 360, angles - 360, angles)
  # Calculate sector for each point based on angle
  sectors <- as.integer(angles / (360 / num_sectors)) + 1 # as.integer
  # Map sectors to colors using base R rainbow function
  colors <- pals::kovesi.cyclic_mrybm_35_75_c68(num_sectors)
  return(colors[sectors])
}

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
        colour = color_points_in_circle(x, y, num_sectors = 18, rotation = 25)
    ) %>% 
    relocate(FID) %>% 
    rowwise %>% 
    mutate(
        size = radius - max(abs(c(q,r,s))) + 1
    ) %>% 
    ungroup
}

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
