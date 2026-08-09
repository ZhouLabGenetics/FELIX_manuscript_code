# 00_people.R — recolour the BioArt unisex icon by ancestry (solid or striped = admixed).
# source() after 00_theme.R. Returns ggplot annotation layers placing a person at (x, y).
suppressMessages({library(png); library(grid)})
.ICON <- png::readPNG(file.path(ROOT, "figures_R", "png", "BIOART-13", "UnisexIcon0001.png"))
.ASPECT <- dim(.ICON)[1] / dim(.ICON)[2]           # height / width (~2.39)
.hex2 <- function(h) c(strtoi(substr(h,2,3),16), strtoi(substr(h,4,5),16), strtoi(substr(h,6,7),16)) / 255

person_raster <- function(hexes) {
  img <- .ICON; m <- img[,,4] > 0.1                # opaque silhouette mask
  if (length(hexes) == 1) {
    for (k in 1:3) { ch <- img[,,k]; ch[m] <- .hex2(hexes)[k]; img[,,k] <- ch }
  } else {
    xb <- cut(col(m), breaks = length(hexes), labels = FALSE)   # vertical stripes = admixture
    for (k in 1:3) { ch <- img[,,k]
      for (s in seq_along(hexes)) { sel <- m & xb == s; ch[sel] <- .hex2(hexes[s])[k] }
      img[,,k] <- ch }
  }
  grid::rasterGrob(img, interpolate = TRUE)
}
# place a person of pixel-height h (data units) centred at (x,y); width = h/aspect
person_layer <- function(x, y, h, hexes) {
  w <- h / .ASPECT
  annotation_custom(person_raster(hexes), xmin = x - w/2, xmax = x + w/2,
                    ymin = y - h/2, ymax = y + h/2)
}
