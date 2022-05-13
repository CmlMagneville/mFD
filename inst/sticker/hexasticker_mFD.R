#'
#' Create an Hexagonal Sticker for the mFD Package
#'


## Load the image:
p <- png::readPNG(here::here("inst", "sticker", "mFD_hexSticker.png"))
grid_mFD <- grid::rasterGrob(p, width = 0.72, height = 0.81, x = 0.5, y = 0.58,
                         interpolate = TRUE)

gg <- ggplot2::ggplot() +
  ggplot2::annotation_custom(grid_mFD) +
  ggplot2::theme_void()



## Export Sticker ----

hexSticker::sticker(
  
  subplot   = gg,
  package   = "mFD",
  filename  = here::here("man", "figures", "hexasticker_mFD.png"),
  dpi       = 2400,
  
  p_size    = 180.0,         # Title
  u_size    =  35.0,         # URL
  p_family  = "Aller_Rg",
  u_family = "Aller_Rg",
  
  p_color   = "#ffffff",
  h_fill    = "#000000",   # Background
  h_color   = "#666666ff",   # Border
  u_color   = "#ffffff",   # URL

  p_x       = 1.00,        # Title
  p_y       = 1.60,        # Title
  s_x       = 1.00,        # Subplot
  s_y       = 0.80,        # Subplot
  
  s_width   = 2.45,         # Subplot 
  s_height  = 2.5,         # Subplot
  
  url       = "multifaceted Functional Diversity",
  
  spotlight = FALSE,
  l_alpha   = 0.10,
  l_width   = 4,
  l_height  = 4
)

