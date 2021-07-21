# install.packages('sysfonts')
# install.packages('hexSticker')
library(hexSticker)
library(showtext)

font_add_google("Ubuntu", "gf")
font_add_google("Roboto mono", "gf")
font_add_google("IBM Plex sans", "gf")
font_add_google("Abril Fatface", "gf")
font_add_google(name = "Goblin One", bold.wt = 800, family = "gf")
font_add_google(name = "Cinzel", bold.wt = 900, family = "gf")
## Automatically use showtext to render text for future devices
showtext_auto()

# library(showtext)
# ## Loading Google fonts (http://www.google.com/fonts)
# font_add_google("Gochi Hand", "gochi")
# ## Automatically use showtext to render text for future devices
# showtext_auto()

imgurl <- "./writing_and_slides/ritter3.png"
s <- sticker(
  subplot = imgurl,
  package = "Ritter", 
  p_family = 'gf',
  p_fontface =  c("bold"),
  p_size = 10.5,
  p_x = 1.025,
  p_y = 1.2,
  p_color = '#a61726', 
  spotlight = T,
  l_height = 5.5,
  l_width = 5.6,
  l_alpha = 0.24,
  l_x = 1, 
  l_y = 1.2,
  h_size = 0.8,
  h_fill = "#125c75",
  h_color = "#999999",
  s_x = 1,
  s_y = 1,
  s_width = 1.2,
  url = 'jmoggridge/ritter',
  filename = "writing_and_slides/ritter_hex_sticker.png"
)
plot(s)

