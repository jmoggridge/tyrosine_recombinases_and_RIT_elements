install.packages('sysfonts')
install.packages('hexSticker')
library(hexSticker)


# library(showtext)
# ## Loading Google fonts (http://www.google.com/fonts)
# font_add_google("Gochi Hand", "gochi")
# ## Automatically use showtext to render text for future devices
# showtext_auto()

imgurl <- "./writing_and_slides/fritter.png"
s <- sticker(subplot = imgurl, package="fritter",
             p_size=10, s_x=.95, s_y=.75,
             h_fill="#125c75", h_color="#000000",
             s_width=.65,filename="writing_and_slides/sticker.png")
plot(s)
