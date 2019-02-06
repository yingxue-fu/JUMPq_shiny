## For information, https://github.com/talgalili/heatmaply

rm(list = ls())
library(heatmaply)

heatmaply(mtcars, k_row = 3, k_col = 2)

## It requires several development packages and it may not work well in clusters