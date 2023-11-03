# Pre - Knit
knitr::knit("vignettes/articles/ABS-Parameter-Recovery.Rmd.config", "vignettes/ABS-Parameter-Recovery.Rmd")
# Move figures
dir.create("vignettes/figure/", showWarnings = F)
file.copy(from="figure/rt-dist-1.png", to = "vignettes/figure/rt-dist-1.png")
