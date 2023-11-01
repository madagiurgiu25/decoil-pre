install.packages('rmarkdown', dependencies = TRUE,quietly = TRUE, update = FALSE)
install.packages("devtools", dependencies = TRUE,quietly = TRUE, update = FALSE)
install.packages("ggplot2", dependencies = TRUE, quietly = TRUE, update = FALSE)
install.packages("dplyr", dependencies = TRUE,quietly = TRUE, update = FALSE)
install.packages("assertthat", dependencies = TRUE,quietly = TRUE, update = FALSE)
install.packages("DT", dependencies = TRUE,quietly = TRUE, update = FALSE)
install.packages("pafr", dependencies = TRUE,quietly = TRUE, update = FALSE)
install.packages("reshape2", dependencies = TRUE,quietly = TRUE, update = FALSE)
install.packages('igraph', repos='http://cran.rstudio.com/', quietly = TRUE)
install.packages("kableExtra", dependencies = TRUE, update = FALSE, quietly = TRUE)

BiocManager::install('GenomicRanges', dependencies = TRUE, update = FALSE)
BiocManager::install('rtracklayer', dependencies = TRUE, update = FALSE)
BiocManager::install('plyranges', dependencies = TRUE, update = FALSE)
BiocManager::install('Gviz', dependencies = TRUE, update = FALSE)

devtools::install_github('mskilab/gUtils', dependencies = TRUE, update = FALSE)
# gTrack_0.1.0.tar.gz
devtools::install_github('mskilab/gTrack', dependencies = FALSE, update = FALSE)
# fishHook_0.1.tar.gz
devtools::install_github('mskilab/fishHook', dependencies = FALSE, update = FALSE)
# bamUtils_0.0.0.9000.tar.gz
devtools::install_github('mskilab/bamUtils', dependencies = FALSE, update = FALSE)
# gGnome_0.1.tar.gz (uses all the above libraries and gChain_0.2.0.tar.gz, igraph, reshape2)
devtools::install_github('mskilab/gGnome', dependencies = TRUE, update = FALSE)
