DEL /F prefmod_*.tar.gz
R CMD build --compact-vignettes=gs+qpdf prefmod
R CMD check --as-cran prefmod_*.tar.gz
PAUSE