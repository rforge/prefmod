DEL /F prefmod_*.tar.gz *check.log
R CMD build --compact-vignettes=gs+qpdf prefmod
R CMD check --as-cran prefmod_*.tar.gz
COPY prefmod.Rcheck\*check.log .
START prefmod.Rcheck\prefmod-manual.pdf
PAUSE