DEL /F prefmod_*.tar.gz
C:\R\R-devel\bin\x64\R CMD build --compact-vignettes=gs+qpdf prefmod
START C:\R\R-devel\bin\x64\R CMD INSTALL --build prefmod_*.tar.gz
C:\R\R-devel\bin\x64\R CMD check --as-cran prefmod_*.tar.gz
PAUSE
