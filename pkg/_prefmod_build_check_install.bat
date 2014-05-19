DEL /F prefmod_*.tar.gz
R CMD build prefmod
START R CMD INSTALL prefmod_*.tar.gz
R CMD check --as-cran prefmod_*.tar.gz
PAUSE
