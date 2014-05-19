DEL /F prefmod_*.tar.gz
R CMD build prefmod
R CMD check --as-cran prefmod_*.tar.gz
PAUSE
