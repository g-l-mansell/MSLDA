#To rebuild package:
system("R CMD build")
system("R CMD INSTALL --preclean --no-multiarch --with-keep.source ../MSLDA")
devtools::document(roclets = c('rd', 'collate', 'namespace'))
system("rm man/documentation.pdf")
system("R CMD Rd2pdf --pdf --title='MS LDA Package Documentation' -o man/documentation.pdf man/*.Rd")
#devtools::test()
#devtools::install_github("g-l-mansell/RcppRidge")
#library(MSLDA)
