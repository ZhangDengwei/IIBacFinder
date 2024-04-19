if (!require("ampir", quietly = TRUE)){
  install.packages("ampir", repos = "http://cran.us.r-project.org")
  require("ampir", quietly = TRUE)
}
  

args <- commandArgs(T)

my_protein_df <- ampir::read_faa(args[1])
my_prediction <- ampir::predict_amps(my_protein_df, model = "mature")

write.table(file=args[2], my_prediction, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
