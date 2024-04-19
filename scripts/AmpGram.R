library("AmpGram", quietly = TRUE)
library("AmpGramModel", quietly = TRUE)
library("itertools", quietly = TRUE)
library("foreach", quietly = TRUE)
library("doParallel", quietly = TRUE)


######################## Package installation
# Note:
# If the R was installed under conda environment, some dependent packages might
# not be compatible, try to install those packages using conda or mamba

#### AmpGram
# source("https://raw.githubusercontent.com/r-lib/remotes/master/install-github.R")$value("michbur/AmpGram")
#### AmpGramModel
# install_AmpGramModel()
#### itertools
# install.packages("itertools", dependencies = TRUE, repos = "http://cran.us.r-project.org")
#### foreach
# install.packages("foreach", dependencies = TRUE, repos = "http://cran.us.r-project.org")
#### doParallel
# install.packages("doParallel", dependencies = TRUE, repos = "http://cran.us.r-project.org")
#############################################


args <- commandArgs(T)


# Set up parallel processing
num_cores <- 5
registerDoParallel(cores = num_cores)

sequences <- read_txt(args[1])

predictions <-  foreach(d=isplitVector(sequences, chunks=10),
                  .combine=c, .packages=c("stats")) %dopar% {
                  # Load libraries in each parallel run
                  library(AmpGram, quietly = TRUE)
                  library(AmpGramModel, quietly = TRUE)
                  predict(AmpGram_model, newdata=d, silent = TRUE)
          }

# Stop parallel processing
stopImplicitCluster()

df_prediction <- data.frame()

for (seq in names(predictions)) {
  seq_id <- seq
  prediction_prop <- as.vector(predictions[[seq]]$single_prot_pred)
  df_temp <- data.frame(seq_name = seq_id,
                        prob_AMP = prediction_prop)
  df_prediction <- rbind(df_prediction, df_temp)
}

write.table(file=args[2], df_prediction, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
