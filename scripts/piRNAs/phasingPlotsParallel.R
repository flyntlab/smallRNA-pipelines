#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
phas <- read.csv(input_file, sep = "\t")
colnames(phas) <- c("nt","a","c","g","t")
phas[2:5] <- scale(phas[2:5])

write.table(phas, file = output_file, row.names = FALSE)

#pdf("phasingTrace.pdf")

#ggplot() +
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
phas <- read.csv(input_file, sep = "\t")
colnames(phas) <- c("nt","a","c","g","t")
phas[2:5] <- scale(phas[2:5])

write.table(phas, file = output_file, row.names = FALSE)
