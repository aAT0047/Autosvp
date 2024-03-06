library(ConsensusClusterPlus)
file_sig <- './gsva_score.txt'
data_sig <- read.csv(file_sig, sep = '\t', header = T, check.names = F, stringsAsFactors = F)
data_sig <- data_sig[grep('boston', row.names(data_sig), ignore.case = T), ]
data_sig <- as.matrix(data_sig)
res <-  ConsensusClusterPlus(d = data_sig, maxK = 6, reps = 50, pItem = 0.8, pFeature = 1, clusterAlg = 'km', title = 'km6_gsva', )
res <- res[[4]]

res <- res$consensusClass
res <- as.data.frame(res)
res$Sample <- row.names(res)

outfile <- './km_cluster.txt'
write.table(res, outfile, sep = '\t', row.names = F, col.names = T, quote = F)
