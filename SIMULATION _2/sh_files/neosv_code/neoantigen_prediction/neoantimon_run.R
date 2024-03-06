library(Neoantimon)
args <- commandArgs(trailingOnly = T)
file_sv <- args[1]
sample <- args[2]
outfile <- args[3]

file_hla <- 'lustre1/shiyang/neosv/neoantimon_hla.txt'
refflat <- '/lustre1/shiyang/reference_genomes/neoantimon/refFlat.grch37.txt'
refmrna <- '/lustre1/shiyang/reference_genomes/neoantimon/refMrna.grch37.fa'
netmhc <- '/lustre1/shiyang/software/netMHCpan-4.0/netMHCpan'

setwd(dirname(outfile))
Result_HLA1_SV <- MainSVFUSIONClass1(input_file = file_sv,
                                     file_name_in_hla_table = sample,
                                     hla_file = file_hla,
                                     refflat_file  = refflat,
                                     refmrna_file = refmrna,
                                     export_dir = paste0(sample, '_tmp'),
                                     netMHCpan_dir = netmhc,
                                     mutation_alt_bnd_column = 5,
                                     gene_symbol_column = 7,
                                     mate_id_column = 8)


write.table(Result_HLA1_SV, outfile, sep = '\t', row.names = F, col.names = T, quote = F)