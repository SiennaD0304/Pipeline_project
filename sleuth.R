library(sleuth)
library(dplyr)


# used snakemake@input for R to recognive input for sleuth rule 
stab = read.table(snakemake@input[["table"]],header = TRUE)
so = sleuth_prep(stab)

#got code structure from slides

so = sleuth_fit(so,~condition, 'full') # fits a model to compare conditions 
so = sleuth_fit(so, ~1, 'reduced') # compares likelihood test 
so = sleuth_lrt(so, 'reduced', 'full') #runs the likelihood test for with different contions 

#getting test results from so 
sleuth_table = sleuth_results(so,"reduced:full", "lrt", show_all = FALSE)

#filtering to only show most significant results of kallisto 
sleuth_significant = dplyr::filter(sleuth_table,qval<=0.05)|>dplyr::arrange(qval)

# taking the information needed for output data
sleuth_important= dplyr::select(sleuth_significant,target_id,test_stat,pval,qval )


write.table(sleuth_important, snakemake@output[["sleuth_out"]],quote = FALSE,row.names = TRUE)
