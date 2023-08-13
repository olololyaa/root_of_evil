#Compare GO term overrepresentation analysis in 2 groups of bacteria
library(stats)
library(stringr)
library(graphics)
my.t.test.p.value <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
    if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}
#a function that shows GO term prevalence in the 2 groups
get_GO_distr <- function(GO_term,to_plot = F, main_ = GO_term,
                         cex.main_ = 4){
    RowNum = which(table1[,1]==GO_term)
    groupnames = c("Health-prevalent and neutral","Health-scarce")
    print("Health-prevalent and neutral")
    print(table(unlist(table1[RowNum,grepl('^MH|^Neutral',colnames(table1))])))
    print("Health-scarce")
    print(table(unlist(table1[RowNum,grepl('^MN',colnames(table1))])))
    if(to_plot){
        occurency_control = as.numeric(unlist(table1[RowNum,grepl('^MH|^Neutral',colnames(table1))]))
        occurency_nonhealthy = as.numeric(unlist(table1[RowNum,grepl('^MN',colnames(table1))]))
        jpeg(filename = str_c("/media/lev-genetik/OS/Yerevan_Root_of_Evil/Results/GO_terms/boxplots/",
                              GO_term,".jpg"), width = 600, height = 400)
        print(boxplot(occurency_control,occurency_nonhealthy,varwidth = T,#names = groupnames,
                main = main_, cex.main = cex.main_, cex.axis = 3, 
                border = c("cadetblue","firebrick"),lwd=3))
        dev.off()
    }
}

table1 <- read.csv("/media/lev-genetik/OS/Yerevan_Root_of_Evil/Results/occurence_table_MN_MH_Neutral.tsv",sep = "\t")

vec1 = c(0,0)
result <- data.frame("Column1" = numeric(0), "Column2" = numeric(0))
colnames(result) = c('MH+Neutral','MN')
for (i in (1:dim(table1)[1])){
    print(i)
    vec1[1] = sum(table1[i,grepl('^MH|^Neutral',colnames(table1))]>0)
    vec1[2] = sum(table1[i,grepl('^MN',colnames(table1))]>0)
    result = rbind(result,vec1)
}
colnames(result) = c('MH+Neutral','MN')
chisq.test(result)


#Wilcoxon Rank sum test - done instead of the T test
Wilcox_test_results_withNeutral = c()
for (i in (1:dim(table1)[1])){
    print(i)
    groupH = as.vector(table1[i,grepl('^MH|^Neutral',colnames(table1))])
    groupH = as.numeric(groupH)
    groupN = as.vector(table1[i,grepl('^MN',colnames(table1))])
    groupN = as.numeric(groupN)
    Wilcox_test_results_withNeutral[i] = wilcox.test(groupH,groupN,exact=F)$p.value
    names(Wilcox_test_results_withNeutral)[i] = as.character(table1[i,1])
}
Wilcox_test_results_withNeutral = Wilcox_test_results_withNeutral[
    !is.na(Wilcox_test_results_withNeutral)] #removing 9 NA values
which.min(Wilcox_test_results_withNeutral)
Wilcox_test_results_withNeutral[which.min(Wilcox_test_results_withNeutral)]

#multiple testing correction
Wilcox_test_results_withNeutral_corr = p.adjust(Wilcox_test_results_withNeutral,
                                           method = 'hochberg')
Wilcox_test_results_withNeutral_ord = Wilcox_test_results_withNeutral_corr[
    order(Wilcox_test_results_withNeutral_corr)]
Wilcox_test_results_withNeutral_ord[1:5]

#get GO terms enriched with corrected p-value < 0.01
GO_champions_001 <- names(Wilcox_test_results_withNeutral_ord)[
    Wilcox_test_results_withNeutral_ord<0.01]
writeLines(GO_champions_001,
           con = "/media/lev-genetik/OS/Yerevan_Root_of_Evil/Results/GO_terms/GO_champions_0_01.txt",
           sep = "\n"
)

#Making boxplots for the 5 top hits
get_GO_distr("GO:0016811",main_ = "Hydrolase activity, acting on carbon-nitrogen bonds", cex.main_ = 2, to_plot = T)
get_GO_distr("GO:0015658",main_ = "Br.-chain AA transmembrane transporter activity", cex.main_ = 2, to_plot = T)
get_GO_distr("GO:0015803",main_ = "Branched-chain amino acid transport", cex.main_ = 2, to_plot = T)
get_GO_distr("GO:0015205",main_ = "Nucleobase transmembrane transporter activity", cex.main_ = 2, to_plot = T)
get_GO_distr("GO:0006790",main_ = "Sulfur compound metabolic process", cex.main_ = 2, to_plot = T)


