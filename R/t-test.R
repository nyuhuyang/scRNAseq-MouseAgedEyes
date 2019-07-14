library(readxl)
library(dplyr)
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

df <- read.csv(paste0("output/20190713/Cell_type_number.csv"),header = T, row.names = 1)
rownames(df) = df[,1]
df = df[,-1]

p_value = c()
ColSum <- colSums(df)

for(i in 1:nrow(df)){
        conting <- rbind(df[i,],ColSum-df[i,])
        FISH <- fisher.test(conting,conf.int = T)
        p_value[i] = FISH$p.value
        #CHI = chisq.test(conting, correct = T)
        #chisq_p_value[i] = CHI$p.value             
}

df$p_value = p_value
df$p_val_adj = p.adjust(p = df$p_value, method = "bonferroni", 
                               n = nrow(df))
df %>% kable %>% kable_styling()
write.csv(df,paste0(path,"Cell_type_number_pval.csv"))
