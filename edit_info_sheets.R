library(stringi)
library(stringr)
library(rebus)
library(lubridate)

HD02_sample_info<-read.delim("HD02_sample_info.txt",stringsAsFactors = F)
HD03_sample_info<-read.delim("HD03_sample_info.txt",stringsAsFactors = F)
HD02_sample_info<-HD02_sample_info[,-2]
HD02_sample_info<-HD02_sample_info[,-24]
HD03_sample_info<-HD03_sample_info[!(HD03_sample_info$pool == "No"),]
names(HD02_sample_info)
names(HD03_sample_info)
HD02_sample_info[HD02_sample_info$batch == 1,"batch"]<-"HD01"
HD02_sample_info[HD02_sample_info$batch == 2,"batch"]<-"HD02"
HD03_sample_info[HD03_sample_info$pool == 1,"pool"]<-"HD03_pool_1"
HD03_sample_info[HD03_sample_info$pool == 2,"pool"]<-"HD03_pool_2"
HD03_sample_info[HD03_sample_info$pool == 3,"pool"]<-"HD03_pool_3"
head(HD02_sample_info)
head(HD03_sample_info)
names(HD03_sample_info) <- names(HD02_sample_info)

HD02_sample_info$sort_date<-ymd(HD02_sample_info$sort_date)
HD02_sample_info$RT_date<-ymd(HD02_sample_info$RT_date)
HD02_sample_info$prep_date<-ymd(HD02_sample_info$prep_date)
HD03_sample_info$sort_date<-dmy(HD03_sample_info$sort_date)
HD03_sample_info$RT_date<-dmy(HD03_sample_info$RT_date)

HD03_sample_info$prep_date<-paste0(HD03_sample_info$prep_date,"-17")
HD03_sample_info$prep_date<-dmy(HD03_sample_info$prep_date)

HD03_sample_info$project<-"HD"

HD03_sample_info$age<-str_replace(HD03_sample_info$age, "e","E")
HD02_sample_info$age<-str_replace(HD02_sample_info$age, capture("12.5"),"E" %R% REF1)

combined_sample_info<-rbind(HD02_sample_info,HD03_sample_info)

?write.table

write.table(combined_sample_info, "sample_info.txt",quote=F,sep="\t")
