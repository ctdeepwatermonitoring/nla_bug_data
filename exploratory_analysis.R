library(vegan) #community ecology package with handy functions
library(ggplot2)

#read in data
spp <- read.csv('data/spp.csv', header=TRUE) #count data
smp <- read.csv('data/samples.csv', header=TRUE) #sample data
mst <- read.csv('data/master_taxa.csv', header=TRUE) #master taxa list

#transform count data from long to wide format with reshape function
cnt_spp <- spp[,c(2,9,10)]
colnames(cnt_spp)[2:3] <- c('OTU','CNT')
cnt_spp <- reshape(cnt_spp, idvar='UID', timevar = 'OTU', direction='wide')
cnt_spp[is.na(cnt_spp)] <- 0 #NAs to 0

#rename rows as UID in spp data
rownames(cnt_spp) <- cnt_spp[,1]
cnt_spp <- cnt_spp[,-1]

#rename rows as UID in smp data and add in year column
rownames(smp) <- smp[,5]
smp$smp_yr <- paste0('20',substr(smp$DATE_COL,8,9))


#rename columns as taxa name
for(i in 1:length(colnames(cnt_spp))){
  colnames(cnt_spp)[i] <- substr(colnames(cnt_spp)[i],5,
                                 nchar(colnames(cnt_spp)[i]))
}

#calculate relative abundance of taxa
rab_spp <- decostand(cnt_spp, 'total')

#calculate presence absence of taxa
pa_spp  <- decostand(cnt_spp, 'pa')

#custom function to sum taxa by column. pa = dataframe of taxa data
#you can reuse a function with another dataframe
#for example use the function with a single year instead of all sample years
sum_spp_ocr <- function(pa){
                  sum_spp <- as.data.frame(colSums(pa))
                  colnames(sum_spp)[1] <- 'spp_sum'
                  sum_spp$taxa <- rownames(sum_spp) 
                  sum_spp <- sum_spp[order(sum_spp$spp_sum, decreasing=TRUE),]
                  
                }

#run function for presence absence data 
sum_spp <- sum_spp_ocr(pa_spp)

#add ordering for ggplot
sum_spp$spp_ord <- ''
for(i in 1:dim(sum_spp)[1]){sum_spp[i,3] <- i}

#plot top 5 most frequently occurring taxa
ggplot(sum_spp[1:5,], aes(spp_ord,spp_sum))+
  geom_col()+
  scale_x_discrete(labels = sum_spp[1:5,2])+
  labs(y = 'Taxa Count', 
       title = 'Top 5 Most Frequently Occurring Taxa Across Samples') +
  theme(axis.title.x = element_blank())

#subset pa by year and call sum spp function to get taxa sum for 2007
smp_07    <- smp[smp$smp_yr == '2007',]
pa_spp_07 <- pa_spp[rownames(smp_07),]
sum_spp_07<- sum_spp_ocr(pa_spp_07) #call function 


