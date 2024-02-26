#' Output for tests 
load("fixtures.RData")


# Set up spider data
X <- as.data.frame(spider$x)
abund <- spider$abund 

pa <-(abund>0)*1


library(labdsv)

# site data
data(brycesite)
brycesite$plotcode=substr(brycesite$plotcode,3,5)
#species data
data(bryceveg)
bryceveg<- bryceveg[,-which(colSums(bryceveg>0) <= 20)]  #most abundant species


#recode data to integer categories
old <- c(0,0.2,0.5,1,2,3,4,5,6) #existing categories
bryceord = bryceveg

for(i in 1:length(old)){
  bryceord[bryceveg==old[i]]=i
}
