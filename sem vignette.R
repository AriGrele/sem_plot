load=function(pkgs,miss=c()){
  for(i in pkgs){if(!suppressWarnings(require(i,character.only=T,quietly=T))){miss=c(miss,i)}}
  if(length(miss)>0){
    if(readline(prompt=paste('install packages: ',paste(miss,collapse=', '),'? (y/n): ',sep=''))=='y'){
      install.packages(miss);load(miss)}}}
load(c('ggplot2','ggdist','ggforce','stringr','lavaan','lavaanPlot'))

devtools::source_url("https://github.com/AriGrele/sem_plot/blob/master/sem.R?raw=TRUE")

data=data.frame('plant'=rnorm(500,10,1))
data$herb=data$plant/10+rnorm(500,0,.1)
data$omni=data$herb/10+data$plant/10+rnorm(500,0,.01)
data$pred=data$herb/10+rnorm(500,0,.01)
plot(data)

model='herb~plant
       omni~herb+plant
       pred~herb'

fit=sem(model,scale(data,center=F))
lavaanPlot(name="example",fit,labels=colnames(data),coefs=T)

format='plant [pos=c(0,5),lab="plant\nbiomass"]
herb [pos=c(5,2.5),lab="herbivore\nbiomass"]
omni [pos=c(10,7.5),lab="omnivore\nbiomass"]
pred [pos=c(10,1),lab="carnivore\nbiomass"]
plant>herb [nudge=c(0,0,0,0),txt="1"]
plant>omni [nudge=c(0,0,0,0),txt="2"]
herb>omni [nudge=c(0,0,0,0),txt="3"]
herb>pred [txtnudge=c(0,0),txt="4"]'

eclair(format)
