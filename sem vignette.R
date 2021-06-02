#######################################################################################################
#### sem plots vignette, 2021/5/27 ####
#######################################################################################################

#### load or install required packages #### 
load=function(pkgs,miss=c()){
  for(i in pkgs){if(!suppressWarnings(require(i,character.only=T,quietly=T))){miss=c(miss,i)}}
  if(length(miss)>0){
    if(readline(prompt=paste('install packages: ',paste(miss,collapse=', '),'? (y/n): ',sep=''))=='y'){
      install.packages(miss);load(miss)}}}
load(c('ggplot2','ggdist','ggforce','stringr','lavaan','lavaanPlot'))

#### install plotting functions from github ####
devtools::source_url("https://github.com/AriGrele/sem_plot/blob/master/sem.R?raw=TRUE")

#### generate example data ####
set.seed(11)
data=data.frame('plant'=rnorm(500,10,1))
data$herb=data$plant/10+rnorm(500,0,.1)
data$omni=data$herb/10+data$plant/10+rnorm(500,0,.01)
data$pred=data$herb/10-data$omni/2+rnorm(500,0,.1)
plot(data)

#### lavaan sem ####
model='herb~plant
       omni~herb+plant
       pred~herb+omni'

fit=sem(model,scale(data,center=F))
lavaanPlot(name="example",fit,labels=colnames(data),coefs=T)

#######################################################################################################
#### create path model figures ####
# path() takes formating string, optional Lavaan or Jags model output
# formating string is split into boxes, paths between boxes
# boxes take arguments of position, size, label
# paths take arguments of x/y axis nudges from default positions, x/y axis nudges for label position, labels
# default arguments:
#   pos      = c(0,0)
#   size     = c(1,1)
#   lab      = ''
#   nudge    = c(0,0,0,0) # source x, source y, destination x, destination y. This moves control points for the bezier curves
#   txtnudge = c(0,0)     # x, y
#   txt      = ''
# missing values are substituted with defaults

format='plant [pos=c(0,5),size=c(1,6),lab="Plant\\nbiomass"]
herb [pos=c(5,2.5),lab="Herbivore\\nbiomass"]
omni [pos=c(10,7.5),lab="Omnivore\\nbiomass"]
pred [pos=c(10,2.5),lab="Carnivore\\nbiomass"]
plant>herb [nudge=c(0,0,0,0),txtnudge=c(0,0),txt="1"]
plant>omni [nudge=c(0,0,0,.33),txt="2"]
herb>omni [nudge=c(0,0,0,-.33),txt="3"]
herb>pred [txt="4"]
omni>pred [txt="5"]'

#### example meta-model ####
# default arguments:
#   model     = no default  # formatting string 
#   fit       = NULL        # lavaan or jags model output
#   scal      = 5           # scaling factor for paths
#   ascal     = c(1,1)*scal # scaling factor for arrowheads
#   size      = 10          # scaling factor for text
#   alpha     = 10          # alpha threshold for dotted paths based on p-value for lavaan, pd for jags
#   outline   = TRUE        # white outline around paths when true
#   mask      = null        # mask if names in formatting string don't match names in lavaan or jags output
#   txtoff    = FALSE       # replaces path labels with parameter estimates when true
#   autonudge = TRUE        # nudges path labels off of paths when true

e1=path(format,scal=3,size=5)
png('meta_example.png',1000,500,type='cairo');e1;dev.off();system2('open','meta_example.png')

#### example fitted model ####
e2=path(format,fit,scal=2,size=5,txtoff=T)
png('fitted_example.png',1000,500,type='cairo');e2;dev.off();system2('open','fitted_example.png')

#### possible to define formatting models where variable names don't match the data using a mask ####
format2='a [pos=c(0,5),lab="Plant\\nbiomass"]
b [pos=c(5,2.5),lab="Herbivore\\nbiomass"]
c [pos=c(10,7.5),lab="Omnivore\\nbiomass"]
d [pos=c(10,2.5),lab="Carnivore\\nbiomass"]
a>b [txt="1"]
a>c [nudge=c(0,0,0,.33),txt="2"]
b>c [nudge=c(0,0,0,-.33),txt="3"]
b>d [txt="4"]
c>d [txt="5"]'

m=c('plant'='a','herb'='b','omni'='c','pred'='d')
e3=path(format2,fit,scal=2,size=5,autonudge=F,alpha=0.05,mask=m,coloralpha=.05)
png('mask_example.png',1000,500,type='cairo');e3;dev.off();system2('open','mask_example.png')

#######################################################################################################
#### these functions and their outputs may change dramatically in the future ####
#### merge parameter estimates from two models ####
(merged = semmerge(fit,sem('herb~plant\npred~plant',scale(data,center=F))))

#### plot parameter estimates and CIs ####
sembars(fit,s = 5,flip=F)
sembars(merged,s = 10,flip=T,group='groups')
merged
