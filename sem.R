#### Checks if vector between two other vectors, used to determine side path comes from ####
inside=function(v1,v2,v3){
  a1=atan2(v1[2],v1[1]);a2=atan2(v2[2],v2[1]);a3=atan2(v3[2],v3[1])
  if(!(a1>=0&0>a2)){if(sign(a1)==-1){a1=a1+2*3.1415926535};if(sign(a2)==-1){a2=a2+2*3.1415926535};if(sign(a3)==-1){a3=a3+2*3.1415926535}}
  return(a1>=a3&a3>a2)}

#### merges multiple sem outputs into single data.frame ####
semmerge=function(f1,f2,...){
  l=c(f1,f2,list(...))
  out=NULL
  for(i in 1:length(l)){
    e=tryCatch(parameterEstimates(l[[i]]),error=function(x){
      cat("Could not extract parameters from model\n")})
    e$groups=(1:length(l))[i]
    if(is.null(out)){out=e}
    else{out=rbind(out,e)}}
  return(out)
}

#### creates ggplot object of paramter estimates and CIs ####
sembars=function(fit=NULL,s=1,mask=NULL,groups=NULL,flip=F,group='groups'){
  if(group!='groups'){group='name';facet='groups'}
  else{facet='name'}
  if(!is.null(fit)){
    e=tryCatch(est_choice(fit),error=function(x){
      if(is.data.frame(fit)){return(fit)}
      cat("Could not extract parameters from model\n")})
    e=e[e$op=='~',]
    e$name=paste(e$rhs,e$lhs)
    if(!is.null(mask)){for(n in names(mask)){e$name=gsub(n,mask[n],e$name)}}
    if(!is.null(groups)){for(n in names(groups)){e$group=gsub(n,groups[n],e$group)}}
    if('groups' %in% names(e)){
      g=ggplot(e,aes_string(x=group))+
        geom_pointinterval(aes(y=est,ymin=ci.lower,ymax=ci.upper),size=s)+
        theme_classic()+
        geom_hline(yintercept = 0)+
        theme(panel.border = element_rect(fill=NA),
              panel.spacing=unit(0,'lines'))+
        guides(color=F,fill=F)
      if(flip){g=g+facet_grid(reformulate('.',facet))}
      else{g=g+facet_grid(.~name)}
    }
    else{
      g=ggplot(e,aes(x=paste(rhs,lhs)))+
        geom_pointinterval(aes(y=est,ymin=ci.lower,ymax=ci.upper),size=s)+
        theme_classic()+
        geom_hline(yintercept = 0)}
    if(flip){g=g+coord_flip()}
    return(g)}
  cat("Model needed\n")
  return(NULL)}

#### bayesian probability of direction ####
pd=function(x){
  side=sign(median(x))
  return(sum(side*x>0)/length(x))}

#### extracts parameter esitmates from jags sem output###
#### formats to match lavaan output, but includes pd, credibility intervals
bparam=function(x){
  x$sims.list[['deviance']]=NULL
  sims=x$sims.list
  out=setNames(as.data.frame(matrix(ncol=9)),c('lhs','op','rhs','est','se','z','pvalue','ci.lower','ci.upper'))
  for(c in 1:length(sims)){
    for(i in 2:ncol(sims[[c]])){
      co=sims[[c]][,i]
      out=rbind(out,data.frame('lhs'=names(x$sims.list)[c],'op'='~','rhs'=i,'est'=mean(co),'se'=sd(co)/sqrt(length(co)),'z'=NA,'pvalue'=pd(co),'ci.lower'=quantile(co,0.025),'ci.upper'=quantile(co,1-0.025)))}}
  return(out[-1,])}

#### chooses which parameter function to use ####
est_choice=function(x){return(switch(class(x),'lavaan'=parameterestimates,'jagsUI'=bparam)(x))}

#### Creates ggplot object of specified path model ####
#### accepts model string, Lavaan or Jags model output, scale factor for paths/bullets/text, scale factor for arrows, alpha value to compare p value, toggels for outlines, text, filter to rename lavaan output to match model string ####
path=function(model,fit=NULL,...){ 
  default=list('pos'=c(0,0),                                                                    #Default values for vars missing from model string
               'size'=c(1,1),
               'lab'='',
               'nudge'=c(0,0,0,0),
               'txtnudge'=c(0,0),
               'txt'='')
  parm=list(...)                                                                                #default values for vars missing from function call
  if(is.null(parm$scal))   {parm$scal=5}
  if(is.null(parm$ascal))  {parm$ascal=c(1,1)*parm$scal}
  else{parm$ascal=parm$ascal*parm$scal*c(1,1)}
  if(is.null(parm$size))   {parm$size=10}
  if(is.null(parm$alpha))  {parm$alpha=10}
  if(is.null(parm$outline)){parm$outline=T}
  if(is.null(parm$mask))   {parm$mask=NULL}
  if(is.null(parm$txtoff)) {parm$txtoff=F}
  if(is.null(parm$autonudge)){parm$autonudge=T}
  if(!is.null(fit)){                                                                             #if model in call:
    e=tryCatch(est_choice(fit),error=function(e){
      cat("Could not extract parameters from model\n")})                                         #get parameter estimates
    e=e[e$op=='~',]
    e$name=paste(e$rhs,'>',e$lhs,sep='')                                                         #convert them to match model string format
    e$pvalue=e$pvalue>parm$alpha                                                                 #add p value 
    e$nest=ceiling(abs(e$est/max(e$est)*parm$scal))*sign(e$est)                                  #add modified estimates
    if(!is.null(parm$mask)){for(n in names(parm$mask)){e$name=gsub(n,parm$mask[n],e$name)}}}     #rename vars
  
  lines=str_split(model,'\\]\n')[[1]]                                                            #split model string into lines
  lines=gsub('\\]$',')',gsub('\\[','list(',lines))                                               #convert to list format
  lines[-length(lines)]=paste(lines[-length(lines)],')',sep='')
  info=list();arrow=list()
  for(l in 1:length(lines)){                                                                     #evaluate lines as lists, extract data into info list for boxes and arrow list for paths
    if('>'%in%strsplit(lines[l],'')[[1]]){arrow[[str_match(lines[l],'^(.+) .+$')[,2]]]=eval(parse(text=str_match(lines[l],'list\\(.+\\)$')))}
    else{info[[str_match(lines[l],'^(.+) list')[,2]]]=eval(parse(text=str_match(lines[l],'list.+')))}}

  for(i in names(info)){                                                                         #fill in default values for info and arrow if missing from names
    for(d in names(default)){
      if(!d%in%names(info[[i]])){info[[i]][[d]]=default[[d]]}}}
  for(i in names(arrow)){
    for(d in names(default)){
      if(!d%in%names(arrow[[i]])){arrow[[i]][[d]]=default[[d]]}}}
  
  box=setNames(as.data.frame(matrix(ncol=5)),c('n','x','y','w','h'))                             #create dataframe of boxes for plot
  for(n in names(info)){                                                                         #populate with data from info
    box=rbind(box,data.frame('n'=info[[n]]$lab,'x'=info[[n]]$pos[1],'y'=info[[n]]$pos[2],'w'=info[[n]]$size[1],'h'=info[[n]]$size[2]))}
  box=na.omit(box)
  arrows=list();points=list();bpoints=list()                                                     #create lists for data about arrow paths, arrow heads
  for(n in names(arrow)){                                                                        #for each path:
    if(is.null(fit)){row=data.frame('est'=1,'nest'=1,'pvalue'=0>parm$alpha)}                     #default values
    else{row=e[e$name==n,]}
    ns=str_split(n,'>')[[1]]                                                                     #split path into source, destination
    v1=info[[ns[1]]];v2=info[[ns[2]]]                                                            #info for source and destination
    sizes=c(v1$size,v2$size)                                                                     #extract box sizes
    v3=(v2$pos-v1$pos)                                                                           #vector between box centers
    ins=c()
    for(i in 1:4){                                                                               #for a vector from the center of the source box to each corner of the source box, find if vector towards destination between two corners
      corners=list(c(-1,1)*v1$size,
                   c(1,1)*v1$size,
                   c(1,-1)*v1$size,
                   c(-1,-1)*v1$size,
                   c(-1,1)*v1$size)
      ins=c(ins,inside(corners[[i]],corners[[i+1]],v3))
    }
    side=match(T,ins)                                                                            #find side based on which angles vector to destination is between
    ascal=list(c(parm$ascal),c(parm$ascal[2],parm$ascal[1]))[[(side%in%c(2,4))+1]]               #rotate arrow scale based on side
    if(sign(row$est)>0){points[[n]]=list(data.frame('x'=c(0,.5,0,-.5),'y'=c(0,2,1.5,2)*.75),     #select arrowhead based on side
                                         data.frame('x'=c(0,2,1.5,2)*.75,'y'=c(0,.5,0,-.5)),
                                         data.frame('x'=c(0,.5,0,-.5),'y'=((2-c(0,2,1.5,2))-2)*.75),
                                         data.frame('x'=((2-c(0,2,1.5,2))-2)*.75,'y'=c(0,.5,0,-.5)))[[c(3,4,1,2)[side]]]/18
    ss=max(abs(c(points[[n]]$x[3],points[[n]]$y[3])))*parm$scal*abs(row$nest)                    #calculate distance factor based on size of arrow & path
    v2$size=v2$size*1.125+ss*2}                                                                  #rescale size to add border
    else{points[[n]]='point';ss=1.5*parm$scal*abs(row$nest)/2/20;v2$size=v2$size*1.25+ss}        #bullets if est negative
                                                                                                 #select control points for bezier and nudge based on data
    nu=arrow[[n]]$nudge*list(c(sizes[1],1,sizes[3],1),c(1,sizes[2],1,sizes[4]))[[(side%in%c(2,4))+1]]
    arrows[[n]]=list(list('x'=c(v1$pos[1],v1$pos[1],v2$pos[1],v2$pos[1])+c(nu[1],nu[1],nu[3],nu[3]),
                          'y'=c(v1$pos[2],v2$pos[2],v1$pos[2],v2$pos[2])+c(0,nu[2],nu[4],0)),
                     list('x'=c(v1$pos[1],v2$pos[1],v1$pos[1],v2$pos[1])+c(0,nu[1],nu[3],0),
                          'y'=c(v1$pos[2],v1$pos[2],v2$pos[2],v2$pos[2])+c(nu[2],nu[2],nu[4],nu[4])))[[(side%in%c(2,4))+1]]
    
    arrows[[n]]$x=                                                                               #adjust paths to point to edge of box, not center
      list(arrows[[n]]$x+c(0,0,0,0),
           arrows[[n]]$x+c(v1$size[1]/2,-v2$size[1]/2,v1$size[1]/2,-v2$size[1]/2),
           arrows[[n]]$x+c(0,0,0,0),
           arrows[[n]]$x+c(-v1$size[1]/2,v2$size[1]/2,-v1$size[1]/2,v2$size[1]/2))[[side]]
    arrows[[n]]$y=
      list(arrows[[n]]$y+c(v1$size[2]/2,-v2$size[2]/2,v1$size[2]/2,-v2$size[2]/2),
           arrows[[n]]$y+c(0,0,0,0),
           arrows[[n]]$y+c(-v1$size[2]/2,v2$size[2]/2,-v1$size[2]/2,v2$size[2]/2),
           arrows[[n]]$y+c(0,0,0,0))[[side]]
    arrows[[n]]$est=row$nest;arrows[[n]]$p=row$pvalue;                                           #extract fit measures
    if(!parm$txtoff){arrows[[n]]$lab=arrow[[n]]$txt}                                             #extract labels
    else{arrows[[n]]$lab=round(row$est,2)}
    as=1.1;                                                                                      #scale values for labels, paths
    ax=c(0,ss,0,-ss)[side];ay=c(ss,0,-ss,0)[side]
    
                                                                                                 #properly scale arrows
    if(length(points[[n]])==1){points[[n]]=list('x'=arrows[[n]]$x[4]+ax*0,'y'=arrows[[n]]$y[4]+ay*0)}
    else{bpoints[[n]]$x=points[[n]]$x*ascal[1]*abs(arrows[[n]]$est)*as+arrows[[n]]$x[4]+ax;bpoints[[n]]$y=points[[n]]$y*ascal[2]*abs(arrows[[n]]$est)*as+arrows[[n]]$y[4]+ay
    points[[n]]$x=points[[n]]$x*ascal[1]*abs(arrows[[n]]$est)+arrows[[n]]$x[4]+ax;points[[n]]$y=points[[n]]$y*ascal[2]*abs(arrows[[n]]$est)+arrows[[n]]$y[4]+ay}
    space=ifelse(parm$autonudge,parm$scal*abs(arrows[[n]]$est)/10,0)                             #space at end of path
    arrows[[n]]$midx=(arrows[[n]]$x[4]+arrows[[n]]$x[1])/2+arrow[[n]]$txtnudge[1]-space          #find midpoints of paths
    arrows[[n]]$midy=(arrows[[n]]$y[4]+arrows[[n]]$y[1])/2+arrow[[n]]$txtnudge[2]+space
  }
  
  pal=c("#CC0C00FF","#5C88DAFF")                                                                 #default palette
  s=ggplot()                                                                                     #new ggplot object
  nl=ifelse(rep(is.null(fit),length(names(arrows))),names(arrows),e$name[order(-abs(e$est))])    #for each path:
  for(n in nl){
    if(parm$outline){s=s+geom_bezier(data=as.data.frame(arrows[[n]]),                            #add white beziers
                    aes(x=x,y=y),
                    size=parm$scal*abs(arrows[[n]]$est)*c(1.5,1)[(arrows[[n]]$p)+1],
                    color='white',
                    linetype=c('solid','dotted')[(arrows[[n]]$p)+1])}
    s=s+geom_bezier(data=as.data.frame(arrows[[n]]),                                             #add main beziers
                    aes(x=x,y=y,color=as.factor(sign(est))),
                    size=parm$scal*abs(arrows[[n]]$est),
                    linetype=c('solid','dotted')[(arrows[[n]]$p)+1])
    if(length(points[[n]]$x)<4){                                                                 #bullets
      s=s+geom_point(data=as.data.frame(points[[n]]),aes(x=x,y=y,color=as.character(-1)),size=1.5*parm$scal*abs(arrows[[n]]$est))}
    else{                                                                                        #arrows
      s=s+geom_polygon(data=as.data.frame(points[[n]]),aes(x=x,y=y,fill=as.character(1),color=as.character(1)),size=1)}
    s=s+geom_text(data=as.data.frame(arrows[[n]]),                                               #labels
                  aes(x=midx,y=midy,label=lab),
                  color='black',
                  size=parm$size)}
                                                                                                 #add boxes and text
  s=s+geom_rect(data=box,aes(xmin=x-w/2,xmax=x+w/2,ymin=y+h/2,ymax=y-h/2),fill='white', colour = "black")+
      geom_text(data=box,aes(x=x,y=y,label=n),size=parm$size,color='black')
  s=s+guides(fill=F,color=F)+                                                                    #remove key, axes
      theme_void()
  if(is.null(fit)){s=s+scale_color_manual(values='black')+                                       #black for metamodels
                      scale_fill_manual(values='black')}
  else{s=s+scale_color_manual(values=c('-1'=pal[1],'1'=pal[2]))+                                 #pal for fit models
           scale_fill_manual(values=c('-1'=pal[1],'1'=pal[2]))}
  return(s)}
