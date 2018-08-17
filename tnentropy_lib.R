# author: Marcin Kulisiewicz
# 
# This script contains methods for computations of entropy for temporal network represented as event 
# sequence. Event sequence basic format is 3 column table. Columns:
# x - sender, initiator of event (if direction taking into account)
# y - receiver
# t - time of event

# load require libraries
require(abind)

# computing set of entropies for given temporal network
allEntropies <- function(data, directed=T){
  o = tnEntropy(data,"o",directed)
  ob = tnEntropy(data,"oboth",directed)
  e = tnEntropy(data,"e",directed)
  s = tnEntropy(data,"s",directed)
  om = tnEntropy(data,"om",directed)
  obm = tnEntropy(data,"obothm",directed)
  em = tnEntropy(data,"em",directed)
  sm = tnEntropy(data,"sm",directed)
  or = o/om
  obr = ob/obm
  er = e/em
  sr = s/sm
  return(list(o=o,ob=ob,e=e,s=s,om=om,obm=obm,em=em,sm=sm,or=or,obr=obr,er=er,sr=sr))
}

# computing entropy for given temporal network - v1.2
tnEntropy <- function(data, mode = "e", directed=T, nodeID=NA){
  if(!directed){
    data[data[,1]>data[,2],] = data[data[,1]>data[,2],c(2,1,3)]
  }
  if(mode == "e"){ # 2nd order entropy
    if(nrow(data)==1){
      return(0)
    }else if(nrow(data)==0){
      return(0)
    }else{
      events = do.call(paste,data[,1:2]) # events
      uevents = unique(events) # unique events
      d=cbind(data,ids=match(events, uevents)) # give ids to unique events
      ejn <- table(d[,4]) # count each event
      en <- sum(ejn) # number of events
      return(-sum(sapply(ejn,function(x){(x/en)*log(x/en)}),na.rm = TRUE))
    }
  }else if(mode == "em"){ # max second order entropy
    if(nrow(data)==1){
      return(0)
    }else if(nrow(data)==0){
      return(0)
    }else{
      nn <- length(unique(c(data[,1],data[,2])))
      if(directed){
        return(log(nn*(nn-1)))
      }else{
        return(log(nn*(nn-1))-log(2))
      }
    }
  }else if(mode == "s"){ # 3rd order entropy
    if(nrow(data)==2){
      return(0)
    }else if(nrow(data)==0 | nrow(data)==1){
      return(0)
    }else{
      events = do.call(paste,data[,1:2]) # events
      uevents = unique(events) # unique events
      data=cbind(data,ids=match(events, uevents)) # give ids to unique events
      succ=as.data.frame(cbind(data[1:(nrow(data)-1),4],data[2:nrow(data),4]))
      successions = do.call(paste,succ) # successions
      usuccessions <- unique(successions) # unique successions
      s=cbind(succ,ids=match(successions,usuccessions)) # give ids to unique successions
      si <- table(s[,3]) # count each succession
      sn <- sum(si) # number of successions
      return(-sum(sapply(si,function(x){(x/sn)*log(x/sn)}),na.rm = TRUE)) # return entropy!
    }
  }else if(mode == "sm"){ # max 3rd order entropy
    if(nrow(data)==2){
      return(1)
    }else if(nrow(data)==0 | nrow(data)==1){
      return(0)
    }else{
      nn <- length(unique(c(data[,1],data[,2])))
      if(directed){
        return(2*log(nn*(nn-1)))
      }else{
        return(2*log(nn*(nn-1)/2))
      }
    }
  }else if(mode == "t"){ # entropy of transition
    if(nrow(data)==2){
      return(0)
    }else if(nrow(data)==0 | nrow(data)==1){
      return(0)
    }else{
      x <- sapply(seq(1,length(data[,1])-1),
                  function(i){
                    if(data[i,2]==data[i+1,1] & data[i,1]!=data[i+1,2]){
                      return(c(data[i,1],data[i,2],data[i+1,2]))
                    }else{return(NA)}}) # find transitions
      if(!all(is.na(x))){ # no transitions at all
        t <- do.call("rbind",x)
        t <- as.data.frame(t[!is.na(t[,1]),])
        if(ncol(t)==1){
          return(-2)
        }
        tu <- unique(t) # unique transitions
        
        temp <- lapply(seq(1:length(tu[,1])),function(i){ # give id for each unique transition
          t[t[,1]==tu[i,1]&t[,2]==tu[i,2]&t[,3]==tu[i,3],4]<-i
          return(t[t[,1]==tu[i,1]&t[,2]==tu[i,2]&t[,3]==tu[i,3],])
        })
        temp <- abind(temp,along=1)
        ti <- table(temp[,4])
        tn <- sum(ti) # number of transitions
        return(-sum(sapply(ti,function(x){(x/tn)*log(x/tn)}),na.rm = TRUE))
        
      }else{
        return(-1)
      }
    }
  }else if(mode == "tm"){ # max entropy of transition
    if(nrow(data)==2){
      return(0)
    }else if(nrow(data)==0 | nrow(data)==1){
      return(0)
    }else{
      nn <- length(unique(c(data[,1],data[,2])))
      return(-log(1/(nn*(nn-1)*(nn-2))))
    }
  }else if(mode == "oboth"){ # 1st order entropy of both senders and receivers
    if(nrow(data)==1){
      return(1)
    }else if(nrow(data)==0){
      return(0)
    }else{
      nodes <- unique(c(data[,1],data[,2])) # unique nodes
      nn = table(nodes)
      nnsum = sum(nn) 
      return(-sum(sapply(nn,function(x){(x/nnsum)*log(x/nnsum)}),na.rm=T))
    }
  }else if(mode == "obothm"){ # max 1st order entropy of both senders and receivers
    if(nrow(data)==1){
      return(1)
    }else if(nrow(data)==0){
      return(0)
    }else{
      nn <- length(unique(c(data[,1],data[,2])))
      return(log(nn))
    }
  }else if(mode == "orec"){ # 1st order entropy of receivers
    if(nrow(data)==1){
      return(1)
    }else if(nrow(data)==0){
      return(0)
    }else{
      nodes <- unique(c(data[,1],data[,2]))
      nn <- length(nodes)
      entr <- 0
      for(i in 1:nn){
        n <- nodes[i]
        prob <- nrow(data[data[,2]==n,])/nrow(data)
        if(prob!=0) entr <- entr + prob*log(prob) 
      }
      return(-entr)
    }
  }else if(mode == "orecm"){ # max 1st order entropy of receivers
    if(nrow(data)==1){
      return(1)
    }else if(nrow(data)==0){
      return(0)
    }else{
      nn <- length(unique(c(data[,1],data[,2])))
      return(log(nn))
    }
  }else if(mode == "o"){ # 1st order entropy of senders
    if(nrow(data)==1){
      return(1)
    }else if(nrow(data)==0){
      return(0)
    }else{
      nodes <- unique(c(data[,1])) # unique nodes
      nn = table(nodes)
      nnsum = sum(nn) 
      return(-sum(sapply(nn,function(x){(x/nnsum)*log(x/nnsum)}),na.rm=T))
    }
  }else if(mode == "om"){ # max 1st order entropy of senders
    if(nrow(data)==1){
      return(1)
    }else if(nrow(data)==0){
      return(0)
    }else{
      nn <- length(unique(c(data[,1],data[,2])))
      return(log(nn))
    }
  }else if(mode == "nodes"){ # entropy of nodes interlocutors
    cNodes=unique(c(data$x,data$y))
    dfEntrs = c()
    for(node in cNodes){
      dfEntrs <- rbind(dfEntrs,c(node,tnEntropy(data,mode="snode",nodeID=node)))
    }
    dfEntrs=as.data.frame(dfEntrs)
    names(dfEntrs) <- c("node","entropy")
    return(dfEntrs)
  }else if(mode == "snode"){ # entropy of interlocutors for given node
    if(is.na(nodeID)){
      return(NA)
    }else{
      dfNodeData=data[data$x==nodeID | data$y==nodeID,]
      freq=as.numeric(table(paste0(dfNodeData$x,"_",dfNodeData$y)))
      probs=freq/sum(freq)
      return(-sum(probs*log(probs)))
    }
  }
}


# computing static measures for first n events of temporal network (incrementaly or not)
allStaticMeasures <- function(data){
  g=simplify(static_net(data))
  lStaticMeasures=staticNetMeasures(g)
}

# computing set of static measures for given graph
staticNetMeasures <- function(g){
  # plot <- ggnet2(g,node.size = 1)
  avg_degree = mean(as.numeric(igraph::degree(g,mode="total")))
  closeness = avg_closeness(g)
  diameter = diameter(g)
  density = igraph::edge_density(g)
  transitivity = transitivity(g,type="global")
  mean_distance = mean_distance(g)
  return(list(avg_degree=avg_degree,closeness=closeness,diameter=diameter,density=density,transitivity=transitivity,mean_distance=mean_distance))
}

# coputing average closeness of given graph
avg_closeness <- function(g){
  mean(igraph::closeness(g,mode="all"))
}

# computing average degree of temporal network (directed/indirected/multiple edges)
avg_degree <- function(data,mode="both",multiple_edges=T){
  if(multiple_edges){
    if(mode=="in"){
      mean(table(data[,2]))
    }else if(mode=="out"){
      mean(table(data[,1]))
    }else if(mode=="both"){
      mean(table(c(data[,1],data[,2])))
    }else{
      print("Unknown mode.")
    }
  }else{
    udata <- data[rownames(data) %in% rownames(unique(data[,1:2])),]
    if(mode=="in"){
      mean(table(udata[,2]))
    }else if(mode=="out"){
      mean(table(udata[,1]))
    }else if(mode=="both"){
      mean(table(c(udata[,1],udata[,2])))
    }else{
      print("Unknown mode.")
    }
  }
}

count_nodes <- function(data){
  return(length(unique(c(data$x,data$y))))
}

# converting temporal network into static network
static_net <- function(data){
  # add vertices
  g <- findAllNodes(data)
  # add edges
  g <- g + edges(as.character(as.vector(rbind(data[,1],data[,2]))))
  g
}

# converting temporal network into graph with equivalent vertices and no edges
findAllNodes <- function(data){
  g <- graph.empty()
  g <- g + vertices(as.character(unique(c(data[,1],data[,2]))))
  V(g)$infected <- rep_len(F,length.out = length(V(g)))
  g
}

# convert data frame columns into numeric
colsToNumeric <- function(df){
  for(i in 1:ncol(df)){
    df[,i] <- as.numeric(as.character(df[,i]))
  }
  df
}

# function for measuring execution time
stoper <- function(mode, sound=F, msg="",stopSound=4,lapSound=2, bLog=F, sFile="stoper_log.csv"){
  if(sound){require(beepr)}
  if(mode=="start"){
    if(!msg==""){print(msg)}
    timeStoper <<- Sys.time()
    timeLast <<- timeStoper
    nLap <<- 0
    if(bLog){
       dfTimes <<- c()
    }
  }else if(mode=="lap"){
    if(!msg==""){print(msg)}
    t <- Sys.time()
    nLap <<- nLap+1
    elapsedTime = as.numeric(difftime(t,timeLast,units = "secs"))
    print(paste(nLap,elapsedTime))
    timeLast <<- t
    if(bLog){
      dfTimes <<- rbind(dfTimes,c(elapsedTime,"lap",msg))
    }
    if(sound){beep(lapSound)}
  }else if(mode=="stop"){
    if(!msg==""){print(msg)}
    elapsedTime=as.numeric(difftime(Sys.time(),timeStoper,units = "secs"))
    if(bLog){
       dfTimes <<- rbind(dfTimes,c(elapsedTime,"stop",msg))
       write.csv2(x = dfTimes, file=sFile)
    }
    print(paste(nLap,elapsedTime))
    if(sound){beep(stopSound)}
    # garbage collector
    rm("timeLast",pos=".GlobalEnv")
    rm("timeStoper",pos=".GlobalEnv")
    rm("nLap",pos=".GlobalEnv")
    if(bLog){
      rm("dfTimes",pos=".GlobalEnv")
    }
  }else{
    print("wrong mode")
  }
}

# generate artificial temporal network (event equence)
generate_tn <- function(N=10,E=200,timeSpan=7,unit="days"){
  # N - number of nodes, E - number of events
  dfData <- c()
  # generate nodes id
  cN = floor(rnorm(N,500,100))  # could be done better
  while(!length(unique(cN))==N){
    cN = floor(rnorm(N,500,100))
  }
  # generate random events
  for(i in 1:E){
    sender = sample.int(N,size=1)
    receiver = sample.int(N,size=1)
    while(sender==receiver){
      receiver = sample.int(N,size=1)
    }
    dfData <- rbind(dfData,c(cN[sender],cN[receiver]))
  }
  firstDate = as.POSIXct("2000-01-01 00:00:00")
  if(unit=="secs"){
    lastDate = firstDate+timeSpan
  }else if(unit=="mins"){
    lastDate = firstDate+timeSpan*60
  }else if(unit=="hours"){
    lastDate = firstDate+timeSpan*60*60
  }else if(unit=="days"){
    lastDate = firstDate+timeSpan*60*60*24
  }
  times = sort(sample(as.numeric(firstDate):as.numeric(lastDate),E),decreasing = F)
  times[1] = firstDate
  times[length(times)] = lastDate
  times = as.POSIXct(times,origin="1970-01-01")
  
  
  dfData <- as.data.frame(dfData)
  dfData$EventTime = times
  names(dfData) <- c("SenderID","ReceiverID","EventTime")
  
  dfData
}

# generate artificial temporal network (event equence) with given distribution of nodes
genTN <- function(n,e, time.mean=10, time.sd=5, timespan=10, timespan_unit="year", sender_distr="norm", receiver_distr="norm", ...){
  if(is.null(timespan) | is.null(timespan_unit)){
    set.seed(as.numeric(Sys.time()))
    times = rnorm(n = e, mean=time.mean, sd=time.sd)
    times[times<0] <- 0
    times=cumsum(times)
  }else{
    times=floor(runif(n=e,min=0,max=timespan))
    times[1]=0
    times[length(times)]=timespan
    times=sort(times)
  }
  
  senders <- c()
  receivers <- c()
  if(sender_distr=="norm"){
    probs = dnorm(0:(n-1),mean=((n-1)/2),sd = ((n-1)/4))/sum(dnorm(0:(n-1),mean=((n-1)/2),
                                                                   sd = ((n-1)/4)))
  }else if(sender_distr=="exp"){
    probs = (dexp(0:(n-1),rate = 10/(n-1))/sum(dexp(0:(n-1),rate = 10/(n-1))))
  }else if(sender_distr=="unif"){
    probs = dunif(0:(n-1),min=0,max=n)
  }
  sdistr = probs
  if(receiver_distr=="norm"){
    probs = dnorm(0:(n-1),mean=((n-1)/2),sd = ((n-1)/4))/sum(dnorm(0:(n-1),mean=((n-1)/2),
                                                                   sd = ((n-1)/4)))
  }else if(receiver_distr=="exp"){
    probs = (dexp(0:(n-1),rate = 10/(n-1))/sum(dexp(0:(n-1),rate = 10/(n-1))))
  }else if(receiver_distr=="unif"){
    probs = dunif(0:(n-1),min=0,max=n)
  }
  rdistr = probs
  for(i in 1:e){
    send = sample(1:n,1,prob=sdistr)
    rec = sample(c(1:n)[-send],1,prob=rdistr[-send])
    senders <- c(senders,send)
    receivers <- c(receivers,rec)
  }
  data <- as.data.frame(cbind(x=senders,y=receivers,t=times),row.names=NA)
  return(data)
}
