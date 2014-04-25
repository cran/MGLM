

##===========================================##
## Make the pairwise scatter plot
##===========================================##

scatter.plot <- function(Y, main=NULL, free=TRUE, ...){
   "ID" <- "variable1"<-"value1" <- "value" <- NULL
  if(is.matrix(Y)){
    Y <- as.data.frame(Y)
    Y$ID <- rownames(Y)
  }else if(is.data.frame(Y)){
    Y$ID <- row.names(Y)
  }else{
    stop("Please provide a matrix or a data frame.")
  }
  meltDf <- melt(Y, id='ID')
  names(meltDf) <- c("ID", "variable1", "value1")
  out <- ddply(meltDf, .(variable1), function(df, df2){
    excl <- unique(df$variable)
    outDf <- merge(df, df2, by="ID")
  }, Y)
  meltDf2 <- melt(out, id=c("ID", "variable1", "value1"))
  meltDf2$value <- ifelse(meltDf2$variable1==meltDf2$variable, NA, 
                          meltDf2$value)
  meltDf2 <- subset(meltDf2, !is.na(value))
  if(!free){
  plot <- ggplot(meltDf2, aes(x=value1, y=value)) + 
    geom_point(...) + 
    facet_grid(variable1~variable,as.table=FALSE) +
    xlab("") + ylab("") + ggtitle(main) + 
    theme_bw()+
    theme(axis.text.x=element_text(size=6),
          axis.text.y=element_text(size=6))
  }else{
    plot <- ggplot(meltDf2, aes(x=value1, y=value)) + 
      geom_point(...) + 
      facet_wrap(variable1~variable,as.table=FALSE, scales="free") +
      xlab("") + ylab("") + ggtitle(main) + 
      theme_bw()+
      theme(axis.text.x=element_text(size=6),
            axis.text.y=element_text(size=6))
  }
  

  return(plot)
}

##===========================================##
## Make the pairwise correlation plot
##===========================================##

corr.plot <- function(Y, main=NULL){
  "Var1" <- "Var2" <- "value" <- NULL
  df <- melt(cor(Y))
  names(df) <- c(c("Var1", "Var2", "value"))
  plot <- qplot(x=Var1, y=Var2, data=subset(df, Var1!=Var2), 
    fill=factor(as.character(sign(value)), levels=c("-1", "1")), 
    alpha=value, geom="tile", drop=FALSE)+
    theme_bw()+
    scale_alpha(guide=FALSE)+
    scale_fill_discrete("", drop=FALSE) + 
    ggtitle(main)+
    xlab("") + ylab("")

  return(plot)
}


##===========================================##
## Make plot for 
##===========================================##

setMethod("plot", signature="MGLMreg", 
          function(x, facet=TRUE, free=TRUE, ...){
            
            colnames(x$fitted) <- colnames(x$data$Y)
            fitDf <- melt(x$fitted*rowSums(x$data$Y))
            oriDf <- melt(x$data$Y)
            names(fitDf)[3]<- "fitted"
            names(oriDf)[3]<- "Y"
            fitDf <- merge(oriDf,fitDf, by=c("Var1", "Var2"), all.x=FALSE)
			      fitDf <- subset(fitDf, !is.na(fitDf$fitted))
			      names(fitDf)[names(fitDf)=="Var2"] <- "Category"
            plot1 <- ggplot(fitDf, aes(x=Y, y=fitted)) + theme_bw()
            if(free) ff <- "free" else ff <- "fixed"
            if(facet==TRUE) {
              plot1 <- plot1+ geom_point(...) + facet_wrap(~Category, scales=ff)
            }else{
              plot1 <- plot1 + geom_point(aes(color=Category, shape=Category), ...)
            }
            return(plot1) 
})



setMethod("plot", signature="MGLMtune", 
          function(x, facet=FALSE, ...){
            
            plot1 <- ggplot(x$path, aes(x=log(Lambda),y=BIC)) + 
              geom_point(...) + geom_line()+
              geom_point(aes(x=log(Lambda[which.min(BIC)]), y=min(BIC)), color="red")+
              theme_bw()+xlab(expression(log(lambda)))
            
            return(plot1)
}
)

