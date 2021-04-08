library(shiny)
library(shinydashboard)
library(rhandsontable)
library(tidyverse)
library(plotly)
library(profvis)


raw_data <- read_csv("data/fGroups_filt_df.csv") %>%
  slice_head(n = 5000)


findhalo_rt <- function(mz,
                     intensity,
                     rt,
                     sf = 79/78.917789,
                     step = 0.001,
                     stepsd1=0.003,
                     stepsd2=0.005,
                     mzc=700,
                     cutoffint = 1000,
                     cutoffr=0.4,
                     clustercf =5){
  mzr <- round(mz)
  sm <- mz*sf
  sd <- ceiling(sm)-sm
  smsd <- ifelse(mz<=mzc,stepsd1,stepsd2)
  smstep <- seq(0,1,step)
  rt <- rt
  
  data <- cbind.data.frame(mz=mz,
                           mzr = mzr,
                           sm = sm,
                           sd =sd,
                           intensity=intensity,
                           rt = rt)
  data2 <<- data
  
  result <- NULL
  for(i in 1:length(smstep)){
    maxi = smstep[i]+smsd
    mini = smstep[i]-smsd
    index = sd<maxi & sd>mini
    
    li <- data[index & intensity > cutoffint,]
    mzt <- mzr[index & intensity > cutoffint]
    rtt <- rt[index & intensity > cutoffint]
    
    if(length(mzt)>=2){
      c <- cutree(hclust(dist(mzt)),h=clustercf)
      t <- cutree(hclust(dist(rtt)), h = clustercf)
      u <- paste0(c, t)
      cn <- length(unique(u))
      lit <- cbind.data.frame(li,u,i)
      for (j in 1:cn){
        li2 <- lit[lit[,7]==j,]
        mzt2 <- lit$mzr[lit[,7]==j]
        
        if(length(mzt2)>=2){
          if(length(unique(li2$intensity))>1){
            ratio <- max(li2$intensity[li2$intensity != max(li2$intensity)]) / max(li2$intensity)
            diff <- abs(li2$mzr[round(li2$intensity) == round(max(li2$intensity[li2$intensity != max(li2$intensity)]))] - li2$mzr[which.max(li2$intensity)])
          }else{
            ratio <- 1
            diff <- abs(li2$mzr[1]-li2$mzr[2])
          }
          
          if(ratio>cutoffr&round(diff)==2){
            result <- rbind.data.frame(result,li2)
          }
        }
      }
    }
  }
  return(result[!duplicated(result$mz), ])
}





profvis(
t <- findhalo_rt(mz = raw_data$mz,
                 intensity = raw_data$intensity,
                 rt = raw_data$rt,
                 sf = 79/78.917789,
                 step = 0.001,
                 stepsd1=0.003,
                 stepsd2=0.005,
                 mzc=700,
                 cutoffint = 1000,
                 cutoffr=0.4,
                 clustercf =5
                 )
)



table1 <- renderRHandsontable({
      rhandsontable(result(), width = 700, height = 400) %>%
        hot_cols(fixedColumnsLeft = 1) %>%
        hot_rows(fixedRowsTop = 1)
      
    })
    
    ###Plot###
    output$plot1 <- renderPlotly({
      
      
      if(input$check_rt){
        df <- result()
        
        if (input$check_int) {
          p <- df %>%
            plot_ly(
              x = df$sm,
              y = df$sd,
              type = "scatter",
              size = df$intensity,
              color = df$i,
              hovertext = paste("RT: ", df$rt)
            )
        } else {
          p <- df %>%
            plot_ly(
              x = df$sm,
              y = df$sd,
              type = "scatter",
              color = df$i,
              hovertext = paste("RT: ", df$rt))
        }
        
        if (input$check_overlay) {
          p <- p %>%
            add_markers(x = ~data2$sm, y = ~data2$sd, 
                        color = I("grey"), opacity = 0.2, hovertext = paste("RT: ", data2$rt))
        } else {
          p <- p
        }
        
        if (input$check_overlay & input$check_int) {
          p <- df %>%
            plot_ly(
              x = df$sm,
              y = df$sd,
              type = "scatter",
              size = df$intensity,
              color = df$i,
              hovertext = paste("RT: ", df$rt)
            ) %>%
            add_markers(x = ~data2$sm, y = ~data2$sd, 
                        color = I("grey"), opacity = 0.2, size = ~data2$intensity, hovertext = paste("RT: ", data2$rt))
        } else {
          p <- p
        }
        
        
      }else{
        
        df <- result()
        
        if (input$check_int) {
          p <- df %>%
            plot_ly(
              x = df$sm,
              y = df$sd,
              type = "scatter",
              size = df$intensity,
              color = df$i)
        } else {
          p <- df %>%
            plot_ly(
              x = df$sm,
              y = df$sd,
              type = "scatter",
              color = df$i)
        }
        
        if (input$check_overlay) {
          p <- p %>%
            add_markers(x = ~data2$sm, y = ~data2$sd, 
                        color = I("grey"), opacity = 0.2)
        } else {
          p <- p
        }
        
        if (input$check_overlay & input$check_int) {
          p <- df %>%
            plot_ly(
              x = df$sm,
              y = df$sd,
              type = "scatter",
              size = df$intensity,
              color = df$i
            ) %>%
            add_markers(x = ~data2$sm, y = ~data2$sd, 
                        color = I("grey"), opacity = 0.2, size = ~data2$intensity)
        } else {
          p <- p
        }
        
        
      }
      
      
      
      
    })
    
  
    