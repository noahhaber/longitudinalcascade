#' This package generates a longitudinal casade, including a graphical representation. This takes a long-formatted list of stage-by-stage events and transforms it into a longitudinal cascade, correcting the orders of events.
#' @title Longitudinal cascade generator
#' @keywords cascade longitudinal survival
#' @param df (required) The main dataframe input parameter. The data frame needs at least the following fields:
#' "ID": (required) Either a numerical or string-based individual identifier, indicating every person in the dataset
#' "date": (required) Date-formatted date son which the event / stage occurred
#' "stage": (required) String indicating the stage achieved by the individual on the specified date. Stages must match the string in the stages.order parameter. Additonal events may be included in the "stage" category, including death, loss to follow up, and interstage events defined in the other parameters.
#' "group": (optional) Strings indicating any relevent groups of data.
#' @param stages.order (required) stages.order is the parameter which defines the events to be considered in the main cascade and their order. This is a vector of strings matching items in the "Stage" column of the main data frame, e.g. c("Stage 1","Stage 2","Stage 3"). 
#' @param groups.order (optional) This is a vector of groups, matching the "group" column of the main data frame
#' @import survival
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @import zoo
#' @import scales
#' @export
#' @examples
#' long.cascade()
long.cascade <- function(events.long,stages.order,groups.order=NA,
                               death.indicator=NA,retention.indicator=NA,censorship.indicator=NA,
                               allow.sub.lines=FALSE,allow.skips=FALSE,
                               x.axis.max=365) {
  # Initialize function
  {
    # Packages
    {
      library(survival)
      #library(survsim)
      library(ggplot2)
      library(tidyr)
      library(dplyr)
      library(zoo)
      library(scales)
    }
    # Settings values
      event.long.orig <- event.long
      x.axis.range <- c(0,x.axis.max)
    # Generate a single "group" if group is unchanged
      if (anyNA(groups.order)) {
        event.long$group <- "All observations"
      }
      else {}
  }
  # Data manipulation
  {
    # Transform data into wide form
    {
      # Generate stage names
        stages <- as.data.frame(matrix(c(stages.order,paste(rep("stage"),as.character(seq(1:length(stages.order))),sep="."),1:length(stages.order)),ncol=3),stringsAsFactors=FALSE)
        colnames(stages) <- c("stage","stage.number","stage.index")
      # Replace stages with an index for stages
        event.long <- merge(event.long,stages,by="stage")
        event.long <- subset( event.long, select = -stage )
        event.long <- event.long[order(event.long$stage.number),]
        event.long$stage.index <- NULL
      # Generate separated wide list of stage events
        events_wide <- event.long %>%
          spread(stage.number, date)
      # Replace names with "date"
        names(events_wide) <- gsub("stage.", "date.stage.",names(events_wide))
    }
    # Fix ordering, so that completion of a later stage indicates completion of earlier stage, and mark where this happens
    {
      for (stage.index in seq(length(stages.order)-1,1,-1)){
        # Determine origin of stage events, and create variable for its origin
        former.date.name <- paste("date.stage.",stage.index,sep="")
        latter.date.name <- paste("date.stage.",stage.index+1,sep="")
        former.date.origin.index.name <- paste("stage.",stage.index,".origin.index",sep="")
        latter.date.origin.index.name <- paste("stage.",stage.index+1,".origin.index",sep="")
        
        # Determine whether or not current stage is eligible for change
          events_wide$temp_flag <- ifelse(events_wide[[latter.date.name]] < events_wide[[former.date.name]] | 
                                          is.na(events_wide[[former.date.name]])==TRUE,1,0)
          events_wide$temp_flag <- ifelse(is.na(events_wide$temp_flag)==TRUE,0,events_wide$temp_flag)
        # Check if this is the last event (determines whether there is an available following stage index)
          if (stage.index==(length(stages.order)-1)){
            # Generate on origin flag index, which = last stage if there is a change, and the current stage if there is none
            events_wide[[former.date.origin.index.name]] <- ifelse(
              events_wide$temp_flag == 1,
              stage.index + 1, stage.index
            )
            events_wide[[former.date.origin.index.name]] <- ifelse(
              is.na(events_wide[[former.date.name]])==TRUE & is.na(events_wide[[latter.date.name]])==TRUE,
              NA,events_wide[[former.date.origin.index.name]])
          } else {
          # Generate on origin flag index, which = the corresponding flag for the next stage if there is a change, and the current stage if there is none
            events_wide[[former.date.origin.index.name]] <- ifelse(
              events_wide$temp_flag == 1 & !is.na(events_wide[[latter.date.origin.index.name]]),
              events_wide[[latter.date.origin.index.name]],NA
            )
            events_wide[[former.date.origin.index.name]] <- ifelse(
              events_wide$temp_flag == 0,
              stage.index,events_wide[[former.date.origin.index.name]]
            )
        }
        # Replace value of date if a change is made
          events_wide[[paste("date.stage.",stage.index,sep="")]] <- ifelse(
            events_wide[[paste("stage.",stage.index,".origin.index",sep="")]] == stage.index,
            events_wide[[paste("date.stage.",stage.index,sep="")]],events_wide[[paste("date.stage.",stage.index +1,sep="")]]
          )
          events_wide[[paste("date.stage.",stage.index,sep="")]] <- as.Date.numeric(events_wide[[paste("date.stage.",stage.index,sep="")]])
      }
      events_wide$temp_flag <- NULL
    }
    # Fill in missing events backward through the cascade, so that subsequent events mark previous events if otherwise missing
    {
      for (i in seq(length(stages.order)-1,1,-1)){
        events_wide[[paste("date.stage.",i,sep="")]] <- ifelse(is.na(events_wide[[paste("date.stage.",i,sep="")]])==TRUE & is.na(events_wide[[paste("date.stage.",(i+1),sep="")]])==FALSE,
                                                               events_wide[[paste("date.stage.",(i+1),sep="")]],
                                                               events_wide[[paste("date.stage.",i,sep="")]])
      }
    }
    # Add a dummy censorship date for any missing dates which is well after the maximum timeline
    {
      for (i in 1:(length(stages.order))){
        events_wide[[paste("date.stage.",i,sep="")]] <- ifelse(is.na(events_wide[[paste("date.stage.",i,sep="")]]),
                                                               events_wide$date.stage.1+length(stages.order)*(x.axis.max+1),
                                                               events_wide[[paste("date.stage.",i,sep="")]])
      }
    }
    # Determine dates of death and censorship, and add to dataset
    {
      # Merge in death events (if any)
        if (is.na(death.indicator)==FALSE){
          events.death <- subset(event.long.orig,stage==death.indicator)
          events.death <- events.death[c("ID","date")]
          colnames(events.death) <- c("ID","date.death")
          events_wide <- merge(events_wide,events.death,by="ID",all.x = TRUE)
        }
      # Merge in censorship events (if any)
        if (is.na(censorship.indicator)==FALSE){
          events.censorship <- subset(event.long.orig,stage==censorship.indicator)
          events.censorship <- events.censorship[c("ID","date")]
          colnames(events.censorship) <- c("ID","date.censorship")
          events_wide <- merge(events_wide,events.censorship,by="ID",all.x = TRUE)
        }
      # Assume that if date of death exists and occurs after censorship event, censorship is false and replace with NA
        if ((is.na(censorship.indicator)==FALSE) & (is.na(death.indicator)==FALSE)){
          events_wide$date.censorship <- ifelse((events_wide$date.censorship < events_wide$date.death) & (is.na(events_wide$date.death)==FALSE) & (is.na(events_wide$date.censorship)==FALSE),
                                                NA,events_wide$date.censorship)
          events_wide$date.censorship <- as.Date.numeric(events_wide$date.censorship)
        }
      # Make dummy date for missing data a censorship data
        for (i in 1:(length(stages.order))){
          events_wide$date.censorship <- ifelse(events_wide[[paste("date.stage.",i,sep="")]]==events_wide$date.stage.1+length(stages.order)*(x.axis.max+1),
                                                events_wide$date.stage.1+length(stages.order)*(x.axis.max+1),
                                                events_wide$date.censorship)
        }
      # Generate time from each stage to death and censorship
        for (i in 1:(length(stages.order)-1)){
          if (is.na(death.indicator)==FALSE){
            events_wide[[paste("time.stage.",i,".to.death",sep="")]] <- 
              events_wide[["date.death"]] - events_wide[[paste("date.stage.",i,sep="")]]
          }
          if (is.na(censorship.indicator)==FALSE){
            events_wide[[paste("time.stage.",i,".to.censorship",sep="")]] <- 
              events_wide[["date.censorship"]] - events_wide[[paste("date.stage.",i,sep="")]]
          }
        }
    }
    # Generate times between all possible events
    {
      for (i in 1:(length(stages.order)-1)){
        for (j in (i+1):(length(stages.order))){
          # Generate time to actual stage event if event occurs
            events_wide[[paste("time.stage.",i,".to.",j,sep="")]] <-
              events_wide[[paste("date.stage.",j,sep="")]] - events_wide[[paste("date.stage.",i,sep="")]]
          # Replace with time to censorship if event never happens and generate a censorship indicator flag
            if (is.na(censorship.indicator)==FALSE){
              events_wide[[paste("censorship.stage.",i,".to.",j,sep="")]] <- ifelse(
                is.na(events_wide[[paste("time.stage.",i,".to.",j,sep="")]])==TRUE & is.na(events_wide[[paste("time.stage.",i,".to.censorship",sep="")]])==FALSE,
                1,0
              )
              events_wide[[paste("time.stage.",i,".to.",j,sep="")]] <- ifelse(
                is.na(events_wide[[paste("time.stage.",i,".to.",j,sep="")]])==TRUE,
                events_wide[[paste("time.stage.",i,".to.censorship",sep="")]],events_wide[[paste("time.stage.",i,".to.",j,sep="")]]
              )
            } else{
              events_wide[[paste("censorship.stage.",i,".to.",j,sep="")]] <- 0
            }
            events_wide[[paste("censorship.stage.",i,".to.",j,sep="")]] <- ifelse(
              is.na(events_wide[[paste("time.stage.",i,".to.",j,sep="")]])==TRUE,
              NA,events_wide[[paste("censorship.stage.",i,".to.",j,sep="")]]
            )
          # Generate transition-specific mortality (i.e. mortality occuring only between each stage)
            if (is.na(death.indicator)==FALSE){
              events_wide[[paste("time.death.stage.",i,".to.",j,sep="")]] <- ifelse(
                (is.na(events_wide[[paste("time.stage.",i,".to.death",sep="")]]) == FALSE) & (is.na(events_wide[[paste("time.stage.",i,".to.",j,sep="")]])==TRUE),
                events_wide[[paste("time.stage.",i,".to.death",sep="")]],NA
              )
              events_wide[[paste("time.death.stage.",i,".to.",j,sep="")]] <- ifelse(
                (is.na(events_wide[[paste("date.stage.",i,sep="")]]) == FALSE) & (is.na(events_wide[[paste("time.death.stage.",i,".to.",j,sep="")]])==TRUE),
                x.axis.max,events_wide[[paste("time.death.stage.",i,".to.",j,sep="")]]
              )
              events_wide[[paste("time.death.stage.",i,".to.",j,sep="")]] <- ifelse(
                events_wide[[paste("time.death.stage.",i,".to.",j,sep="")]] > x.axis.max,
                x.axis.max,events_wide[[paste("time.death.stage.",i,".to.",j,sep="")]]
              )
              events_wide[[paste("censorship.death.stage.",i,".to.",j,sep="")]] <- ifelse(
                events_wide[[paste("time.death.stage.",i,".to.",j,sep="")]] >= x.axis.max,
                1,0
              )
            }
          # Replace transition times where initializing event occurs but finishing does not, with censorship after x axis limit
            events_wide[[paste("time.stage.",i,".to.",j,sep="")]] <- ifelse(
              is.na(events_wide[[paste("time.stage.",i,".to.",j,sep="")]]) == TRUE & is.na(events_wide[[paste("date.stage.",i,sep="")]]) == FALSE,
              x.axis.max, events_wide[[paste("time.stage.",i,".to.",j,sep="")]]
            )
          # Replace transition times greater than axis limit, with censorship after x axis limit
            events_wide[[paste("time.stage.",i,".to.",j,sep="")]] <- ifelse(
              events_wide[[paste("time.stage.",i,".to.",j,sep="")]] > x.axis.max,
              x.axis.max, events_wide[[paste("time.stage.",i,".to.",j,sep="")]]
            )
          # Change to censorship event if censored after axis limit
            events_wide[[paste("censorship.stage.",i,".to.",j,sep="")]] <- ifelse(
              events_wide[[paste("time.stage.",i,".to.",j,sep="")]] > x.axis.max,
              1, events_wide[[paste("censorship.stage.",i,".to.",j,sep="")]]
            )
        }
      }
    }
    # Generate full survival dataset
    {
      # Main events
      {
        # Function to generate survival data for a single curve
        gen.survival.data <- function(start.stage.index,end.stage.index,group.index){
          
          # Select data
            if (allow.skips==TRUE){
              events_wide$selection <- is.na(events_wide[[paste("time.stage.",start.stage.index,".to.",end.stage.index,sep="")]])==FALSE &
                events_wide$group == groups.order[group.index]
            } else {
              events_wide$selection <- is.na(events_wide[[paste("time.stage.",start.stage.index,".to.",end.stage.index,sep="")]])==FALSE &
                events_wide[[paste("time.stage.",start.stage.index,".to.",end.stage.index,sep="")]]!=0 &
                events_wide$group == groups.order[group.index]
            }
          # Remove those who are disqualified from the "0+1" stage events from "0+n" events
            events_wide$selection <- ifelse(events_wide[[paste("time.stage.",start.stage.index,".to.",(start.stage.index+1),sep="")]]==0,
                                            FALSE,
                                            events_wide$selection)
          
          # Generate list of events
            chart.time <- as.integer(events_wide[[paste("time.stage.",start.stage.index,".to.",end.stage.index,sep="")]][events_wide$selection==TRUE])
            chart.event <- as.integer(!events_wide[[paste("censorship.stage.",start.stage.index,".to.",end.stage.index,sep="")]][events_wide$selection==TRUE])
          # Add dummy event to allow graphical continuation of data past last known event
            chart.time <- c(chart.time, x.axis.max)
            chart.event <- c(chart.event,1)
            chart.time <- c(chart.time, x.axis.max + 1)
            chart.event <- c(chart.event,1)
          # Generate survival data
            chart.surv <- survfit(Surv(time = chart.time, event = chart.event) ~ 1)
            surv.time <- chart.surv$time
            surv.surv <- 1-chart.surv$surv
            surv.data <- data.frame(surv.time,surv.surv)
            surv.data$start.stage.index <- start.stage.index
            surv.data$end.stage.index <- end.stage.index
            surv.data$group.index <- group.index
          # Return
            return(surv.data)
        }
        # Generate wrapper functions
          wrapper.sub.curves <- function(start.stage.index,group.index){
            surv.combined <- do.call("rbind",
                                     lapply((start.stage.index+1):(length(stages.order)),function(x) gen.survival.data(start.stage.index,x,group.index)))
            return(surv.combined)
          }
          wrapper.transitions <- function(group.index){
            surv.combined <- do.call("rbind",
                                     lapply(1:(length(stages.order)-1),function(x) wrapper.sub.curves(x,group.index)))
            return(surv.combined)
          }
        # Call wrapped functions to generate full long dataset of survival events
          surv.combined <- do.call("rbind",
                                   lapply(1:length(groups.order),function(x) wrapper.transitions(x)))
        # Change stages and groups to factors for easier charting
          surv.combined$start.stage.factor <- ordered(surv.combined$start.stage.index,levels=c(1:(length(stages.order)-1)),
                                                      labels = stages.order[1:(length(stages.order)-1)])
          surv.combined$end.stage.factor <- ordered(surv.combined$end.stage.index,levels=c(2:(length(stages.order))),
                                                    labels = stages.order[2:(length(stages.order))])
          surv.combined$group.factor <- ordered(surv.combined$group.index,levels=c(1:(length(groups.order))),
                                                labels = groups.order[1:(length(groups.order))])
        # Revise data based on options
          if (allow.sub.lines == FALSE) {
            surv.combined <- subset(surv.combined,surv.combined$start.stage.index==(surv.combined$end.stage.index)-1)
          }
      }
      # Death events
        if (is.na(death.indicator)==FALSE){
          # Function to generate survival data for a single curve
            gen.death.data <- function(start.stage.index,end.stage.index,group.index){
            
              # Select data
                events_wide$selection <- is.na(events_wide[[paste("time.death.stage.",start.stage.index,".to.",end.stage.index,sep="")]])==FALSE &
                  events_wide$group == groups.order[group.index]
              # Generate list of events
                chart.time <- as.integer(events_wide[[paste("time.death.stage.",start.stage.index,".to.",end.stage.index,sep="")]][events_wide$selection==TRUE])
                chart.event <- as.integer(!events_wide[[paste("censorship.death.stage.",start.stage.index,".to.",end.stage.index,sep="")]][events_wide$selection==TRUE])
              # Add dummy event to allow graphical continuation of data past last known event
                chart.time <- c(chart.time, (x.axis.max+1))
                chart.event <- c(chart.event,1)
              # Generate survival data
                chart.surv <- survfit(Surv(time = chart.time, event = chart.event) ~ 1)
                surv.time <- chart.surv$time
                surv.surv <- 1-chart.surv$surv
                surv.data <- data.frame(surv.time,surv.surv)
                surv.data$start.stage.index <- start.stage.index
                surv.data$end.stage.index <- end.stage.index
                surv.data$group.index <- group.index
              # Return
              return(surv.data)
          }
          # Generate wrapper functions
            wrapper.death.transitions <- function(group.index){
              surv.death.combined <- do.call("rbind",
                                             lapply(1:(length(stages.order)-1),function(x) gen.death.data(x,x+1,group.index)))
              return(surv.death.combined)
            }
          # Call wrapped functions to generate full long dataset of survival events
            surv.death.combined <- do.call("rbind",
                                           lapply(1:length(groups.order),function(x) wrapper.death.transitions(x)))
          # Change stages and groups to factors for easier charting
            surv.death.combined$start.stage.factor <- ordered(surv.death.combined$start.stage.index,levels=c(1:(length(stages.order)-1)),
                                                              labels = stages.order[1:(length(stages.order)-1)])
            surv.death.combined$end.stage.factor <- ordered(surv.death.combined$end.stage.index,levels=c(2:(length(stages.order))),
                                                            labels = stages.order[2:(length(stages.order))])
            surv.death.combined$group.factor <- ordered(surv.death.combined$group.index,levels=c(1:(length(groups.order))),
                                                        labels = groups.order[1:(length(groups.order))])
        }
    }
  }
  # Generate chart graphics
  {
    # Generate transition charts
    # Color gradient creator
      color.gradient <- function(original.color,number.divisions,darken=0){
        # Generate hsv color
          hsv.color.orig <- rgb2hsv(col2rgb(original.color))
        # Manipulate if set to lighten
          if (darken == 0) {
            hsv.color.new <- t(matrix(c(rep(hsv.color.orig[1,1],number.divisions),
                                        hsv.color.orig[2,1]*seq(1,1/number.divisions,-1/number.divisions),
                                        (1-hsv.color.orig[3,1])*seq(0,(number.divisions-1)/(number.divisions),1/(number.divisions))+hsv.color.orig[3,1]),
                                      nrow=number.divisions))
          } else {
            hsv.color.new <- t(matrix(c(rep(hsv.color.orig[1,1],number.divisions),
                                        (1-hsv.color.orig[2,1])*seq(0,(number.divisions-1)/(number.divisions),1/(number.divisions))+hsv.color.orig[2,1],
                                        hsv.color.orig[3,1]*seq(1,1/number.divisions,-1/number.divisions)),
                                      nrow=number.divisions))
            
          }
          rownames(hsv.color.new) <- c("h","s","v")
        # Convert back to regular color
          return(unname(sapply(data.frame(hsv.color.new), function(x) do.call(hsv, as.list(x)))))
      }
    # Generate colors
      main.line.colors <- color.gradient("#4472C4",(length(stages.order)-1))
    # Stacked chart
      # Temporary for putting in years
        surv.combined.chart <- surv.combined
        surv.combined.chart$surv.time = surv.combined.chart$surv.time/365
        
        chart <- ggplot(data=surv.combined.chart) +
          geom_rect(aes(xmin=surv.time,xmax=lead(surv.time),ymin=0,ymax=surv.surv,fill=end.stage.factor),alpha=1) +
          scale_fill_manual(values=main.line.colors) +
          geom_step(aes(x=surv.time,y=surv.surv,color=end.stage.factor)) +
          theme_bw() %+replace%
          theme(
            panel.grid = element_blank(),
            plot.margin = unit(c(.1,.1,.1,.1), "cm"),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            legend.position="bottom",
            legend.title=element_blank(),
            axis.text = element_text(colour="black",size=10),
            strip.text = element_text(size = 12,hjust=0),
            panel.spacing = unit(1, "lines")
          ) +
          scale_x_continuous(limits = (x.axis.range/365),expand = c(0, 0)) +
          scale_y_continuous(limits = c(0, 1),expand = c(0, 0),labels=percent) +
          scale_color_manual(values=c(rep("black",length(stages.order)-1))) +
          facet_grid(group.factor ~ start.stage.factor,
                     switch="y") +
          theme(strip.background = element_blank(),
                strip.placement = "outside")
      # Add death event if present
        if (is.na(death.indicator)==FALSE){
          chart <- chart +
            geom_step(data=surv.death.combined,aes(x=surv.time,y=1-surv.surv)) + 
            geom_rect(data=surv.death.combined,aes(xmin=surv.time,xmax=lead(surv.time),ymin=1-surv.surv,ymax=1),alpha=1,fill="indianred1")
        }
  }
  # Generate between-group statistical tests
  {
    # Test for heterogeneity of survival curves for main between-stage events between groups
    if (length(groups.order)>1){
      # Generate tests function
        generate.surv.diff <- function(init.stage){
          chart.time <- as.integer(events_wide[[paste("time.stage.",(init.stage),".to.",(init.stage+1),sep="")]])
          chart.event <- as.integer(!events_wide[[paste("censorship.stage.",(init.stage),".to.",(init.stage+1),sep="")]])
          groups <- ordered(events_wide$group,labels=groups.order,levels=groups.order)
          test.name <- paste0("curve.het.stage.",(init.stage),"to",(init.stage+1))
          surv.diff <- list(survdiff(Surv(time = chart.time, event = chart.event) ~ groups))
          names(surv.diff)[1] <- test.name
          surv.diff
        }
      # Package tests
        surv.diffs.combined <- do.call("list",
                                 lapply(1:(length(stages.order)-1),function(x) generate.surv.diff(x)))
        
    }
    
  }
  # Prepare export data
  {
    if (length(groups.order)>1){
      output.df <- list("chart" = chart,"surv.dataset" = surv.combined,"events_wide" = events_wide,"surv.diffs" = surv.diffs.combined)
    } else {
      output.df <- list("chart" = chart,"surv.dataset" = surv.combined,"events_wide" = events_wide)
    }
    
  }
}