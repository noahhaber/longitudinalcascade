#' This package generates a longitudinal casade, including a graphical representation. This takes a long-formatted list of stage-by-stage events and transforms it into a longitudinal cascade, correcting the orders of events.
#' @title Longitudinal cascade statistics and charts
#' @name longitudinalcascade
#' @keywords cascade longitudinal survival
#' @param events.long (required) The main dataframe input parameter. The data frame needs at least the following fields:
#' "ID": (required) A string-based individual identifier, indicating every person in the dataset.
#' "date": (required) Date-formatted date on which the event / stage occurred
#' "stage": (required) String indicating the stage achieved by the individual on the specified date. Stages must match the string in the stages.order parameter. Additonal events may be included in the "stage" category, including death, loss to follow up, and interstage events defined in the other parameters.
#' "group": (optional) String indicating any relevent groups of data.
#' @param stages.order (required) stages.order is the parameter which defines the events to be considered in the main cascade and their order. This is a vector of strings matching items in the "Stage" column of the main data frame, e.g. c("Stage 1","Stage 2","Stage 3").
#' @param groups.order (optional) This is a vector of groups, matching the "group" column of the main data frame. If left blank, no group comparisons will be performed. For the chart, each group will have its own row.
#' @param groups.date.breaks (optional) If groups.date.breaks is filled in, the grouping will be defined by date range of entry event for each transition, rather than groups of individuals. Each transition will independently determine its own groups, based on the time in which the entrance event occurs. Times are determined by a vector of date breaks. Each group is defined as starting from a given date break value and continuing until it reaches the subsequent date break, not including data from that ending break value. For example, setting the break values to be January 1, 2011, January 1, 2012, and January 1, 2013 will create two groups. The first group will take individuals who entered each stage from January 1, 2011 to Dec 31, 2011, and the second will take individuals who entered into the stage from January 1, 2012 to Dec 31, 2012.
#' @param death.indicator (optional) This parameter is the string which indicates a death event in the dataset. If specified, between-stage mortality will be estimated and shown as a KM curve on the top of the chart(s). If left blank, death events will not be estimated.
#' @param censorship.indicator (optional) This parameter is the string which indicates a right-censorship event. Most commonly, this will indicate permanent loss to follow up and/or end of data collection.
#' @param allow.sub.lines Sub-lines indicate subsequent transitions across the cascade. If TRUE, the main chart will show transitions to all possible subsequent events. For example, if there are 4 stages (1-4), the leftmost chart will show each transition from 1-2, 1-3, and 1-4, while the next chart will show 2-3 and 2-4, and the last chart will show only 3-4. If FALSE, the charts will only show transition to the subsequent stage.
#' @param allow.skips (To be replaced with additional options) This option shows "skips" across the cascade in each chart, as indicated by the y intercept. If FALSE, each stage will start only with people who have not moved on to a subsequent stage, i.e. the y intercept will always be 0. If TRUE, an individual can enter into a stage even if they have "skipped" it. For example, an individual may go straight from stage 1 to stage 3, skipping 2. If this indicator is FALSE, the stage transition chart from 2-3 will not contain this individual in the denomenator. If TRUE, this individual will be counted in the denomenator for this transition, but will be counted as having transitioned into stage 3 immediately upon entering stage 2.
#' @param x.axis.max This option shows the maximum range of the x axis in days. Defaults to 365 days (1 year).
#' @param nochart Setting this to TRUE prevents the function from generating the main chart.
#' @param risk.pool.size.line Setting to TRUE adds an indicator of risk pool remaining to the main charts as a line reflected beneath the main chart, showing the proportion of the original risk pool remaining at each time point, starting from 100%. Defaults to FALSE.
#' @param main.fill.colors (optional) This defines the color scheme of the stage transition graphs, as a string indicator for color or a c() list of colors. If the colors contain only one color, the color scheme will automatically generate progressively faded versions of the initial color provided for the remaining stage transitions. Otherwise, a list which is exactly one fewer than the # of stages must be provided, in the order of stage trasitions.
#' @param death.fill.color (optional) This defines the color scheme for the death stage transition, as a string indicator for color.
#' @param risk.pool.fill.color (optional) This defines the color scheme for the risk pool graphic, as a string indicator for color.
#' @import survival ggplot2 dplyr tidyr zoo scales grDevices
#' @importFrom stats relevel
#' @importFrom rlang .data
#' @export
#' @references Haber et al. (2017) Lancet HIV 4(5):e223-e230
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/28153470}{PubMed})
#' @examples
#' # Pull in data from example simulated dataset
#' library(longitudinalcascade)
#' data(events_long)
#'
#' # Set up options
#' stages.order <- c("First tested positive","Knows status","Linked to care","Eligible for ART",
#' "Initiated ART","Therapeutic response")
#' groups.order <- c("Group 1","Group 2","Group 3")
#' death.indicator <- "Death"
#' retention.indicator <- "Clinic visit"
#' censorship.indicator <- "LTFU"
#' allow.sub.lines <- TRUE
#' allow.skips <- TRUE
#'
#' # Create cascade object
#' longitudinalcascade.sim <- longitudinalcascade(events_long,stages.order=stages.order,
#' groups.order=groups.order,death.indicator=death.indicator,
#' censorship.indicator=censorship.indicator,
#' allow.sub.lines=allow.sub.lines,allow.skips=allow.skips)
#'
#' # Print/output main chart
#' longitudinalcascade.sim$chart
#' # Output full survival dataset generated, as a data frame
#' df.longitudainalcascade.survival <- longitudinalcascade.sim$surv.dataset
#' # Output heterogeneity test
#' longitudinalcascade.sim$surv.diffs
#' # Output original long-formatted list of events
#' df.events.long <- longitudinalcascade.sim$events.long
#' # Output generated wide-formatted list of events
#' df.events.wide <- longitudinalcascade.sim$events.wide

longitudinalcascade <- function(events.long,stages.order,groups.order=NA,
                         death.indicator=NA,censorship.indicator=NA,
                         allow.sub.lines=FALSE,allow.skips=FALSE,
                         groups.date.breaks=NA,
                         x.axis.max=365,
                         main.fill.colors = "#4472C4",death.fill.color = "#FF6A6A",
                         nochart=FALSE,risk.pool.size.line=FALSE,
                         risk.pool.fill.color = "#90dbb2") {
  # Data manipulation
  {
    # Transform data into wide form
    {
      # Preserve data
        events.long.orig <- events.long
      # Generate stage names
        stages <- as.data.frame(matrix(c(stages.order,paste(rep("stage"),as.character(seq(1:length(stages.order))),sep="."),1:length(stages.order)),ncol=3),stringsAsFactors=FALSE)
        colnames(stages) <- c("stage","stage.number","stage.index")
      # Replace stages with an index for stages
        events.long <- merge(events.long,stages,by="stage")
        events.long <- events.long[,!(names(events.long) %in% c("stage"))]
        #events.long <- subset( events.long, select = -stage )
        events.long <- events.long[order(events.long$stage.number),]
        events.long$stage.index <- NULL
      # Generate a single "group" if groups are not specified
        if (anyNA(groups.order)) {
          events.long$group <- "All observations"
          groups.order <- c("All observations")
        } else{}
      # Generate separated wide list of stage events
        events.wide <- events.long %>%
          tidyr::spread(.data$stage.number, date)
      # Replace names with "date"
        names(events.wide) <- gsub("stage.", "date.stage.",names(events.wide))
    }
    # Generate time-based groups if time breaks are specified
    if (anyNA(groups.date.breaks)==FALSE){
      # Generate group names from breaks
        groups.order <- c(paste0(as.character(groups.date.breaks[1])," to ",as.character(groups.date.breaks[2]-1)))
        for (i in 2:(length(groups.date.breaks)-1)){
          groups.order <- c(groups.order,paste0(as.character(groups.date.breaks[i])," to ",as.character(groups.date.breaks[i+1]-1)))
        }
      # Generate grouping for each event except the last to determine whether the start event is within time breaks
        for (stage.index in 1:(length(stages.order)-1)){
          events.wide[[paste0("date.stage.",stage.index,".group")]] <- NA
          for (break.index in 1:(length(groups.date.breaks)-1)){
            events.wide[[paste0("date.stage.",stage.index,".group")]] <- ifelse(
              (events.wide[[paste0("date.stage.",stage.index)]] >= groups.date.breaks[break.index]) & (events.wide[[paste0("date.stage.",stage.index)]] < groups.date.breaks[break.index+1]),
              groups.order[break.index],
              events.wide[[paste0("date.stage.",stage.index,".group")]])
              #events.wide[[paste0("date.stage.group.",group.index)]])
          }
          events.wide[[paste0("date.stage.",stage.index,".group")]] <- ifelse(
            is.na(events.wide[[paste0("date.stage.",stage.index,".group")]]),
            "No group membership",
            events.wide[[paste0("date.stage.",stage.index,".group")]])
        }
    } else{}
    # Fix ordering, so that completion of a later stage indicates completion of earlier stage, and mark where this happens
    {
      for (stage.index in seq(length(stages.order)-1,1,-1)){
        # Determine origin of stage events, and create variable for its origin
          former.date.name <- paste("date.stage.",stage.index,sep="")
          latter.date.name <- paste("date.stage.",stage.index+1,sep="")
          former.date.origin.index.name <- paste("stage.",stage.index,".origin.index",sep="")
          latter.date.origin.index.name <- paste("stage.",stage.index+1,".origin.index",sep="")

        # Determine whether or not current stage is eligible for change
          events.wide$temp_flag <- ifelse(events.wide[[latter.date.name]] < events.wide[[former.date.name]] |
                                            is.na(events.wide[[former.date.name]])==TRUE,1,0)
          events.wide$temp_flag <- ifelse(is.na(events.wide$temp_flag)==TRUE,0,events.wide$temp_flag)
        # Check if this is the last event (determines whether there is an available following stage index)
          if (stage.index==(length(stages.order)-1)){
            # Generate on origin flag index, which = last stage if there is a change, and the current stage if there is none
              events.wide[[former.date.origin.index.name]] <- ifelse(
                events.wide$temp_flag == 1,
                stage.index + 1, stage.index
              )
              events.wide[[former.date.origin.index.name]] <- ifelse(
                is.na(events.wide[[former.date.name]])==TRUE & is.na(events.wide[[latter.date.name]])==TRUE,
                NA,events.wide[[former.date.origin.index.name]])
          } else {
            # Generate on origin flag index, which = the corresponding flag for the next stage if there is a change, and the current stage if there is none
              events.wide[[former.date.origin.index.name]] <- ifelse(
                events.wide$temp_flag == 1 & !is.na(events.wide[[latter.date.origin.index.name]]),
                events.wide[[latter.date.origin.index.name]],NA
              )
              events.wide[[former.date.origin.index.name]] <- ifelse(
                events.wide$temp_flag == 0,
                stage.index,events.wide[[former.date.origin.index.name]]
              )
          }
        # Replace value of date if a change is made
          events.wide[[paste("date.stage.",stage.index,sep="")]] <- ifelse(
            events.wide[[paste("stage.",stage.index,".origin.index",sep="")]] == stage.index,
            events.wide[[paste("date.stage.",stage.index,sep="")]],events.wide[[paste("date.stage.",stage.index +1,sep="")]]
          )
          events.wide[[paste("date.stage.",stage.index,sep="")]] <- as.Date.numeric(events.wide[[paste("date.stage.",stage.index,sep="")]])
      }
      events.wide$temp_flag <- NULL
    }
    # Determine dates of death and censorship, and add to dataset
    {
      # Merge in death events (if any)
        if (is.na(death.indicator)==FALSE){
          events.death <- subset(events.long.orig,events.long.orig$stage==death.indicator)
          events.death <- events.death[c("ID","date")]
          colnames(events.death) <- c("ID","date.death")
          events.wide <- merge(events.wide,events.death,by="ID",all.x = TRUE)
        }
      # Merge in censorship events (if any)
        if (is.na(censorship.indicator)==FALSE){
          events.censorship <- subset(events.long.orig,events.long.orig$stage==censorship.indicator)
          events.censorship <- events.censorship[c("ID","date")]
          colnames(events.censorship) <- c("ID","date.censorship")
          events.wide <- merge(events.wide,events.censorship,by="ID",all.x = TRUE)
        }
      # Assume that if date of death exists and occurs after censorship event, censorship is false and replace with NA
        if ((is.na(censorship.indicator)==FALSE) & (is.na(death.indicator)==FALSE)){
          events.wide$date.censorship <- ifelse((events.wide$date.censorship < events.wide$date.death) & (is.na(events.wide$date.death)==FALSE) & (is.na(events.wide$date.censorship)==FALSE),
                                                NA,events.wide$date.censorship)
          events.wide$date.censorship <- as.Date.numeric(events.wide$date.censorship)
        }
      # Generate time from each stage to death, censorship, and end of data collection
        for (i in 1:(length(stages.order)-1)){
          if (is.na(death.indicator)==FALSE){
            events.wide[[paste("time.stage.",i,".to.death",sep="")]] <-
              events.wide[["date.death"]] - events.wide[[paste("date.stage.",i,sep="")]]
          }
          if (is.na(censorship.indicator)==FALSE){
            events.wide[[paste("time.stage.",i,".to.censorship",sep="")]] <-
              events.wide[["date.censorship"]] - events.wide[[paste("date.stage.",i,sep="")]]
          }
        }
    }
    # Generate times between all possible events
    {
      for (i in 1:(length(stages.order)-1)){
        for (j in (i+1):(length(stages.order))){
          # Generate time to actual stage event if event occurs
            events.wide[[paste("time.stage.",i,".to.",j,sep="")]] <-
              events.wide[[paste("date.stage.",j,sep="")]] - events.wide[[paste("date.stage.",i,sep="")]]
          # Generate maximum data collection time
            max.time <- max(events.wide[[paste("time.stage.",i,".to.",j,sep="")]])
          # Generate last date collected
            max.date <- max(events.wide[[paste("date.stage.",j,sep="")]],na.rm=TRUE)
          # Generate maximum possible time for transition (date of event start to last data collection)
            events.wide[[paste("time.stage.",i,".to.",j,".last.date",sep="")]] <- max.date - events.wide[[paste("date.stage.",i,sep="")]]
          # Replace with time to censorship if event never happens and generate a censorship indicator flag
            if (is.na(censorship.indicator)==FALSE){
              events.wide[[paste("censorship.stage.",i,".to.",j,sep="")]] <- ifelse(
                is.na(events.wide[[paste("time.stage.",i,".to.",j,sep="")]])==TRUE & is.na(events.wide[[paste("time.stage.",i,".to.censorship",sep="")]])==FALSE,
                1,0
              )
              events.wide[[paste("time.stage.",i,".to.",j,sep="")]] <- ifelse(
                is.na(events.wide[[paste("time.stage.",i,".to.",j,sep="")]])==TRUE,
                events.wide[[paste("time.stage.",i,".to.censorship",sep="")]],events.wide[[paste("time.stage.",i,".to.",j,sep="")]]
              )
            } else{
              events.wide[[paste("censorship.stage.",i,".to.",j,sep="")]] <- 0
            }
            events.wide[[paste("censorship.stage.",i,".to.",j,sep="")]] <- ifelse(
              is.na(events.wide[[paste("time.stage.",i,".to.",j,sep="")]])==TRUE,
              NA,events.wide[[paste("censorship.stage.",i,".to.",j,sep="")]]
            )
          # Generate transition-specific mortality (i.e. mortality occuring only between each stage)
            if (is.na(death.indicator)==FALSE){
              events.wide[[paste("time.death.stage.",i,".to.",j,sep="")]] <- ifelse(
                (is.na(events.wide[[paste("time.stage.",i,".to.death",sep="")]]) == FALSE) & (is.na(events.wide[[paste("time.stage.",i,".to.",j,sep="")]])==TRUE),
                events.wide[[paste("time.stage.",i,".to.death",sep="")]],NA
              )
              events.wide[[paste("time.death.stage.",i,".to.",j,sep="")]] <- ifelse(
                (is.na(events.wide[[paste("date.stage.",i,sep="")]]) == FALSE) & (is.na(events.wide[[paste("time.death.stage.",i,".to.",j,sep="")]])==TRUE),
                events.wide[[paste("time.stage.",i,".to.",j,".last.date",sep="")]],events.wide[[paste("time.death.stage.",i,".to.",j,sep="")]]
              )
              events.wide[[paste("time.death.stage.",i,".to.",j,sep="")]] <- ifelse(
                events.wide[[paste("time.death.stage.",i,".to.",j,sep="")]] > events.wide[[paste("time.stage.",i,".to.",j,".last.date",sep="")]],
                events.wide[[paste("time.stage.",i,".to.",j,".last.date",sep="")]],events.wide[[paste("time.death.stage.",i,".to.",j,sep="")]]
              )
              events.wide[[paste("censorship.death.stage.",i,".to.",j,sep="")]] <- ifelse(
                events.wide[[paste("time.death.stage.",i,".to.",j,sep="")]] >= events.wide[[paste("time.stage.",i,".to.",j,".last.date",sep="")]],
                1,0
              )
            }
          # Replace transition times where initializing event occurs but finishing does not, with censorship after maximum data collection time
            events.wide[[paste("time.stage.",i,".to.",j,sep="")]] <- ifelse(
              is.na(events.wide[[paste("time.stage.",i,".to.",j,sep="")]]) == TRUE & is.na(events.wide[[paste("date.stage.",i,sep="")]]) == FALSE,
              events.wide[[paste("time.stage.",i,".to.",j,".last.date",sep="")]], events.wide[[paste("time.stage.",i,".to.",j,sep="")]]
            )
          # Replace transition times greater than axis limit, with censorship after maximum data collection time
            events.wide[[paste("time.stage.",i,".to.",j,sep="")]] <- ifelse(
              events.wide[[paste("time.stage.",i,".to.",j,sep="")]] > events.wide[[paste("time.stage.",i,".to.",j,".last.date",sep="")]],
              events.wide[[paste("time.stage.",i,".to.",j,".last.date",sep="")]], events.wide[[paste("time.stage.",i,".to.",j,sep="")]]
            )
          # Change to censorship event if censored after maximum data collection time
            events.wide[[paste("censorship.stage.",i,".to.",j,sep="")]] <- ifelse(
              events.wide[[paste("time.stage.",i,".to.",j,sep="")]] >= events.wide[[paste("time.stage.",i,".to.",j,".last.date",sep="")]],
              1, events.wide[[paste("censorship.stage.",i,".to.",j,sep="")]]
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
              if (anyNA(groups.date.breaks)==TRUE){
                events.wide$selection <- is.na(events.wide[[paste("time.stage.",start.stage.index,".to.",end.stage.index,sep="")]])==FALSE &
                events.wide$group == groups.order[group.index]
              } else {
                events.wide$selection <- is.na(events.wide[[paste("time.stage.",start.stage.index,".to.",end.stage.index,sep="")]])==FALSE &
                events.wide[[paste0("date.stage.",start.stage.index,".group")]] == groups.order[group.index]
              }
            } else {
              if (anyNA(groups.date.breaks)==TRUE){
                events.wide$selection <- is.na(events.wide[[paste("time.stage.",start.stage.index,".to.",end.stage.index,sep="")]])==FALSE &
                  events.wide[[paste("time.stage.",start.stage.index,".to.",end.stage.index,sep="")]]!=0 &
                  events.wide$group == groups.order[group.index]
              } else {
                events.wide$selection <- is.na(events.wide[[paste("time.stage.",start.stage.index,".to.",end.stage.index,sep="")]])==FALSE &
                  events.wide[[paste("time.stage.",start.stage.index,".to.",end.stage.index,sep="")]]!=0 &
                  events.wide[[paste0("date.stage.",start.stage.index,".group")]] == groups.order[group.index]
              }
              # Remove those who are disqualified from the "0+1" stage events from "0+n"
                events.wide$selection <- ifelse(events.wide[[paste("time.stage.",start.stage.index,".to.",(start.stage.index+1),sep="")]]==0,
                                                FALSE,
                                                events.wide$selection)
            }
          # Generate list of events
            chart.time <- as.integer(events.wide[[paste("time.stage.",start.stage.index,".to.",end.stage.index,sep="")]][events.wide$selection==TRUE])
            chart.event <- as.integer(!events.wide[[paste("censorship.stage.",start.stage.index,".to.",end.stage.index,sep="")]][events.wide$selection==TRUE])
          # Generate survival data
            chart.surv <- survival::survfit(survival::Surv(time = chart.time, event = chart.event) ~ 1)
            surv.time <- chart.surv$time
            surv.p <- 1-chart.surv$surv
            surv.p.UB <- 1-chart.surv$lower
            surv.p.LB <- 1-chart.surv$upper
            surv.n.t0 <- chart.surv$n
            surv.n.atrisk <- chart.surv$n.risk
            surv.p.atrisk <- chart.surv$n.risk/surv.n.t0
            surv.data <- data.frame(surv.time,surv.p,surv.p.UB,surv.p.LB,surv.n.t0,surv.n.atrisk,surv.p.atrisk)
            #surv.data <- data.frame(surv.time,surv.p)
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
          if (allow.sub.lines == FALSE) {surv.combined <- subset(surv.combined,surv.combined$start.stage.index==(surv.combined$end.stage.index)-1)}
      }
      # Death events
        if (is.na(death.indicator)==FALSE){
        # Function to generate survival data for a single curve
          gen.death.data <- function(start.stage.index,end.stage.index,group.index){
          # Select data
            if (anyNA(groups.date.breaks)==TRUE){
              events.wide$selection <- is.na(events.wide[[paste("time.death.stage.",start.stage.index,".to.",end.stage.index,sep="")]])==FALSE &
                events.wide$group == groups.order[group.index]
            } else {
              events.wide$selection <- is.na(events.wide[[paste("time.death.stage.",start.stage.index,".to.",end.stage.index,sep="")]])==FALSE &
              events.wide[[paste0("date.stage.",start.stage.index,".group")]] == groups.order[group.index]
            }
          # Generate list of events
            chart.time <- as.integer(events.wide[[paste("time.death.stage.",start.stage.index,".to.",end.stage.index,sep="")]][events.wide$selection==TRUE])
            chart.event <- as.integer(!events.wide[[paste("censorship.death.stage.",start.stage.index,".to.",end.stage.index,sep="")]][events.wide$selection==TRUE])
          # Generate survival data
            chart.surv <- survival::survfit(survival::Surv(time = chart.time, event = chart.event) ~ 1)
            surv.time <- chart.surv$time
            surv.p <- 1-chart.surv$surv
            surv.data <- data.frame(surv.time,surv.p)
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
    if (nochart==FALSE){
      # Settings
      {
        # Color gradient creator
          color.gradient <- function(original.color,number.divisions,darken=0){
            # Generate hsv color
              hsv.color.orig <- grDevices::rgb2hsv(grDevices::col2rgb(original.color))
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
        # X scale creator functions
          round2 <- function(x, n=0) {scale<-10^n; trunc(x*scale+sign(x)*0.5)/scale}
          x.scale.function <- function(x) round2(x,0)
          x.axis.range <- c(0,x.axis.max)
        # Generate / change colors
          if (length(main.fill.colors)==1){
            main.fill.colors <- color.gradient(main.fill.colors,(length(stages.order)-1))
          } else {}

      }
      # Graphical functions
      {
        # Step ribbon function. Note: Original author of code was Triad sou from the RcmdrPlugin.KMggplot2 package
        {
          geom_stepribbon <- function(
            mapping = NULL, data = NULL, stat = "identity", position = "identity",
            na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, kmplot = FALSE, ...) {
            layer(
              data = data,
              mapping = mapping,
              stat = stat,
              geom = GeomStepribbon,
              position = position,
              show.legend = show.legend,
              inherit.aes = inherit.aes,
              params = list(
                na.rm = na.rm,
                kmplot = kmplot,
                ...
              )
            )
          }
          GeomStepribbon <- ggproto(
            "GeomStepribbon", GeomRibbon,
            extra_params = c("na.rm", "kmplot"),
            draw_group = function(data, panel_scales, coord, na.rm = FALSE) {
              if (na.rm) data <- data[complete.cases(data[c("x", "ymin", "ymax")]), ]
              data <- rbind(data, data)
              data <- data[order(data$x), ]
              data$x <- c(data$x[2:nrow(data)], NA)
              data <- data[complete.cases(data["x"]), ]
              GeomRibbon$draw_group(data, panel_scales, coord, na.rm = FALSE)
            },
            setup_data = function(data, params) {
              if (params$kmplot) {
                data <- data[order(data$PANEL, data$group, data$x), ]
                tmpmin <- tmpmax <- NA
                for (i in 1:nrow(data)) {
                  if (is.na(data$ymin[i])) {
                    data$ymin[i] <- tmpmin
                  }
                  if (is.na(data$ymax[i])) {
                    data$ymax[i] <- tmpmax
                  }
                  tmpmin <- data$ymin[i]
                  tmpmax <- data$ymax[i]
                }
              }
              data
            }
          )
        }
      }
      # Data manipulation
      {
        # Main events
        {
          # Generate data
            surv.combined.chart <- surv.combined
          # Drop any events past x axis range to speed up calcs
            surv.combined.chart <- subset(surv.combined.chart,surv.combined.chart$surv.time<=x.axis.max)
          # Drop out duplicate boxes (due to censoring) to reduce drawing time
            surv.combined.chart <- surv.combined.chart %>%
              dplyr::arrange(.data$group.factor,.data$start.stage.factor,.data$end.stage.factor,.data$surv.time) %>%
              dplyr::group_by(.data$surv.p,.data$start.stage.factor,.data$end.stage.factor,.data$group.factor) %>%
              dplyr::slice(c(n()))
          # Generate beginning axis events to keep fill graphics going from start of chart at 0
            surv.combined.chart.beginning <- surv.combined.chart
            surv.combined.chart.beginning$surv.p <- 0
            surv.combined.chart.beginning$surv.time <- -1
            surv.combined.chart.beginning$surv.p.LB <- 0
            surv.combined.chart.beginning$surv.p.UB <- 0
            surv.combined.chart.beginning$surv.n <- 0
            surv.combined.chart.beginning$surv.p.atrisk <- 1
            surv.combined.chart.beginning <- unique(surv.combined.chart.beginning)
            surv.combined.chart.beginning <- surv.combined.chart.beginning[, colnames(surv.combined.chart)]
            surv.combined.chart <- rbind(surv.combined.chart,surv.combined.chart.beginning)
            # Rearrange and sort for drawing
              surv.combined.chart <- surv.combined.chart %>%
                dplyr::arrange(.data$group.factor,.data$start.stage.factor,.data$end.stage.factor,.data$surv.time)
          # Generate end of axis events to keep fill graphics going to end of chart
            surv.combined.chart.end <- surv.combined.chart %>%
              dplyr::group_by(.data$start.stage.index,.data$end.stage.index,.data$group.index) %>%
              dplyr::slice(which.max(.data$surv.p))
            surv.combined.chart.end$surv.time <- x.axis.max + 1
            surv.combined.chart <- rbind(surv.combined.chart,surv.combined.chart.end)
            surv.combined.chart.end$surv.p <- 0
            surv.combined.chart.end$surv.p.atrisk <- 1
            surv.combined.chart <- rbind(surv.combined.chart,surv.combined.chart.end)
          # Rearrange and sort for drawing
            surv.combined.chart <- surv.combined.chart %>%
              dplyr::arrange(.data$group.factor,.data$start.stage.factor,.data$end.stage.factor,.data$surv.time)
          # Temporary for putting in years
            surv.combined.chart$surv.time = surv.combined.chart$surv.time/365
          # Generate risk pool dataset
            surv.combined.chart.risk.pool <- surv.combined.chart[surv.combined.chart$end.stage.index == (surv.combined.chart$start.stage.index +1) ,]
          # Add a dummy factor for legend
            surv.combined.chart.risk.pool$event.factor <- factor("Risk pool",labels = c("Risk pool"),levels = c("Risk pool"))
          # Generate version which includes end steps for stepwise fill
        }
        # Death events
        {
          if (is.na(death.indicator)==FALSE){
            # Generate data
              surv.death.combined.chart <- surv.death.combined
            # Drop any events past x axis range
              surv.death.combined.chart <- subset(surv.death.combined.chart,surv.death.combined.chart$surv.time<=x.axis.max)
            # Drop out duplicate boxes (due to censoring) to reduce drawing time
              surv.death.combined.chart <- surv.death.combined.chart %>%
                dplyr::arrange(.data$group.factor,.data$start.stage.factor,.data$end.stage.factor,.data$surv.time) %>%
                dplyr::group_by(.data$surv.p,.data$start.stage.factor,.data$end.stage.factor,.data$group.factor) %>%
                dplyr::slice(c(n()))
            # Generate beginning of axis events to keep graphics going until the end of chart
              surv.death.combined.chart.extra <- surv.death.combined.chart
              surv.death.combined.chart.extra$surv.p <- 0
              surv.death.combined.chart.extra$surv.time <- 0
              surv.death.combined.chart.extra <- unique(surv.death.combined.chart.extra)
              surv.death.combined.chart.extra <- surv.death.combined.chart.extra[, colnames(surv.death.combined.chart)]
              surv.death.combined.chart <- rbind(surv.death.combined.chart,surv.death.combined.chart.extra)
            # Rearrange and sort for drawing
              surv.death.combined.chart <- surv.death.combined.chart %>%
                dplyr::arrange(.data$group.factor,.data$start.stage.factor,.data$end.stage.factor,.data$surv.time)
            # Generate end of axis events to keep graphics going until the end of chart
              surv.death.combined.chart.extra <- surv.death.combined.chart
              surv.death.combined.chart.extra$surv.p <- 0
              surv.death.combined.chart.extra$surv.time <- x.axis.max
              surv.death.combined.chart.extra <- unique(surv.death.combined.chart.extra)
              surv.death.combined.chart.extra <- surv.death.combined.chart.extra[, colnames(surv.death.combined.chart)]
              surv.death.combined.chart <- rbind(surv.death.combined.chart,surv.death.combined.chart.extra)
            # Rearrange and sort for drawing
              surv.death.combined.chart <- surv.death.combined.chart %>%
                dplyr::arrange(.data$group.factor,.data$start.stage.factor,.data$end.stage.factor,.data$surv.time)
            # Temporary for putting in years
              surv.death.combined.chart$surv.time = surv.death.combined.chart$surv.time/365
            # Add dummy event factor for legends
              surv.death.combined.chart$event.factor <- factor(death.indicator,labels = c(death.indicator),levels = c(death.indicator))
          }
        }
      }
      # Generate chart
      {
        chart <- ggplot2::ggplot() +
          geom_stepribbon( data=surv.combined.chart,aes(x=.data$surv.time,ymax=.data$surv.p,fill=.data$end.stage.factor,ymin=0),
                                                 alpha=1) +
          ggplot2::geom_step(data=surv.combined.chart,aes(x=.data$surv.time,y=.data$surv.p,color=.data$end.stage.factor),
                             show.legend = FALSE) +
          ggplot2::theme_bw() %+replace%
          ggplot2::theme(
            panel.grid = element_blank(),
            plot.margin = unit(c(.1,.1,.1,.1), "cm"),
            axis.title.y=element_blank(),
            legend.position="bottom",
            legend.title=element_blank(),
            axis.text = element_text(colour="black",size=10),
            strip.text = element_text(size = 12),
            strip.text.x = element_text(hjust=0),
            panel.spacing = unit(1, "lines")
          ) +
          ggplot2::guides(color = guide_legend(override.aes = list(linetype = 0))) +
          ggplot2::scale_x_continuous(expand = c(0, 0),
                            labels=x.scale.function,
                            breaks = c(0,round2(x.axis.range/365,0))) +
          ggplot2::xlab("Time (years) from start of stage") +
          ggplot2::scale_y_continuous(expand = c(0, 0),labels=percent) +
          ggplot2::scale_color_manual(values=c(rep("black",length(stages.order)-1))) +
          ggplot2::facet_grid(group.factor ~ start.stage.factor,
                     switch="y") +
          ggplot2::theme(strip.background = element_blank(),
                strip.placement = "outside")
      # Add risk pool proportion indicator if indicated
        if (risk.pool.size.line==TRUE){
          chart <- chart +
            ggplot2::coord_cartesian(xlim=c(0,(x.axis.range/365)),ylim = c(-.2, 1)) +
            geom_stepribbon(data = surv.combined.chart.risk.pool,
                                  aes(x=.data$surv.time,ymin=((.data$surv.p.atrisk-1)/5),ymax=0),
                                  alpha=1,fill=risk.pool.fill.color) +
            ggplot2::geom_step(data = surv.combined.chart.risk.pool,aes(x=.data$surv.time,y=(.data$surv.p.atrisk-1)/5))
        } else {
          chart <- chart +
            ggplot2::coord_cartesian(xlim=c(0,(x.axis.range/365)),ylim = c(0, 1))
        }
      # Add death event if present
        if (is.na(death.indicator)==FALSE){
          chart <- chart +
            geom_stepribbon(data=surv.death.combined.chart,aes(x=.data$surv.time,ymin=1-.data$surv.p,ymax=1),
                                  alpha=1,fill=death.fill.color) +
            #ggplot2::geom_polygon(data=surv.death.combined.chart,aes(x=.data$surv.time,y=1-.data$surv.p,fill=.data$event.factor),
            #                      alpha=1) +
            ggplot2::geom_step(data=surv.death.combined.chart,aes(x=.data$surv.time,y=1-.data$surv.p))
        } else {}
      # Add touchup graphics
        chart <- chart +
          ggplot2::geom_hline(yintercept=0)
      # Set fill colors
        chart <- chart +
          ggplot2::scale_fill_manual(values = main.fill.colors)
      }
    }
    else {
      chart = FALSE
    }
  }
  # Generate between-group statistical tests
  {
    # Test for heterogeneity of survival curves for main between-stage events between groups
    if (length(groups.order)>1){
      # Generate tests function
        generate.surv.diff <- function(init.stage){
          # Generate internal data
            if (anyNA(groups.date.breaks)==FALSE){
              #events.wide.test <- events.wide[0,]
              generate.date.break.data <- function(stage.index,group.index){
                df.output <- events.wide[events.wide[[paste0("date.stage.",stage.index,".group")]]==groups.order[group.index],]
                df.output$group <- groups.order[group.index]
                return(df.output)
              }
              events.wide.test <- do.call("rbind",
                                   lapply(1:length(groups.order),function(x) generate.date.break.data(init.stage,x)))
            } else {
              events.wide.test <- events.wide
            }
          # Perform test
            chart.time <- as.integer(events.wide.test[[paste("time.stage.",(init.stage),".to.",(init.stage+1),sep="")]])
            chart.event <- as.integer(!events.wide.test[[paste("censorship.stage.",(init.stage),".to.",(init.stage+1),sep="")]])
            groups <- ordered(events.wide.test$group,labels=groups.order,levels=groups.order)
            test.name <- paste0("curve.het.stage.",(init.stage),"to",(init.stage+1))
            surv.diff <- list(survdiff(Surv(time = chart.time, event = chart.event) ~ groups))
            names(surv.diff)[1] <- test.name
          # Return results
            return(surv.diff)
        }
      # Package tests
        surv.diffs.combined <- do.call("list",
                                 lapply(1:(length(stages.order)-1),function(x) generate.surv.diff(x)))
        names.list <- paste0("init.stage.",1:(length(stages.order)-1))
        names(surv.diffs.combined) <- names.list
    }
    # Test for between-group differences, using Cox proportional hazards at each stage transition
    if (length(groups.order)>1){
      # Function to take two groups ad output Cox PH regression stats
        gen.cox.ph <- function(init.stage,reference.group){
          # Generate internal data
            if (anyNA(groups.date.breaks)==FALSE){
              #events.wide.test <- events.wide[0,]
              generate.date.break.data <- function(stage.index,group.index){
                df.output <- events.wide[events.wide[[paste0("date.stage.",stage.index,".group")]]==groups.order[group.index],]
                df.output$group <- groups.order[group.index]
                return(df.output)
              }
              events.wide.test <- do.call("rbind",
                                   lapply(1:length(groups.order),function(x) generate.date.break.data(init.stage,x)))
            } else {
              events.wide.test <- events.wide
            }
          # Keep only relevent variables for internal manipulation
            time.var <- paste("time.stage.",(init.stage),".to.",(init.stage+1),sep="")
            events.var <- paste("censorship.stage.",(init.stage),".to.",(init.stage+1),sep="")
            events.wide.internal <- events.wide.test[c("group",time.var,events.var)]
            colnames(events.wide.internal) <- c("group","time","event")
          # Designate reference vs. comparator groups
            #group.order.internal <- c(reference.group,groups.order[groups.order!=reference.group])
            events.wide.internal$group <- factor(events.wide.internal$group,levels=groups.order,labels=groups.order)
            events.wide.internal$group <- stats::relevel(events.wide.internal$group, ref = reference.group)
          # Run Cox Proportional hazards
            cox.ph <- coxph(Surv(time = time, event = event) ~ group,data=events.wide.internal)
            # name <- paste0("ref.group.",reference.group)
            # list.form <- list(cox.ph)
            # names(list.form) <- name
            return(cox.ph)
        }
      # Generate tests for each reference group, within a beginning stage
        gen.cox.ph.groups <- function(init.stage){
          output <- do.call("list",
                  lapply(1:(length(groups.order)),function(x) gen.cox.ph(init.stage,x)))
          names.list <- paste0("ref.group.",1:(length(groups.order)))
          names(output) <- names.list
          return(output)
        }
      # Generate set of tests for each beginning stage
        gen.cox.ph.stages <- function(){
          output <- do.call("list",
                            lapply(1:(length(stages.order)-1),function(x) gen.cox.ph.groups(x)))
          names.list <- paste0("init.stage.",1:(length(stages.order)-1))
          names(output) <- names.list
          return(output)
        }
        cox.ph.combined <- gen.cox.ph.stages()
    }

  }
  # Prepare export data
  {
    if (length(groups.order)>1){
      if (is.na(death.indicator)){
        output.df <- list("chart" = chart,"surv.dataset" = surv.combined,
                          "events.long" = events.long,"events.wide" = events.wide,
                          "surv.cox.ph" = cox.ph.combined,"surv.diffs" = surv.diffs.combined)
      } else{
        output.df <- list("chart" = chart,"surv.dataset" = surv.combined,"surv.death.combined" = surv.death.combined,
                          "events.long" = events.long,"events.wide" = events.wide,
                          "surv.cox.ph" = cox.ph.combined,"surv.diffs" = surv.diffs.combined)
      }
    } else {
      if (is.na(death.indicator)){
        output.df <- list("chart" = chart,"surv.dataset" = surv.combined,
                          "events.long" = events.long,"events.wide" = events.wide)
      } else {
        output.df <- list("chart" = chart,"surv.dataset" = surv.combined,"surv.death.combined" = surv.death.combined,
                          "events.long" = events.long,"events.wide" = events.wide)
      }
    }
  }
}
