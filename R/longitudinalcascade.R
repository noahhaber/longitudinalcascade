#' This package generates a longitudinal cascade, including a graphical representation. This takes a long-formatted list of stage-by-stage events and transforms it into a longitudinal cascade, correcting the orders of events.
#' @title Longitudinal cascade statistics and charts
#' @name longitudinalcascade
#' @keywords cascade longitudinal survival
#' @param events.long (required) The main dataframe input parameter. The data frame needs at least the following fields:
#' "ID": (required) A string-based individual identifier, indicating every individual in the dataset.
#' "date": (required) Date-formatted date on which the event / stage occurred
#' "stage": (required) String indicating the stage achieved by the individual on the specified date. Stages must match the string in the stages.order parameter. Additonal events may be included in the "stage" category, including death, loss to follow up, and interstage events defined in the other parameters.
#' "group": (optional) String indicating any relevant groups of data.
#' @param stages.order (required) stages.order is the parameter which defines the events to be considered in the main cascade and their order. This is a vector of strings matching items in the "Stage" column of the main data frame, e.g. c("Stage 1","Stage 2","Stage 3").
#' @param groups.order (optional) This is a vector of groups, matching the "group" column of the main data frame. If left blank, no group comparisons will be performed. For the chart, each group will have its own row.
#' @param groups.date.breaks (optional) If groups.date.breaks is filled in, the grouping will be defined by date range of entry event for each transition, rather than groups of individuals. Each transition will independently determine its own groups, based on the time in which the entrance event occurs. Times are determined by a vector of date breaks. Each group is defined as starting from a given date break value and continuing until it reaches the subsequent date break, not including data from that ending break value. For example, setting the break values to be January 1, 2011, January 1, 2012, and January 1, 2013 will create two groups. The first group will take individuals who entered each stage from January 1, 2011 to Dec 31, 2011, and the second will take individuals who entered into the stage from January 1, 2012 to Dec 31, 2012.
#' @param groups.date.breaks.labels (optional) Changes the default labelling of the groups when you are using date break groupings. Entered as a vector of strings, and must be the same length as the number of groups.
#' @param death.indicator (optional) This parameter is the string which indicates a death event in the dataset. If specified, between-stage mortality will be estimated and shown as a KM curve on the top of the chart(s). If left blank, death events will not be estimated.
#' @param censorship.indicator (optional) This parameter is the string which indicates a right-censorship event. Most commonly, this will indicate permanent loss to follow up and/or end of data collection.
#' @param censorship.date (optional) By default, censorship is set to the last date of data collection. If you would prefer to set a different date than that, enter it into censorship.date argument as a date.
#' @param allow.sub.stages Sub-lines indicate subsequent transitions across the cascade. If TRUE, the main chart will show transitions to all possible subsequent events. For example, if there are 4 stages (1-4), the leftmost chart will show each transition from 1-2, 1-3, and 1-4, while the next chart will show 2-3 and 2-4, and the last chart will show only 3-4. If FALSE, the charts will only show transition to the subsequent stage.
#' @param allow.sub.stage.mortality Sub-stage-mortality indicate subsequent mortality transitions across the cascade. If TRUE, the main chart will show transitions to all possible subsequent events. For example, if there are 4 stages (1-4), the leftmost chart will show each transition from 1-2, 1-3, and 1-4, while the next chart will show 2-3 and 2-4, and the last chart will show only 3-4. If FALSE, the charts will only show transition to the subsequent stage.
#' @param sub.stage.mortality.mode By default, sub-stage mortality is shown as a transition underneath the main transition ("standard"). If set to "shifted" substage mortality will be shifted to the top of the chart, and all substages will be shifted downward accordingly
#' @param skip.mode This option shows "skips" across the cascade in each chart, as indicated by the y intercept. If "none" (default) each stage will start only with people who have not moved on to a subsequent stage, i.e. the y intercept will always be 0. If set to "internal" an individual can enter into a stage even if they have "skipped" through it. For example, an individual may go straight from stage 1 to stage 3, skipping 2. If this indicator is FALSE, the stage transition chart from 2-3 will not contain this individual in the denomenator. If TRUE, this individual will be counted in the denomenator for this transition, but will be counted as having transitioned into stage 3 immediately upon entering stage 2. If "external" individuals contribute person-time and are in the y-axis of transitions even prior to their first recorded stage date.
#' @param time.horizon This option shows the maximum range of each stage in days. Defaults to 365 days (1 year).
#' @param TTE.quantiles This option sets the quantiles measured for the quantile time to event outputs, using a c() list. By default, this is set to 0.2, 0.5 (i.e. the median), and 0.75.
#' @param chart.mode By default, the chart is set to a stage-by-stage panel view ("stage panels"). Alternatively, it may be desirable to have only the first panel showing the overall experience from the first entry condition, as indicated by the "first transition" option.
#' @param ts.indicator (experimental) (optional) This indicates the name of the events indicating the transient stage of interest. Can be either a single indicator, or a c() vector if there are multiple transient stages defined. Defaults to NA, disabling this feature.
#' @param ts.gap.time (experimental) (optional) This indicates the time between events required until an indidual is considered "off" for the given transient stage. Can be either a single indicator, or a c() vector if there are multiple transient stages defined.
#' @param ts.start.stage (experimental) (optional) This indicates the stage at which the transient stage of interest starts. Can be either a single indicator, or a c() vector if there are multiple transient stages defined.
#' @param ts.end.stage (experimental) (optional) This indicates the stage at which the transient stage of interest ends Can be either a single indicator, or a c() vector if there are multiple transient stages defined.
#' @param ts.color (experimental) (optional) Indicates the color used for transient stages
#' @param risk.pool.size.line Setting to TRUE adds an indicator of risk pool remaining to the main charts as a line reflected beneath the main chart, showing the proportion of the original risk pool remaining at each time point. Defaults to FALSE.
#' @param main.fill.colors (optional) This defines the color scheme of the stage transition graphs, as a string indicator for color or a c() vector of colors. If the colors contain only one color, the color scheme will automatically generate progressively faded versions of the initial color provided for the remaining stage transitions. Otherwise, a list which is exactly one fewer than the # of stages must be provided, in the order of stage trasitions.
#' @param death.fill.color (optional) This defines the color scheme for the death stage transition, as a string indicator for color.
#' @param risk.pool.fill.color (optional) This defines the color scheme for the risk pool graphic, as a string indicator for color.
#' @param background.prior.event (optional) This changes the background of the faceted chart to be the color for the prior event.
#' @param suppress.messages (optional) Suppresses tips and messages about the dataset
#' @param x.axis.title (optional) Changes the x axis label
#' @param direct.label (optional) Adds direct labeling of the survival lines via the geom_textpath command
#' @param legend.position (optional) Changes legend position (passing on to ggplot theme)
#' @return description All data are output to an object containing a the main chart ($chart), a survival-formatted dataset ($surv.dataset), the data underlying the main chart ($surv.dataset.chart), the underlying original dataset in long ($events.long) and wide ($events.wide), individual time to event data ($TTE.ind), TTE data by quantiles ($quantile.TTE), the equivalent functions for deat and transient events, and group difference tests
#' @import survival ggplot2 dplyr tidyr zoo scales grDevices
#' @importFrom stats relevel
#' @importFrom rlang .data
#' @importFrom lubridate is.Date
#' @importFrom stats na.omit
#' @export
#' @references Haber et al. (2017) Lancet HIV 4(5):e223-e230
#' (\href{https://pubmed.ncbi.nlm.nih.gov/28153470/}{PubMed})
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
#' allow.sub.stages <- TRUE
#'
#' # Create cascade object
#' longitudinalcascade.sim <- longitudinalcascade(events_long,stages.order=stages.order,
#' groups.order=groups.order,death.indicator=death.indicator,
#' censorship.indicator=censorship.indicator,
#' allow.sub.stages=allow.sub.stages)
#'
#' # Print/output main multipanel chart
#' longitudinalcascade.sim$chart
#' # Output full survival dataset generated, as a data frame
#' df.longitudainalcascade.survival <- longitudinalcascade.sim$surv.dataset
#' # Output heterogeneity test
#' longitudinalcascade.sim$surv.diffs
#' # Output original long-formatted list of events
#' df.events.long <- longitudinalcascade.sim$events.long
#' # Output generated wide-formatted list of events
#' df.events.wide <- longitudinalcascade.sim$events.wide

longitudinalcascade <- function(events.long,stages.order,
                         groups.order=NA,groups.date.breaks=NA,groups.date.breaks.labels=NA,
                         death.indicator=NA,censorship.indicator=NA,censorship.date = "lastdate",
                         allow.sub.stages=FALSE,allow.sub.stage.mortality=FALSE,sub.stage.mortality.mode="standard",
                         skip.mode="none",
                         time.horizon=365,
                         TTE.quantiles=c(0.25,0.50,0.75),
                         main.fill.colors = "#4472C4",death.fill.color = "#FF6A6A",
                         chart.mode = "stage panels",
                         ts.indicator = NA,ts.gap.time = 90,ts.start.stage=1,ts.end.stage=NA,ts.color="#a2f2da",
                         risk.pool.size.line=FALSE,
                         risk.pool.fill.color = "#90dbb2",
                         background.prior.event=TRUE,
                         suppress.messages = FALSE,
                         x.axis.title = "Time (years) from start of stage",
                         direct.label = TRUE,
                         legend.position = "bottom") {

  # Functions for general use
  {
    # Function to handle warnings and warning suppression
      warning.f <- function(warning.message){
        if (suppress.messages==FALSE){
          message(warning.message)
        } else {}
      }
    # Function for adding curves on top of/beneath other curves
      relative.curve <- function(reference,input,over.under="over"){
        # Sort inputs
          input <- input[order(input$surv.time),]
          reference <- reference[order(reference$surv.time),]
        # Generate indexes for merging
          input$index <- 1:nrow(input)
          reference$index <- 1:nrow(reference)
        # Generate mutual times for each curve
          surv.time.combined <- sort(unique(c(reference$surv.time,input$surv.time)))
          relative <- data.frame(surv.time.combined)
        # Find the index for which time is equal to or just under the time of interest for reference curve
          relative$index.obs.reference <- sapply(relative$surv.time.combined,function(x) sum(reference$surv.time <= x))
          relative$index.obs.input <- sapply(relative$surv.time.combined,function(x) sum(input$surv.time <= x))
        # Pull in the survival levels at each curve
          relative <- merge(relative,input,by.x="index.obs.input",by.y="index",all.x=TRUE)
          relative$surv.p.input <- relative$surv.p
          relative$surv.p <- NULL
          relative$surv.time.input <- relative$surv.time
          relative$surv.time <- NULL
          relative <- merge(relative,reference,by.x="index.obs.reference",by.y="index",all.x=TRUE)
          relative$surv.p.reference <- relative$surv.p
          relative$surv.p <- NULL
          relative$surv.time.reference <- relative$surv.time
          relative$surv.time <- NULL
          relative$surv.time <- relative$surv.time.combined
        # Generate relative curve
          if (over.under == "over") {
            relative$surv.p <- relative$surv.p.reference + relative$surv.p.input
          } else {
            relative$surv.p <- relative$surv.p.reference - relative$surv.p.input
          }
        # Trim and output new curves
          relative <- relative[c("surv.p","surv.time","surv.p.reference","surv.p.input")]
          relative <- relative[!is.na(relative$surv.p),]
          relative <- relative[order(relative$surv.time),]
          return(relative)
      }
    # Function for selecting which stats to output from survival curves
      surv.data.stats <- function(time.to.event,censorship,weights=NA,...){
        # Main stats calculations
          event <- 1-censorship
          if (!anyNA(weights)){
            surv <- survival::survfit(survival::Surv(time = time.to.event, event = event) ~ 1,weights=weights)
          } else {
            surv <- survival::survfit(survival::Surv(time = time.to.event, event = event) ~ 1)
          }
          surv.time <- surv$time
          surv.p <- 1-surv$surv
          surv.p.UB <- 1-surv$lower
          surv.p.LB <- 1-surv$upper
          surv.n.t0 <- surv$n
          surv.n.atrisk <- surv$n.risk
          surv.p.atrisk <- surv$n.risk/surv.n.t0
          surv.data <- data.frame(surv.time,surv.p,surv.p.UB,surv.p.LB,surv.n.t0,surv.n.atrisk,surv.p.atrisk)
        # Add a time-horizon event if one does not already exist
          if (max(surv.data$surv.time)<time.horizon){
            surv.data <- rbind(surv.data,surv.data[nrow(surv.data),])
            surv.data[nrow(surv.data),]$surv.time <- time.horizon
          } else {}
        # Add a day 0 event if one does not already exist
          if (surv.data$surv.time[1]>0){
            surv.data <- rbind(surv.data[1,],surv.data)
            surv.data[1,]$surv.time <- 0
            surv.data[1,]$surv.p <- 0
            surv.data[1,]$surv.p.UB <- 0
            surv.data[1,]$surv.p.LB <- 0
          } else if (surv.data$surv.time[1]<0){
            message("Error, negative time disallowed")
          } else {}
        # Generate percentile TTE (i.e. median, etc)
          TTE.quantile <- quantile.survfit(surv, TTE.quantiles,conf.int=TRUE)
        # Prepare export data
          output <- list("surv.curve" = surv.data,"survfit.object" = surv,"TTE.quantile" = TTE.quantile)
        # Export data
          return(output)
      }
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
          return(rev(unname(sapply(data.frame(hsv.color.new), function(x) do.call(hsv, as.list(x))))))
      }
    # X scale creator functions
      round2 <- function(x, n=0) {scale<-10^n; trunc(x*scale+sign(x)*0.5)/scale}
      x.scale.function <- function(x) round2(x,0)
    # Step ribbon function. Note: Original author of code was Triad Sou from the RcmdrPlugin.KMggplot2 package
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
          #if (na.rm) data <- data[complete.cases(data[c("x", "ymin", "ymax")]), ]
          data <- rbind(data, data)
          data <- data[order(data$x), ]
          data$x <- c(data$x[2:nrow(data)], NA)
          #data <- data[complete.cases(data["x"]), ]
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
    # Quantile survfit, Note: Original author of this code is Terry M Therneau
      {
        findq <- function(x, y, p, tol) {
            # This case occurs for a survival curve whose upper limit never drops below 1
            if (max(y, na.rm=T) < min(p)) return(rep(NA, length(p)))

            # Remove duplicate y values, i.e., the censors, since dups cause
            #  issues for approx
            xmax <- x[length(x)]
            dups <- duplicated(y)
            if (any(dups)) {
                x <- x[!dups]
                y <- y[!dups]
            }
            n <- length(y)
            indx1 <- stats::approx(y+tol, 1:n, p, method="constant", f=1)$y
            indx2 <- stats::approx(y-tol, 1:n, p, method="constant", f=1)$y
            quant <- (x[indx1] + x[indx2])/2
            quant[p==0] <- x[1]
            if (!is.na(y[n])) {
                lastpt <- (abs(p- y[n]) < tol)  # end of the curve
                if (any(lastpt)) quant[lastpt] <- (x[indx1[lastpt]] + xmax)/2
            }
            quant
            }
        doquant <- function(p, time, surv, upper, lower, firstx, tol) {
            qq <- findq(c(firstx,time), c(0, 1-surv), p, tol)
            if (missing(upper)) qq
            else rbind(qq, findq(c(firstx, time), c(0, 1-lower), p, tol),
                           findq(c(firstx, time), c(0, 1-upper), p, tol))
            }
        quantile.survfit <- function(x, probs=c(.25, .5, .75), conf.int=TRUE,
                                     tolerance= sqrt(.Machine$double.eps), ...) {
            if (!inherits(x, "survfit")) stop("Must be a survfit object")
            if (any(!is.numeric(probs)) || any(is.na(probs)))
                stop("invalid probability")
            if (any(probs <0 | probs >1)) stop("Invalid probability")
            if (is.null(x$lower)) conf.int <- FALSE
            nprob <- length(probs)
            pname <- format(probs*100)

            # What do we report for p=0?  Use x$start.time if it exists, 0 otherwise
            xmin <- if (is.null(x$start.time)) 0 else x$start.time

            # There are 8 cases: strata yes/no
            #                    ncol(x$surv) =1 or >1
            #                    conf.int = T/F
            if (is.null(x$strata)) {
                if (is.matrix(x$surv) && ncol(x$surv) >1) {
                    qmat <- matrix(0., ncol=nprob, nrow=ncol(x$surv))
                    dimnames(qmat) <- list(dimnames(x$surv)[[2]], pname)
                    if (conf.int) {
                        qupper <- qlower <- qmat
                        for (i in 1:ncol(x$surv)) {
                            temp <- doquant(probs, x$time, x$surv[,i], x$upper[,i],
                                            x$lower[,i], xmin, tolerance)
                            qmat[i,] <- temp[1,]
                            qupper[i,] <- temp[3,]
                            qlower[i,] <- temp[2,]
                        }
                        return(list(quantile=qmat, lower=qlower, upper=qupper))
                    }
                    else {
                        for (i in 1:ncol(x$surv))
                            qmat[i,] <- doquant(probs, x$time, x$surv[,i], firstx=xmin,
                                                tol=tolerance)
                        return(qmat)
                    }
                }
                else {
                    # No strata and no matrix
                    if (conf.int) {
                        temp <- doquant(probs, x$time, x$surv, x$upper, x$lower, xmin,
                                        tolerance)
                        dimnames(temp) <- list(NULL, pname)
                        return(list(quantile=temp[1,], lower=temp[2,], upper=temp[3,]))
                    }
                    else {
                        temp <- doquant(probs, x$time, x$surv, firstx=xmin,
                                        tol =tolerance)
                        names(temp) <- pname
                        return(temp)
                    }
                }
            }

            else {
                nstrat <- length(x$strata)
                if (is.matrix(x$surv) && ncol(x$surv) >1) {
                    # uncommon case, e.g., predicted survivals from a Cox model
                    # return an array with strata as the first dimension, and
                    #   the probabilites as the third.
                    qmat <- array(0., dim=c(nstrat, ncol(x$surv), nprob))
                    dimnames(qmat) <-list(names(x$strata), dimnames(x$surv)[[2]], pname)
                    if (conf.int) {
                        qupper <- qlower <- qmat
                        for (strat in 1:nstrat) {
                            z <- x[strat,]
                            for (i in 1:ncol(z$surv)) {
                                temp <- doquant(probs, z$time, z$surv[,i],
                                                z$upper[,i], z$lower[,i], xmin,tolerance)
                                qmat[strat,i,] <- temp[1,]
                                qupper[strat,i,] <- temp[3,]
                                qlower[strat,i,] <- temp[2,]
                            }
                        }
                        return(list(quantile=qmat, lower=qlower, upper=qupper))
                    }
                    else {
                        for (strat in 1:nstrat) {
                            z <- x[strat]
                            for (i in 1:ncol(z$surv))
                                qmat[strat,i,] <- doquant(probs, z$time, z$surv[,i],
                                                          firstx=xmin, tol=tolerance)
                        }
                        return(qmat)
                    }
                 }
                else {
                    # Only a strata, the most common case
                    qmat <- matrix(0., nstrat, nprob)
                    dimnames(qmat) <- list(names(x$strata), pname)
                    if (conf.int) {
                        qupper <- qlower <- qmat
                        for (i in 1:nstrat) {
                            z <- x[i]
                            temp <- doquant(probs, z$time, z$surv, z$upper, z$lower,
                                            xmin, tolerance)
                            qmat[i,] <- temp[1,]
                            qupper[i,] <- temp[3,]
                            qlower[i,] <- temp[2,]
                        }
                        return(list(quantile=qmat, lower=qlower, upper=qupper))
                    }
                    else {
                        for (i in 1:nstrat) {
                            z <- x[i]
                            qmat[i,] <- doquant(probs, z$time, z$surv, firstx=xmin,
                                                tol = tolerance)
                        }
                        return(qmat)
                    }
                }
            }
        }
      }
    # Function for rbinding ALL of the outputs of a function
      rbind.lists <- function(a, b) {
        lapply(1:(length(a)),function(x) rbind(a[[x]],b[[x]]))
      }
  }
  # Check settings to make sure compatible, give an error, warning, and/or fix settings if something is wrong
  {
    # Background prior events
      if ((background.prior.event == TRUE) & length(stages.order) == 2 ){
        #warning.f("No background color manipulation if only two stages.")
        background.prior.event=FALSE
      } else {}
    # Sub-mortality required mortality to be specified
      if ((allow.sub.stage.mortality==TRUE) & is.na(death.indicator)){
        stop("Sub-stage mortality requires specifying an indicator for death.")
      } else {}
  }
  # Generate main wide dataset
  {
    # Dataset cleanup
    {
      # Preserve original data
        events.long.orig <- events.long
      # Write list of columns to keep
        cols.to.keep <- c("ID","stage","date")
        if (!anyNA(groups.order)) {cols.to.keep <- c(cols.to.keep,"group")}
      # Strip out any columns that are not relevant
        events.long <- events.long[cols.to.keep]
      # Strip out any rows with missing stage, group, or other data
        events.long <- na.omit(events.long)
    }
    # Transform data into wide form
    {
      # Generate stage names
        stages <- as.data.frame(matrix(c(stages.order,paste(rep("stage"),as.character(seq(1:length(stages.order))),sep="."),1:length(stages.order)),ncol=3),stringsAsFactors=FALSE)
        colnames(stages) <- c("stage","stage.number","stage.index")
      # Replace stages with an index for stages
        events.long <- merge(events.long,stages,by="stage")
        events.long <- events.long[,!(names(events.long) %in% c("stage"))]
        events.long <- events.long[order(events.long$stage.number),]
        events.long$stage.index <- NULL
      # Generate a single "group" if groups are not specified
        groups.order.orig <- groups.order
        if (anyNA(groups.order)) {
          events.long$group <- "All observations"
          groups.order <- c("All observations")
        } else{}
      # Generate separated wide list of stage events
        events.wide <- events.long %>%
          tidyr::spread(.data$stage.number, date)
      # Preserve original dates and record which is the first actual recorded event in the data
        events.wide$first.recorded.stage <- NA
        for (i in 1:(nrow(stages))){
          events.wide[[paste0("orig.stage.",i)]] <- events.wide[[paste0("stage.",i)]]
          if(length(events.wide[!is.na(events.wide[[paste0("stage.",i)]]) & is.na(events.wide$first.recorded.stage),]$first.recorded.stage)!=0){
            events.wide[!is.na(events.wide[[paste0("stage.",i)]]) & is.na(events.wide$first.recorded.stage),]$first.recorded.stage  <- i
          } else {}
        }
      # Replace names with "date"
        names(events.wide) <- gsub("stage.", "date.stage.",names(events.wide))
    }
    # Fix ordering, so that completion of a later stage indicates completion of earlier stage (and skips), and mark where this happens
    for (k in 1:(length(stages.order))) {
        for (stage.index in seq(length(stages.order)-1,1,-1)){
          # Determine origin of stage events, and create variable for its origin
            former.date.name <- paste("date.stage.",stage.index,sep="")
            latter.date.name <- paste("date.stage.",stage.index+1,sep="")
            former.date.origin.index.name <- paste("stage.",stage.index,".origin.index",sep="")
            latter.date.origin.index.name <- paste("stage.",stage.index+1,".origin.index",sep="")
          # Determine whether or not current stage is eligible for change due to there existing a subsequent stage
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
            events.wide[[paste("date.stage.",stage.index,sep="")]] <- zoo::as.Date.numeric(events.wide[[paste("date.stage.",stage.index,sep="")]])
        }
        # Cleanup
          events.wide$temp_flag <- NULL
          for (i in 1:(length(stages.order)-1)){ events.wide[[paste("stage.",i,".origin.index",sep="")]] <- NULL}
      }
    # Generate time-based groups if time breaks are specified
    if (!anyNA(groups.date.breaks)){
      # Generate group names from breaks
        if (!anyNA(groups.date.breaks.labels)){
          groups.order <- groups.date.breaks.labels
        } else {
          groups.order <- c(paste0(as.character(groups.date.breaks[1])," to ",as.character(groups.date.breaks[2]-1)))
          for (i in 2:(length(groups.date.breaks)-1)){
            groups.order <- c(groups.order,paste0(as.character(groups.date.breaks[i])," to ",as.character(groups.date.breaks[i+1]-1)))
          }
        }
      # Generate grouping for each event except the last to determine whether the start event is within time breaks
        for (stage.index in 1:(length(stages.order)-1)){
          events.wide[[paste0("date.stage.",stage.index,".group")]] <- NA
          for (break.index in 1:(length(groups.date.breaks)-1)){
            events.wide[[paste0("date.stage.",stage.index,".group")]] <- ifelse(
              (events.wide[[paste0("date.stage.",stage.index)]] >= groups.date.breaks[break.index]) & (events.wide[[paste0("date.stage.",stage.index)]] < groups.date.breaks[break.index+1]),
              groups.order[break.index],
              events.wide[[paste0("date.stage.",stage.index,".group")]])
          }
          events.wide[[paste0("date.stage.",stage.index,".group")]] <- ifelse(
            is.na(events.wide[[paste0("date.stage.",stage.index,".group")]]),
            "No group membership",
            events.wide[[paste0("date.stage.",stage.index,".group")]])
        }
    } else{}
    # Determine dates of death and censorship, and add to dataset
    {
      # Merge in death events (if any)
      if (!is.na(death.indicator)){
        events.death <- subset(events.long.orig,events.long.orig$stage==death.indicator)
        events.death <- events.death[c("ID","date")]
        colnames(events.death) <- c("ID","date.death")
        events.wide <- merge(events.wide,events.death,by="ID",all.x = TRUE)
      }
      # Merge in censorship events (if any)
      if (!is.na(censorship.indicator)){
        events.censorship <- subset(events.long.orig,events.long.orig$stage==censorship.indicator)
        events.censorship <- events.censorship[c("ID","date")]
        colnames(events.censorship) <- c("ID","date.censorship")
        events.wide <- merge(events.wide,events.censorship,by="ID",all.x = TRUE)
      } else {
        events.wide$date.censorship <- as.Date(NA)
      }
      # Add in date of censorship based on end of data collection if option checked
        last.date <- max(events.long.orig$date)
        if (lubridate::is.Date(censorship.date)){
          events.wide$date.censorship <- ifelse(is.na(events.wide$date.censorship),
                                                censorship.date,events.wide$date.censorship)
        } else if (censorship.date == "lastdate") {
          events.wide$date.censorship <- ifelse(is.na(events.wide$date.censorship),
                                                last.date,events.wide$date.censorship)
        } else {
          stop("Incompatible date for censorship. Use a date format or specify 'lastdate'")
        }
        events.wide$date.censorship <- zoo::as.Date.numeric(events.wide$date.censorship)
      # Assume that if date of death exists and occurs after censorship event, censorship is false and replace with last date recorded
      if (!is.na(censorship.indicator) & (!is.na(death.indicator))){
        # Check to see if any individuals fall into this category, and warn user
        if (any((events.wide$date.censorship < events.wide$date.death) & (!is.na(events.wide$date.death)) & (!is.na(events.wide$date.censorship)))){
          warning.f(paste0("At least one censorship event takes place after death date. Changed censorship date(s) to last date recorded for those cases."))
          events.wide$date.censorship <- ifelse((events.wide$date.censorship < events.wide$date.death) & (!is.na(events.wide$date.death)) & (!is.na(events.wide$date.censorship)),
                                      last.date,events.wide$date.censorship)
          events.wide$date.censorship <- zoo::as.Date.numeric(events.wide$date.censorship)
        } else {}
      }
      # Deal with scenarios where censorship event occurs on or before stage event, by assuming censorship is false, and sets it to the last date
        for (i in 1:(length(stages.order))){
          #warning.f(paste0("At least one censorship event occurs before stage ",stages.order[i],". Changed censorship date(s) to last date recorded for those cases."))
          if (any((events.wide$date.censorship < events.wide[[paste0("date.stage.",i)]]) & !is.na(events.wide[[paste0("date.stage.",i)]]) & (!is.na(events.wide$date.censorship)))){
            events.wide$date.censorship <- ifelse((events.wide$date.censorship < events.wide[[paste0("date.stage.",i)]]) & !is.na(events.wide[[paste0("date.stage.",i)]]) & (!is.na(events.wide$date.censorship)),
                                    last.date,events.wide$date.censorship)
            events.wide$date.censorship <- zoo::as.Date.numeric(events.wide$date.censorship)
          }
        }
    }
    # Generate dataset for generating individual TTE data
    {
      TTE <- events.wide %>%
        select("ID",ends_with("group"))
    }
  }
  # Create stage-by-stage events and survival
  {
    # Function for generating survival curves
      stage.survival.run <- function(start.stage.index,group.index,run.type="main") {
        # Determine if individual is eligible for contributing person-time to transition frame
        {
          events.wide$eligible <- 1
          # Remove if individual has no date for the starting event
            events.wide$eligible <- ifelse(is.na(events.wide[[paste0("date.stage.",start.stage.index)]]),
                                           0,events.wide$eligible)
          # Determine correct group column to apply given settings
            if (!anyNA(groups.date.breaks)){
              reference.group <- paste0("date.stage.",start.stage.index,".group")
            } else {reference.group <- "group" }
          # Remove if in wrong group
            events.wide$eligible <- ifelse(events.wide$eligible == 1 & events.wide[[reference.group]]!=groups.order[group.index],
                                           0,events.wide$eligible)
          # Apply skip rules
            # No skipping (remove anyone who has already moved to a subsequent phase)
              if (skip.mode=="none"){
                events.wide$eligible <- ifelse(events.wide$eligible == 1 & !is.na(events.wide[[paste0("date.stage.",start.stage.index + 1)]]) & (events.wide[[paste0("date.stage.",start.stage.index + 1)]] <= events.wide[[paste0("date.stage.",start.stage.index)]]),
                                               0,events.wide$eligible)
              } else {}
            # Internal skips only (remove people for whom the recorded event was not the first recorded event)
              if (skip.mode=="internal"){
                events.wide$eligible <- ifelse(events.wide$first.recorded.stage > start.stage.index,
                                               0,events.wide$eligible)
              } else {}
        }
        # Apply survival functions for events
        {
          # Main stage transitions
          if (run.type == "main") {
            # Function for stats for all stage transitions
              surv.stats.main.stages <- function(df.surv,start.stage.index,end.stage.index,group.index){
              # Assign start date
                df.surv$start.date <- df.surv[[paste0("date.stage.",start.stage.index)]]
              # Determine time to event and censorship
                  df.surv$endpoint <- 0
                  df.surv$end.date <- as.Date(NA)
                  df.surv$censorship <- NA
                # Option 1: Main next stage event occurs
                  df.surv$endpoint <- ifelse(!is.na(df.surv[[paste0("date.stage.",end.stage.index)]]),
                                             1,df.surv$endpoint)
                  df.surv[df.surv$endpoint==1,]$end.date <- df.surv[df.surv$endpoint==1,][[paste0("date.stage.",end.stage.index)]]
                  df.surv[df.surv$endpoint==1,]$censorship <- rep(0,sum(df.surv$endpoint==1))
                # Option 2: If endpoint is never reached, check if death occurs
                  # if (!is.na(death.indicator)){
                  #   df.surv$endpoint <- ifelse(df.surv$endpoint==0 & !is.na(df.surv$date.death),
                  #                                                           2,df.surv$endpoint)
                  # }
                  # df.surv[df.surv$endpoint==2,]$end.date <- df.surv[df.surv$endpoint==2,]$start.date + time.horizon
                  # df.surv[df.surv$endpoint==2,]$censorship <- rep(1,sum(df.surv$endpoint==2))
                # Option 3: Otherwise, use the censorship date
                  df.surv$endpoint <- ifelse(df.surv$endpoint==0,3,df.surv$endpoint)
                  df.surv[df.surv$endpoint==3,]$end.date <- df.surv[df.surv$endpoint==3,]$date.censorship
                  df.surv[df.surv$endpoint==3,]$censorship <- rep(1,sum(df.surv$endpoint==3))
                # Revise events which occur after the time horizon of interest to the time horizon, and change to censorship
                  df.surv$endpoint <- ifelse((df.surv$end.date - df.surv$start.date) > time.horizon,
                                             4,df.surv$endpoint)
                  df.surv[df.surv$endpoint==4,]$end.date <- df.surv[df.surv$endpoint==4,]$start.date + time.horizon
                  df.surv[df.surv$endpoint==4,]$censorship <- rep(1,sum(df.surv$endpoint==4))
                # Get time to event
                  df.surv$time.to.event <- df.surv$end.date-df.surv$start.date
                # Generate TTE dataset output and merge with larger dataset
                  colname.TTE <- paste0("time.stage.",start.stage.index,".to.",end.stage.index)
                  colname.censorship <- paste0("censorship.stage.",start.stage.index,".to.",end.stage.index)
                  df.surv[[colname.TTE]] <- df.surv$time.to.event
                  df.surv[[colname.censorship]] <- df.surv$censorship
                  TTE <- df.surv[c("ID","time.to.event","censorship")]
                  TTE$transition.name <- colname.TTE
                  TTE$censorship.name <- colname.censorship
                # Pull survival curves
                  surv.data <- surv.data.stats(df.surv$time.to.event,df.surv$censorship)
                # Generate quantile TTEs
                  quantile.TTE.q <- t(unlist(surv.data$TTE.quantile[1]))
                  quantile.TTE.q.lb <- t(unlist(surv.data$TTE.quantile[2]))
                  quantile.TTE.q.ub <- t(unlist(surv.data$TTE.quantile[3]))
                  colnames(quantile.TTE.q) <- paste0("q.",colnames(quantile.TTE.q))
                  colnames(quantile.TTE.q.lb) <- paste0("q.",colnames(quantile.TTE.q.lb))
                  colnames(quantile.TTE.q.ub) <- paste0("q.",colnames(quantile.TTE.q.ub))
                  quantile.TTE <- as.data.frame(cbind(start.stage.index,end.stage.index,group.index,quantile.TTE.q,quantile.TTE.q.lb,quantile.TTE.q.ub))
                # Generate main survival curve outputs
                  surv.curve <- surv.data$surv.curve
                  surv.curve$start.stage.index <- start.stage.index
                  surv.curve$end.stage.index <- end.stage.index
                  surv.curve$group.index <- group.index
                # Return output
                  return(list("surv.curve" = surv.curve,"TTE.ind" = TTE,"quantile.TTE" = quantile.TTE))
              }
            # Determine which stages to run
              if (allow.sub.stages==TRUE) {
                end.stage.index.list <- (start.stage.index+1):(length(stages.order))
              } else {
                end.stage.index.list <- (start.stage.index+1):(start.stage.index+1)
              }
            # Run main stages
              # surv <- do.call("rbind",lapply(end.stage.index.list,function(x)
              #   surv.stats.main.stages(events.wide[events.wide$eligible==1,],start.stage.index,x,group.index)))
              outputs <- Reduce(rbind.lists, lapply(end.stage.index.list,function(x)
                surv.stats.main.stages(events.wide[events.wide$eligible==1,],start.stage.index,x,group.index)))
              surv <- outputs[[1]]
              TTE.ind <- outputs[[2]]
              quantile.TTE <- outputs[[3]]

          } else {}
          # Death stage transitions
          if (run.type=="death") {
            # Function for mortality stages and substages
              surv.stats.mortality <- function(df.surv,start.stage.index,end.stage.index,group.index){
            # Assign start date
              df.surv$start.date <- df.surv[[paste0("date.stage.",start.stage.index)]]
              df.surv$end.stage.date <- df.surv[[paste0("date.stage.",end.stage.index)]]
              df.surv$prev.stage.date <- df.surv[[paste0("date.stage.",end.stage.index-1)]]
            # Determine time to event and censorship
                df.surv$endpoint <- 0
                df.surv$end.date <- as.Date(NA)
                df.surv$censorship <- NA
              # Option 1: Mortality occurs after event before end stage, and end date never occurs (died before end date)
                df.surv$endpoint <- ifelse(!is.na(df.surv$date.death) & is.na(df.surv$end.stage.date) & !is.na(df.surv$prev.stage.date) &
                                             (df.surv$date.death > df.surv$start.date),
                                           1,df.surv$endpoint)
                df.surv[df.surv$endpoint==1,]$end.date <- df.surv[df.surv$endpoint==1,]$date.death
                df.surv[df.surv$endpoint==1,]$censorship <- rep(0,sum(df.surv$endpoint==1))
              # Option 2: Censorship occurs before time horizon reached
                df.surv$endpoint <- ifelse(df.surv$endpoint==0 & (df.surv$date.censorship - df.surv$start.date) < time.horizon,
                                           2,df.surv$endpoint)
                df.surv[df.surv$endpoint==2,]$end.date <- df.surv[df.surv$endpoint==2,]$date.censorship
                df.surv[df.surv$endpoint==2,]$censorship <- rep(1,sum(df.surv$endpoint==2))
              # Option 3: Otherwise, use the censorship date
                df.surv$endpoint <- ifelse(df.surv$endpoint==0,3,df.surv$endpoint)
                df.surv[df.surv$endpoint==3,]$end.date <- df.surv[df.surv$endpoint==3,]$date.censorship
                df.surv[df.surv$endpoint==3,]$censorship <- rep(1,sum(df.surv$endpoint==3))
              # Revise events which occur after the time horizon of interest to the time horizon, and change to censorship
                df.surv$endpoint <- ifelse((df.surv$end.date - df.surv$start.date) > time.horizon,
                                           4,df.surv$endpoint)
                df.surv[df.surv$endpoint==4,]$end.date <- df.surv[df.surv$endpoint==4,]$start.date + time.horizon
                df.surv[df.surv$endpoint==4,]$censorship <- rep(1,sum(df.surv$endpoint==4))
              # Get time to event
                df.surv$time.to.event <- df.surv$end.date-df.surv$start.date
              # Generate TTE dataset output and merge with larger dataset
                colname.TTE <- paste0("time.stage.",start.stage.index,".to.death.",end.stage.index-1,".to.",end.stage.index)
                colname.censorship <- paste0("censorship.stage.",start.stage.index,".to.death.",end.stage.index-1,".to.",end.stage.index)
                df.surv[[colname.TTE]] <- df.surv$time.to.event
                df.surv[[colname.censorship]] <- df.surv$censorship
                TTE <- df.surv[c("ID","time.to.event","censorship")]
                TTE$transition.name <- colname.TTE
                TTE$censorship.name <- colname.censorship
              # Pull survival curves
                output <- surv.data.stats(df.surv$time.to.event,df.surv$censorship)$surv.curve
                output$start.stage.index <- start.stage.index
                output$end.stage.index <- end.stage.index
                output$reference.stage.index <- end.stage.index - 1
                output$group.index <- group.index
              # Check if this is the last stage, and if so, run one additional mortality check
                if (end.stage.index == length(stages.order) & allow.sub.stage.mortality==TRUE){
                  # Assign start date
                    df.surv$start.date <- df.surv[[paste0("date.stage.",start.stage.index)]]
                    df.surv$end.stage.date <- df.surv[[paste0("date.stage.",end.stage.index)]]
                  # Determine time to event and censorship
                    df.surv$endpoint <- 0
                    df.surv$end.date <- as.Date(NA)
                    df.surv$censorship <- NA
                  # Option 1: Mortality occurs after the end stage date
                    df.surv$endpoint <- ifelse(!is.na(df.surv$date.death) & !is.na(df.surv$end.stage.date) & (df.surv$date.death > df.surv$end.stage.date),
                                               1,df.surv$endpoint)
                    df.surv[df.surv$endpoint==1,]$end.date <- df.surv[df.surv$endpoint==1,]$date.death
                    df.surv[df.surv$endpoint==1,]$censorship <- rep(0,sum(df.surv$endpoint==1))
                  # Option 2: Censorship occurs before time horizon reached
                    df.surv$endpoint <- ifelse(df.surv$endpoint==0 & (df.surv$date.censorship - df.surv$start.date) < time.horizon,
                                               2,df.surv$endpoint)
                    df.surv[df.surv$endpoint==2,]$end.date <- df.surv[df.surv$endpoint==2,]$date.censorship
                    df.surv[df.surv$endpoint==2,]$censorship <- rep(1,sum(df.surv$endpoint==2))
                  # Option 3: Otherwise, use the censorship date
                    df.surv$endpoint <- ifelse(df.surv$endpoint==0,3,df.surv$endpoint)
                    df.surv[df.surv$endpoint==3,]$end.date <- df.surv[df.surv$endpoint==3,]$date.censorship
                    df.surv[df.surv$endpoint==3,]$censorship <- rep(1,sum(df.surv$endpoint==3))
                  # Revise events which occur after the time horizon of interest to the time horizon, and change to censorship
                    df.surv$endpoint <- ifelse((df.surv$end.date - df.surv$start.date) > time.horizon,
                                               4,df.surv$endpoint)
                    df.surv[df.surv$endpoint==4,]$end.date <- df.surv[df.surv$endpoint==4,]$start.date + time.horizon
                    df.surv[df.surv$endpoint==4,]$censorship <- rep(1,sum(df.surv$endpoint==4))
                  # Get time to event
                    df.surv$time.to.event <- df.surv$end.date-df.surv$start.date
                  # Generate TTE dataset output and merge with larger dataset
                    colname.TTE <- paste0("time.stage.",start.stage.index,".to.death.",end.stage.index,".to.",end.stage.index+1)
                    colname.censorship <- paste0("censorship.stage.",start.stage.index,".to.death.",end.stage.index,".to.",end.stage.index+1)
                    df.surv[[colname.TTE]] <- df.surv$time.to.event
                    df.surv[[colname.censorship]] <- df.surv$censorship
                    TTE <- df.surv[c("ID","time.to.event","censorship")]
                    TTE$transition.name <- colname.TTE
                    TTE$censorship.name <- colname.censorship
                  # Pull survival curves
                    output2 <- surv.data.stats(df.surv$time.to.event,df.surv$censorship)$surv.curve
                    output2$start.stage.index <- start.stage.index
                    output2$end.stage.index <- end.stage.index + 1
                    output2$reference.stage.index <- end.stage.index
                    output2$group.index <- group.index
                  # Merge back in with main mortality output
                    output <- rbind(output,output2)
                }
              # Return output
                return(list("surv.curve" = output, "TTE.ind" = TTE))
            }
            # Determine which stages to run
              if (allow.sub.stage.mortality==TRUE) {
                end.stage.index.list <- (start.stage.index+1):(length(stages.order))
              } else {
                end.stage.index.list <- (start.stage.index+1):(start.stage.index+1)
              }
            # Run mortality stages
              outputs <- Reduce(rbind.lists, lapply(end.stage.index.list,function(x)
                surv.stats.mortality(events.wide[events.wide$eligible==1,],start.stage.index,x,group.index)))
              surv <- outputs[[1]]
              TTE.ind <- outputs[[2]]
            } else {}
          # Transient status transitions
          if (run.type=="transient") {
            # Function for transient stages and substages
              surv.stats.transient <- function(events.wide.ts, ts.indicator,ts.gap.time,
                                                  start.stage.index = 1,end.stage.index = NA,group.index=NA) {
              # Generate initial data
              {
                # Get events wide and long versions for all events, with the given set of eligible IDs
                  df.ts <- events.long.orig[events.long.orig$stage==ts.indicator & (events.long.orig$ID %in% events.wide.ts$ID),][c("ID","date","stage")]
                # Add stage start events to df.ts (so that time is measured from start of stage)
                  df.ts.date.stage <- stats::na.omit(events.wide.ts[c("ID",paste0("date.stage.",start.stage.index))])
                  df.ts.date.stage$date <- df.ts.date.stage[[paste0("date.stage.",start.stage.index)]]
                  df.ts.date.stage[[paste0("date.stage.",start.stage.index)]] <- NULL
                  df.ts.date.stage$stage <- "Stage start"
                  df.ts <- rbind(df.ts,df.ts.date.stage)
                  rm(df.ts.date.stage)
                  df.ts <- df.ts %>%
                    dplyr::arrange(.data$ID,date)
                # Merge in relevant data from events.wide
                  cols.to.keep <- c("ID",paste0("date.stage.",start.stage.index),"date.censorship")
                  if (!is.na(end.stage.index)){
                    cols.to.keep <- c(cols.to.keep,paste0("date.stage.",end.stage.index))
                  } else {}
                  if(!is.na(death.indicator)){
                    cols.to.keep <- c(cols.to.keep,"date.death")
                  } else {}
                  df.ts <- merge(df.ts,events.wide.ts[cols.to.keep],
                                 by="ID",all.x=TRUE,all.y=FALSE)
                # Trim to date bounds
                  df.ts <- df.ts[df.ts$date>=df.ts[[paste0("date.stage.",start.stage.index)]],]
                  df.ts <- df.ts[df.ts$date<=df.ts$date.censorship,]
                  if (!is.na(end.stage.index)){
                    df.ts <- df.ts[df.ts$date<df.ts[[paste0("date.stage.",end.stage.index)]] | is.na(df.ts[[paste0("date.stage.",end.stage.index)]]),]
                  } else {}
                  if (!is.na(death.indicator)){
                    df.ts <- df.ts[df.ts$date<df.ts$date.death | is.na(df.ts$date.death),]
                    # Keep death events only if it occured within staging bounds
                      if (is.na(end.stage.index)){
                        df.ts$date.death <- ifelse(df.ts$date.death>df.ts[[paste0("date.stage.",start.stage.index)]],
                                                 NA,df.ts$date.death)
                      }
                  } else {}
                # Determine final "off" date and type (whichever happens first: ts.gap, next stage, death, or censorship)
                    df.ts$date.end <- as.Date(NA)
                  # If both death and next stage events are available
                    if (!is.na(death.indicator) & !is.na(end.stage.index)){
                      df.ts$date.end <- pmin(df.ts$date.death,df.ts[[paste0("date.stage.",end.stage.index)]],df.ts$date+ts.gap.time, na.rm = TRUE)
                    } else {}
                  # If only next stage events are available
                    if (!is.na(end.stage.index)){
                      df.ts$date.end <- pmin(df.ts[[paste0("date.stage.",end.stage.index)]],df.ts$date+ts.gap.time,na.rm=TRUE)
                    } else {}
                  # If only death is available
                    if (!is.na(death.indicator) & is.na(end.stage.index)){
                      df.ts$date.end <- pmin(df.ts$date.death,df.ts$date+ts.gap.time,na.rm=TRUE)
                    } else {}
                  # Revise so that death
                  # If nothing else available, end at censorship date
                    df.ts$date.end <- zoo::as.Date.numeric(ifelse(is.na(df.ts$date.end),df.ts$date.censorship,df.ts$date.end))
              }
              # Identify gaps and consolidate periods
              {
                # Determine if there is a gap of ts.gap.time days or more
                  df.ts <- df.ts %>%
                    dplyr::arrange(.data$ID, date) %>%
                    dplyr::group_by(.data$ID) %>%
                    dplyr::mutate(gap = dplyr::lead(date) - date>ts.gap.time)
                # Create an "off" event when there is a gap
                  df.ts$date.off <- NA
                  df.ts$date.off <- zoo::as.Date.numeric(ifelse(df.ts$gap,df.ts$date+ts.gap.time,df.ts$date.off))
                  df.ts$gap <- NULL
                # Replace dates backward
                  # Find the number of times needed to run
                    df.ts <- df.ts %>% dplyr::group_by(.data$ID) %>% dplyr::mutate(count = row_number())
                    n.run <- max(df.ts$count)
                    df.ts$count <- NULL
                  # Give an end date if it's the last one
                    df.ts$date.off <- zoo::as.Date.numeric(ifelse(dplyr::lead(df.ts$ID)!=df.ts$ID | is.na(dplyr::lead(df.ts$ID)),
                                             df.ts$date.end,df.ts$date.off))
                    if (!is.na(death.indicator)){
                      df.ts$date.off <- zoo::as.Date.numeric(ifelse((dplyr::lead(df.ts$ID)!=df.ts$ID | is.na(dplyr::lead(df.ts$ID))) & df.ts$date.death < df.ts$date.off  & !is.na(df.ts$date.death),
                                               df.ts$date.death,df.ts$date.off))
                    } else {}
                  # Trace backward until beginning of period
                    for (i in 1:n.run){
                      df.ts$date <- zoo::as.Date.numeric(ifelse(!is.na(df.ts$date.off) & !is.na(dplyr::lag(df.ts$ID)) & dplyr::lag(df.ts$ID)==df.ts$ID & is.na(dplyr::lag(df.ts$date.off)),
                                           dplyr::lag(df.ts$date),df.ts$date))
                      df.ts <- df.ts[!(dplyr::lead(df.ts$ID)==df.ts$ID & dplyr::lead(df.ts$date)==df.ts$date & !is.na(dplyr::lead(df.ts$date.off))),]
                    }
              }
              # Generate repeat #s and weights
                # Create count repeat #s
                  df.ts <- df.ts %>%
                    dplyr::arrange(.data$ID,date)%>%
                    dplyr::group_by(.data$ID) %>% dplyr::mutate(count = row_number())
                  max.count <- max(df.ts$count)
              # Merge back into the wide dataset under appropriate labels
              {
                for (i in 1:max.count){
                  df.ts.temp <- df.ts[df.ts$count==i,][c("ID","date","date.off")]
                  df.ts.temp[[paste0("date.on.",i)]] <- df.ts.temp$date
                  df.ts.temp[[paste0("date.off.",i)]] <- df.ts.temp$date.off
                  df.ts.temp$date <- NULL
                  df.ts.temp$date.off <- NULL
                  events.wide.ts <- merge(events.wide.ts,df.ts.temp,by="ID",all.x = TRUE)
                  rm(df.ts.temp)
                }
              }
              # Function for survival stats
                surv.stats.transient <- function(df.surv,start.stage.index,end.date.colname){
                # Assign start date
                  df.surv$start.date <- df.surv[[paste0("date.stage.",start.stage.index)]]
                  df.surv$event.date <- df.surv[[end.date.colname]]
                # Determine time to event and censorship
                    df.surv$endpoint <- 0
                    df.surv$end.date <- as.Date(NA)
                    df.surv$censorship <- NA
                  # Option 1: Main on/off event occurs
                    df.surv$endpoint <- ifelse(!is.na(df.surv$event.date),
                                               1,df.surv$endpoint)
                    df.surv[df.surv$endpoint==1,]$end.date <- df.surv[df.surv$endpoint==1,]$event.date
                    df.surv[df.surv$endpoint==1,]$censorship <- rep(0,sum(df.surv$endpoint==1))
                  # Option 2: If end date is the censorship date, censor on that date
                    df.surv$endpoint <- ifelse(df.surv$event.date >= df.surv$date.censorship & !is.na(df.surv$event.date),
                                               2,df.surv$endpoint)
                    df.surv[df.surv$endpoint==2,]$end.date <- df.surv[df.surv$endpoint==2,]$date.censorship
                    df.surv[df.surv$endpoint==2,]$censorship <- rep(1,sum(df.surv$endpoint==2))
                  # Option 3: Otherwise (if missing), use the censorship date
                    df.surv$endpoint <- ifelse(df.surv$endpoint==0,3,df.surv$endpoint)
                    df.surv[df.surv$endpoint==3,]$end.date <- df.surv[df.surv$endpoint==3,]$date.censorship
                    df.surv[df.surv$endpoint==3,]$censorship <- rep(1,sum(df.surv$endpoint==3))
                  # Revise events which occur after the time horizon of interest to the time horizon, and change to censorship
                    df.surv$endpoint <- ifelse((df.surv$end.date - df.surv$start.date) > time.horizon,
                                               4,df.surv$endpoint)
                    df.surv[df.surv$endpoint==4,]$end.date <- df.surv[df.surv$endpoint==4,]$start.date + time.horizon
                    df.surv[df.surv$endpoint==4,]$censorship <- rep(1,sum(df.surv$endpoint==4))
                  # Get time to event
                    df.surv$time.to.event <- df.surv$end.date-df.surv$start.date
                  # Pull survival curves
                    #output <- surv.data.stats(df.surv$time.to.event,df.surv$censorship,df.surv$weight)$surv.curve
                    output <- surv.data.stats(df.surv$time.to.event,df.surv$censorship)$surv.curve
                  # Return output
                    return(output)
                }
              # Generate on/off curves and combine
              {
                # Start off with the first event
                  trans.surv.on <- surv.stats.transient(events.wide.ts,start.stage.index,paste0("date.on.",1))[c("surv.p","surv.time")]
                  trans.surv.off <- surv.stats.transient(events.wide.ts,start.stage.index,paste0("date.off.",1))[c("surv.p","surv.time")]
                # Add additional events as necessary
                  if (max(df.ts$count)>1){
                    for (i in 2:max(df.ts$count)){
                      trans.surv.on <- relative.curve(trans.surv.on,surv.stats.transient(events.wide.ts,start.stage.index,paste0("date.on.",i))[c("surv.p","surv.time")])[c("surv.p","surv.time")]
                      trans.surv.off <- relative.curve(trans.surv.off,surv.stats.transient(events.wide.ts,start.stage.index,paste0("date.off.",i))[c("surv.p","surv.time")])[c("surv.p","surv.time")]
                    }
                  } else {}
                # Combine into one curve
                  output <- relative.curve(trans.surv.on,trans.surv.off,"under")[c("surv.p","surv.time")]
              }
              # Output data
              {
                output$start.stage.index <- start.stage.index
                output$end.stage.index <- end.stage.index
                output$group.index <- group.index
                output$ts.indicator <- ts.indicator
                output$ts.gap.time <- ts.gap.time
              # Return output
                return(output)
              }
            }
            # Find each combination of transient statuses which start with the given start stage
              ts.settings.rows <- which(ts.start.stage %in% start.stage.index)
            # Run transient stages
              surv <- do.call("rbind",lapply(ts.settings.rows,function(x)
                surv.stats.transient(events.wide.ts = events.wide[events.wide$eligible==1,],
                                     ts.indicator = ts.indicator[x],ts.gap.time = ts.gap.time[x],
                                     start.stage.index = start.stage.index,end.stage.index = ts.end.stage[x],group.index=group.index)))
            } else {}
        }
        # Change stages and groups to factors for easier charting
          surv$start.stage.factor <- ordered(surv$start.stage.index,levels=c(1:(length(stages.order)-1)),
                                                      labels = stages.order[1:(length(stages.order)-1)])
          surv$end.stage.factor <- ordered(surv$end.stage.index,levels=c(2:(length(stages.order))),
                                                    labels = stages.order[2:(length(stages.order))])
          surv$group.factor <- ordered(surv$group.index,levels=c(1:(length(groups.order))),
                                                labels = groups.order[1:(length(groups.order))])
        # Return data
          if (run.type=="main"){
            quantile.TTE$start.stage.factor <- ordered(quantile.TTE$start.stage.index,levels=c(1:(length(stages.order)-1)),
                                                      labels = stages.order[1:(length(stages.order)-1)])
            quantile.TTE$end.stage.factor <- ordered(quantile.TTE$end.stage.index,levels=c(2:(length(stages.order))),
                                                      labels = stages.order[2:(length(stages.order))])
            quantile.TTE$group.factor <- ordered(quantile.TTE$group.index,levels=c(1:(length(groups.order))),
                                                  labels = groups.order[1:(length(groups.order))])
            return(list("surv" = surv,"TTE.ind" = TTE.ind,"quantile.TTE" = quantile.TTE))
          } else if (run.type=="death"){
            return(list("surv" = surv,"TTE.ind" = TTE.ind))
          } else { return(surv)}
      }
    # Generate survival functions for all combinatorials of events
      comb <- expand.grid(stages = 1:(length(stages.order)-1),groups = 1:(length(groups.order)))
    # Run all specified transitions to generate main datasets
      #surv.main <- do.call("rbind",mapply(stage.survival.run,start.stage.index = comb$stages,group.index = comb$groups,SIMPLIFY = FALSE))
      outputs <- Reduce(rbind.lists, mapply(stage.survival.run,start.stage.index = comb$stages,group.index = comb$groups,SIMPLIFY = FALSE))
      surv.main <- outputs[[1]]
      TTE.ind <- outputs[[2]]
      quantile.TTE <- outputs[[3]]
      if (!is.na(death.indicator)){
        comb$run.type <- "death"
        outputs <- Reduce(rbind.lists, mapply(stage.survival.run,start.stage.index = comb$stages,group.index = comb$groups,run.type = comb$run.type,SIMPLIFY = FALSE))
        surv.death <- outputs[[1]]
        TTE.ind.death <- outputs[[2]]
      } else {}
      if (!anyNA(ts.indicator)){
        comb$run.type <- "transient"
        surv.transient <- do.call("rbind",mapply(stage.survival.run,start.stage.index = comb$stages,group.index = comb$groups,run.type = comb$run.type,SIMPLIFY = FALSE))
      } else {}
    # Modify TTE outputs to be in wide form
      # Main events
        TTE.times <- TTE.ind[c("ID","time.to.event","transition.name")] %>%
        tidyr::spread(.data$transition.name,c(.data$time.to.event))
        TTE.censorship <- TTE.ind[c("ID","censorship","censorship.name")] %>%
        tidyr::spread(.data$censorship.name,c(.data$censorship))
        TTE.ind <- merge(TTE.times,TTE.censorship,by="ID",all.x=TRUE)
        TTE.ind <- merge(TTE.ind,events.wide[c("ID","group")],by="ID")
        if (!anyNA(groups.date.breaks)){
          TTE.ind <- merge(TTE.ind,events.wide[c("ID",paste0("date.stage.",1:(length(stages.order)-1),".group"))])
        }
      # Death events
      if (!is.na(death.indicator)){
        TTE.times <- TTE.ind.death[c("ID","time.to.event","transition.name")] %>%
        tidyr::spread(.data$transition.name,c(.data$time.to.event))
        TTE.censorship <- TTE.ind.death[c("ID","censorship","censorship.name")] %>%
        tidyr::spread(.data$censorship.name,c(.data$censorship))
        TTE.ind.death <- merge(TTE.times,TTE.censorship,by="ID",all.x=TRUE)
        TTE.ind.death <- merge(TTE.ind.death,events.wide[c("ID","group")],by="ID")
        if (!anyNA(groups.date.breaks)){
          TTE.ind <- merge(TTE.ind,events.wide[c("ID",paste0("date.stage.",1:(length(stages.order)-1),".group"))])
        }
      }
  }
  # Generate chart graphics
  {
    # Data manipulation
    {
      # Main events
      {
        # Generate data
          surv.main.chart <- surv.main
        # Add background area in for foundation event
          if (background.prior.event==TRUE){
            surv.main.chart.temp <- surv.main.chart
            surv.main.chart.temp$surv.p <- 0
            surv.main.chart.temp$surv.time <- 0
            surv.main.chart.temp$surv.p.LB <- 0
            surv.main.chart.temp$surv.p.UB <- 0
            surv.main.chart.temp$surv.p.atrisk <- 1
            surv.main.chart.temp$surv.n.atrisk <- 0
            surv.main.chart.temp$start.stage.factor <- as.character(surv.main.chart.temp$end.stage.factor)
            surv.main.chart.temp$start.stage.index <- surv.main.chart.temp$end.stage.index
            surv.main.chart.temp <- unique(surv.main.chart.temp)
            surv.main.chart.temp <- surv.main.chart.temp[surv.main.chart.temp$start.stage.index!=length(stages.order),]
            surv.main.chart.temp$start.stage.factor <- ordered(surv.main.chart.temp$start.stage.factor,labels=levels(surv.main.chart$start.stage.factor),levels=levels(surv.main.chart$start.stage.factor))
            surv.main.chart.temp$surv.p <- 1
            surv.main.chart <- rbind(surv.main.chart,surv.main.chart.temp)
            surv.main.chart.temp$surv.time <- time.horizon
            surv.main.chart <- rbind(surv.main.chart,surv.main.chart.temp)
            rm(surv.main.chart.temp)
          } else {}
        # Rearrange and sort for drawing
          surv.main.chart <- surv.main.chart %>%
            dplyr::arrange(.data$group.factor,.data$start.stage.factor,.data$end.stage.factor,.data$surv.time)
        # Change time to years
          surv.main.chart$surv.time = surv.main.chart$surv.time/365
        # Generate risk pool dataset
          surv.main.chart.risk.pool <- surv.main.chart[surv.main.chart$end.stage.index == (surv.main.chart$start.stage.index +1) ,]
        # Add a dummy factor for legend
          surv.main.chart.risk.pool$event.factor <- factor("Risk pool",labels = c("Risk pool"),levels = c("Risk pool"))
      }
      # Death events
      if (!is.na(death.indicator)){
        # Generate data
          surv.death.chart <- surv.death
        # Trim data if only a single panel is wanted
          if (chart.mode == "first transition"){
            surv.death.chart <- surv.death.chart[surv.death.chart$start.stage.index == 1,]
          } else {}
        # Rearrange and sort for drawing
          surv.death.chart <- surv.death.chart %>%
            dplyr::arrange(.data$group.factor,.data$start.stage.factor,.data$end.stage.factor,.data$surv.time)
        # Change time to years
          surv.death.chart$surv.time = surv.death.chart$surv.time/365
        # Revise and add mortality to relevent reference frames
          gen.substage.mort.df <- function(start.stage.index,reference.stage.index,group.index){
            input <- surv.death.chart[surv.death.chart$start.stage.index==start.stage.index & surv.death.chart$reference.stage.index==reference.stage.index & surv.death.chart$group.index==group.index,]
            if (start.stage.index==reference.stage.index){
              reference <- surv.main.chart[surv.main.chart$start.stage.index==start.stage.index & surv.main.chart$end.stage.index==reference.stage.index + 1 & surv.main.chart$group.index==group.index,]
              reference <- reference[(1:2),]
              reference$surv.time <- c(0,time.horizon)
              reference$surv.p <- c(1,1)
            } else {
              reference <- surv.main.chart[surv.main.chart$start.stage.index==start.stage.index & surv.main.chart$end.stage.index==(reference.stage.index) & surv.main.chart$group.index==group.index,]
            }
            df <- relative.curve(reference,input,over.under="under")
            df$start.stage.index <- start.stage.index
            df$start.stage.factor <- reference$start.stage.factor[1]
            df$end.stage.index <- reference.stage.index
            df$end.stage.factor <- reference$end.stage.factor[1]
            df$group.index <- group.index
            df$group.factor <- reference$group.factor[1]
            df$reference.stage.index <- reference.stage.index
            return(df)
          }
          df.substage.stages <- unique(surv.death.chart[c("start.stage.index","reference.stage.index","group.index")])
          surv.substage.death.chart <- do.call("rbind",lapply(1:nrow(df.substage.stages),function(x) gen.substage.mort.df(df.substage.stages$start.stage.index[x],df.substage.stages$reference.stage.index[x],df.substage.stages$group.index[x])))
      } else {}
      # Transient statuses
      if (!anyNA(ts.indicator)){
        # Generate data
          surv.transient.chart <- surv.transient
        # Trim data if only a single panel is wanted
          if (chart.mode == "first transition"){
            surv.transient.chart <- surv.transient.chart[surv.transient.chart$start.stage.index == 1,]
          } else {}
        # Rearrange and sort for drawing
          surv.transient.chart <- surv.transient.chart %>%
            dplyr::arrange(.data$group.factor,.data$start.stage.factor,.data$end.stage.factor,.data$surv.time)
        # Change time to years
          surv.transient.chart$surv.time = surv.transient.chart$surv.time/365
        # Revise and add transient status to relevent reference frames
          gen.substage.transient.df <- function(start.stage.index,reference.stage.index,group.index){
            # Get main curve for transient event
              input <- surv.transient.chart[surv.transient.chart$start.stage.index==start.stage.index & surv.transient.chart$end.stage.index==reference.stage.index & surv.transient.chart$group.index==group.index,]

              reference <- surv.main.chart[surv.main.chart$start.stage.index==start.stage.index & surv.main.chart$end.stage.index==(reference.stage.index) & surv.main.chart$group.index==group.index,]
              df <- relative.curve(reference,input,over.under="over")
              df$start.stage.index <- start.stage.index
              df$start.stage.factor <- reference$start.stage.factor[1]
              df$end.stage.index <- reference.stage.index
              df$end.stage.factor <- reference$end.stage.factor[1]
              df$group.index <- group.index
              df$group.factor <- reference$group.factor[1]
              df$reference.stage.index <- reference.stage.index
              return(df)
          }
          df.substage.stages <- unique(surv.transient.chart[c("start.stage.index","group.index")])
          df.substage.transient <- do.call("rbind",lapply(1:nrow(df.substage.stages),function(x) gen.substage.transient.df(df.substage.stages$start.stage.index[x],df.substage.stages$start.stage.index[x]+1,df.substage.stages$group.index[x])))
      } else {}
      # Data modifications due to chart options
      {
        # Trim data if only a single panel is wanted
          if (chart.mode == "first transition"){
            surv.main.chart <- surv.main.chart[surv.main.chart$start.stage.index == 1,]
            if (!is.na(death.indicator)){
              surv.death.chart <- surv.death.chart[surv.death.chart$start.stage.index == 1,]
            } else {}
            if (!anyNA(ts.indicator)){
              surv.transient.chart <- surv.transient.chart[surv.transient.chart$start.stage.index == 1,]
            } else {}
          } else {}
        # Shift sub-mortality to top, and drop sub-transitions accordingly
          if (sub.stage.mortality.mode=="shifted"){
            if (chart.mode == "first transition"){
              max.start.stage <- 1
            } else {
              max.start.stage <- length(stages.order)-1
            }
            # Generate shifted death
              surv.substage.death.chart.temp <- surv.substage.death.chart[0,]
              for (group.index in 1:length(groups.order)){
                for (start.stage.index in 1:(max.start.stage)){
                #for (start.stage.index in 1:(length(stages.order)-1)){
                    reference <- surv.substage.death.chart[surv.substage.death.chart$start.stage.index==start.stage.index &
                                                             surv.substage.death.chart$end.stage.index==start.stage.index &
                                                             surv.substage.death.chart$group.index==group.index,]
                    surv.substage.death.chart.temp <- rbind(surv.substage.death.chart.temp,reference)
                  for (end.stage.index in (start.stage.index+1):length(stages.order)){
                    reference <- surv.substage.death.chart.temp[surv.substage.death.chart.temp$start.stage.index==start.stage.index &
                                                             surv.substage.death.chart.temp$end.stage.index==end.stage.index-1 &
                                                             surv.substage.death.chart.temp$group.index==group.index,]
                    input <- surv.substage.death.chart[surv.substage.death.chart$start.stage.index==start.stage.index &
                                                             surv.substage.death.chart$end.stage.index==end.stage.index &
                                                             surv.substage.death.chart$group.index==group.index,]
                    input$surv.p <- input$surv.p.reference-input$surv.p
                    combined <- relative.curve(reference[c("surv.p","surv.time")],input[c("surv.p","surv.time")],"under")
                    combined$com <- 1
                    input$com <- 1
                    combined <- merge(combined,input[1,!names(input) %in% c("surv.p","surv.time","surv.p.reference","surv.p.input")],by="com",all.x=TRUE)
                    combined$com <- NULL
                    surv.substage.death.chart.temp <- rbind(surv.substage.death.chart.temp,reference,combined)
                  }
                }
              }
            # Generate shifted main stages
              surv.dataset.chart.temp <- surv.main.chart[0,]
              for (group.index in 1:length(groups.order)){
                # for (start.stage.index in 1:(length(stages.order)-1)){
                for (start.stage.index in 1:(max.start.stage)){
                  for (end.stage.index in (start.stage.index+1):length(stages.order)){
                    reference <- surv.main.chart[surv.main.chart$start.stage.index==start.stage.index &
                                                      surv.main.chart$end.stage.index==end.stage.index &
                                                      surv.main.chart$group.index==group.index,]
                    input.1 <- surv.substage.death.chart.temp[surv.substage.death.chart.temp$start.stage.index==start.stage.index &
                                                             surv.substage.death.chart.temp$end.stage.index==end.stage.index &
                                                             surv.substage.death.chart.temp$group.index==group.index,]
                    input.2 <- surv.substage.death.chart.temp[surv.substage.death.chart.temp$start.stage.index==start.stage.index &
                                                             surv.substage.death.chart.temp$end.stage.index==length(stages.order) &
                                                             surv.substage.death.chart.temp$group.index==group.index,]

                    input.1$surv.p <- input.1$surv.p.reference
                    combined <- relative.curve(input.1[c("surv.p","surv.time")],input.2[c("surv.p","surv.time")],"under")
                    combined <- relative.curve(reference[c("surv.p","surv.time")],combined[c("surv.p","surv.time")],"under")
                    combined$com <- 1
                    reference$com <- 1
                    combined <- merge(combined,reference[1,!names(reference) %in% c("surv.p","surv.time","surv.p.reference","surv.p.input")],by="com",all.x=TRUE)
                    combined$com <- NULL
                    reference$com <- NULL
                    #surv.dataset.chart.temp <- rbind(surv.dataset.chart.temp,reference[,!names(reference) %in% c("surv.p.reference","surv.p.input")],combined[,!names(combined) %in% c("surv.p.reference","surv.p.input")])
                    surv.dataset.chart.temp <- rbind(surv.dataset.chart.temp,combined[,!names(combined) %in% c("surv.p.reference","surv.p.input")])
                  }
                }
              }
            # Replace death and survival datasets with shifted version
              surv.main.chart <- surv.dataset.chart.temp
              surv.substage.death.chart <- surv.substage.death.chart.temp
        } else {}
      }
    }
    # Generate ggplot chart
    {
      # Generate / change colors and legend
      {
        if (length(main.fill.colors)==1){
          main.fill.colors <- color.gradient(main.fill.colors,(length(stages.order)-1))
        } else {}
        # Add colors and generate fill list depending on options
          legend.states <- stages.order
          legend.fill.colors <- c("white",main.fill.colors)
          if (!is.na(death.indicator)){
            legend.states <- c(legend.states,death.indicator)
            legend.fill.colors <- c(legend.fill.colors,death.fill.color)
          } else {}
          if (risk.pool.size.line==TRUE){
            legend.states <- c(legend.states,"Risk pool size")
            legend.fill.colors <- c(legend.fill.colors,risk.pool.fill.color)
          } else {}
          legend.states <- ordered(legend.states,labels=legend.states,levels=legend.states)
        # Generate main chart values with larger factor list
          surv.main.chart$end.stage.factor.legend <- ordered(as.character(surv.main.chart$end.stage.factor),
                                                             labels = legend.states,levels=legend.states)
      }
      # Main plot
      chart <- ggplot2::ggplot() +
        geom_stepribbon(data=surv.main.chart,aes(x=.data$surv.time,ymax=.data$surv.p,fill=.data$end.stage.factor.legend,ymin=0),
                                               alpha=1,show.legend = TRUE,color="black") +
        ggplot2::theme_bw() %+replace%
        ggplot2::theme(
          panel.grid = element_blank(),
          plot.margin = unit(c(.1,.1,.1,.1), "cm"),
          axis.title.y=element_blank(),
          legend.position=legend.position,
          legend.title=element_blank(),
          legend.spacing.x = unit(0.1, 'cm'),
          axis.text = element_text(colour="black",size=10),
          strip.text = element_text(size = 12),
          strip.text.x = element_text(hjust=0),
          panel.spacing = unit(1, "lines")
        ) +
        ggplot2::scale_x_continuous(expand = c(0, 0),
                          labels=x.scale.function,
                          breaks = c(0,round2(time.horizon/365,0))) +
        ggplot2::xlab(x.axis.title) +
        ggplot2::scale_y_continuous(expand = c(0, 0),labels=scales::percent) +
        ggplot2::scale_color_manual(values=c(rep("black",length(legend.fill.colors)))) +
        ggplot2::theme(strip.background = element_blank(),
              strip.placement = "outside") +
        ggplot2::scale_fill_manual(values = legend.fill.colors,limits=levels(legend.states)) +
        ggplot2::guides(colour = guide_legend(override.aes = list(color="black",alpha = 1,size = 0.5))) +
        ggplot2::guides(fill = guide_legend(nrow = 1))

      # Add direct labeling
      if (direct.label==TRUE){
        chart <- chart +
          geomtextpath::geom_textline(data=surv.main.chart,aes(x=.data$surv.time,y=.data$surv.p,group=.data$end.stage.factor.legend,label=.data$end.stage.factor.legend),
                                      color=NA,textcolour="black",text_smoothing = 30,vjust = 1.3,
                                      size=3.5)
      }

    # Specify facet grid types
      if (chart.mode=="first transition"){
        chart <- chart +
          ggplot2::facet_grid(group.factor ~ 1,
                   switch="y") +
          ggplot2::theme(strip.text.x = element_blank())
      } else {
        chart <- chart +
          ggplot2::facet_grid(group.factor ~ start.stage.factor,
                   switch="y")
      }
    # Add risk pool proportion indicator if indicated
      if (risk.pool.size.line==TRUE){
        chart <- chart +
          ggplot2::coord_cartesian(xlim=c(0,(time.horizon/365)),ylim = c(-.2, 1)) +
          geom_stepribbon(data = surv.main.chart.risk.pool,
                                aes(x=.data$surv.time,ymin=((.data$surv.p.atrisk-1)/5),ymax=0),
                                alpha=1,fill=risk.pool.fill.color,show.legend = FALSE,color="black")
      } else {
        chart <- chart +
          ggplot2::coord_cartesian(xlim=c(0,(time.horizon/365)),ylim = c(0, 1))
      }
    # Add transient events if indicated
      if (!anyNA(ts.indicator)){
        chart <- chart +
          geom_stepribbon(data=df.substage.transient,
                          aes(x=.data$surv.time,ymin=.data$surv.p,ymax=.data$surv.p.reference,group=.data$reference.stage.index),
                                alpha=1,fill=ts.color,show.legend = FALSE,color="black")
      } else {}
    # Add death events if indicated
      if (!is.na(death.indicator)){
        chart <- chart +
          geom_stepribbon(data=surv.substage.death.chart,
                          aes(x=.data$surv.time,ymin=.data$surv.p,ymax=.data$surv.p.reference,group=.data$reference.stage.index),
                                alpha=1,fill=death.fill.color,show.legend = FALSE,color="black")
      } else {}
    # Remove y axis facet label if there are no groups defined
      if (anyNA(groups.order.orig)==TRUE && anyNA(groups.date.breaks)==TRUE){
        chart <- chart +
          ggplot2::theme(strip.text.y = element_blank())
      } else {}
    # Add overlays
      chart <- chart +
        ggplot2::geom_hline(yintercept=0)
    }
  }
  # Generate between-group statistical tests
  {
    # Test for heterogeneity of survival curves for main between-stage events between groups
    if (length(groups.order)>1){
      # Generate tests function
        generate.surv.diff <- function(init.stage){
          # Generate internal data
            if (!anyNA(groups.date.breaks)){
              generate.date.break.data <- function(stage.index,group.index){
                df.output <- TTE.ind[TTE.ind[[paste0("date.stage.",stage.index,".group")]]==groups.order[group.index],]
                df.output$group <- groups.order[group.index]
                return(df.output)
              }
              TTE.test <- do.call("rbind",
                                   lapply(1:length(groups.order),function(x) generate.date.break.data(init.stage,x)))
            } else {
              TTE.test <- TTE.ind
            }
          # Perform test
            chart.time <- as.integer(TTE.test[[paste("time.stage.",(init.stage),".to.",(init.stage+1),sep="")]])
            chart.event <- as.integer(!TTE.test[[paste("censorship.stage.",(init.stage),".to.",(init.stage+1),sep="")]])
            groups <- ordered(TTE.test$group,labels=groups.order,levels=groups.order)
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
              #TTE.test <- TTE[0,]
              generate.date.break.data <- function(stage.index,group.index){
                df.output <- TTE.ind[TTE.ind[[paste0("date.stage.",stage.index,".group")]]==groups.order[group.index],]
                df.output$group <- groups.order[group.index]
                return(df.output)
              }
              TTE.test <- do.call("rbind",
                                   lapply(1:length(groups.order),function(x) generate.date.break.data(init.stage,x)))
            } else {
              TTE.test <- TTE.ind
            }
          # Keep only relevent variables for internal manipulation
            time.var <- paste("time.stage.",(init.stage),".to.",(init.stage+1),sep="")
            events.var <- paste("censorship.stage.",(init.stage),".to.",(init.stage+1),sep="")
            TTE.internal <- TTE.test[c("group",time.var,events.var)]
            colnames(TTE.internal) <- c("group","time","event")
          # Designate reference vs. comparator groups
            #group.order.internal <- c(reference.group,groups.order[groups.order!=reference.group])
            TTE.internal$group <- factor(TTE.internal$group,levels=groups.order,labels=groups.order)
            TTE.internal$group <- stats::relevel(TTE.internal$group, ref = reference.group)
          # Run Cox Proportional hazards
            cox.ph <- coxph(Surv(time = time, event = event) ~ group,data=TTE.internal)
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
    output.df <- list("chart" = chart,"surv.dataset" = surv.main,"surv.dataset.chart" = surv.main.chart,
                      "events.long" = events.long,"events.wide" = events.wide,"TTE.ind" = TTE.ind,"TTE.quantile" = quantile.TTE)
    if (!is.na(death.indicator)) { output.df <- c(output.df,
                      list("surv.death" = surv.death,"surv.death.chart" = surv.death.chart,"TTE.ind.death" = TTE.ind.death))
    } else {}
    if (!anyNA(ts.indicator)) { output.df <- c(output.df,
                      list("surv.transient" = surv.transient,"surv.transient.chart" = surv.transient.chart))
    } else {}
    if ((length(groups.order)>1) | (length(groups.date.breaks)>1)) { output.df <- c(output.df,
                      list("surv.cox.ph" = cox.ph.combined,"surv.diffs" = surv.diffs.combined))
    } else {}
    # Temporary for debugging
    if (!is.na(death.indicator) & allow.sub.stage.mortality==TRUE){
      output.df <- c(output.df,
                      list("surv.substage.death.chart" = surv.substage.death.chart))
    }

    return(output.df)
  }
}
