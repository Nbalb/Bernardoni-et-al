require(gplots)
require(survival)


### Classic Survival: split dataset on the median of gene expression
plotclassicsurv <-
  function(oritrack,
           survival,
           mygene,
           title = "Kaplan Meier Curves",
           colorizeSignificant = FALSE,
           plot = TRUE,
           sizelegend = 1) {
    # Define types for survivalfit
    track <- oritrack
    track[] <- "Other"
    track[oritrack > median(oritrack)] <- paste0(mygene, " high")
    track[oritrack <= median(oritrack)] <- paste0(mygene, " low")
    track <- as.factor(track)
    track <- relevel(track, ref = paste0(mygene, " low"))
    
    # Fit model
    sfit <- survfit(survival ~ track)
    
    # Survdiff for all vs. all p-value
    if (length(levels(track)) == 1) {
      return(1)
    }
    sdiff <- survdiff(survival ~ track)
    p <- 1 - pchisq(sdiff$chisq, df = length(sdiff$n) - 1)
    
    # Sign? When I did this I understood it, now not anymore
    events_difference <- sdiff$obs - sdiff$exp
    if (events_difference[1] > events_difference[length(events_difference)]) {
      sign <- "neg"
    } else {
      sign <- "pos"
    }
    
    # OE<-(sdiff$obs-sdiff$exp)^2/sdiff$exp
    # sign<-"neg"
    # if(OE[1]<=OE[2]){sign<-"pos"}
    
    if (plot) {
      # Color matching
      colors <- rep(NA, 2)
      colors[grep("high", names(sfit$strata))] <- "red3"
      colors[grep("low", names(sfit$strata))] <- "cornflowerblue"
      
      # Plot
      plot(
        sfit,
        lwd = 3,
        main = title,
        xlab = "days",
        ylab = "Percent Survival",
        col = colors,
        mark.time = TRUE
      )
      mtext(paste(
        "p =",
        signif(p, 4),
        " (Samples: ",
        nrow(survival),
        " - Deaths: ",
        sum(survival[, 2], na.rm = TRUE),
        ")",
        sep = ""
      ),
      cex = 0.65)
      if (colorizeSignificant) {
        if (p <= 0.05 &
            sign == "pos") {
          rect(par("usr")[1],
               par("usr")[3],
               par("usr")[2],
               par("usr")[4],
               col = "#FF000022")
        }
        if (p <= 0.05 &
            sign == "neg") {
          rect(par("usr")[1],
               par("usr")[3],
               par("usr")[2],
               par("usr")[4],
               col = "#FF000022")
        }
      }
      grid(col = "grey")
      legend(
        "bottomleft",
        col = colors,
        legend = substr(names(sfit$strata), 7, nchar(names(sfit$strata))),
        lty = 1,
        lwd = 3,
        bg = "white",
        cex = sizelegend
      )
    }
    if (sign == "neg") {
      p <- -p
    }
    return(p)
  }


### Step Survival: divide data into overlapping groups, as in the paper by Gonzalo

plotstepsurv <-
  function(oritrack,
           survival,
           ngroups = 10,
           ylab = "Gene expression",
           title = NULL,
           plot = TRUE,
           layout = TRUE,
           cox_survival = FALSE,
           sizelegend = 0.9,
           mark.time = TRUE) {
    if (layout & plot) {
      layout(t(matrix(c(1, 2))), width = c(1, 2))
    }
    
    # Sort by gene exp
    oritrack <- sort(oritrack)
    survival <- survival[names(oritrack), ]
    if (plot) {
      plot(oritrack,
           ylab = ylab,
           pch = 20,
           cex.axis = 0.6)
      #grid()
      mtext(
        title,
        side = 3,
        line = -1.5,
        outer = TRUE,
        cex = 1.5,
        font = 2
      )
    }
    
    # Gene groups
    chunkfun <-
      function(x, n)
        split(x, cut(seq_along(x), n, labels = FALSE))
    track <- chunkfun(oritrack, ngroups)
    survcols <- colorpanel(ngroups - 1, "blue", "gray90", "red")
    
    for (i in 1:(ngroups - 1)) {
      ycut <- max(track[[i]])
      xcut <- which(oritrack == ycut)
      seglength <- floor(length(oritrack) / 10)
      if (plot) {
        segments(xcut - seglength,
                 ycut,
                 xcut + seglength,
                 ycut,
                 col = survcols[i],
                 lwd = 6)
      }
    }
    
    # Define types for survfit
    track <- cut(seq_along(oritrack), ngroups, labels = FALSE)
    samplesingroups <- c()
    samplestotake <- c()
    for (i in 1:(ngroups - 1)) {
      samplesingroups <-
        c(samplesingroups, rep(i, length(which(
          track == i | track == (i + 1)
        ))))
      samplestotake <- c(samplestotake, which(track == i |
                                                track == (i + 1)))
    }
    stepsurv <- survival[samplestotake, ]
    stepexp <- oritrack[samplestotake]
    track <- as.factor(samplesingroups)
    names(track) <- names(stepexp)
    
    # Fit model
    sfit <- survfit(stepsurv ~ track)
    
    # Survdiff for all vs. all p-value
    sdiff <- survdiff(stepsurv ~ track)
    p <- 1 - pchisq(sdiff$chisq, df = 2)
    
    
    # Sign?
    events_difference <- sdiff$obs - sdiff$exp
    if (events_difference[1] > events_difference[length(events_difference)]) {
      sign <- "neg"
    } else {
      sign <- "pos"
    }
    
    #Cox model
    if (cox_survival) {
      res.cox <- summary(coxph(survival ~ oritrack))
      cox_p <- signif(res.cox$coefficients[1, 5], 3)
      cox_coef <- signif(res.cox$coefficients[1, 1], 3)
    }
    
    # OE<-(sdiff$obs-sdiff$exp)^2/sdiff$exp
    # sign<-"neg"
    # if(OE[1]<=OE[2]){sign<-"pos"}
    
    # Plot
    # ?plot.survfit
    if (plot) {
      plot(
        sfit,
        lwd = 5,
        xlab = "Survival time (days)",
        ylab = "Percent Survival",
        col = survcols,
        mark.time = mark.time,
        main = ""
      )
      #grid(col="grey")
      legend(
        "bottomleft",
        legend = c(paste0(mygene, " High"), paste0(mygene, " Low")),
        col = c("red", "blue"),
        lty = 1,
        lwd = 3,
        bg = "white",
        cex = sizelegend
      )
      
      mtext(
        paste(
          "KM model p = ",
          scales::scientific(p, digits = 4),
          " (Samples: ",
          nrow(survival),
          " - Deaths: ",
          sum(survival[, 2], na.rm = TRUE),
          ")",
          sep = ""
        ),
        side = 3,
        line = -2.5,
        outer = TRUE,
        cex = sizelegend
      )
      
      if (cox_survival) {
        mtext(
          paste0(
            "Cox model p = ",
            scales::scientific(cox_p, digits = 4),
            "    ",
            "Cox coefficient = ",
            cox_coef
          ),
          cex = sizelegend,
          line = -3.5,
          outer = TRUE
        )
      }
      
      if (sign == "neg") {
        p <- -p
      }
    }
    return(p)
  }