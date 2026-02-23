#!/usr/bin/Rscript

#----
# Setup environment
#----

options(width=180)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(stringr)
library(fuzzyjoin)
set.seed(117)


##############################

### Define miami function

##############################

miami <- function (
    x, y = NULL,
    chr = "chr", bp = "pos", p = "p", snp = "rsid",
    col1 = c("#e34a33", "#fdbb84"), col2 = c("#2b8cbe","#a6bddb"),
    chrlabs = NULL, suggestiveline = -log10(1e-05),
    genomewideline = -log10(5e-08), ymin = -15, ymax = 15,
    x.name = "x", y.name = "y", highlight = NULL, highlight_col = "red",
    label = NULL, logp = TRUE, ymid = 0,
    special_colouring = FALSE, special_flag = NULL, special_colour = NULL,
    ...
  ) {

  ## Local variable placeholders to avoid NOTE checks
  CHR = BP = P = index = NULL

  # small transform to handle mid offset
  ymid <- abs(ymid)
  ymin = ymin + ymid
  ymax = ymax - ymid

  # ---- helper: column presence/type check ----
  column_check <- function(column, df) {
    if (!(column %in% names(df))) stop(paste0("Column ", column, " not found in ", deparse(substitute(df)), "!"))
    if (!(column == snp) & !is.numeric(df[[column]])) stop(paste(column, "column in", deparse(substitute(df)), "should be numeric."))
  }
  for (column in c(bp, p, chr, snp)) column_check(column, x)
  if (!is.null(y)) for (column in c(bp, p, chr, snp)) column_check(column, y)

  # ---- LOAD D1 (primary/top) ----
  print("##LOAD D1")
  d1 <- data.table::data.table(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]])) d1$SNP <- x[[snp]]
  if (!is.null(highlight)) d1$highlight <- x[[highlight]]

  # --- attach special columns BEFORE any subsetting so they stay aligned ---
  if (special_colouring && !is.null(special_flag) && !is.null(special_colour)) {
    # helper to coerce flexible flag encodings to logical (TRUE/FALSE/NA)
    coerce_flag <- function(sf) {
      if (is.null(sf)) return(rep(NA, nrow(d1)))
      # accept TRUE/FALSE, 1/0, "1"/"0", "TRUE"/"FALSE", "T"/"F"
      v <- ifelse(is.na(sf), NA,
                  sf %in% c(TRUE, 1, "1", "TRUE", "T"))
      return(v)
    }
    if (special_flag %in% names(x)) {
      d1$special_flag <- coerce_flag(x[[special_flag]])
    } else {
      d1$special_flag <- NA
    }
    if (special_colour %in% names(x)) {
      sc <- as.character(x[[special_colour]])
      sc[sc == ""] <- NA_character_
      d1$special_colour <- sc
    } else {
      d1$special_colour <- NA_character_
    }
  }

  # keep only numeric CHR/BP/P rows
  d1 <- subset(d1, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d1 <- d1[order(d1$CHR, d1$BP), ]
  if (logp) {
    d1$logp <- -log10(d1$P) - ymid
    d1 <- d1[d1$logp > 0, ]
  } else {
    d1$logp <- d1$P - ymid
    d1 <- d1[d1$logp > 0, ]
  }

  # ---- LOAD D2 (bottom) only if provided ----
  if (!is.null(y)) {
    print("##LOAD D2")
    d2 <- data.table::data.table(CHR = y[[chr]], BP = y[[bp]], P = y[[p]])
    if (!is.null(y[[snp]])) d2$SNP <- y[[snp]]
    if (!is.null(highlight)) d2$highlight <- y[[highlight]]

    # attach special columns for d2 (same coercion logic)
    if (special_colouring && !is.null(special_flag) && !is.null(special_colour)) {
      coerce_flag2 <- function(sf) {
        if (is.null(sf)) return(rep(NA, nrow(d2)))
        v <- ifelse(is.na(sf), NA,
                    sf %in% c(TRUE, 1, "1", "TRUE", "T"))
        return(v)
      }
      if (special_flag %in% names(y)) {
        d2$special_flag <- coerce_flag2(y[[special_flag]])
      } else {
        d2$special_flag <- NA
      }
      if (special_colour %in% names(y)) {
        sc2 <- as.character(y[[special_colour]])
        sc2[sc2 == ""] <- NA_character_
        d2$special_colour <- sc2
      } else {
        d2$special_colour <- NA_character_
      }
    }

    d2 <- subset(d2, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d2 <- d2[order(d2$CHR, d2$BP), ]
    if (logp) {
      d2$logp <- log10(d2$P) + ymid   # negatives for plotting downwards
      d2 <- d2[d2$logp < 0, ]
    } else {
      d2$logp <- -d2$P + ymid
      d2 <- d2[d2$logp < 0, ]
    }

    # combine for axis calculations etc.
    d <- data.frame(rbind(d1, d2))
  } else {
    d <- data.frame(d1)
  }

  # ---- Define plotting positions / chromosome ticks ----
  d$pos <- NA
  d$index <- NA
  ind <- 0
  for (i in unique(d$CHR)) {
    ind <- ind + 1
    d[d$CHR == i, ]$index <- ind
  }
  nchr <- length(unique(d$CHR))
  if (nchr == 1) {
    d$pos <- d$BP
    ticks <- floor(length(d$pos)) / 2 + 1
    xlabel <- paste("Chromosome", unique(d$CHR), "position")
    labs <- ticks
  } else {
    lastbase <- 0
    ticks <- NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos <- d[d$index == i, ]$BP
      } else {
        lastbase <- lastbase + tail(subset(d, index == i - 1)$BP, 1)
        d[d$index == i, ]$pos <- d[d$index == i, ]$BP + lastbase
      }
      ticks <- c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index == i, ]$pos)) / 2 + 1)
    }
    xlabel <- "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax <- ceiling(max(d$pos) * 1.03)
  xmin <- floor(max(d$pos) * -0.03)

  print("##CALL PLOT")

  # ---- default plotting arguments ----
  def_args <- list(
    xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", yaxt = "n",
    cex.axis = 1, cex.lab = 1.5, las = 1, pch = 20, cex = 1.5,
    xlim = c(xmin, xmax), ylim = c(ymin - 1, ymax + 1),
    xlab = xlabel, ylab = expression(-log[10](italic(p))),
    mar = 5 * c(5.1, 6.1, 4.1, 2.1)
  )
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))

  # ---- chromosome labels (optional) ----
  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      } else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    } else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }

  # ---- y ticks/labels ----
  yticks <- seq(-5 * ceiling(abs(ymin / 5)), 5 * ceiling(abs(ymax / 5)), by = 5)
  yticks <- sort(c(yticks + sign(yticks) * (5 - ymid), 5 - ymid, -5 + ymid))
  if (is.null(y)) yticks <- yticks[yticks >= 0]
  ylabels <- abs(yticks) + ymid

  # ---- draw significance lines ----
  if ((abs(suggestiveline) - ymid) > 0) {
    abline(h = suggestiveline - ymid, col = "blue", lwd = 3)
    if (!is.null(y)) abline(h = -suggestiveline + ymid, col = "blue", lwd = 3)
  }
  if ((abs(genomewideline) - ymid) > 0) {
    abline(h = genomewideline - ymid, col = "red", lwd = 3)
    if (!is.null(y)) abline(h = -genomewideline + ymid, col = "red", lwd = 3)
  }
  if (!is.null(y)) abline(h = 0, col = "black")

  # ---- legend ----
  if (!is.null(y)) legend(0, y = -ymax + ymid + 5 * strheight(1), c(x.name, y.name), fill = c(col1[1], col2[1]), bty = "n", cex = 1.2, title = NULL, yjust = 0.5)

  # ---- axes ----
  if (nchr == 1) {
    axis(1, cex.axis = 1.5, cex.lab = 1.5, las = 1, ...)
    axis(2, at = yticks, labels = ylabels, las = 1, cex.axis = 1.5, cex.lab = 1.5, ...)
  } else {
    axis(1, at = ticks, labels = labs, las = 1, cex.axis = 1.5, cex.lab = 1.5, ...)
    axis(2, at = yticks, labels = ylabels, las = 1, cex.axis = 1.5, cex.lab = 1.5, ...)
  }

  # ---- prepare colours per chromosome ----
  col1 <- rep(col1, max(d$CHR))
  col2 <- rep(col2, max(d$CHR))

  # split positive/negative sets for plotting
  d1 <- d[d$logp > 0, ]
  d2 <- d[d$logp < 0, ]

  # ---- Plotting points ----
  if (nchr == 1) {
    # single-chromosome simpler plotting
    with(d1[d1$logp < ymax, ], points(pos, logp, pch = 20, col = col1[1], ...))
    with(d1[d1$logp > ymax, ], points(pos, rep(ymax + 0.5, length(pos)), pch = 2, col = col1[1], ...))
    with(d2[d2$logp > ymin, ], points(pos, logp, pch = 20, col = col2[1], ...))
    with(d2[d2$logp < ymin, ], points(pos, rep(ymin - 0.5, length(pos)), pch = 6, col = col2[1], ...))
  } else {
    # multi-chromosome plotting: alternate colours per chromosome index
    print("##PLOT X")
    icol <- 1
    for (idx in unique(d1$index)) {

      # If highlight column provided, draw highlight shapes (takes precedence for shapes)
      if (!is.null(highlight)) {
        with(d1[d1$index == idx & d1$logp < ymax & d1$highlight == 0, ], points(pos, logp, col = col1[icol], pch = 20, ...))
        with(d1[d1$index == idx & d1$logp < ymax & d1$highlight == 1, ], points(pos, logp, col = highlight_col, pch = 4, cex = 4, lwd = 5, ...))
        with(d1[d1$index == idx & d1$logp < ymax & d1$highlight == 2, ], points(pos, logp, col = highlight_col, bg = highlight_col, pch = 24, cex = 1.5, lwd = 2, ...))
        with(d1[d1$index == idx & d1$logp > ymax, ], points(pos, rep(ymax + 0.5, length(pos)), col = highlight_col, bg = highlight_col, pch = 24, cex = 1.5, lwd = 2, ...))
      }

      # special colouring: draw points flagged as special using their provided colour (fallback to default)
      if (special_colouring) {
        # non-special (or NA) points first
        other <- d1[d1$index == idx & d1$logp < ymax & !( !is.na(d1$special_flag) & d1$special_flag == TRUE), ]
        if (nrow(other) > 0) points(other$pos, other$logp, col = col1[icol], pch = 20, ...)
        # special points with explicit colours
        tmp <- d1[d1$index == idx & d1$logp < ymax & !is.na(d1$special_flag) & d1$special_flag == TRUE, ]
        if (nrow(tmp) > 0) {
          tmpcols <- ifelse(is.na(tmp$special_colour) | tmp$special_colour == "", col1[icol], tmp$special_colour)
          points(tmp$pos, tmp$logp, col = tmpcols, pch = 20, ...)
        }
        # clipped-up points
        upclip <- d1[d1$index == idx & d1$logp > ymax, ]
        if (nrow(upclip) > 0) points(upclip$pos, rep(ymax + 0.5, nrow(upclip)), pch = 24, col = col1[icol], bg = col1[icol], cex = 1.5, lwd = 2, ...)
      } else if (is.null(highlight)) {
        # neither special nor highlight -> plain plotting
        with(d1[d1$index == idx & d1$logp < ymax, ], points(pos, logp, col = col1[icol], pch = 20, ...))
        with(d1[d1$index == idx & d1$logp > ymax, ], points(pos, rep(ymax + 0.5, length(pos)), pch = 24, col = col1[icol], bg = col1[icol], cex = 1.5, lwd = 2, ...))
      }

      icol <- icol + 1
    }

    print("##PLOT Y")
    if (!is.null(y)) {
      icol <- 1
      for (idx in unique(d2$index)) {

        if (!is.null(highlight)) {
          with(d2[d2$index == idx & d2$logp > ymin & d2$highlight == 0, ], points(pos, logp, col = col2[icol], pch = 20, ...))
          with(d2[d2$index == idx & d2$logp < ymax & d2$highlight == 1, ], points(pos, logp, col = highlight_col, pch = 4, cex = 4, lwd = 5, ...))
          with(d2[d2$index == idx & d2$logp < ymax & d2$highlight == 2, ], points(pos, logp, col = highlight_col, bg = highlight_col, pch = 25, cex = 1.5, lwd = 2, ...))
          with(d2[d2$index == idx & d2$logp < ymin, ], points(pos, rep(ymin - 0.5, length(pos)), col = highlight_col, bg = highlight_col, pch = 25, cex = 1.5, lwd = 2, ...))
        }

        if (special_colouring) {
          other <- d2[d2$index == idx & d2$logp > ymin & !( !is.na(d2$special_flag) & d2$special_flag == TRUE), ]
          if (nrow(other) > 0) points(other$pos, other$logp, col = col2[icol], pch = 20, ...)
          tmp <- d2[d2$index == idx & d2$logp > ymin & !is.na(d2$special_flag) & d2$special_flag == TRUE, ]
          if (nrow(tmp) > 0) {
            tmpcols <- ifelse(is.na(tmp$special_colour) | tmp$special_colour == "", col2[icol], tmp$special_colour)
            points(tmp$pos, tmp$logp, col = tmpcols, pch = 20, ...)
          }
          downclip <- d2[d2$index == idx & d2$logp < ymin, ]
          if (nrow(downclip) > 0) points(downclip$pos, rep(ymin - 0.5, nrow(downclip)), pch = 25, col = col2[icol], bg = col2[icol], cex = 1.5, lwd = 2, ...)
        } else if (is.null(highlight)) {
          with(d2[d2$index == idx & d2$logp > ymin, ], points(pos, logp, col = col2[icol], pch = 20, ...))
          with(d2[d2$index == idx & d2$logp < ymin, ], points(pos, rep(ymin - 0.5, length(pos)), pch = 25, col = col2[icol], bg = col2[icol], cex = 1.5, lwd = 2, ...))
        }

        icol <- icol + 1
      }
    }
  }

  # ---- Labels block (unchanged logic but safer forbidden construction) ----
  if (!is.null(label)) {
    required_cols <- c("LABEL", "CHR", "POS", "P")
    if (is.null(colnames(label))) warning(paste("label must be a data.frame with col.names:", paste0(required_cols, collapse = " ")))
    label <- data.frame(label)

    if (!all(required_cols %in% colnames(label))) {
      warning(paste("label data.frame is missing columns:", paste0(required_cols[!required_cols %in% colnames(label)], collapse = " ")))
    } else {

      # compute label positions in global coordinates
      label$pos <- label$POS
      for (i in unique(label$CHR)[!unique(label$CHR) == 1]) {
        lastbase <- tail(d[d$index == (i - 1), "pos"], n = 1)
        label[label$CHR == i, "pos"] <- label[label$CHR == i, "pos"] + lastbase
      }

      # adjust label P positions and clamp inside ymin/ymax
      label$P <- ifelse(label$P > 0, label$P - ymid + 2 * strheight(label$LABEL), label$P + ymid - 2 * strheight(label$LABEL))
      label$P <- ifelse(label$P > ymax, ymax - strheight(label$LABEL) / 2, label$P)
      label$P <- ifelse(label$P < ymin, ymin + strheight(label$LABEL) / 2, label$P)
      label <- label[order(-sign(label$P), abs(label$P), label$pos), ]
      row.names(label) <- 1:nrow(label)

      print("##POSITION LABELS")

      x.buffer <- strwidth("A")
      y.buffer <- strheight("A")

      label$YMIN <- label$P - strheight(label$LABEL) / 2
      label$YMAX <- label$P + strheight(label$LABEL) / 2

      # safer forbidden points creation: use d$highlight explicitly
      if (!is.null(highlight) && "highlight" %in% names(d)) {
        fh <- which(d$highlight == TRUE)
        if (length(fh) > 0) {
          forbidden <- data.frame(
            pos = d$pos[fh],
            P = d$logp[fh],
            XMIN = d$pos[fh] - strwidth("A") / 2,
            XMAX = d$pos[fh] + strwidth("A") / 2,
            YMIN = ifelse(d$logp[fh] < ymin, ymin, ifelse(d$logp[fh] > ymax, ymax, d$logp[fh] - sign(d$logp[fh]) * strheight("A") / 2)),
            YMAX = ifelse(d$logp[fh] > ymax, ymax, ifelse(d$logp[fh] < ymin, ymin, d$logp[fh] + sign(d$logp[fh]) * strheight("A") / 2))
          )
        } else {
          forbidden <- data.frame(pos = numeric(0), P = numeric(0), XMIN = numeric(0), XMAX = numeric(0), YMIN = numeric(0), YMAX = numeric(0))
        }
      } else {
        forbidden <- data.frame(pos = numeric(0), P = numeric(0), XMIN = numeric(0), XMAX = numeric(0), YMIN = numeric(0), YMAX = numeric(0))
      }

      rect_overlap <- function(df1, df2) {
        return(!(df2$XMAX <= df1$XMIN | df2$XMIN >= df1$XMAX | df2$YMAX <= df1$YMIN | df2$YMIN >= df1$YMAX))
      }

      # iterative moving to avoid overlaps
      for (precision in 1:10) {
        for (i in 1:nrow(label)) {

          n <- 0
          repeat {
            overlap_highlight <- forbidden[rect_overlap(df1 = label[i, ], df2 = forbidden), ]

            if (nrow(overlap_highlight) > 0) {
              if (all(overlap_highlight$P > 0)) {
                label[i, "YMAX"] <- max(overlap_highlight$YMAX) + strheight("A") / 2 * 3 + y.buffer
                label[i, "YMIN"] <- max(overlap_highlight$YMAX) + strheight("A") / 2 + y.buffer
              } else {
                label[i, "YMIN"] <- min(overlap_highlight$YMIN) - strheight("A") / 2 * 3 - y.buffer
                label[i, "YMAX"] <- min(overlap_highlight$YMIN) - strheight("A") / 2 - y.buffer
              }
            }

            overlap <- rbind(label[i, ], label[-i, ][rect_overlap(df1 = label[i, ], df2 = label[-i, ]), ])

            if (nrow(overlap) == 1) break
            n <- n + 1
            if (n == 1) cat("\n", label[i, "LABEL"], "")
            cat(n, "")

            if (all(overlap$P > 0)) {
              overlap[1, "YMAX"] <- max(overlap[-1, "YMAX"]) + strheight("A") / 2 * 3 + y.buffer
              overlap[1, "YMIN"] <- max(overlap[-1, "YMAX"]) + strheight("A") / 2 + y.buffer
            } else {
              overlap[1, "YMIN"] <- min(overlap[-1, "YMIN"]) - strheight("A") / 2 * 3 - y.buffer
              overlap[1, "YMAX"] <- min(overlap[-1, "YMIN"]) - strheight("A") / 2 - y.buffer
            }

            label[match(paste0(overlap$LABEL, overlap$P), paste0(label$LABEL, label$P)), ] <- overlap
          }
        }
      }

      cat("\n")

      label$P <- ifelse(label$P > ymax, ymax - strheight(label$LABEL) / 6, label$P)
      label$P <- ifelse(label$P < ymin, ymin + strheight(label$LABEL) / 6, label$P)
      moved <- data.table::data.table(label)[P != (YMIN + YMAX) / 2, .(x = pos, y1 = P, y2 = (YMIN + YMAX) / 2)]

      for (i in 1:nrow(moved)) lines(x = moved[i, x] + c(0, 0), y = c(moved[i, y1], moved[i, y2]) - sign(moved[i, y1]) * strheight("A"))

      shadowtext <- function(x, y = NULL, labels, col = 'white', bg = 'black',
                             theta = seq(0, 2 * pi, length.out = 50), r = 0.1, ...) {
        xy <- xy.coords(x, y)
        xo <- r * strwidth('A')
        yo <- r * strheight('A')
        for (i in theta) {
          text(xy$x + cos(i) * xo, xy$y + sin(i) * yo, labels, col = bg, ...)
        }
        text(xy$x, xy$y, labels, col = col, ...)
      }

      print("##DRAW LABELS")
      with(label, shadowtext(x = pos, y = (YMIN + YMAX) / 2, labels = LABEL, cex = 1.25, col = "black", bg = "white", font = 4, r = 0.2))
    }
  }

  par(xpd = FALSE)
}


##############################

### Define QQ functions

##############################

add_expected_pvalue_to_gwas_results <- function (gwas_results, p_value_threshold = 1, sample_rate = 1,
                                                    p_col = "p")
{
    ord <- order(gwas_results[,p_col])
    gwas_results <- gwas_results[ord, ]
    iSignificant <- which(gwas_results[,p_col] <= p_value_threshold)
    nSignificant <- length(iSignificant)
    nNotSignificant <- nrow(gwas_results) - nSignificant
    nTotal <- nSignificant + sample_rate * nNotSignificant
    gwas_results$rank <- NA
    if (nSignificant > 0) {
        gwas_results$rank[1:nSignificant] <- 1:nSignificant
    }
    if (nNotSignificant > 0) {
        gwas_results$rank[(nSignificant + 1):nrow(gwas_results)] <- nSignificant +
            sample_rate * (1:nNotSignificant)
    }
    gwas_results$expected_p_value <- gwas_results$rank/nTotal
    return(gwas_results)
}

plot_QQ <- function (gwas_results, p_value_threshold = 1, sample_rate = 1,
    maf_threshold = 0, p_col = "p", freq_col = "freq1")
{
    pal <- brewer.pal(3, name = "Set2")
    gwas_results$maf <- ifelse(gwas_results[,freq_col] <= 0.5, gwas_results[,freq_col],
        1 - gwas_results[,freq_col])
    if (maf_threshold > 0) {
        i <- which(gwas_results$maf > maf_threshold)
        gwas_results <- gwas_results[i, ]
        cat("MAF threshold of", maf_threshold, "applied to QQ plot\n")
    }
    i <- which(gwas_results$maf > 0.05)
    gwas_results_0.05_0.5 <- gwas_results[i, ]
    gwas_results_0.05_0.5 <- add_expected_pvalue_to_gwas_results(gwas_results_0.05_0.5,
        p_value_threshold, sample_rate, p_col = p_col)
    plot(-log10(gwas_results_0.05_0.5$expected_p_value), -log10(gwas_results_0.05_0.5[, p_col]),
        xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
        col = pal[1], pch = 20, ylim = c(0, ceiling(max(-log10(gwas_results[, p_col]),
            na.rm = T))))
    i <- which((gwas_results$maf <= 0.05) & (gwas_results$maf >
        0.01))
    gwas_results_0.01_0.05 <- gwas_results[i, ]
    gwas_results_0.01_0.05 <- add_expected_pvalue_to_gwas_results(gwas_results_0.01_0.05,
        p_value_threshold, sample_rate, p_col = p_col)
    points(-log10(gwas_results_0.01_0.05$expected_p_value), -log10(gwas_results_0.01_0.05[, p_col]),
        xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
        col = pal[2], pch = 20)
    i <- which(gwas_results$maf <= 0.01)
    if(length(i) == 0){
		abline(0, 1, col = "red")
		legend("topleft", c("MAF > 5%", "5% >= MAF > 1%"),
			col = pal, pch = c(20, 20))
    } else {
		gwas_results_0.01 <- gwas_results[i, ]
		gwas_results_0.01 <- add_expected_pvalue_to_gwas_results(gwas_results_0.01,
			p_value_threshold, sample_rate, p_col = p_col)
		points(-log10(gwas_results_0.01$expected_p_value), -log10(gwas_results_0.01[, p_col]),
			xlab = "Expected -log10(p)", ylab = "Observed -log10(p)",
			col = pal[3], pch = 20)
		abline(0, 1, col = "red")
		legend("topleft", c("MAF > 5%", "5% >= MAF > 1%", "MAF <= 1%"),
			col = pal, pch = c(20, 20, 20))
    }
}


##############################

### Organise data for plotting

##############################

all_files = c(Sys.glob("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/*.TBL.gz"),
              Sys.glob("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/*.TBL"))

all_files = all_files[!(grepl("herit", all_files))]

print(all_files)

### COLOUR DEFINITION
### CURRENT STANDARD IS:
### TOTAL/GP = BLACK/GREEN
### HES/PRESCRIPT = BLUE/ORANGE

green_set = c("#33BBBB", "#009988", "#006655") # Primary
blue_set = c("#3399DD", "#0077BB", "#004488") # Inpatient
purple_set = c("#FF5599", "#EE3377", "#AA1144") # INOUT patient
black_set = c("#252525", "#525252", "#737373") # Total
orange_set = c("#FF9955", "#EE7733", "#BB5500") # Prescription

### Column names from METAL

chrom_col = "Chromosome"
pos_col = "Position"
freq_col = "Freq1"
p_col = "P-value"

for (i in 1:length(all_files)){

    ### Reading file and getting variable and column names

    curr_file_name = basename(all_files[i])

    cat(paste0("\nRunning on ", curr_file_name, ". CHECK1.", i, " COMPLETE.\n\n==================\n\n"))

    curr_pheno = fread(all_files[i], data.table = FALSE, tmpdir = "/scratch/project_2007428/projects/prj_001_cost_gwas/tmpdir/")

    cohort = "Meta_v4"
    trait_name = strsplit(curr_file_name, "_")[[1]][1]
    stratif_name = strsplit(curr_file_name, "_")[[1]][2]
    sample_size = as.numeric(max(curr_pheno$Weight))

    nice_name = case_when(trait_name == "ALL" ~ "total",
                            trait_name == "INOUT" ~ "inpatient + outpatient",
                            trait_name == "DRUG" ~ "prescription",
                            trait_name == "PRIM" ~ "primary care",
                            trait_name == "IN" ~ "inpatient",
                            .default = as.character(trait_name))

    nice_stratif_name = case_when(stratif_name == "36_55" ~ "age between 36-55",
                                    stratif_name == "19_35" ~ "age between 19-35",
                                    stratif_name == "5_18" ~ "age between 5-18",
                                    stratif_name == "56_75" ~ "age between 56-75",
                                    stratif_name == "76_95" ~ "age 76+",
                                    stratif_name == "M" ~ "males",
                                    stratif_name == "F" ~ "females",
                                    .default = as.character(stratif_name))

    cat(paste0("Variable names collected.\nCohort: ", cohort,
                "\nTrait: ", nice_name,
                "\nStratification: ", nice_stratif_name,
                "\nFrequency column: ", freq_col,
                "\nCHECK2.", i, " COMPLETE.\n\n==================\n\n"))

    print(nrow(curr_pheno))

    p_val_thresh = 0.01

    if (any(curr_pheno[,p_col] == 0)){

    	curr_pheno[curr_pheno[,p_col] == 0, p_col] = 1e-300

    }

    ### Replace X-chrom with "23"

    curr_pheno[,chrom_col] = gsub("X", 23, curr_pheno[,chrom_col])
    curr_pheno[,chrom_col] = as.numeric(curr_pheno[,chrom_col])

    curr_pheno_manhattan = curr_pheno[curr_pheno[, p_col] < p_val_thresh,]
    curr_pheno_manhattan = curr_pheno_manhattan[curr_pheno_manhattan[, freq_col] >= 0.01,]

    if (any(curr_pheno_manhattan[, freq_col] > 0.99)){

      curr_pheno_manhattan = curr_pheno_manhattan[curr_pheno_manhattan[, freq_col] <= 0.99,]

      curr_pheno = curr_pheno[curr_pheno[, freq_col] <= 0.99,]

    }

    cat(paste0("\nAF and INFO filter applied.\n\n"))

    print(nrow(curr_pheno_manhattan))
    print(summary(curr_pheno_manhattan))

    colour_use = case_when(trait_name == "ALL" ~ black_set,
                            trait_name == "INOUT" ~ purple_set,
                            trait_name == "DRUG" ~ orange_set,
                            trait_name == "PRIM" ~ green_set,
                            trait_name == "IN" ~ blue_set,
                            .default = c("red", "blue", "green"))

    ### SAVE AS BOTH PDF AND PNG

    #plot_out = miami(x = curr_pheno_1,
    #    y = curr_pheno_2,
    #    chr = "CHROM", bp = "GENPOS", snp = "ID", p = "PVAL", 
    #    genomewideline = -log10(5e-8),
    #    ymin = -25,
    #    ymax = 25,
    #    col1 = colour_1,
    #    col2 = colour_2,
    #    x.name = paste0(cohort, " ", trait_1, " costs"),
    #    y.name = paste0(cohort, " ", trait_2, " costs")
    #    )

    manhattan_figDir = paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/figures/manhattan_plots/", cohort, "/")
    qq_figDir = paste0("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/figures/qq_plots/", cohort, "/")

    system(paste0("mkdir -p ", manhattan_figDir))
    system(paste0("mkdir -p ", qq_figDir))

    fig_name = ifelse(stratif_name == "ALL", paste0(trait_name, "_costs_", cohort, "_n_", sample_size), paste0(trait_name, "_costs_", stratif_name, "_", cohort, "_n_", sample_size))
    legend_name = ifelse(stratif_name == "ALL", paste0(cohort, " ", nice_name, " costs, N = ", sample_size), paste0(cohort, " ", nice_name, " costs (stratified by ", nice_stratif_name, "), N = ", sample_size))

    max_p = 5 * (ceiling(-log10(min(curr_pheno_manhattan[, p_col])) / 5))

    if (max_p < 10){

      max_p = 10

    }

    cat(paste0("All processing complete. Plotting.\nLegend name: ", legend_name, ". CHECK3.", i, " COMPLETE.\n\n==================\n\n"))

    #pdf(paste0(manhattan_figDir, fig_name, ".pdf"), width = 30, height = 20)
    #par(mar=c(5.1,5.1,4.1,2.1))

    #miami(x = curr_pheno_manhattan,
    #    chr = chrom_col, bp = pos_col, snp = "SNPID", p = p_col, 
    #    genomewideline = -log10(5e-8),
    #    ymax = max_p,
    #    ymin = 0,
    #    col1 = colour_use,
    #    x.name = legend_name
    #    )

    #dev.off()

    png(paste0(manhattan_figDir, fig_name, ".png"), width = 2800, height = 1000)
    #par(mar=c(5.1,5.1,4.1,2.1))

    miami(x = curr_pheno_manhattan,
        chr = chrom_col, bp = pos_col, snp = "MarkerName", p = p_col,
        genomewideline = -log10(5e-8),
        ymax = max_p,
        ymin = 0,
        col1 = colour_use,
        x.name = legend_name
        )

    dev.off()

    cat(paste0("Manhattans made. Doing QQ plot. CHECK4.", i, " COMPLETE.\n\n==================\n\n"))

    curr_pheno = curr_pheno[sample(nrow(curr_pheno), round(nrow(curr_pheno) / 10)),]

    print(summary(curr_pheno))

    png(paste0(qq_figDir, fig_name, ".png"), width = 1000, height = 1000)

      plot_QQ(gwas_results = curr_pheno, p_col = p_col, freq_col = freq_col)

    dev.off()

    cat(paste0("QQs made. Doing QQ plot. CHECK5.", i, " COMPLETE.\n\n==================\n\n"))

}

### Make the Miami plots

cat(paste0("MAIN LOOP COMPLETE.\n\n==================\n\n"))

mainDir = "/scratch/project_2007428/projects/prj_001_cost_gwas/"
inDir = paste0(mainDir, "outputs/METAL_v4/")
snpDir = paste0(mainDir, "outputs/gcta_cojo_v4_MAF_0_001/jma_combined/")
manhattan_figDir = paste0(mainDir, "outputs/figures/manhattan_plots/Meta_v4/")
qq_figDir = paste0(mainDir, "outputs/figures/qq_plots/Meta_v4/")

### Miami for IN + DRUG

primary_pheno_1 = fread(paste0(inDir, "IN_ALL_metal_output_1.TBL.gz"), data.table = FALSE)
primary_pheno_2 = fread(paste0(inDir, "DRUG_ALL_metal_output_1.TBL.gz"), data.table = FALSE)

cat(paste0("\nREADING COMPLETE.\n\n==================\n\n"))

### P-value threshold to not plot all variants

p_val_thresh = 0.05

### Any P-values of 0 converted to minimum value

if (any(primary_pheno_1[,p_col] == 0)){

	primary_pheno_1[primary_pheno_1[,p_col] == 0, p_col] = 1e-300

}

### Apply threshold and MAF < 0.01 filter

primary_pheno_1 = primary_pheno_1[primary_pheno_1[, p_col] < p_val_thresh,]
primary_pheno_1 = primary_pheno_1[primary_pheno_1[, freq_col] >= 0.01,]

if (any(primary_pheno_1[, freq_col] > 0.99)){

	primary_pheno_1 = primary_pheno_1[primary_pheno_1[, freq_col] <= 0.99,]

	primary_pheno_1 = primary_pheno_1[primary_pheno_1[, freq_col] <= 0.99,]

}

### Repeat this for the second phenotype
### Any P-values of 0 converted to minimum value

if (any(primary_pheno_2[,p_col] == 0)){

	primary_pheno_2[primary_pheno_2[,p_col] == 0, p_col] = 1e-300

}

### Apply threshold and MAF < 0.01 filter

primary_pheno_2 = primary_pheno_2[primary_pheno_2[, p_col] < p_val_thresh,]
primary_pheno_2 = primary_pheno_2[primary_pheno_2[, freq_col] >= 0.01,]

if (any(primary_pheno_2[, freq_col] > 0.99)){

	primary_pheno_2 = primary_pheno_2[primary_pheno_2[, freq_col] <= 0.99,]

	primary_pheno_2 = primary_pheno_2[primary_pheno_2[, freq_col] <= 0.99,]

}

cat(paste0("\nPRUNING COMPLETE.\n\n==================\n\n"))

### Create the "highlight" column

primary_pheno_1$highlight_col = FALSE
primary_pheno_2$highlight_col = FALSE

### SNPs for highlighting

prim_pheno_1_snps = fread(paste0(snpDir, "IN_ALL_gcta_jma_out_annotated.txt"), data.table = FALSE)
prim_pheno_2_snps = fread(paste0(snpDir, "DRUG_ALL_gcta_jma_out_annotated.txt"), data.table = FALSE)

prim_pheno_1_snps = prim_pheno_1_snps[prim_pheno_1_snps$p < 5e-8,]
prim_pheno_2_snps = prim_pheno_2_snps[prim_pheno_2_snps$p < 5e-8,]

### Switch to TRUE

primary_pheno_1$highlight_col[which(primary_pheno_1$MarkerName %in% prim_pheno_1_snps$SNP)] = TRUE
primary_pheno_2$highlight_col[which(primary_pheno_2$MarkerName %in% prim_pheno_2_snps$SNP)] = TRUE

### Change X to 23

primary_pheno_1[,chrom_col] = gsub("X", 23, primary_pheno_1[,chrom_col])
primary_pheno_1[,chrom_col] = as.numeric(primary_pheno_1[,chrom_col])

primary_pheno_2[,chrom_col] = gsub("X", 23, primary_pheno_2[,chrom_col])
primary_pheno_2[,chrom_col] = as.numeric(primary_pheno_2[,chrom_col])

cat(paste0("\nHIGHLIGHTING COMPLETE.\n\n==================\n\n"))

### Now add special colour column to highlight overlap

locus_table = fread("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/table/top_variants_main_analysis.tsv", data.table = FALSE)

keep_phenos = c("Inpatient", "Prescription drugs")
locus_table = locus_table[locus_table$phenotype %in% keep_phenos,]

locus_table$pheno = case_when(locus_table$phenotype == "Inpatient" ~ "IN",
                              locus_table$phenotype == "Prescription drugs" ~ "DRUG",
                              locus_table$phenotype == "Inpatient + outpatient" ~ "INOUT",
                              locus_table$phenotype == "Primary care" ~ "PRIM",
                              .default = as.character(locus_table$phenotype))

locus_table = locus_table %>%
                    group_by(locus_gene) %>%
                    filter(n_distinct(pheno) == length(unique(locus_table$pheno))) %>%
                    ungroup()

### Locus definition table

locus_range_table = locus_table %>%
                      group_by(locus_gene, chromosome) %>%
                        summarise(
                          start_bp = min(position_38),
                          end_bp   = max(position_38),
                          .groups = "drop"
                        ) %>%
                      arrange(chromosome, start_bp)

locus_range_table = locus_table %>%
                      group_by(locus_gene, chromosome) %>%
                        summarise(
                          start_bp = min(position_38),
                          end_bp   = max(position_38),
                          .groups = "drop"
                        ) %>%
                      arrange(chromosome, start_bp)

locus_range_table$start_bp = locus_range_table$start_bp - 500000
locus_range_table$end_bp = locus_range_table$end_bp + 500000

### Now colour the SNPs based on the locus ranges

primary_pheno_1_flag = primary_pheno_1 %>%
  # Perform an interval join between SNP position and locus ranges
  fuzzy_inner_join(
    locus_range_table,
    by = c("Chromosome" = "chromosome",
           "Position" = "start_bp",
           "Position" = "end_bp"),
    match_fun = list(`==`, `>=`, `<=`)
  ) %>%
  # Add flag + colour columns
  mutate(
    special_flag = TRUE,
    special_colour = "black"
  ) %>%
  # Keep GWAS columns plus new flags (drop extra locus info unless you want it)
  select(names(primary_pheno_1), special_flag, special_colour)

### Extract from original data frame and rbind back

primary_pheno_1 = primary_pheno_1[!(primary_pheno_1$MarkerName %in% primary_pheno_1_flag$MarkerName),]

primary_pheno_1$special_flag = FALSE
primary_pheno_1$special_colour = NA

primary_pheno_1 = rbind(primary_pheno_1, primary_pheno_1_flag)

cat(paste0("\nPHENO ONE OVERLAP COMPLETE.\n\n==================\n\n"))

### Repeat for second phenotype
### Now colour the SNPs based on the locus ranges

primary_pheno_2_flag = primary_pheno_2 %>%
  # Perform an interval join between SNP position and locus ranges
  fuzzy_inner_join(
    locus_range_table,
    by = c("Chromosome" = "chromosome",
           "Position" = "start_bp",
           "Position" = "end_bp"),
    match_fun = list(`==`, `>=`, `<=`)
  ) %>%
  # Add flag + colour columns
  mutate(
    special_flag = TRUE,
    special_colour = "black"
  ) %>%
  # Keep GWAS columns plus new flags (drop extra locus info unless you want it)
  select(names(primary_pheno_2), special_flag, special_colour)

### Extract from original data frame and rbind back

primary_pheno_2 = primary_pheno_2[!(primary_pheno_2$MarkerName %in% primary_pheno_2_flag$MarkerName),]

primary_pheno_2$special_flag = FALSE
primary_pheno_2$special_colour = NA

primary_pheno_2 = rbind(primary_pheno_2, primary_pheno_2_flag)

cat(paste0("\nPHENO TWO OVERLAP COMPLETE.\n\n==================\n\n"))

fwrite(primary_pheno_1, "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/pruned_meta/IN_ALL_miami_plot_frame.txt.gz",
        quote = FALSE, sep = "\t", row.names = FALSE, na = "NA", compress = "gzip")
fwrite(primary_pheno_2, "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/pruned_meta/DRUG_ALL_miami_plot_frame.txt.gz",
        quote = FALSE, sep = "\t", row.names = FALSE, na = "NA", compress = "gzip")

### PLOT

png(paste0(manhattan_figDir, "IN_plus_DRUG_miami_plot_highlighted.png"), width = 3500, height = 2200, res = 300)
#par(mar=c(5.1,5.1,4.1,2.1))

miami(x = primary_pheno_1,
    y = primary_pheno_2,
	  chr = chrom_col, bp = pos_col, snp = "MarkerName", p = p_col,
    #highlight = "highlight_col",
    #highlight_col = "black",
    genomewideline = -log10(5e-8),
    ymin = -25,
    ymax = 25,
    col1 = blue_set,
    col2 = orange_set,
    x.name = paste0("Inpatient costs"),
    y.name = paste0("Prescription drug costs"),
    special_colouring = TRUE,
    special_flag = "special_flag",
    special_colour = "special_colour"
    )

dev.off()

### Now making the final MH

cat(paste0("MIAMI ONE created. Making the second.\n\n==================\n\n"))

primary_pheno_1 = fread(paste0(inDir, "INOUT_ALL_metal_output_1.TBL.gz"), data.table = FALSE)
primary_pheno_2 = fread(paste0(inDir, "PRIM_ALL_metal_output_1.TBL.gz"), data.table = FALSE)

cat(paste0("\nREADING COMPLETE.\n\n==================\n\n"))

### P-value threshold to not plot all variants

p_val_thresh = 0.05

### Any P-values of 0 converted to minimum value

if (any(primary_pheno_1[,p_col] == 0)){

	primary_pheno_1[primary_pheno_1[,p_col] == 0, p_col] = 1e-300

}

### Apply threshold and MAF < 0.01 filter

primary_pheno_1 = primary_pheno_1[primary_pheno_1[, p_col] < p_val_thresh,]
primary_pheno_1 = primary_pheno_1[primary_pheno_1[, freq_col] >= 0.01,]

if (any(primary_pheno_1[, freq_col] > 0.99)){

	primary_pheno_1 = primary_pheno_1[primary_pheno_1[, freq_col] <= 0.99,]

	primary_pheno_1 = primary_pheno_1[primary_pheno_1[, freq_col] <= 0.99,]

}

### Repeat this for the second phenotype
### Any P-values of 0 converted to minimum value

if (any(primary_pheno_2[,p_col] == 0)){

	primary_pheno_2[primary_pheno_2[,p_col] == 0, p_col] = 1e-300

}

### Apply threshold and MAF < 0.01 filter

primary_pheno_2 = primary_pheno_2[primary_pheno_2[, p_col] < p_val_thresh,]
primary_pheno_2 = primary_pheno_2[primary_pheno_2[, freq_col] >= 0.01,]

if (any(primary_pheno_2[, freq_col] > 0.99)){

	primary_pheno_2 = primary_pheno_2[primary_pheno_2[, freq_col] <= 0.99,]

	primary_pheno_2 = primary_pheno_2[primary_pheno_2[, freq_col] <= 0.99,]

}

cat(paste0("\nPRUNING COMPLETE.\n\n==================\n\n"))

### Create the "highlight" column

primary_pheno_1$highlight_col = FALSE
primary_pheno_2$highlight_col = FALSE

### SNPs for highlighting

prim_pheno_1_snps = fread(paste0(snpDir, "INOUT_ALL_gcta_jma_out_annotated.txt"), data.table = FALSE)
prim_pheno_2_snps = fread(paste0(snpDir, "PRIM_ALL_gcta_jma_out_annotated.txt"), data.table = FALSE)

prim_pheno_1_snps = prim_pheno_1_snps[prim_pheno_1_snps$p < 5e-8,]
prim_pheno_2_snps = prim_pheno_2_snps[prim_pheno_2_snps$p < 5e-8,]

### Switch to TRUE

primary_pheno_1$highlight_col[which(primary_pheno_1$MarkerName %in% prim_pheno_1_snps$SNP)] = TRUE
primary_pheno_2$highlight_col[which(primary_pheno_2$MarkerName %in% prim_pheno_2_snps$SNP)] = TRUE

### Change X to 23

primary_pheno_1[,chrom_col] = gsub("X", 23, primary_pheno_1[,chrom_col])
primary_pheno_1[,chrom_col] = as.numeric(primary_pheno_1[,chrom_col])

primary_pheno_2[,chrom_col] = gsub("X", 23, primary_pheno_2[,chrom_col])
primary_pheno_2[,chrom_col] = as.numeric(primary_pheno_2[,chrom_col])

cat(paste0("\nHIGHLIGHTING COMPLETE.\n\n==================\n\n"))

### Now add special colour column to highlight overlap

locus_table = fread("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/table/top_variants_main_analysis_MAF_0_001_250kb.tsv", data.table = FALSE)

keep_phenos = c("Inpatient + Outpatient", "Primary Care")
locus_table = locus_table[locus_table$phenotype %in% keep_phenos,]

locus_table$pheno = case_when(locus_table$phenotype == "Inpatient" ~ "IN",
                              locus_table$phenotype == "Prescription drugs" ~ "DRUG",
                              locus_table$phenotype == "Inpatient + outpatient" ~ "INOUT",
                              locus_table$phenotype == "Primary care" ~ "PRIM",
                              .default = as.character(locus_table$phenotype))

locus_table = locus_table %>%
                    group_by(locus_gene) %>%
                    filter(n_distinct(pheno) == length(unique(locus_table$pheno))) %>%
                    ungroup()

### Locus definition table

locus_range_table = locus_table %>%
                      group_by(locus_gene, chromosome) %>%
                        summarise(
                          start_bp = min(position_38),
                          end_bp   = max(position_38),
                          .groups = "drop"
                        ) %>%
                      arrange(chromosome, start_bp)

locus_range_table = locus_table %>%
                      group_by(locus_gene, chromosome) %>%
                        summarise(
                          start_bp = min(position_38),
                          end_bp   = max(position_38),
                          .groups = "drop"
                        ) %>%
                      arrange(chromosome, start_bp)

locus_range_table$start_bp = locus_range_table$start_bp - 500000
locus_range_table$end_bp = locus_range_table$end_bp + 500000

### Now colour the SNPs based on the locus ranges

primary_pheno_1_flag = primary_pheno_1 %>%
  # Perform an interval join between SNP position and locus ranges
  fuzzy_inner_join(
    locus_range_table,
    by = c("Chromosome" = "chromosome",
           "Position" = "start_bp",
           "Position" = "end_bp"),
    match_fun = list(`==`, `>=`, `<=`)
  ) %>%
  # Add flag + colour columns
  mutate(
    special_flag = TRUE,
    special_colour = "black"
  ) %>%
  # Keep GWAS columns plus new flags (drop extra locus info unless you want it)
  select(names(primary_pheno_1), special_flag, special_colour)

### Extract from original data frame and rbind back

primary_pheno_1 = primary_pheno_1[!(primary_pheno_1$MarkerName %in% primary_pheno_1_flag$MarkerName),]

primary_pheno_1$special_flag = FALSE
primary_pheno_1$special_colour = NA

primary_pheno_1 = rbind(primary_pheno_1, primary_pheno_1_flag)

cat(paste0("\nPHENO ONE OVERLAP COMPLETE.\n\n==================\n\n"))

### Repeat for second phenotype
### Now colour the SNPs based on the locus ranges

primary_pheno_2_flag = primary_pheno_2 %>%
  # Perform an interval join between SNP position and locus ranges
  fuzzy_inner_join(
    locus_range_table,
    by = c("Chromosome" = "chromosome",
           "Position" = "start_bp",
           "Position" = "end_bp"),
    match_fun = list(`==`, `>=`, `<=`)
  ) %>%
  # Add flag + colour columns
  mutate(
    special_flag = TRUE,
    special_colour = "black"
  ) %>%
  # Keep GWAS columns plus new flags (drop extra locus info unless you want it)
  select(names(primary_pheno_2), special_flag, special_colour)

### Extract from original data frame and rbind back

primary_pheno_2 = primary_pheno_2[!(primary_pheno_2$MarkerName %in% primary_pheno_2_flag$MarkerName),]

primary_pheno_2$special_flag = FALSE
primary_pheno_2$special_colour = NA

primary_pheno_2 = rbind(primary_pheno_2, primary_pheno_2_flag)

cat(paste0("\nPHENO TWO OVERLAP COMPLETE.\n\n==================\n\n"))

fwrite(primary_pheno_1, "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/pruned_meta/INOUT_ALL_miami_plot_frame.txt.gz",
        quote = FALSE, sep = "\t", row.names = FALSE, na = "NA", compress = "gzip")
fwrite(primary_pheno_2, "/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/METAL_v4/pruned_meta/PRIM_ALL_miami_plot_frame.txt.gz",
        quote = FALSE, sep = "\t", row.names = FALSE, na = "NA", compress = "gzip")

### PLOT

png(paste0(manhattan_figDir, "INOUT_plus_PRIM_miami_plot_highlighted.png"), width = 3500, height = 2200, res = 300)
#par(mar=c(5.1,5.1,4.1,2.1))

miami(x = primary_pheno_1,
    y = primary_pheno_2,
	  chr = chrom_col, bp = pos_col, snp = "MarkerName", p = p_col,
    #highlight = "highlight_col",
    #highlight_col = "black",
    genomewideline = -log10(5e-8),
    ymin = -25,
    ymax = 25,
    col1 = purple_set,
    col2 = green_set,
    x.name = paste0("Inpatient + outpatient costs"),
    y.name = paste0("Primary care costs"),
    special_colouring = TRUE,
    special_flag = "special_flag",
    special_colour = "special_colour"
    )

dev.off()

