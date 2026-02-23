#!/usr/bin/Rscript

#----
# Setup environment
#----

options(width=180)
library(data.table)
library(dplyr)
library(RColorBrewer)


##############################

### Define miami function

##############################

miami <- function (x, y=NULL, chr = "chr", bp = "pos", p = "p", snp = "rsid", col1 = c("#e34a33", "#fdbb84"), col2 = c("#2b8cbe","#a6bddb"), chrlabs = NULL, suggestiveline = -log10(1e-05),
          genomewideline = -log10(5e-08), ymin = -15, ymax = 15, x.name = "x", y.name = "y", highlight = NULL, highlight_col="red", label = NULL, logp = TRUE, ymid = 0, ...)
  {
  CHR = BP = P = index = NULL
  ymid <- abs(ymid)
  ymin = ymin + ymid
  ymax = ymax - ymid

  #----
  # Check variables
  #----

  column_check <- function(column, df) {
    if (!(column %in% names(df)))
      stop(paste0("Column ", column, " not found in ",deparse(substitute(df)),"!"))
    if (!(column == snp) & !is.numeric(df[[column]]))
      stop(paste(column, "column in",deparse(substitute(df)),"should be numeric."))
  }

  for (column in c(bp,p,chr,snp)) column_check(column, x)
  if(!is.null(y)) for (column in c(bp,p,chr,snp)) column_check(column, y)

  #----
  # Create data frame
  #----
  print("##LOAD D1")
  
  d1 = data.table(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
  if (!is.null(x[[snp]])) d1$SNP <- x[[snp]]
  if (!is.null(highlight)) d1$highlight <- x[[highlight]]
  d1 <- subset(d1, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
  d1 <- d1[order(d1$CHR, d1$BP), ]
  if (logp) {
    d1$logp <- -log10(d1$P) - ymid
    d1 <- d1[d1$logp > 0,]
  } else {
    d1$logp <- d1$P - ymid
    d1 <- d1[d1$logp > 0,]
  }
  
  
  if(!is.null(y)) {
    print("##LOAD D2")
    d2 = data.table(CHR = y[[chr]], BP = y[[bp]], P = y[[p]])
    if (!is.null(y[[snp]])) d2$SNP <- y[[snp]]
    if (!is.null(highlight)) d2$highlight <- y[[highlight]]
    d2 <- subset(d2, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d2 <- d2[order(d2$CHR, d2$BP), ]
    if (logp) {
      d2$logp <- log10(d2$P) + ymid
      d2 <- d2[d2$logp < 0,]
    } else {
      d2$logp <- -d2$P + ymid
      d2 <- d2[d2$logp < 0,]
    }

  d <- data.frame(rbind(d1,d2))
  
  } else {
    d <- data.frame(d1)
  }

  #d <- d[abs(d$logp) > ymid,]
  
  #----
  # Define plotting parameters
  #----

  d$pos = NA
  d$index = NA
  ind = 0
  for (i in unique(d$CHR)) {
    ind = ind + 1
    d[d$CHR == i, ]$index = ind
  }
  nchr = length(unique(d$CHR))
  if (nchr == 1) {
    d$pos = d$BP
    ticks = floor(length(d$pos))/2 + 1
    xlabel = paste("Chromosome", unique(d$CHR), "position")
    labs = ticks
  } else {
    lastbase = 0
    ticks = NULL
    for (i in unique(d$index)) {
      if (i == 1) {
        d[d$index == i, ]$pos = d[d$index == i, ]$BP
      }
      else {
        lastbase = lastbase + tail(subset(d, index ==
                                            i - 1)$BP, 1)
        d[d$index == i, ]$pos = d[d$index == i, ]$BP +
          lastbase
      }
      ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index ==
                                                             i, ]$pos))/2 + 1)
    }
    xlabel = "Chromosome"
    labs <- unique(d$CHR)
  }
  xmax = ceiling(max(d$pos) * 1.03)
  xmin = floor(max(d$pos) * -0.03)

  print("##CALL PLOT")
  
  def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", yaxt = "n", cex.axis = 2, cex.lab = 2,
                   las = 1, pch = 20, xlim = c(xmin, xmax), ylim = c(ymin-1,ymax+1), xlab = xlabel, ylab = expression(-log[10](italic(p))),mar=5*c(5.1,6.1,4.1,2.1))
  dotargs <- list(...)
  do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))

  if (!is.null(chrlabs)) {
    if (is.character(chrlabs)) {
      if (length(chrlabs) == length(labs)) {
        labs <- chrlabs
      }
      else {
        warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
      }
    }
    else {
      warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
    }
  }
  yticks <- seq(-5*ceiling(abs(ymin/5)), 5*ceiling(abs(ymax/5)), by = 5)
  yticks <- sort(c(yticks+sign(yticks)*(5-ymid),5-ymid,-5+ymid))
  if(is.null(y)) yticks <- yticks[yticks>=0]
  ylabels <- abs(yticks)+ymid
  
  if ((abs(suggestiveline) - ymid) > 0) {
    abline(h = suggestiveline - ymid, col = "blue")
    if(!is.null(y)) abline(h = -suggestiveline + ymid, col = "blue")
    }
  if ((abs(genomewideline) - ymid) > 0) {
    abline(h = genomewideline - ymid, col = "red")
    if(!is.null(y)) abline(h = -genomewideline + ymid, col = "red")
    }
  if(!is.null(y)) abline(h = 0, col = "black")
  #if(!is.null(y)) legend(0, y = ymax - ymid - 5*strheight(1), x.name, fill = col1[1], bty="n", cex = 3, title=NULL, yjust=0.5)
  #if(!is.null(y)) legend(0, y = -ymax + ymid + 5*strheight(1), y.name, fill = col2[1], bty="n", cex = 3, title=NULL, yjust=0.5)
  if(!is.null(y)) legend(0, y = -ymax + ymid + 5*strheight(1), c(x.name,y.name), fill = c(col1[1],col2[1]), bty="n", cex = 3, title=NULL, yjust=0.5)

  if (nchr == 1) {
    axis(1, cex.axis = 1.5, cex.lab = 1.5, las=1, ...)
    axis(2, at = yticks, labels = ylabels, las=1, cex.axis=1.5, cex.lab = 1.5, ...)
  } else {
    axis(1, at = ticks, labels = labs, las=1, cex.axis = 1.5, cex.lab = 1.5, ...)
    axis(2, at = yticks, labels = ylabels, las=1, cex.axis = 1.5, cex.lab = 1.5, ...)
  }
  col1 = rep(col1, max(d$CHR))
  col2 = rep(col2, max(d$CHR))
  
  d1 <- d[d$logp>0,]
  d2 <- d[d$logp<0,]

  if (nchr == 1) {
    with(d1[d1$logp<ymax,], points(pos, logp, pch = 20, col = col1[1], ...))
    with(d1[d1$logp>ymax,], points(pos, rep(ymax+0.5, length(pos)), pch = 2, col = col1[1], ...))
    with(d2[d$logp>ymin,], points(pos, logp, pch = 20, col = col2[1], ...))
    with(d2[d2$logp<ymin,], points(pos, rep(ymin-0.5, length(pos)), pch = 6, col = col2[1], ...))
  } else {
    print("##PLOT X")
    icol = 1
    for (i in unique(d1$index)) {
      
      if(!is.null(highlight)) {
        with(d1[d1$index == unique(d1$index)[i] & d1$logp<ymax & d1$highlight == 0, ], points(pos, logp, col = col1[icol], pch = 20, ...))  
        with(d1[d1$index == unique(d1$index)[i] & d1$logp<ymax & d1$highlight == 1, ], points(pos, logp, col = highlight_col, pch = 4, cex= 1.5, lwd= 3, ...))
        with(d1[d1$index == unique(d1$index)[i] & d1$logp<ymax & d1$highlight == 2, ], points(pos, logp, col = highlight_col, bg=highlight_col, pch = 24, cex= 1.5, lwd= 2, ...))
        with(d1[d1$index == unique(d1$index)[i] & d1$logp>ymax,], points(pos, rep(ymax+0.5, length(pos)), col = highlight_col, bg=highlight_col, pch = 24, cex=1.5, lwd = 2, ...))
      } else {
        with(d1[d1$index == unique(d1$index)[i] & d1$logp<ymax, ], points(pos, logp, col = col1[icol], pch = 20, ...))  
        with(d1[d1$index == unique(d1$index)[i] & d1$logp>ymax,], points(pos, rep(ymax+0.5, length(pos)), pch = 24, col = col1[icol], bg=col1[icol], cex=1.5, lwd = 2, ...))
      }
      
      
      icol = icol + 1
    }
    
    print("##PLOT Y")
    if(!is.null(y)) {

      icol = 1
      for (i in unique(d2$index)) {
      
      if(!is.null(highlight)) {
        with(d2[d2$index == unique(d2$index)[i] & d2$logp>ymin & d2$highlight == 0, ], points(pos, logp, col = col2[icol], pch = 20, ...))
        with(d2[d2$index == unique(d2$index)[i] & d2$logp<ymax & d2$highlight == 1, ], points(pos, logp, col = highlight_col, pch = 4, cex = 1.5, lwd=3, ...))
        with(d2[d2$index == unique(d2$index)[i] & d2$logp<ymax & d2$highlight == 2, ], points(pos, logp, col = highlight_col, bg=highlight_col, pch = 25, cex = 1.5, lwd=2, ...))
        with(d2[d2$index == unique(d2$index)[i] & d2$logp<ymin,], points(pos, rep(ymin-0.5, length(pos)), col = highlight_col, bg=highlight_col, pch = 25, cex=1.5, lwd=2, ...))
      } else {
        with(d2[d2$index == unique(d2$index)[i] & d2$logp>ymin, ], points(pos, logp, col = col2[icol], pch = 20, ...))
        with(d2[d2$index == unique(d2$index)[i] & d2$logp<ymin,], points(pos, rep(ymin-0.5, length(pos)), pch = 25, col = col2[icol], bg=col2[icol], cex=1.5, lwd=2, ...))
      }
      
      
      icol = icol + 1
      }
    }
  }



  #----
  # Plot names
  #----
  
  if (!is.null(label)) {
    required_cols <- c("LABEL","CHR","POS","P")
    if (is.null(colnames(label))) warning(paste("label must be a data.frame with col.names:",paste0(required_cols,collapse=" ")))
    label <- data.frame(label)
    
    if (all(required_cols %in% colnames(label))) {  
      label$pos <- label$POS
      for (i in unique(label$CHR)[!unique(label$CHR)==1]) {
        lastbase <- tail(d[d$index==(i-1),"pos"],n=1)
        label[label$CHR==i,"pos"] <- label[label$CHR==i,"pos"]+lastbase
      }
      
      #with(label, points(pos, P, col = "#974d26", pch = 4, cex=2, ...))

      label$P <- ifelse(label$P>0, label$P-ymid+2*strheight(label$LABEL), label$P+ymid-2*strheight(label$LABEL))
      label$P <- ifelse(label$P>ymax,ymax-strheight(label$LABEL)/2,label$P)
      label$P <- ifelse(label$P<ymin,ymin+strheight(label$LABEL)/2,label$P)
      label <- label[order(-sign(label$P),abs(label$P),label$pos),]
      row.names(label) <- 1:nrow(label)

      print("##POSITION LABELS")
      
      x.buffer <- strwidth("A")
      y.buffer <- strheight("A")

      label$YMIN <- label$P - strheight(label$LABEL)/2
      label$YMAX <- label$P + strheight(label$LABEL)/2

      label$XMIN <- label$pos - strwidth(label$LABEL)/2 - x.buffer
      label$XMAX <- label$pos + strwidth(label$LABEL)/2 + x.buffer
      
      forbidden <- data.frame(data.table(d)[highlight==TRUE,.(pos,P=logp,XMIN=pos-strwidth("A")/2, XMAX=pos+strwidth("A")/2, YMIN=ifelse(logp<ymin,ymin,ifelse(logp>ymax,ymax,logp - sign(logp) * strheight("A")/2)), YMAX=ifelse(logp>ymax, ymax, ifelse(logp<ymin,ymin,logp + sign(logp) * strheight("A")/2)))])

      #ggplot(label, aes(x=pos, xmin=XMIN, xmax=XMAX, y=(YMIN+YMAX)/2, ymin=YMIN, ymax=YMAX)) + geom_point(data=forbidden, aes(x=pos, y=YMAX), colour="red") + geom_rect() + geom_text(aes(label=LABEL), size=4, color="white") 

      rect_overlap <- function(df1, df2) {
        return(!(df2$XMAX <= df1$XMIN | df2$XMIN >= df1$XMAX | df2$YMAX <= df1$YMIN | df2$YMIN >= df1$YMAX))
      }


      for (precision in 1:10) {
      
      for (i in 1:nrow(label)) {
      

      n <- 0
      repeat{
      overlap_highlight <- forbidden[rect_overlap(df1=label[i,], df2=forbidden),]
      
      if(nrow(overlap_highlight) > 0) {
        if(all(overlap_highlight$P > 0)) {
          label[i,"YMAX"] <- max(overlap_highlight$YMAX) + strheight("A")/2*3 + y.buffer
          label[i,"YMIN"] <- max(overlap_highlight$YMAX) + strheight("A")/2 + y.buffer
        } else {
          label[i,"YMIN"] <- min(overlap_highlight$YMIN) - strheight("A")/2*3 - y.buffer
          label[i,"YMAX"] <- min(overlap_highlight$YMIN) - strheight("A")/2 - y.buffer
        }
      }

      overlap <- rbind(label[i,], label[-i,][rect_overlap(df1=label[i,], df2=label[-i,]),])
      

      if(nrow(overlap) == 1 ) break

      n <- n+1
      if(n == 1) cat("\n",label[i,"LABEL"],"")
      cat(n,"")
      

      if(all(overlap$P > 0)) {
        overlap[1,"YMAX"] <- max(overlap[-1,"YMAX"]) + strheight("A")/2*3 + y.buffer
        overlap[1,"YMIN"] <- max(overlap[-1,"YMAX"]) + strheight("A")/2 + y.buffer
      } else {
        overlap[1,"YMIN"] <- min(overlap[-1,"YMIN"]) - strheight("A")/2*3 - y.buffer
        overlap[1,"YMAX"] <- min(overlap[-1,"YMIN"]) - strheight("A")/2 - y.buffer
      }

      label[match(paste0(overlap$LABEL,overlap$P),paste0(label$LABEL,label$P)),] <- overlap
      }
      }
    }
    cat("\n")

      label$P <- ifelse(label$P>ymax,ymax-strheight(label$LABEL)/6,label$P)
      label$P <- ifelse(label$P<ymin,ymin+strheight(label$LABEL)/6,label$P)
      moved <- data.table(label)[P!=(YMIN+YMAX)/2, .(x=pos, y1=P, y2=(YMIN+YMAX)/2)]

      for(i in 1:nrow(moved)) lines(x=moved[i,x] + c(0,0), y=c(moved[i,y1], moved[i,y2]) - sign(moved[i,y1]) * strheight("A"))


      shadowtext <- function(x, y=NULL, labels, col='white', bg='black', 
                       theta= seq(0, 2*pi, length.out=50), r=0.1, ... ) {
        xy <- xy.coords(x,y)
        xo <- r*strwidth('A')
        yo <- r*strheight('A')
      # draw background text with small shift in x and y in background colour
      for (i in theta) {
        text( xy$x + cos(i)*xo, xy$y + sin(i)*yo, labels, col=bg, ... )
      }
        # draw actual text in exact xy position in foreground colour
        text(xy$x, xy$y, labels, col=col, ... )
      }
      
      print("##DRAW LABELS")
      with(label, shadowtext(x=pos,y=(YMIN+YMAX)/2,labels=LABEL,cex=1.25,col="black",bg="white",font=4,r=0.2))

     } else {
      warning(paste("label data.frame is missing columns:",paste0(required_cols[!required_cols%in%colnames(label)],collapse=" ")))
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

all_files = Sys.glob("/scratch/project_2007428/projects/prj_001_cost_gwas/outputs/cohort_sumstats/*/*.gz")
all_files = all_files[!(grepl("aligned", all_files))]

#for (i in 1:length(all_files)){ # FOR ALL
for (i in 147:length(all_files)){ # FOR SUBSET

  ### Reading file and getting variable and column names

  curr_file_name = basename(all_files[i])

  cat(paste0("\nRunning on ", curr_file_name, ". CHECK1.", i, " COMPLETE.\n\n==================\n\n"))

  curr_pheno = fread(all_files[i], data.table = FALSE, tmpdir = "/scratch/project_2007428/projects/prj_001_cost_gwas/tmpdir/")

  cohort = strsplit(curr_file_name, "\\.")[[1]][1]
  trait_name = strsplit(curr_file_name, "\\.")[[1]][3]
  stratif_name_1 = strsplit(curr_file_name, "\\.")[[1]][4]
  stratif_name_2 = strsplit(curr_file_name, "\\.")[[1]][5]
  sample_size = as.numeric(strsplit(curr_file_name, "\\.")[[1]][7])

  population = strsplit(curr_file_name, "\\.")[[1]][6]

  if (cohort == "QGP"){

      realColumns = colnames(curr_pheno)[-(which(colnames(curr_pheno) == "SNPID"))]

      curr_pheno[,ncol(curr_pheno)] = NULL

      colnames(curr_pheno) = realColumns

      print(head(curr_pheno))

  }

  if (is.na(sample_size)){

    sample_size = as.numeric(strsplit(curr_file_name, "\\.")[[1]][8])

  }

  chrom_col = case_when(cohort == "UKB" ~ "CHROM",
                          cohort == "GNH" ~ "CHROM",
                          .default = "CHR")

  pos_col = case_when(cohort == "UKB" ~ "GENPOS",
                        cohort == "QGP" ~ "POS",
                        cohort == "EstBB" ~ "POS",
                        cohort == "GNH" ~ "GENPOS",
                        .default = "BP")

  if (cohort == "UKB" | cohort == "GNH"){

    curr_pheno$Pval = 10 ^ (-curr_pheno$LOG10P)

  }

  p_col = case_when(cohort == "CHB" ~ "P",
                      cohort == "QGP" ~ "p.value",
                      cohort == "EstBB" ~ "p.value",
                      cohort == "GE" ~ "Pvalue",
                      .default = "Pval")

  freq_col = case_when(cohort == "CHB" ~ "A1FREQ",
                        cohort == "UKB" ~ "A1FREQ",
                        cohort == "QGP" ~ "AF_Allele2",
                        cohort == "EstBB" ~ "A1FREQ",
                        cohort == "GNH" ~ "A1FREQ",
                        .default = "EAF")

  info_col = case_when(cohort == "QGP" ~ "imputationInfo",
                        .default = "INFO")

  if (cohort == "EstBB"){

    colnames(curr_pheno)[3] = "rsid"

  }

  curr_pheno$SNPID = paste0(curr_pheno[, chrom_col], "_", curr_pheno[, pos_col])

  trait_name = gsub("INPATIENTS", "IN", trait_name)
  trait_name = gsub("INPATIENT", "IN", trait_name)
  trait_name = gsub("TOTAL", "ALL", trait_name)
  
  stratif_name_1 = gsub("75plus", "76_95", stratif_name_1)

  nice_name = case_when(trait_name == "ALL" ~ "total",
                          trait_name == "INOUT" ~ "inpatient + outpatient",
                          trait_name == "DRUG" ~ "prescription",
                          trait_name == "PRIM" ~ "primary care",
                          trait_name == "IN" ~ "inpatient",
                          .default = as.character(trait_name))

  stratif_name = case_when(stratif_name_1 == "ALL" & stratif_name_2 == "ALL" ~ "ALL",
                              stratif_name_1 != "ALL" & stratif_name_2 == "ALL" ~ stratif_name_1,
                              stratif_name_1 == "ALL" & stratif_name_2 != "ALL" ~ stratif_name_2,
                              .default = as.character(stratif_name_1))

  nice_stratif_name = case_when(stratif_name == "36_55" ~ "age between 36-55",
                                  stratif_name == "19_35" ~ "age between 19-35",
                                  stratif_name == "5_18" ~ "age between 5-18",
                                  stratif_name == "56_75" ~ "age between 56-75",
                                  stratif_name == "76_95" ~ "age 76+",
                                  stratif_name == "M" ~ "males",
                                  stratif_name == "F" ~ "females",
                                  .default = as.character(stratif_name_1))

  cat(paste0("Variable names collected.\nCohort: ", cohort,
              "\nTrait: ", nice_name,
              "\nStratification: ", nice_stratif_name,
              "\nFrequency column: ", freq_col,
              "\nCHECK2.", i, " COMPLETE.\n\n==================\n\n"))

  print(nrow(curr_pheno))
  print(summary(curr_pheno))

  curr_pheno[, chrom_col] = gsub("X", 23, curr_pheno[, chrom_col])

  p_val_thresh = 0.01

  curr_pheno[, p_col] = as.numeric(curr_pheno[, p_col])
  curr_pheno[, freq_col] = as.numeric(curr_pheno[, freq_col])
  curr_pheno[, pos_col] = as.numeric(curr_pheno[, pos_col])
  curr_pheno[, chrom_col] = as.numeric(curr_pheno[, chrom_col])

  if (any(is.na(curr_pheno[, freq_col])) | any(is.na(curr_pheno[, p_col])) | any(is.na(curr_pheno[, pos_col])) | any(is.na(curr_pheno[, chrom_col]))){

    curr_pheno = curr_pheno[complete.cases(curr_pheno[,c(p_col, freq_col, pos_col, chrom_col)]),]

  }

  cat(paste0("\nColumns switched to numeric. NAs removed\n\n"))

  print(nrow(curr_pheno))
  print(summary(curr_pheno))

  if (any(curr_pheno[,p_col] == 0)){

    curr_pheno[curr_pheno[,p_col] == 0, p_col] = 1e-300

  }

  curr_pheno_manhattan = curr_pheno[curr_pheno[, p_col] < p_val_thresh,]
  curr_pheno_manhattan = curr_pheno_manhattan[curr_pheno_manhattan[, freq_col] >= 0.01,]

  ### NTR has no INFO column

  if (cohort != "NTR" & cohort != "AOU" & cohort != "GNH"){

    curr_pheno[, info_col] = as.numeric(curr_pheno[, info_col])
    curr_pheno_manhattan = curr_pheno_manhattan[curr_pheno_manhattan[, info_col] >= 0.6,]

  }

  if (any(curr_pheno_manhattan[, freq_col] > 0.99)){

    curr_pheno_manhattan = curr_pheno_manhattan[curr_pheno_manhattan[, freq_col] <= 0.99,]

  }

  if (any(curr_pheno[, freq_col] > 0.99)){

    curr_pheno = curr_pheno[curr_pheno[, freq_col] <= 0.99,]

  }

  cat(paste0("\nAF and INFO filter applied.\n\n"))

  print(nrow(curr_pheno_manhattan))
  print(summary(curr_pheno_manhattan))

  if (any(is.na(curr_pheno_manhattan[,chrom_col])) | any(is.na(curr_pheno_manhattan[,pos_col])) | any(is.na(curr_pheno_manhattan[,p_col]))){

    curr_pheno_manhattan = curr_pheno_manhattan[complete.cases(curr_pheno_manhattan[,c(chrom_col, pos_col, p_col)]),]

  }

  ### COLOUR DEFINITION
  ### CURRENT STANDARD IS:
  ### TOTAL/GP = BLACK/GREEN
  ### HES/PRESCRIPT = BLUE/ORANGE

  green_set = c("#006d2c", "#238b45", "#41ab5d") # Primary
  blue_set = c("#2166ac","#2171b5","#4393c3") # Inpatient
  purple_set = c("#542788","#6a51a3","#8073ac") # INOUT patient
  black_set = c("#252525", "#525252", "#737373") # Total
  orange_set = c("#ff7b00", "#ffa500", "#ffba00") # Prescription

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

  if (cohort == "UKB"){

    cohort = paste0(cohort, "_", population)

    nice_cohort_name = paste0(cohort, " ", population)

  } else {

    nice_cohort_name = cohort

  }

  fig_name = ifelse(stratif_name == "ALL", paste0(trait_name, "_costs_", cohort, "_n_", sample_size), paste0(trait_name, "_costs_", stratif_name, "_", cohort, "_n_", sample_size))
  legend_name = ifelse(stratif_name == "ALL", paste0(nice_cohort_name, " ", nice_name, " costs, N = ", sample_size), paste0(nice_cohort_name, " ", nice_name, " costs (stratified by ", nice_stratif_name, "), N = ", sample_size))

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
      chr = chrom_col, bp = pos_col, snp = "SNPID", p = p_col,  
      genomewideline = -log10(5e-8),
      ymax = max_p,
      ymin = 0,
      col1 = colour_use,
      x.name = legend_name
      )

  dev.off()

  cat(paste0("Manhattans made. Doing QQ plot. CHECK4.", i, " COMPLETE.\n\n==================\n\n"))

  png(paste0(qq_figDir, fig_name, ".png"), width = 1000, height = 1000)

  plot_QQ(gwas_results = curr_pheno, p_col = p_col, freq_col = freq_col)

  dev.off()

  cat(paste0("QQs made. Doing QQ plot. CHECK5.", i, " COMPLETE.\n\n==================\n\n"))

}


