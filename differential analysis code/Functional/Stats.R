
vkey2 <- function (map, title = NA, side = 2, stretch = 1.4, x, y, 
		wh)
{
	opar <- par(xpd = NA)
	on.exit(par(opar))
	n <- length(map$breaks) + 1
	dy <- strheight("A")
	aspect <- diff(grconvertX(1:2, from = "inches"))/diff(grconvertY(1:2,
					from = "inches"))
	dx <- dy * aspect
	if (missing(wh)) {
		
		wh <- 1:(n-1)
		
	}
	labs <- format(map$breaks[wh])
	maxlabwidth <- max(strwidth(labs))
	if (missing(x)) {
		x <- grconvertX(1, from = "nfc") - (2 * dx)
		if (side == 4)
			x <- x - maxlabwidth - dx
	}
	else {
		if (is.list(x)) {
			y <- x$y
			x <- x$x
		}
	}
	if (missing(y))
		y <- par("usr")[3] + dy
	ybord <- y + ((0:(n - 1)) * dy * stretch)
	rect(x, ybord[-n], x + dx, ybord[-1], col = map$colors, border = NA)
	if (side == 4) {
		xtext <- x + dx
		text(x = x, y = ybord[n] + (1.5 * dy), title, adj = c(0,
						0))
	}
	if (side == 2) {
		xtext <- x
		text(x = x + dx, y = ybord[n] + (1.5 * dy), title, adj = c(1,
						0))
	}
	text(x = xtext, y = ybord[wh] + 0.5 * dy, labels = labs, pos = side)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	require(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
				ncol = cols, nrow = ceiling(numPlots/cols), byrow=TRUE)
	}
	
	if (numPlots==1) {
		print(plots[[1]])
		
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
							layout.pos.col = matchidx$col))
		}
	}
}


heatmap.3 <- function(x,
		Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
		distfun = dist,
		hclustfun = hclust,
		dendrogram = c("both","row", "column", "none"),
		symm = FALSE,
		scale = c("none","row", "column"),
		na.rm = TRUE,
		revC = identical(Colv,"Rowv"),
		add.expr,
		breaks,
		symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
		col = "heat.colors",
		colsep,
		rowsep,
		sepcolor = "white",
		sepwidth = c(0.05, 0.05),
		cellnote,
		notecex = 1,
		notecol = "cyan",
		na.color = par("bg"),
		trace = c("none", "column","row", "both"),
		tracecol = "cyan",
		hline = median(breaks),
		vline = median(breaks),
		linecol = tracecol,
		margins = c(5,5),
		ColSideColors = NULL,
		RowSideColors = NULL,
		side.height.fraction=0.3,
		cexRow = 0.2 + 1/log10(nr),
		cexCol = 0.2 + 1/log10(nc),
		labRow = NULL,
		labCol = NULL,
		key = TRUE,
		keysize = 1.5,
		density.info = c("none", "histogram", "density"),
		denscol = tracecol,
		symkey = max(x < 0, na.rm = TRUE) || symbreaks,
		densadj = 0.25,
		main = NULL,
		xlab = NULL,
		ylab = NULL,
		lmat = NULL,
		lhei = NULL,
		lwid = NULL,
		NumColSideColors = 1,
		NumRowSideColors = 1,
		KeyValueName="Value",...){
	
	invalid <- function (x) {
		if (missing(x) || is.null(x) || length(x) == 0)
			return(TRUE)
		if (is.list(x))
			return(all(sapply(x, invalid)))
		else if (is.vector(x))
			return(all(is.na(x)))
		else return(FALSE)
	}
	
	x <- as.matrix(x)
	scale01 <- function(x, low = min(x), high = max(x)) {
		x <- (x - low)/(high - low)
		x
	}
	retval <- list()
	scale <- if (symm && missing(scale))
				"none"
			else match.arg(scale)
	dendrogram <- match.arg(dendrogram)
	trace <- match.arg(trace)
	density.info <- match.arg(density.info)
	if (length(col) == 1 && is.character(col))
		col <- get(col, mode = "function")
	if (!missing(breaks) && (scale != "none"))
		warning("Using scale=\"row\" or scale=\"column\" when breaks are",
				"specified can produce unpredictable results.", "Please consider using only one or the other.")
	if (is.null(Rowv) || is.na(Rowv))
		Rowv <- FALSE
	if (is.null(Colv) || is.na(Colv))
		Colv <- FALSE
	else if (Colv == "Rowv" && !isTRUE(Rowv))
		Colv <- FALSE
	if (length(di <- dim(x)) != 2 || !is.numeric(x))
		stop("`x' must be a numeric matrix")
	nr <- di[1]
	nc <- di[2]
	if (nr <= 1 || nc <= 1)
		stop("`x' must have at least 2 rows and 2 columns")
	if (!is.numeric(margins) || length(margins) != 2)
		stop("`margins' must be a numeric vector of length 2")
	if (missing(cellnote))
		cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
	if (!inherits(Rowv, "dendrogram")) {
		if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
					c("both", "row"))) {
			if (is.logical(Colv) && (Colv))
				dendrogram <- "column"
			else dedrogram <- "none"
			warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
					dendrogram, "'. Omitting row dendogram.")
		}
	}
	if (!inherits(Colv, "dendrogram")) {
		if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
					c("both", "column"))) {
			if (is.logical(Rowv) && (Rowv))
				dendrogram <- "row"
			else dendrogram <- "none"
			warning("Discrepancy: Colv is FALSE, while dendrogram is `",
					dendrogram, "'. Omitting column dendogram.")
		}
	}
	if (inherits(Rowv, "dendrogram")) {
		ddr <- Rowv
		rowInd <- order.dendrogram(ddr)
	}
	else if (is.integer(Rowv)) {
		hcr <- hclustfun(distfun(x))
		ddr <- as.dendrogram(hcr)
		ddr <- reorder(ddr, Rowv)
		rowInd <- order.dendrogram(ddr)
		if (nr != length(rowInd))
			stop("row dendrogram ordering gave index of wrong length")
	}
	else if (isTRUE(Rowv)) {
		Rowv <- rowMeans(x, na.rm = na.rm)
		hcr <- hclustfun(distfun(x))
		ddr <- as.dendrogram(hcr)
		ddr <- reorder(ddr, Rowv)
		rowInd <- order.dendrogram(ddr)
		if (nr != length(rowInd))
			stop("row dendrogram ordering gave index of wrong length")
	}
	else {
		rowInd <- nr:1
	}
	if (inherits(Colv, "dendrogram")) {
		ddc <- Colv
		colInd <- order.dendrogram(ddc)
	}
	else if (identical(Colv, "Rowv")) {
		if (nr != nc)
			stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
		if (exists("ddr")) {
			ddc <- ddr
			colInd <- order.dendrogram(ddc)
		}
		else colInd <- rowInd
	}
	else if (is.integer(Colv)) {
		hcc <- hclustfun(distfun(if (symm)
									x
								else t(x)))
		ddc <- as.dendrogram(hcc)
		ddc <- reorder(ddc, Colv)
		colInd <- order.dendrogram(ddc)
		if (nc != length(colInd))
			stop("column dendrogram ordering gave index of wrong length")
	}
	else if (isTRUE(Colv)) {
		Colv <- colMeans(x, na.rm = na.rm)
		hcc <- hclustfun(distfun(if (symm)
									x
								else t(x)))
		ddc <- as.dendrogram(hcc)
		ddc <- reorder(ddc, Colv)
		colInd <- order.dendrogram(ddc)
		if (nc != length(colInd))
			stop("column dendrogram ordering gave index of wrong length")
	}
	else {
		colInd <- 1:nc
	}
	retval$rowInd <- rowInd
	retval$colInd <- colInd
	retval$call <- match.call()
	x <- x[rowInd, colInd]
	x.unscaled <- x
	cellnote <- cellnote[rowInd, colInd]
	if (is.null(labRow))
		labRow <- if (is.null(rownames(x)))
					(1:nr)[rowInd]
				else rownames(x)
	else labRow <- labRow[rowInd]
	if (is.null(labCol))
		labCol <- if (is.null(colnames(x)))
					(1:nc)[colInd]
				else colnames(x)
	else labCol <- labCol[colInd]
	if (scale == "row") {
		retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
		x <- sweep(x, 1, rm)
		retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
		x <- sweep(x, 1, sx, "/")
	}
	else if (scale == "column") {
		retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
		x <- sweep(x, 2, rm)
		retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
		x <- sweep(x, 2, sx, "/")
	}
	if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
		if (missing(col) || is.function(col))
			breaks <- 16
		else breaks <- length(col) + 1
	}
	if (length(breaks) == 1) {
		if (!symbreaks)
			breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
					length = breaks)
		else {
			extreme <- max(abs(x), na.rm = TRUE)
			breaks <- seq(-extreme, extreme, length = breaks)
		}
	}
	nbr <- length(breaks)
	ncol <- length(breaks) - 1
	if (class(col) == "function")
		col <- col(ncol)
	min.breaks <- min(breaks)
	max.breaks <- max(breaks)
	x[x < min.breaks] <- min.breaks
	x[x > max.breaks] <- max.breaks
	if (missing(lhei) || is.null(lhei))
		lhei <- c(keysize, 4)
	if (missing(lwid) || is.null(lwid))
		lwid <- c(keysize, 4)
	if (missing(lmat) || is.null(lmat)) {
		lmat <- rbind(4:3, 2:1)
		
		if (!is.null(ColSideColors)) {
			#if (!is.matrix(ColSideColors))
			#stop("'ColSideColors' must be a matrix")
			if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
				stop("'ColSideColors' must be a matrix of nrow(x) rows")
			lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
			#lhei <- c(lhei[1], 0.2, lhei[2])
			lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
		}
		
		if (!is.null(RowSideColors)) {
			#if (!is.matrix(RowSideColors))
			#stop("'RowSideColors' must be a matrix")
			if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
				stop("'RowSideColors' must be a matrix of ncol(x) columns")
			lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
			#lwid <- c(lwid[1], 0.2, lwid[2])
			lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
		}
		lmat[is.na(lmat)] <- 0
	}
	
	if (length(lhei) != nrow(lmat))
		stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
	if (length(lwid) != ncol(lmat))
		stop("lwid must have length = ncol(lmat) =", ncol(lmat))
	op <- par(no.readonly = TRUE)
	on.exit(par(op))
	
	layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
	
	if (!is.null(RowSideColors)) {
		if (!is.matrix(RowSideColors)){
			par(mar = c(margins[1], 0, 0, 0.5))
			image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
			box()
		} else {
			par(mar = c(margins[1], 0, 0, 0.5))
			rsc = t(RowSideColors[,rowInd, drop=F])
			rsc.colors = matrix()
			rsc.names = names(table(rsc))
			rsc.i = 1
			for (rsc.name in rsc.names) {
				rsc.colors[rsc.i] = rsc.name
				rsc[rsc == rsc.name] = rsc.i
				rsc.i = rsc.i + 1
			}
			rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
			image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
			
			if (length(rownames(RowSideColors)) > 0) {
				axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), rownames(RowSideColors), las = 2, tick = FALSE)
			}
			box()
		}
	}
	
	if (!is.null(ColSideColors)) {
		
		if (!is.matrix(ColSideColors)){
			par(mar = c(0.5, 0, 0, margins[2]))
			image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
			box()
		} else {
			par(mar = c(0.5, 0, 0, margins[2]))
			csc = ColSideColors[colInd, , drop=F]
			csc.colors = matrix()
			csc.names = names(table(csc))
			csc.i = 1
			for (csc.name in csc.names) {
				csc.colors[csc.i] = csc.name
				csc[csc == csc.name] = csc.i
				csc.i = csc.i + 1
			}
			csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
			image(csc, col = as.vector(csc.colors), axes = FALSE)
			if (length(colnames(ColSideColors)) > 0) {
				axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
			}
			box()
		}
	}
	
	par(mar = c(margins[1], 0, 0, margins[2]))
	x <- t(x)
	cellnote <- t(cellnote)
	if (revC) {
		iy <- nr:1
		if (exists("ddr"))
			ddr <- rev(ddr)
		x <- x[, iy]
		cellnote <- cellnote[, iy]
	}
	else iy <- 1:nr
	image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
	retval$carpet <- x
	if (exists("ddr"))
		retval$rowDendrogram <- ddr
	if (exists("ddc"))
		retval$colDendrogram <- ddc
	retval$breaks <- breaks
	retval$col <- col
	if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
		mmat <- ifelse(is.na(x), 1, NA)
		image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
				col = na.color, add = TRUE)
	}
	axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
			cex.axis = cexCol)
	if (!is.null(xlab))
		mtext(xlab, side = 1, line = margins[1] - 1.25)
	axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
			cex.axis = cexRow)
	if (!is.null(ylab))
		mtext(ylab, side = 4, line = margins[2] - 1.25)
	if (!missing(add.expr))
		eval(substitute(add.expr))
	if (!missing(colsep))
		for (csep in colsep) rect(xleft = csep + 0.5, ybottom = 0, xright = csep + 0.5 + sepwidth[1], ytop = ncol(x) + 1, lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
	if (!missing(rowsep))
		for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
	min.scale <- min(breaks)
	max.scale <- max(breaks)
	x.scaled <- scale01(t(x), min.scale, max.scale)
	if (trace %in% c("both", "column")) {
		retval$vline <- vline
		vline.vals <- scale01(vline, min.scale, max.scale)
		for (i in colInd) {
			if (!is.null(vline)) {
				abline(v = i - 0.5 + vline.vals, col = linecol,
						lty = 2)
			}
			xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
			xv <- c(xv[1], xv)
			yv <- 1:length(xv) - 0.5
			lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
		}
	}
	if (trace %in% c("both", "row")) {
		retval$hline <- hline
		hline.vals <- scale01(hline, min.scale, max.scale)
		for (i in rowInd) {
			if (!is.null(hline)) {
				abline(h = i + hline, col = linecol, lty = 2)
			}
			yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
			yv <- rev(c(yv[1], yv))
			xv <- length(yv):1 - 0.5
			lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
		}
	}
	if (!missing(cellnote))
		text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
				col = notecol, cex = notecex)
	par(mar = c(margins[1], 0, 0, 0))
	if (dendrogram %in% c("both", "row")) {
		plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
	}
	else plot.new()
	par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
	if (dendrogram %in% c("both", "column")) {
		plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
	}
	else plot.new()
	if (!is.null(main))
		title(main, cex.main = 1.5 * op[["cex.main"]])
	if (key) {
		par(mar = c(5, 4, 2, 1), cex = 0.75)
		tmpbreaks <- breaks
		if (symkey) {
			max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
			min.raw <- -max.raw
			tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
			tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
		}
		else {
			min.raw <- min(x, na.rm = TRUE)
			max.raw <- max(x, na.rm = TRUE)
		}
		
		z <- seq(min.raw, max.raw, length = length(col))
		image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
				xaxt = "n", yaxt = "n")
		par(usr = c(0, 1, 0, 1))
		lv <- pretty(breaks)
		xv <- scale01(as.numeric(lv), min.raw, max.raw)
		axis(1, at = xv, labels = lv)
		if (scale == "row")
			mtext(side = 1, "Row Z-Score", line = 2)
		else if (scale == "column")
			mtext(side = 1, "Column Z-Score", line = 2)
		else mtext(side = 1, KeyValueName, line = 2)
		if (density.info == "density") {
			dens <- density(x, adjust = densadj, na.rm = TRUE)
			omit <- dens$x < min(breaks) | dens$x > max(breaks)
			dens$x <- dens$x[-omit]
			dens$y <- dens$y[-omit]
			dens$x <- scale01(dens$x, min.raw, max.raw)
			lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
					lwd = 1)
			axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
			title("Color Key\nand Density Plot")
			par(cex = 0.5)
			mtext(side = 2, "Density", line = 2)
		}
		else if (density.info == "histogram") {
			h <- hist(x, plot = FALSE, breaks = breaks)
			hx <- scale01(breaks, min.raw, max.raw)
			hy <- c(h$counts, h$counts[length(h$counts)])
			lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
					col = denscol)
			axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
			title("Color Key\nand Histogram")
			par(cex = 0.5)
			mtext(side = 2, "Count", line = 2)
		}
		else title("Color Key")
	}
	else plot.new()
	retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
			high = retval$breaks[-1], color = retval$col)
	invisible(retval)
}

load_package <- function (package.names=c('MASS', 'ggbiplot', 'ape', 'vegan', 'pscl', 'glmmADMB', 'aod', 'nlme', 'MiRKAT', 'matrixStats', 
				'gplots', 'scales', 'ggplot2', 'GUniFrac', 'rpart', 'qvalue', 'DESeq2',
				'phangorn', 'phyloseq', 'RColorBrewer', 'squash', 'rhdf5', 'biom', 'reshape', 'randomForest', 'Boruta', 'ade4', 'Daim', 'rlocal')) {
	for (package.name in package.names) {
		require(package.name, character.only=T)
	}
}

getPermuteMatrixBlock <- function(permutations, strata) {
	strata <- factor(strata)
	strata <- factor(as.numeric(strata))
	res <- sapply(1:permutations, function(i) {
				strata1 <- strata
				levels(strata1) <- sample(levels(strata1))
				order(as.numeric(as.character(strata1)))
			})
	t(res)
}

# adonis2: permutate the covariate instead of data matrix - validated
# The variable of interest should be put at the last; Only p values of variable of interest will be reported
# Currently it can only be applied to subject-specific variable (nonvarying within subject) but covariates are allowed to vary
# Need to improve the speed
adonis2 <- function (formula, data = NULL, permutations = 999, method = "bray", 
		strata = NULL, block.perm = TRUE, contr.unordered = "contr.sum", contr.ordered = "contr.poly", 
		...) 
{
	TOL <- 1e-07
	Terms <- terms(formula, data = data, keep.order=TRUE)
	lhs <- formula[[2]]
	lhs <- eval(lhs, data, parent.frame())
	formula[[2]] <- NULL
	rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
	op.c <- options()$contrasts
	options(contrasts = c(contr.unordered, contr.ordered))
	rhs <- model.matrix(formula, rhs.frame)
	options(contrasts = op.c)
	grps <- attr(rhs, "assign")
	qrhs <- qr(rhs)
	rhs <- rhs[, qrhs$pivot, drop = FALSE]
	rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
	grps <- grps[qrhs$pivot][1:qrhs$rank]
	u.grps <- unique(grps)
	nterms <- length(u.grps) - 1
	st <- ifelse(nterms == 1, 2, nterms)
	H.s <- lapply(st:(nterms+1), function(j) {
				Xj <- rhs[, grps %in% u.grps[1:j]]
				qrX <- qr(Xj, tol = TOL)
				Q <- qr.Q(qrX)
				tcrossprod(Q[, 1:qrX$rank])
			})
	if (inherits(lhs, "dist")) {
		if (any(lhs < -TOL)) 
			stop("dissimilarities must be non-negative")
		dmat <- as.matrix(lhs^2)
	}
	else {
		dist.lhs <- as.matrix(vegdist(lhs, method = method, ...))
		dmat <- dist.lhs^2
	}
	n <- nrow(dmat)
	G <- -sweep(dmat, 1, rowMeans(dmat))/2
	SS.Exp.comb <- sapply(H.s, function(hat) sum(G * t(hat)))
	if (nterms == 1) {
		SS.Exp.each <- SS.Exp.comb
	} else {
		SS.Exp.each <- SS.Exp.comb[2] - SS.Exp.comb[1]
	}
	
	H.snterm <- H.s[[length(H.s)]]
	tIH.snterm <- t(diag(n) - H.snterm)
	
#	if (length(H.s) > 1) {
#		H.s[[2]] <- H.s[[2]] - H.s[[1]]
#	}
	
	SS.Res <- sum(G * tIH.snterm)
	df.Exp <- sum(grps == u.grps[length(u.grps)])
	df.Res <- n - qrhs$rank
	
	if (inherits(lhs, "dist")) {
		beta.sites <- qr.coef(qrhs, as.matrix(lhs))
		beta.spp <- NULL
	}
	else {
		beta.sites <- qr.coef(qrhs, dist.lhs)
		beta.spp <- qr.coef(qrhs, as.matrix(lhs))
	}
	colnames(beta.spp) <- colnames(lhs)
	colnames(beta.sites) <- rownames(lhs)
	
	F.Mod <- (SS.Exp.each/df.Exp)/(SS.Res/df.Res)
	
	f.test <- function(tH, G, df.Exp, df.Res, tIH.snterm) {
		(sum(G * tH)/df.Exp)/(sum(G * tIH.snterm)/df.Res)
	}
	
	rhs.1 <- rhs[, grps %in% u.grps[1:nterms], drop=F]
	rhs.2 <- rhs[, grps %in% u.grps[(nterms+1)], drop=F]	
	
	if (missing(strata)) strata <- NULL
	
	if (block.perm == FALSE) {
		permute.ind <- vegan:::getPermuteMatrix(permutations, n, strata = strata)
	} else {
		if (is.null(strata)) stop('Block permutation requires strata!\n')
		strata.u <- unique(strata)
		reorder.ind <- unlist(lapply(strata.u, function(x) which(strata == x)))
		expand.ind <- rep(1:length(strata.u), sapply(strata.u, function(x) sum(strata == x)))
		rhs.2.u <- rhs.2[sapply(strata.u, function(x) which(strata == x)[1]), , drop=F]
	}
	
	# To be checked and validated. 
	f.perms <- 
			sapply(1:permutations, function(i) {
						if (block.perm == FALSE) {
							rhs.2 <- rhs[permute.ind[i, ], grps %in% u.grps[(nterms+1)], drop=F]	
						} else {
							rhs.2[reorder.ind, ] <- rhs.2.u[sample(nrow(rhs.2.u)), , drop=F][expand.ind, , drop=F]
						}		
						Xj <- cbind(rhs.1, rhs.2)
						qrX <- qr(Xj, tol = TOL)
						Q <- qr.Q(qrX)
						tH.snterm <- t(tcrossprod(Q[, 1:qrX$rank]))
						tIH.snterm <- diag(n) - tH.snterm
						
						if (nterms > 1) {
							tH.snterm <- tH.snterm - t(H.s[[1]])
						}							
						f.test(tH.snterm, G, df.Exp, df.Res, 
								tIH.snterm)
					}
			)
	
	f.perms <- round(f.perms, 12)
	F.Mod <- round(F.Mod, 12)
	SumsOfSqs = c(SS.Exp.each, SS.Res, SS.Exp.comb[length(SS.Exp.comb)] + SS.Res)
	tab <- data.frame(Df = c(df.Exp, df.Res, n - 1), SumsOfSqs = SumsOfSqs, 
			MeanSqs = c(SS.Exp.each/df.Exp, SS.Res/df.Res, NA), F.Model = c(F.Mod, 
					NA, NA), R2 = SumsOfSqs/SumsOfSqs[length(SumsOfSqs)], 
			P = c((sum(f.perms >= F.Mod) + 1)/(permutations + 
								1), NA, NA))
	rownames(tab) <- c(attr(attr(rhs.frame, "terms"), "term.labels")[u.grps[length(u.grps)]], 
			"Residuals", "Total")
	colnames(tab)[ncol(tab)] <- "Pr(>F)"
	attr(tab, "heading") <- "Terms added sequentially (first to last)\n"
	class(tab) <- c("anova", class(tab))
	out <- list(aov.tab = tab, call = match.call(), coefficients = beta.spp, 
			coef.sites = beta.sites, f.perms = as.matrix(f.perms), model.matrix = rhs, 
			terms = Terms)
	class(out) <- "adonis"
	out
}


########################################
# This generates the matrix columns-wise
# From JnPaulson
generate_matrix <- function(x){
	indptr  = x$sample$matrix$indptr+1
	indices = x$sample$matrix$indices+1
	data    = x$sample$matrix$data
	nr = length(x$observation$ids)
	
	counts = sapply(2:length(indptr),function(i){
				x = rep(0,nr)
				seq = indptr[i-1]:(indptr[i]-1)
				x[indices[seq]] = data[seq]
				x
			})
	rownames(counts) = x$observation$ids
	colnames(counts) = x$sample$ids
	# I wish this next line wasn't necessary
	lapply(1:nrow(counts),function(i){
				counts[i,]
			})
}
generate_metadata <- function(x){
	metadata = x$metadata
	metadata = lapply(1:length(x$ids),function(i){
				id_metadata = lapply(metadata,function(j){
							if(length(dim(j))>1){ as.vector(j[,i,drop=FALSE]) }
							else{ j[i] }
						})
				list(id = x$ids[i],metadata=id_metadata)
			})
	return(metadata)
}
namedList <- function(...) {
	L <- list(...)
	snm <- sapply(substitute(list(...)),deparse)[-1]
	if (is.null(nm <- names(L))) nm <- snm
	if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
	setNames(L,nm)
}

read_hdf5_biom <- function(file_input){
	x = h5read(file_input,"/",read.attributes = TRUE)
	data = generate_matrix(x)
	rows = generate_metadata(x$observation)
	columns = generate_metadata(x$sample)
	shape = c(length(data),length(data[[1]])) # dim(data)
	
	# Experimental -- need to actually load these from file
	id = attr(x,"id")
	vs = attr(x,"format-version")
	format = sprintf("Biological Observation Matrix %s.%s",vs[1],vs[2])
	format_url = attr(x,"format-url")
	type = "OTU table"
	#type=attr(x,"type")
	generated_by = attr(x,"generated-by")
	date = attr(x,"creation-date")
	matrix_type = "dense"
	matrix_element_type = "int"
	
	namedList(id,format,format_url,type,generated_by,date,matrix_type,matrix_element_type,
			rows,columns,shape,data)
}
############################################
# comm are otu counts, row: otus, column samples
# intersect.no: Pairwise ratio calculated on pairs with at least 'intersect.no' common taxa; 
# Rev: 2016_12_12 add prevalence filter, which is not necessary
# Rev: 2017_02_07 
GMPR.old <- function (comm, intersect.no=4, prev.filter=0.0) {
	comm <- comm[rowMeans(comm != 0) >= prev.filter, ]
	
	ind.vec <- numeric(ncol(comm))
	res <- sapply(1:ncol(comm),  function(i) {
				x <- comm[, i]
				pr <- sapply(1:ncol(comm),  function(j) {
							y <- comm[, j]
							ind <- x != 0 & y != 0
							if (sum(ind) >= intersect.no) {
								res <- median(x[ind] / y[ind])
							} else {
								res <- NA
							}
						})
				if (sum(is.na(pr)) != 0) ind.vec[i] <<- 1
				exp(mean(log(pr[!is.na(pr)])))
			}
	
	)
	if (sum(ind.vec)) warnings(paste0('Sample ', paste(which(ind.vec!=0), collapse=' '), ' do not have at least ', 
						intersect.no, ' common taxa with some of the other samples. '))
	if (sum(is.nan(res))) {
		ind <- is.nan(res)
		res[ind] <- 1
		warning(paste0('Sample ', paste(which(ind), collapse=' '), ' do not have at least ', intersect.no, ' common taxa for any of the other samples. ',
						'For these samples, we force the size factors to be 1! ', 
						'You may consider removing these samples since they are very different from other samples!\n'))
	}
	cat('Finished!')
	res
}

# New: 2017_02_07
GMPR <- function (comm, intersect.no=4, ct.min=2) {
	# Computes the GMPR size factor
	#
	# Args:
	#   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
	#   intersect.no: the minimum number of shared features between sample pair, where the ratio is calculated
	#   ct.min: the minimum number of counts required to calculate ratios ct.min = 5 has better results
	
	#
	# Returns:
	#   a list that contains:
	#      gmpr： the GMPR size factors for all samples; Samples with distinct sets of features will be output as NA.
	#      nss:   number of samples with significant sharing (> intersect.no) including itself
	
	# mask counts < ct.min
	comm[comm < ct.min] <- 0
	
	if (is.null(colnames(comm))) {
		colnames(comm) <- paste0('S', 1:ncol(comm))
	}
	
	cat('Begin GMPR size factor calculation ...\n')
	
	comm.no <- numeric(ncol(comm))
	gmpr <- sapply(1:ncol(comm),  function(i) {		
				if (i %% 50 == 0) {
					cat(i, '\n')
				}
				x <- comm[, i]
				# Compute the pairwise ratio
				pr <- x / comm
				# Handling of the NA, NaN, Inf
				pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
				# Counting the number of non-NA, NaN, Inf
				incl.no <- colSums(!is.na(pr))		
				# Calculate the median of PR
				pr.median <- colMedians(pr, na.rm=TRUE)
				# Record the number of samples used for calculating the GMPR
				comm.no[i] <<- sum(incl.no >= intersect.no)
				# Geometric mean of PR median
				if (comm.no[i] > 1) {
					return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
				} else {
					return(NA)
				}
			}
	)
	
	if (sum(is.na(gmpr))) {
		warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'), 
						'\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
						'For these samples, their size factors are set to be NA! \n', 
						'You may consider removing these samples since they are potentially outliers or negative controls!\n',
						'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
	}
	
	cat('Completed!\n')
	cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
	attr(gmpr, 'NSS') <- comm.no
	# Rev: 2017_09_07
    gmpr <- gmpr * median(colSums(comm))
	names(gmpr) <- colnames(comm)
	return(gmpr)
}


# Rev: 2016_10_25
uniquefy_taxa_names <- function (data.obj) {
	for (level in names(data.obj$abund.list)) {
		obj <- data.obj$abund.list[[level]]
#		rownames(obj) <- gsub('unclassified', paste0('Unclassified', substr(level, 1, 1)), rownames(obj))
#		rownames(obj) <- gsub('unassigned', paste0('Unclassified', substr(level, 1, 1)), rownames(obj))
		rownames(obj) <- paste0(rownames(obj), substr(level, 1, 1))
		data.obj$abund.list[[level]] <- obj
	}
#	suffix <- c('K', 'P', 'C', 'O', 'F', 'G', 'S')
#	for (i in 1:ncol(data.obj$otu.name)) {
#		data.obj$otu.name[, i] <- paste0(data.obj$otu.name[, i], suffix[i])
#	}
	return(data.obj)
}


# Rev: 2016_09_22 Add load.map
# Rev: 2016_12_12 Reogranize rarefy, normalize, winsorize
load_data <- function (otu.file, map.file, tree.file=NULL,  load.map=TRUE, parseFunction=parse_taxonomy_greengenes, version='Old', 
		 species=TRUE, filter.no=1, rep.seq=NULL,
		 rff=FALSE, dep=NULL,
	     norm='TSS', level='OTU', intersect.no=4,
		 winsor=FALSE, winsor.qt=0.97,
		 ko.file=NULL, cog.file=NULL, ko.ann.file=NULL,
		 meta.sep='\t', quote="\"", comment="",
		 read.gg=FALSE,  seed=1234, ...) {
	# ko and cog file are not rarefied	
	# filter.no: filter the OTUs with read support less than filter.no (default is filtering singleton); singleton will not be filtered after rarefaction
	# winsorization and GMPR should be further studied. Current default is false and GMPR is on the genus level
	act.seq <- NULL
	set.seed(seed)
	if (load.map == TRUE) {
		cat("Load meta file...\n")
		if (grepl("csv$", map.file)) {
			meta.dat <- read.csv(map.file, header=T, check.names=F, row.names=1, comment=comment, quote=quote, ...)
		} else {
			meta.dat <- read.table(map.file, header=T, check.names=F, row.names=1, comment=comment, sep=meta.sep, quote=quote, ...)
		}
	} else {
		meta.dat <- NULL
	}
	
	# Load Tree
	if (!is.null(tree.file)) {
		cat("Load tree file ...\n")
		if (read.gg == F) {
			tree.12 <- read.tree(tree.file)
		} else {
			tree.12 <- read_tree_greengenes(tree.file)
		}
		
		if (is.rooted(tree.12) == F) {
			tree.12 <- midpoint(tree.12)
		}
		
	} else {
		tree.12 <- NULL
	}
	
	cat("Load OTU file...\n")  # Rewrite load new biom file, rev:2016-06-20
	if (version != 'New') {
		biom.obj <-  import_biom(otu.file, parseFunction = parseFunction)  
		
		otu.tab.12 <- otu_table(biom.obj)@.Data
		otu.ind <- rowSums(otu.tab.12) > filter.no  # change otu.tab.12 != 0, rev:2016-06-20
		otu.tab.12 <- otu.tab.12[otu.ind, ]
		# OTU names
		otu.name.full <- as.matrix(biom.obj@tax_table[, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')])
		otu.name.full <- otu.name.full[otu.ind, ]
		
		otu.name.12 <- otu.name.full[, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')]
		otu.name.12[is.na(otu.name.12)] <- 'unclassified'
		otu.name.12 <- otu.name.12@.Data
		
		otu.name.12[, ] <- gsub('\\[', '', otu.name.12)
		otu.name.12[, ] <- gsub('\\]', '', otu.name.12)
		otu.name.full <- otu.name.full@.Data
	} else {
		temp <-  read_hdf5_biom(otu.file)
#		otu.file <- paste0(otu.file, '.old')
#		write_biom(temp, otu.file) 
#		biom.obj <- import_biom(otu.file, parseFunction = parseFunction)  
		otu.tab.12 <- matrix(unlist(temp$data), byrow=T, nrow=temp$shape[1], ncol=temp$shape[2])
	
		otu.ids <- sapply(temp$rows, function(x) x[['id']])
		sam.ids <- sapply(temp$columns, function(x) x[['id']])
		
#		res <- t(sapply(temp$rows, function(i) {
#					i$metadata$taxonomy
#				}))
#		rownames(res) <- otu.ids 
	
		# From phyloseq: to be double checked!! Checked!
		if (all(sapply(sapply(temp$rows, function(i) {
									i$metadata
								}), is.null))) {
			otu.name.full <- NULL
		} else {
			taxlist = lapply(temp$rows, function(i) {
						parseFunction(i$metadata$taxonomy)
					})
			names(taxlist) = sapply(temp$rows, function(i) {
						i$id
					})
			otu.name.full = build_tax_table(taxlist)
		}

		# Rev: 2017_04_18
		otu.name.full <- otu.name.full@.Data
		
		rownames(otu.tab.12) <- otu.ids
		colnames(otu.tab.12) <- sam.ids
		
		otu.ind <- rowSums(otu.tab.12) > filter.no  # change otu.tab.12 != 0, rev:2016-06-20
		otu.tab.12 <- otu.tab.12[otu.ind, ]	
		otu.name.full <- otu.name.full[otu.ind, ]
		otu.name.12 <- otu.name.full[, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')]
		otu.name.12[is.na(otu.name.12)] <- 'unclassified'
		otu.name.12[, ] <- gsub('\\[', '', otu.name.12)
		otu.name.12[, ] <- gsub('\\]', '', otu.name.12)
	}
	
	if (rff == TRUE) {
		cat('Rarefaction ...\n')
		if (is.null(dep)) {
			otu.tab.12 <- t(Rarefy(t(otu.tab.12))$otu.tab.rff)
		} else {
			otu.tab.12 <- t(Rarefy(t(otu.tab.12), dep)$otu.tab.rff)
		}	
		dep <- colSums(otu.tab.12)[1]
		cat("Depth ", dep, '\n')
		# Remove empty OTUs
		otu.ind <- rowSums(otu.tab.12) > 0  # rev:2016-06-28
		otu.tab.12 <- otu.tab.12[otu.ind, ]	
		otu.name.12 <- otu.name.12[otu.ind, ]
		otu.name.full <- otu.name.full[otu.ind, ]
		act.seq <- paste0(act.seq, 'R')
	} else {
		rff <- FALSE
		dep <- NULL
	}

	# Load mapping file
	if (load.map == TRUE) {
		samIDs <- intersect(rownames(meta.dat), colnames(otu.tab.12))
		if (length(samIDs) == 0) {
			stop('Sample names in the meta file and biom are completely different?\n')
		} 
		if (length(samIDs) < length(colnames(otu.tab.12)) | length(samIDs) < length(rownames(meta.dat))) {
			warning('Sample names in the meta file and biom differ! May be due to rarefaction?\n')
		}
		meta.dat <- meta.dat[samIDs, ]
		otu.tab.12 <- otu.tab.12[, samIDs]
	} else {
		samIDs <- colnames(otu.tab.12)
	}

   # Create abundance list
    cat("Create taxa abundance list ...\n")
	abund.list.12 <- list()
	hierachs <- c('Phylum', 'Class', 'Order', 'Family', 'Genus')
	for (hierach in hierachs) {	
		if (hierach != 'Phylum') {
			single.names <- otu.name.12[, hierach]
		#	single.names[grepl('unclassified', single.names, ignore.case=T)] <- paste0('Unclassified',substr(hierach, 1, 1))
			tax.family <- paste(otu.name.12[, 'Phylum'], single.names, sep=";")
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- paste0('Unclassified_', hierach)
		} else {
			tax.family <- otu.name.12[, 'Phylum']
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- 'Unclassified_Phylum'
		}
		family <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
		rownames(family) <- family[, 1]
		family <- as.matrix(family[, -1])
		abund.list.12[[hierach]] <- family
	}
	
	if (species) {
		abund.list.12[['Species']] <- otu.tab.12
		rownames(abund.list.12[['Species']]) <- paste0("OTU", rownames(otu.tab.12), ":", otu.name.12[, 'Phylum'], ";", otu.name.12[, 'Genus'])
	}
	
    cat('Normalize (size factor) ...\n')
	if (rff == TRUE) {
		cat('For rarefied data, the size factor for samples can still be different!\n')
	}
		
	if (level == 'OTU') {
		data <- otu.tab.12
	} else {
		if (level %in% names(abund.list.12)) {
			data <- abund.list.12[[level]]
		} else {
			data <- otu.tab.12
			level <-'OTU'
			warning('No or wrong level specified! OTU level will be used!\n')
		}
	}
	
	if (norm == 'GMPR') {
		sf <- GMPR(data, intersect.no=intersect.no)
		warning('GMPR is only suitable for samples from the same body location!\n')
		norm <- 'GMPR'
		names(sf) <- colnames(data)
		act.seq <- paste0(act.seq, 'N')
	} else {
		if (norm == 'TSS') {
			sf <- colSums(data)
			norm <- 'TSS'
			act.seq <- paste0(act.seq, 'N')
		} else {
			warning('Normalization method not specified or unknown! TSS is used!\n')
			sf <- colSums(data)
			norm <- 'TSS'
			act.seq <- paste0(act.seq, 'N')
		}
		
	}

	if (winsor == TRUE) {
		act.seq <- paste0(act.seq, 'W')
		if (rff == TRUE) {
			warning('Winsorization after rarefaction will make the data have different total numbers!\n')
		}
		cat('Winsorize ...\n')
		if (is.null(winsor.qt)) {
			winsor.qt <- 0.97
		}
		# Addressing the outlier (97% percent) or at least one outlier
		abund.list.12 <- sapply(abund.list.12, function(genus) {
					genus.p <- t(t(genus) / sf)
					genus.p <- apply(genus.p, 1, function(x) {
								cutoff <- quantile(x, winsor.qt)
								x[x >= cutoff] <- cutoff
								x
							}
					)
					# column/row switch
					genus.w <- t(round(genus.p * sf))
				})
		# OTU table
		otu.tab.12.p <- t(t(otu.tab.12) / sf)
		otu.tab.12.p <- apply(otu.tab.12.p, 1, function(x) {
					cutoff <- quantile(x, winsor.qt)
					x[x >= cutoff] <- cutoff
					x
				}
		)
		# column/row switch
		otu.tab.12 <- t(round(otu.tab.12.p * sf))
	} else {
		winsor <- FALSE
		winsor.qt <- NULL
	}

	# Rarefaction/Normalizing factors are not calculated for functional data
    # Rev: 2017_01_19 the sample IDs for functional data are not ordered! Potential Danger! augment with NA
	if (!is.null(ko.file)) {
		
		cat("Load kegg file...\n")
		ko <- read_biom(ko.file)
		ko.dat <- as.matrix(biom_data(ko))
		
		if (sum(!(samIDs %in% colnames(ko.dat))) != 0) {
			missingIDs <- setdiff(samIDs, colnames(ko.dat))
			aug.mat <- matrix(NA, nrow(ko.dat), length(missingIDs))
			colnames(aug.mat) <- missingIDs
			rownames(aug.mat) <- rownames(ko.dat)
			ko.dat <- cbind(ko.dat, aug.mat) 
		}	
		
		ko.dat <- ko.dat[, samIDs]
		# Rarefaction?
	  
	    if (is.null(ko.ann.file)) {
			# Old - back compatability
			ko.ann <- observation_metadata(ko)
			ko.ann <- cbind(KEGG_Pathways1=sapply(ko.ann, function(x) x['KEGG_Pathways1']), 
					KEGG_Pathways2=sapply(ko.ann, function(x) x['KEGG_Pathways2']), 
					KEGG_Pathways3=sapply(ko.ann, function(x) x['KEGG_Pathways3']))
			rownames(ko.ann) <- rownames(ko.dat)		
			ko.ann[is.na(ko.ann)] <- 'Unclassified'
			
			hierachs <- c("KEGG_Pathways1", "KEGG_Pathways2", "KEGG_Pathways3")
			for (hierach in hierachs) {	
				tax.family <- ko.ann[, hierach]
				family <- aggregate(ko.dat, by=list(tax.family), FUN=sum)
				rownames(family) <- family[, 1]
				family <- as.matrix(family[, -1])
				abund.list.12[[hierach]] <- family
			}
		} else {
			# New
			load(ko.ann.file)
			#
			kos <- rownames(ko.dat)
			abund.list.12[["KEGG_Pathways3"]] <- NULL
			kos.id <- NULL
			for (ko.item in names(kegg.map)) {
				kos.common <- intersect(kos, kegg.map[[ko.item]])
				if (length(kos.common) != 0) {
					abund.list.12[["KEGG_Pathways3"]] <- rbind(abund.list.12[["KEGG_Pathways3"]], colSums(ko.dat[kos.common, , drop=F]))
					kos.id <- c(kos.id, ko.item)
				}
			}
			rownames(abund.list.12[["KEGG_Pathways3"]]) <- kos.id
			
			abund.list.12[["KEGG_Metabolism"]] <- abund.list.12[["KEGG_Pathways3"]][intersect(kos.id, unlist(kegg.ann[['Metabolism']])), ]
			rownames(abund.list.12[["KEGG_Metabolism"]]) <- paste0('M', rownames(abund.list.12[["KEGG_Metabolism"]])) 
			
			abund.list.12[["KEGG_Defense"]] <- NULL
			kos.id <- NULL
			for (ko.item in names(defense.map)) {
				kos.common <- intersect(kos, defense.map[[ko.item]])
				if (length(kos.common) != 0) {
					abund.list.12[["KEGG_Defense"]] <- rbind(abund.list.12[["KEGG_Defense"]], colSums(ko.dat[kos.common, , drop=F]))
					kos.id <- c(kos.id, ko.item)
				}
			}
			rownames(abund.list.12[["KEGG_Defense"]]) <- kos.id
			
			abund.list.12[["KEGG_Toxin"]] <- NULL
			kos.id <- NULL
			for (ko.item in names(toxin.map)) {
				kos.common <- intersect(kos, toxin.map[[ko.item]])
				if (length(kos.common) != 0) {
					abund.list.12[["KEGG_Toxin"]] <- rbind(abund.list.12[["KEGG_Toxin"]], colSums(ko.dat[kos.common, , drop=F]))
					kos.id <- c(kos.id, ko.item)
				}
			}
			rownames(abund.list.12[["KEGG_Toxin"]]) <- kos.id
		}
	
	}
	
	if (!is.null(cog.file)) {
		cat("Load cog file...\n")
		cog <- read_biom(cog.file)
		cog.dat <- as.matrix(biom_data(cog))
		
		if (sum(!(samIDs %in% colnames(cog.dat ))) != 0) {
			missingIDs <- setdiff(samIDs, colnames(cog.dat ))
			aug.mat <- matrix(NA, nrow(cog.dat), length(missingIDs))
			colnames(aug.mat) <- missingIDs
			rownames(aug.mat) <- rownames(cog.dat)
			cog.dat  <- cbind(cog.dat, aug.mat) 
		}	
		
		cog.dat <- cog.dat[, samIDs]
		
		# rarefaction?
		cog.ann <- observation_metadata(cog)
		hierachs <- c("COG_Category1", "COG_Category2")
		for (hierach in hierachs) {	
			tax.family <- sapply(cog.ann, function(x) x[hierach])
			family <- aggregate(cog.dat, by=list(tax.family), FUN=sum)
			rownames(family) <- family[, 1]
			family <- as.matrix(family[, -1])
			abund.list.12[[hierach]] <- family
		}
	}
		
	# Drop tree tips
    if (!is.null(tree.12)) {
		absent <- tree.12$tip.label[!(tree.12$tip.label %in% rownames(otu.tab.12))]
		if (length(absent) != 0) {
			tree.12 <- drop.tip(tree.12, absent)
			warning("The tree has OTUs not in the OTU table!")
		}
	}

	data.obj <- list(otu.tab=otu.tab.12, abund.list=abund.list.12, meta.dat=meta.dat, tree=tree.12,
			otu.name=otu.name.12, otu.name.full=otu.name.full, 
			size.factor=sf, norm.method=norm, norm.level=level, 
			winsor=winsor, winsor.qt=winsor.qt,
			rff=rff, rff.dep=dep, act.seq=act.seq,
			call=match.call())
}

# New: 2016_12_12, separate rarafaction and load data
# Functional data will not be rarefied
rarefy_data <- function (data.obj, dep=NULL) {
	otu.tab.12 <- data.obj$otu.tab
	otu.name.12 <- data.obj$otu.name
	otu.name.full <- data.obj$otu.name.full
	abund.list.12 <- data.obj$abund.list
	siz.factor <- data.obj$size.factor
	meta.dat <- data.obj$meta.dat
	
	cat('Rarefaction ...\n')
	if (is.null(dep)) {
		otu.tab.12 <- t(Rarefy(t(otu.tab.12))$otu.tab.rff)
	} else {
		otu.tab.12 <- t(Rarefy(t(otu.tab.12), dep)$otu.tab.rff)
	}	
	dep <- colSums(otu.tab.12)[1]
	cat("Depth ", dep, '\n')
	# Remove empty OTUs
	otu.ind <- rowSums(otu.tab.12) > 0  # rev:2016-06-28
	otu.tab.12 <- otu.tab.12[otu.ind, ]	
	otu.name.12 <- otu.name.12[otu.ind, ]
	otu.name.full <- otu.name.full[otu.ind, ]
	
	samIDs <- intersect(rownames(meta.dat), colnames(otu.tab.12))

	if (length(samIDs) < nrow(meta.dat)) {
		warning('Some samples were lost during rarefaction!\n')
	}
	
	meta.dat <- meta.dat[samIDs, ]
	otu.tab.12 <- otu.tab.12[, samIDs]
	
	# Create abundance list
	cat("Recreate taxa abundance list ...\n")
	hierachs <- intersect(c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), names(abund.list.12))
	for (hierach in hierachs) {	

		if (hierach != 'Phylum') {
			single.names <- otu.name.12[, hierach]
			#	single.names[grepl('unclassified', single.names, ignore.case=T)] <- paste0('Unclassified',substr(hierach, 1, 1))
			tax.family <- paste(otu.name.12[, 'Phylum'], single.names, sep=";")
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- paste0('Unclassified_', hierach)
		} else {
			tax.family <- otu.name.12[, 'Phylum']
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- 'Unclassified_Phylum'
		}
		family <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
		rownames(family) <- family[, 1]
		family <- as.matrix(family[, -1])
		abund.list.12[[hierach]] <- family
	}
	
	if ('Species' %in% hierachs) {
		abund.list.12[['Species']] <- otu.tab.12
		rownames(abund.list.12[['Species']]) <- paste0("OTU", rownames(otu.tab.12), ":", otu.name.12[, 'Phylum'], ";", otu.name.12[, 'Genus'])
	}
	
	data.obj$meta.dat <- meta.dat  # Rev: 2017_01_19
	data.obj$otu.tab <- otu.tab.12
	data.obj$otu.name <- otu.name.12
	data.obj$otu.name.full <- otu.name.full
	data.obj$abund.list <- abund.list.12
	data.obj$rff <- TRUE
	data.obj$rff.depth <- dep
	cat('After rarefaction, the size factor is automatically set to be the rarefaction depth!\n')
	data.obj$norm.method <- 'TSS'
	data.obj$norm.level <- 'OTU'
	data.obj$size.factor <- colSums(otu.tab.12)
	
	# Futher subset if it contains functional data
    data.obj <- subset_data(data.obj, samIDs)
	data.obj$act.seq <- paste0(data.obj$act.seq, 'R')
#	data.obj$act.seq <- paste0(data.obj$act.seq, 'N')
	return(data.obj)
}

# New: 2017_04_18
remove_otu <- function (data.obj, ind) {
	data.obj$otu.tab <- data.obj$otu.tab[ind, ]
	data.obj$otu.name <- data.obj$otu.name[ind, ]
	data.obj$otu.name.full <- data.obj$otu.name.full[ind, ]
	
# Re-create Abundance 
	otu.name <- data.obj$otu.name
	otu.tab <- data.obj$otu.tab
	abund.list <- list()
	hierachs <- c('Phylum', 'Class', 'Order', 'Family', 'Genus')
	for (hierach in hierachs) {	
		if (hierach != 'Phylum') {
			single.names <- otu.name[, hierach]
			#	single.names[grepl('unclassified', single.names, ignore.case=T)] <- paste0('Unclassified',substr(hierach, 1, 1))
			tax.family <- paste(otu.name[, 'Phylum'], single.names, sep=";")
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- paste0('Unclassified_', hierach)
		} else {
			tax.family <- otu.name[, 'Phylum']
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- 'Unclassified_Phylum'
		}
		family <- aggregate(otu.tab, by=list(tax.family), FUN=sum)
		rownames(family) <- family[, 1]
		family <- as.matrix(family[, -1])
		abund.list[[hierach]] <- family
	}
	
	data.obj$abund.list <- abund.list
	return(data.obj)
}

# New: 2017_10_17
update_name <- function (data.obj, new.name) {
	
	if (length(new.name) != nrow(data.obj$meta.dat)) stop('The number of sample do not agree!\n')
	if (length(new.name) != length(unique(new.name))) stop ('The new names have duplicates!\n')
	
	rownames(data.obj$meta.dat) <- new.name
	data.obj$abund.list <- lapply(data.obj$abund.list, function (x) {colnames(x) <- new.name; x})
	colnames(data.obj$otu.tab) <- new.name
    if(!is.null(data.obj$size.factor))  names(data.obj$size.factor) <- new.name
	
	return(data.obj)
}

# Rev: 2016_09_22
update_data <- function (data.obj, map.file, meta.sep='\t', quote="\"", comment="", ...) {
	cat("Load meta file...\n")
	if (is.character(map.file)) {
		if (grepl("csv$", map.file)) {
			meta.dat <- read.csv(map.file, header=T, check.names=F, row.names=1, comment=comment, quote=quote, ...)
		} else {
			meta.dat <- read.table(map.file, header=T, check.names=F, row.names=1, comment=comment, sep=meta.sep, quote=quote, ...)
		}
	} else {
		meta.dat <- map.file
	}
    
	data.obj$meta.dat <- meta.dat
	
	samIDs <- intersect(rownames(meta.dat), colnames(data.obj$otu.tab))
	if (length(samIDs) == 0)  stop('Sample names in the meta file and biom file differ?\n')
    data.obj <- subset_data(data.obj, samIDs)
	return(data.obj)
}

# New: 2017_02_16  aggregate data
# New: 2017_11_02  add more variable
aggregate_data <- function (data.obj, subject, strata = NULL) {
	cat('Aggregate data by a factor ...\n')
	subject <- (data.obj$meta.dat[,  subject])
	
	if (!is.null(strata)) {
		grp <- (data.obj$meta.dat[,  strata])
	}

	data.obj.new <- list()
	abund.list<- list()
	for (name in names(data.obj$abund.list)) {
		abund <- data.obj$abund.list[[name]]
		if (is.null(strata)) {
			obj <- aggregate(t(abund), by=list(subject), sum)
			abund.list[[name]] <- t(as.matrix(obj[, -1]))
			colnames(abund.list[[name]]) <- obj[, 1]
		} else {
			obj <- aggregate(t(abund), by=list(subject, grp), sum)
			abund.list[[name]] <- t(as.matrix(obj[, -(1:2)]))
			colnames(abund.list[[name]]) <- paste(obj[, 1], obj[, 2], sep="_")
		}
	}
	data.obj.new$abund.list <- abund.list
	
	abund <- data.obj$otu.tab
	if (is.null(strata)) {
		obj <- aggregate(t(abund), by=list(subject), sum)
		otu.tab <- t(as.matrix(obj[, -1]))
		colnames(otu.tab) <- obj[, 1]
	} else {
		obj <- aggregate(t(abund), by=list(subject, grp), sum)
		otu.tab <- t(as.matrix(obj[, -(1:2)]))
		colnames(otu.tab) <- paste(obj[, 1], obj[, 2], sep="_")
	}

	data.obj.new$otu.tab <- otu.tab
	
	data.obj.new$tree <- data.obj$tree
	data.obj.new$otu.name <- data.obj$otu.name
	data.obj.new$otu.name.full <- data.obj$otu.name.full
	
	if (is.null(strata)) {
		unique.ct <- aggregate(data.obj$meta.dat, by=list(subject), function(x) {
					length(unique(x))
				})
		
		unique.ct <- as.matrix(unique.ct[, -1])
		ind <- colSums(unique.ct) == nlevels(factor(subject))
		
		meta.dat <- aggregate(data.obj$meta.dat, by=list(subject), function(x) {
					x[1]
				})
		rownames(meta.dat) <- meta.dat[, 1]
		meta.dat <- meta.dat[, -1]
		meta.dat <- meta.dat[, ind]
	} else {
		unique.ct <- aggregate(data.obj$meta.dat, by=list(subject, grp), function(x) {
					length(unique(x))
				})
		
		unique.ct <- as.matrix(unique.ct[, -(1:2)])
		ind <- colSums(unique.ct) == nrow(unique.ct)
		
		meta.dat <- aggregate(data.obj$meta.dat, by=list(subject, grp), function(x) {
					x[1]
				})
		rownames(meta.dat) <- paste(meta.dat[, 1], meta.dat[, 2], sep="_")
		meta.dat <- meta.dat[, -(1:2)]
		meta.dat <- meta.dat[, ind]
	}

	data.obj.new$meta.dat <- meta.dat
	data.obj.new$winsor <- data.obj$winsor
	data.obj.new$winsor.qt <- data.obj$winsor.qt
	data.obj.new$rff <- data.obj$rff
	data.obj.new$act.seq <- paste0(data.obj$act.seq, 'A')

	cat('Finished!\n')
	return(data.obj.new)
}

# New: 2017_06_01 Remove reads belonging to specific taxa
remove_taxa <- function (data.obj, taxa.level='Genus', taxa.names) {
	# Remove Methanothermococcus
	if (taxa.level == 'OTU') {
		if (is.numeric(taxa.names)) warning('taxa.names should be characters!\n')
		ind <- !(rownames(data.obj$otu.tab) %in% paste(taxa.names))
		if (sum(!ind) == 0) {
			warning('Do not detect the given OTUs! Check the OTU names!\n')
			return(data.obj)
		}
	} else {
		if (!taxa.level %in% colnames(data.obj$otu.name)) {
			stop("Can't find the specified taxa levels!\n")
		}
		ind <- !(data.obj$otu.name[, taxa.level] %in% taxa.names)
		
		if (sum(!ind) == 0) {
			warning('Do not detect the given taxa! Check the taxa level!\n')
			return(data.obj)
		}
	}

	cat(sum(!ind), ' OTUs will be removed!\n')
	
	data.obj$otu.tab <- data.obj$otu.tab[ind, ]
	data.obj$otu.name <- data.obj$otu.name[ind, ]
	data.obj$otu.name.full <- data.obj$otu.name.full[ind, ]
	
# Re-create Abundance 
	otu.name <- data.obj$otu.name
	otu.tab <- data.obj$otu.tab
	abund.list <- list()
	hierachs <- c('Phylum', 'Class', 'Order', 'Family', 'Genus')
	for (hierach in hierachs) {	
		if (hierach != 'Phylum') {
			single.names <- otu.name[, hierach]
			#	single.names[grepl('unclassified', single.names, ignore.case=T)] <- paste0('Unclassified',substr(hierach, 1, 1))
			tax.family <- paste(otu.name[, 'Phylum'], single.names, sep=";")
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- paste0('Unclassified_', hierach)
		} else {
			tax.family <- otu.name[, 'Phylum']
			tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- 'Unclassified_Phylum'
		}
		family <- aggregate(otu.tab, by=list(tax.family), FUN=sum)
		rownames(family) <- family[, 1]
		family <- as.matrix(family[, -1])
		abund.list[[hierach]] <- family
	}
	if ('Species' %in% names(data.obj$abund.list)) {
		abund.list[['Species']] <- data.obj$abund.list[['Species']][ind, ]
	}
	data.obj$abund.list <- abund.list
	# Erase size factor
	data.obj$size.factor <- NULL
	data.obj$norm.method <- NULL
	data.obj$norm.level <- NULL
	if (!is.na.null(data.obj$rff)) {
		if (data.obj$rff == TRUE) {
			warning('The function is not intended to be applied to rarefied object! Rarefaction should be performed again!\n')
			data.obj$rff <- FALSE
			data.obj$rff.dep <- NULL
		}
	}
    data.obj$act.seq <- NULL
	cat('Remove taxa suceeded!\n')
	return(data.obj)
}
# New: 2016_12_12
# intersect.no=4, winsor.qt=0.97
calculate_size <- function (data.obj, level='OTU', norm='GMPR', ...) {
	if (!is.null(data.obj$rff)) {
		if (data.obj$rff)
		warning("The data has been rarefied! Better calculate the size factor for unrarefied data!\n")
	}
	if (!is.null(data.obj$winsor)) {
		if (data.obj$winsor) {
			warning("The data has been winsorized! Better calculate the size factor before winsorization!\n")
		}
	}

	if (level == 'OTU') {
		data <- data.obj$otu.tab
	} else {
		if (level %in% names(data.obj$abund.list)) {
			data <- data.obj$abund.list[[level]]
		} else {
			data <- data.obj$otu.tab
			level <-'OTU'
			warning('No or wrong level specified! OTU level will be used!\n')
		}
	}
	if (norm == 'TSS') {
		sf <- colSums(data)
		norm <- 'TSS'
	}
	if (norm == 'GMPR') {
		sf <- GMPR(data, ...)
#		warning('GMPR is only suitable for samples from a same body location!\n')
		norm <- 'GMPR'
	}
	names(sf) <- colnames(data)
	data.obj$size.factor <- sf
	data.obj$norm.method <- norm
	data.obj$norm.level <- level
	data.obj$act.seq <- paste0(data.obj$act.seq, 'N')
	return(data.obj)
}

# Rev: 2016_12_08
winsor_data <- function (data.obj, winsor.qt=0.97) {
	if (data.obj$rff) {
		warning("The data has been rarefied! Better winsorize without rarefaction!\n")
	}
    cat('Winsorize ...\n')
	if (is.null(winsor.qt)) {
		winsor.qt <- 0.97
	}
	otu.tab.12 <- data.obj$otu.tab
	abund.list.12 <- data.obj$abund.list
	
    # Rev: 2016_12_08
	if (is.null(data.obj$size.factor)) {
		sf <- colSums(otu.tab.12)
		data.obj$norm.method <- 'TSS'
		data.obj$norm.level <- 'OTU'
		data.obj$size.factor <- sf
		cat('Size factor not available! Calculating using total sum! You may consider using better method such as GMPR!\n')
		data.obj$act.seq <- paste0(data.obj$act.seq, 'N')
	} else {
		sf <- data.obj$size.factor
	}
	
	# Addressing the outlier (97% percent) or at least one outlier
	abund.list.12 <- lapply(abund.list.12, function(genus) {
				genus.p <- t(t(genus) / sf)
				genus.p <- apply(genus.p, 1, function(x) {
							cutoff <- quantile(x, winsor.qt)
							x[x >= cutoff] <- cutoff
							x
						}
				)
				# column/row switch
				genus.w <- t(round(genus.p * sf))
			})
	# OTU table
	otu.tab.12.p <- t(t(otu.tab.12) / sf)
	otu.tab.12.p <- apply(otu.tab.12.p, 1, function(x) {
				cutoff <- quantile(x, winsor.qt)
				x[x >= cutoff] <- cutoff
				x
			}
	)
	# column/row switch
	otu.tab.12 <- t(round(otu.tab.12.p * sf))
	data.obj$winsor <- TRUE
	data.obj$winsor.qt <- winsor.qt
    data.obj$otu.tab <- otu.tab.12
	data.obj$abund.list <- abund.list.12
	data.obj$act.seq <- paste0(data.obj$act.seq, 'W')
	#data.obj$size.factor <- sf
	return(data.obj)
}

# Rev: 2016_11_28, Bray-curtis use normalized data.
construct_distance <- function (data.obj, unifrac.file=NULL,  Phylum='All', dist.RData=NULL, save.RData=NULL, 
		filter.no=0, rff=FALSE, dep=NULL, seed=1234) {
	set.seed(seed)
	
	if (!is.null(dist.RData)) {
		load(dist.RData, envir=.GlobalEnv )
	} else {
		dist.list.12 <- list()
		cat("Generalized UniFrac ...\n")
		
		otu.tab <- t(data.obj$otu.tab)
		
		if (rff == TRUE) {
			if (is.null(dep)) {
				otu.tab <- Rarefy(otu.tab)$otu.tab.rff
			} else {
				otu.tab <- Rarefy(otu.tab, dep)$otu.tab.rff
			}	
		}
				
		if (Phylum != 'All') {
			ind <- data.obj$otu.name[, 'Phylum'] == Phylum
			otu.tab <- otu.tab[, ind]
		}
		
		# Filter otus with reads <= filter.no
		otu.tab <- otu.tab[, colSums(otu.tab) > filter.no]
        
		# Remove samples with no reads
		if (sum(rowSums(otu.tab) == 0) >= 1) {
			otu.tab <- otu.tab[rowSums(otu.tab) != 0, ]
			warning('Some samples do not have reads after rarefaction! Please be careful!\n')
		}
		
		# To make sure the OTUs in otu.tab are in the tree (rev:2016-06-19)
		
		if (sum(!(colnames(otu.tab) %in% data.obj$tree$tip.label))) {
			warning('Some OTU names are not in the tree! An intersection set will be used!\n')	
		}
		common.otus <- intersect(colnames(otu.tab), data.obj$tree$tip.label)
		unifrac12 <- GUniFrac(otu.tab[, common.otus], data.obj$tree)$unifracs
	
		dist.list.12[['WUniFrac']] <- unifrac12[, , 'd_1']
		dist.list.12[['GUniFrac']] <- unifrac12[, , 'd_0.5']
		if (is.null(unifrac.file)) {
			dist.list.12[['UniFrac']] <- unifrac12[, , 'd_UW']
		} else {
			# The orders may be different
			dist.list.12[['UniFrac']] <- as.matrix(read.table(unifrac.file, row.names=1, header=T)) # Rarefaction
		}

		# Need speed up
		# Suggest using rarefied counts 
	    # If case/control has different sequencing depth, it will result in false clustering! 
        # Rev: 2016_11_28 Use normalized data for BC distance calculation to reduce noise in case of variable library size
		dist.list.12[['BC']] <-as.matrix(vegdist(otu.tab / rowSums(otu.tab)))
		
		genus <- t(data.obj$abund.list[['Genus']])
		genus <- genus / rowSums(genus)
		dist.list.12[['Euc']] <- as.matrix(dist(genus))
		
		genus <- sqrt(genus)
		dist.list.12[['Hel']] <-as.matrix(dist(genus))
#		dist.list.12[['JS']] <- as.matrix(distance(otu_table(data.obj$abund.list[['Genus']], taxa_are_rows=T), method='jsd'))
		
		if (!is.null(save.RData)) {
			save(dist.list.12, file=save.RData)
		}
	}

    return(dist.list.12)
}


outlier_detect <- function (data.obj, dist.obj, min.dep=2000) {
	# Future development
	samIDs <- colnames(data.obj$otu.tab)[colSums(data.obj$otu.tab) >= 2000]
	return(samIDs)
}

is.na.null <- function (x) {
	if (is.null(x)) {
		return(TRUE)
	} else {
		if (is.na(x)[1]) {
			return(TRUE)
		}  else {
			return(FALSE)
		}
	}
	
}

# Rev: 2016_09_26 remove empty OTUs/taxa
# Rev: 2016_12_01 add more logical controls
subset_data <- function (data.obj, samIDs) {
	
	# Rev: 2016_1_19 to add error protection
	# Transform logical samIDs into characer samIDs
	if (is.logical(samIDs) | is.numeric(samIDs)) {
		samIDs <- rownames(data.obj$meta.dat)[samIDs]
	}
	
	data.obj$meta.dat <- data.obj$meta.dat[samIDs, , drop=FALSE]

	if (!is.na.null(data.obj$otu.tab)) {

		data.obj$otu.tab <- data.obj$otu.tab[, samIDs, drop=FALSE]
		data.obj$otu.tab <- data.obj$otu.tab[rowSums(data.obj$otu.tab) != 0, , drop=FALSE]
		data.obj$otu.name <- data.obj$otu.name[rownames(data.obj$otu.tab), , drop=FALSE]

		if (!is.na.null(data.obj$otu.name.full)) {
			data.obj$otu.name.full <- data.obj$otu.name.full[rownames(data.obj$otu.tab), , drop=FALSE]
		}
	}

	if (!is.na.null(data.obj$abund.list)) {
		data.obj$abund.list <- lapply(data.obj$abund.list, function(x) {
					xx <- x[, samIDs, drop=FALSE]
					xx <- xx[rowSums(xx) != 0, , drop=FALSE]
				})
	}

	if (!is.na.null(data.obj$size.factor)) {
		data.obj$size.factor <- data.obj$size.factor[samIDs]
	}

	if (!is.na.null(data.obj$ko.list)) {
		data.obj$ko.list <- lapply(data.obj$ko.list, function(x) {
					xx <- x[, samIDs, drop=FALSE]
					xx <- xx[rowSums(xx) != 0, , drop=FALSE]
				})
	}

	if (!is.na.null(data.obj$cog.list)) {
		data.obj$cog.list <- lapply(data.obj$cog.list, function(x) {
					xx <- x[, samIDs, drop=FALSE]
					xx <- xx[rowSums(xx) != 0, , drop=FALSE]
				})
	}
	return(data.obj)
}

# Rev: 2016_12_01 add more logical controls
# Rev: 2016_02_01 fix one error
subset_dist <- function (dist.obj, samIDs) {
	
	# Rev: 2016_1_19 to add error protection
	# Transform logical samIDs into character samIDs
	if (is.logical(samIDs) | is.numeric(samIDs)) {
		samIDs <- rownames(dist.obj[[1]])[samIDs]
	}
	
	lapply(dist.obj, function(x) {
				if(!is.na.null(x)){
					x <- x[samIDs, samIDs]
				} else {
					x
				}
				x
			})
}


perform_sequence_stat_analysis <- function (data.obj, ann='') {
	sink(paste0('Sequence_Analysis_Statistics_', ann, '.txt'))
	otu.tab <- data.obj$otu.tab
	
	# Sequencing depth
	otu.abund <- rowSums(otu.tab)
	sam.abund <- colSums(otu.tab)
	otu.prev <- rowSums(otu.tab!=0)/ncol(otu.tab)
	
	otu.abund <- otu.abund[otu.abund >= 1]
	sam.abund <- sam.abund[sam.abund >= 1]
	cat('This data set contains ', length(sam.abund), ' samples after quality controls.\n')
	cat('16S rDNA targeted sequencing yields ', median(sam.abund), 'reads/sample on average (range:', min(sam.abund), '-', max(sam.abund), ').\n')
	cat('Clustering of these 16S sequence tags produces ', sum(otu.abund > 0), ' OTUs at 97% similarity level.\n')
	
	pdf(paste0('Sequence_Analysis_Statistics_', ann, '.pdf'), height=5, width=5)
	obj <- ggplot2::ggplot(data=data.frame(x=otu.abund), aes(x=x)) + geom_histogram(col='black', fill='gray') + ylab('Frequency') + xlab('Abundance(Total counts)') +
			scale_x_log10(breaks=c(1, 10, 100, 1000, 10000, 100000, 100000))
	print(obj)
	obj <- ggplot2::ggplot(data=data.frame(x=sam.abund), aes(x=x)) + geom_histogram(col='black', fill='gray')  + ylab('Frequency') + xlab('Sequencing depth')
	print(obj)
	obj <- ggplot2::ggplot(data=data.frame(x=otu.prev), aes(x=x))  + ylab('Frequency') + xlab('Prevalence(Occurence frequency)') + geom_histogram(col='black', fill='gray')
	print(obj)
	dev.off()
	
	phy.abund <- data.obj$abund.list[['Phylum']]
	fam.abund <- data.obj$abund.list[['Family']]
	gen.abund <- data.obj$abund.list[['Genus']]
	
	phy.prev <- rowSums(phy.abund != 0) / ncol(phy.abund)
	fam.prev <- rowSums(fam.abund != 0) / ncol(phy.abund)
	gen.prev <- rowSums(gen.abund != 0) / ncol(phy.abund)
	
	phy.abund <- rowMeans(t(t(phy.abund) / sam.abund))
	fam.abund <- rowMeans(t(t(fam.abund) / sam.abund))
	gen.abund <- rowMeans(t(t(gen.abund) / sam.abund))
	
	cat('These OTUs belong to ', sum(phy.abund > 0), ' phyla,', sum(fam.abund > 0), ' families and ', sum(gen.abund > 0), 'genera.\n\n')
	
	phy.prev <- sort(phy.prev, decr=T)
	phy.prev <- round(phy.prev[phy.prev >= 0.90] * 100, 2)
	
	fam.prev <- sort(fam.prev, decr=T)
	fam.prev <- round(fam.prev[fam.prev >= 0.90] * 100, 2)
	
	gen.prev <- sort(gen.prev, decr=T)
	gen.prev <- round(gen.prev[gen.prev >= 0.90] * 100, 2)
	
	# Rev: 2017_02_19 ' ' -> '\n'
	cat('\nThe most prevalent phyla are:\n', paste(paste0(names(phy.prev), '(', phy.prev, '%)'), collapse='\n'), '\n')
	cat('\nThe most prevalent families are:\n', paste(paste0(names(fam.prev), '(', fam.prev, '%)'), collapse='\n'), '\n')
	cat('\nand the most prevalent genera are:\n', paste(paste0(names(gen.prev), '(', gen.prev, '%)'), collapse='\n'), '\n\n')
	
	phy.abund <- sort(phy.abund, decr=T)
	phy.abund <- round(phy.abund[phy.abund >= 0.05] * 100, 2)
	
	fam.abund <- sort(fam.abund, decr=T)
	fam.abund <- round(fam.abund[fam.abund >= 0.05] * 100, 2)
	
	gen.abund <- sort(gen.abund, decr=T)
	gen.abund <- round(gen.abund[gen.abund >= 0.05] * 100, 2)
	
	cat('\nThe most abundant phyla are:\n', paste(paste0(names(phy.abund), '(', phy.abund, '%)'), collapse='\n'), '\n')
	cat('\nThe most abundant families are:\n', paste(paste0(names(fam.abund), '(', fam.abund, '%)'), collapse='\n'), '\n')
	cat('\nand the most abundant genera are:\n', paste(paste0(names(gen.abund), '(', gen.abund, '%)'), collapse='\n'), '\n\n')
	sink()
}

perform_demograph_analysis <- function (data.obj, grp.name) {
	obj <- summary(data.obj$meta.dat)
	write.csv(obj, "meta.data.summary.csv", quote=F)

	grp <- data.obj$meta.dat[, grp.name]
	res <- NULL
	pv.vec <- NULL
	if (is.factor(grp)) {
		for (var1 in setdiff(colnames(data.obj$meta.dat), grp.name)) {
			temp <- data.obj$meta.dat[, var1]
			res <- rbind(res, c("", var1, rep("", nlevels(grp) - 1)))
			res <- rbind(res, c("", levels(grp)))
			if (is.factor(temp)) {
				res <- rbind(res, cbind(levels(temp), table(temp, grp)))
				if (nlevels(temp) == 1) {
					pv <- NA
				} else {
					err <- try(
							pv <- formatC(fisher.test(table(temp, grp))$p.value, digit=3)
					)
					if (inherits(err, "try-error")) {
						pv <- NA
					}
				}
				res <- rbind(res, c('Fisher p', pv, rep("", nlevels(grp) - 1)))

			} else {
				res <- rbind(res, c('mean', aggregate(temp, by=list(grp), FUN='mean')[, 2]))
				res <- rbind(res, c('sd', aggregate(temp, by=list(grp), FUN='sd')[, 2]))
				
				err <- try(
						pv <- formatC(summary(aov(temp ~ grp))[[1]][1, 'Pr(>F)'], digit=3)
				)
				if (!inherits(err, "try-error")) {
					res <- rbind(res, c('ANOVA p', pv, rep("", nlevels(grp) - 1)))
				} else {
					res <- rbind(res, c('ANOVA p', 'NA', rep("", nlevels(grp) - 1)))
					pv <- NA
				}

			}
			res <- rbind(res, rep("", nlevels(grp)+1))
			pv.vec <- c(pv.vec, pv)
		}
		
	} else {

		
	}

	write.csv(res, "meta.data.by.grp.csv", row.names=F, quote=F)
	names(pv.vec) <- setdiff(colnames(data.obj$meta.dat), grp.name)
	return(pv.vec)
}


# Rev: 2017_08_23 add automatically create phylo.obj
generate_rarefy_curve <- function (data.obj, phylo.obj=NULL, grp.name, depth=NULL, npoint=10, iter.no=5,
		measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'), ann='', gg.cmd="theme(legend.justification=c(1,0), legend.position=c(1,0))", wid=5, hei=5) {
	cat("Create rarefaction curves!\n")
	# Rev: 2017_08_23
	if (is.null(phylo.obj)) {
		phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows=T), phy_tree(data.obj$tree), 
				tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))
	}
	if (is.null(depth)) {
		depth <- min(sample_sums(phylo.obj))
		phylo.even <- rarefy_even_depth(phylo.obj, rngseed=12345)
	} else {
		if (depth > min(sample_sums(phylo.obj))) {
			ind <- sample_sums(phylo.obj) >= depth
			cat(sum(!ind), " samples do not have sufficient number of reads!\n")
			sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ] 
			data.obj <- subset_data(data.obj, ind)
		}
		phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed=12345)
	}
	
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	# Rev: 2016_12_12
	if (is.character(grp)) {
		grp <- factor(grp)
	}

	if (!is.factor(grp)) {
		stop('Rarefaction curve needs a factor!\n')
	}
	res <- NULL
	incr <- depth %/% npoint
	sink('temp.txt')
	for (dep in c(10, incr*(1:npoint))) {
		x <- 0
		for (i in 1:iter.no) {
			phylo.even <- rarefy_even_depth(phylo.obj, dep, rngseed=12345+i)
			x <- x + estimate_richness(phylo.even, measures=measures)
		}
		
		res <- rbind(res, t(x[, measures, drop=F]/iter.no))
	}
	colnames(res) <- rownames(df)
	sink()
	
	pdf(paste0("Alpha_diversity_Rarefaction_Curve_", ann, ".pdf"), width=wid, height=hei)
	for (i in 1:length(measures)) {
		measure <- measures[i]
		cat("Measure: ", measure, "\n")
		res2 <- res[(0:(npoint))*length(measures)+i, , drop=F]
		m <- t(apply(res2, 1, function(x) tapply(x, grp, mean)))
		se <- t(apply(res2, 1, function(x) tapply(x, grp, function(y) sd(y)/sqrt(length(y)))))
		uci <- m+se
		lci <- m-se
		
		m <- melt(m)
		uci <- melt(uci)
		lci <- melt(lci)
		
		res2 <- cbind(c(10, incr*(1:npoint)), m[, 2:3], uci[, 3], lci[, 3])
		colnames(res2) <- c('Depth', 'Group', 'mean', 'max', 'min')
		
		res2 <- as.data.frame(res2)
		res2$Group <- factor(res2$Group, levels=levels(grp))
		
		#write.table(res2, paste0("Alpha_diversity_Rarefaction_", ann, "_", measure, ".txt"))
		
		obj <- ggplot(res2, aes(x=Depth, y=mean, color=Group, group=Group)) +
				geom_errorbar(aes(ymin=min, ymax=max), alpha=0.5, width=.25, position=position_dodge(.2)) + 
				geom_line() + 
				geom_point(size=3, shape=21, fill="white") +
				labs(y=measure)

		if (!is.null(gg.cmd)) {
			obj <- obj + eval(parse(text=gg.cmd))
		}
		
		print(obj)		
	}
	dev.off()
}

# Rev: 2016_12_25  Add anova, and record the results
perform_alpha_test <- function (data.obj, phylo.obj=NULL, rarefy=TRUE, depth=NULL, iter.no=5, 
		measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'),  model='lm', 
		formula=NULL, grp.name=NULL, adj.name=NULL, ann='', ...) {
	# Implement future random effects model
	# Rev: 2017_08_23
	if (is.null(phylo.obj)) {
		phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows=T), phy_tree(data.obj$tree), 
				tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))
	}
	
	result <- list()
	if (rarefy == TRUE) {
		if (is.null(depth)) {
			depth <- min(sample_sums(phylo.obj))
		} else {
			if (depth > min(sample_sums(phylo.obj))) {
				ind <- sample_sums(phylo.obj) >= depth
				cat(sum(!ind), " samples do not have sufficient number of reads!\n")
				sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ] 
				data.obj <- subset_data(data.obj, ind)
			}
		}
			
		x <- 0
		sink('temp.txt')
		for (i in 1:iter.no) {
			phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed=12345+i)
			x <- x + estimate_richness(phylo.even, measures=measures)
		}
		sink()
		x <- x / iter.no
	} else {
		x <- estimate_richness(phylo.obj, measures=measures)
	}

	result$alpha.diversity <- x
	
	df <- data.obj$meta.dat
	
	fitted.obj <- list()
	
	if (rarefy == T) {
		sink(paste0('Alpha_diversity_test_results_rarefied_', ann, '.txt'))
	} else {
		sink(paste0('Alpha_diversity_test_results_unrarefied_', ann, '.txt'))
	}
	
	date()
	
	# variables to adjust always come first in anova analyses
	if (is.null(formula)) {
		if (is.null(adj.name)) {
			formula <- paste('~', grp.name)
		} else {
			formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
		}
	}
	
	for (measure in measures) {
		cat("Alpha diversity:", measure, "\n")
		xx <- x[, measure]
		if (model == 'lm') {
			cat('Linear model:\n')

			lm.obj <- lm(as.formula(paste('xx ', formula)), df, ...)
			prmatrix(summary(lm.obj)$coefficients)
			cat('\nANOVA:\n')
			print(anova(lm.obj))
			cat('\n')
			fitted.obj[[measure]] <- lm.obj
		}
		if (model == 'lme') {
			df$xx <- xx
			cat('Linear mixed effects model:\n')
			lm.obj <- lme(as.formula(paste('xx ', formula)), df, method='ML', ...)
			prmatrix(summary(lm.obj)$tTable)
			cat('\nANOVA:\n')
			print(anova(lm.obj))
			cat('\n')
			fitted.obj[[measure]] <- lm.obj
		}
		cat("\n")
	}
	sink()
	
	result$fitted.obj <- fitted.obj
	# Rev: 2016_12_25
	return(invisible(result))

}

# Rev: 2016_09_10
# Rev: 2016_11_28
# Rev: 2017_04_18
generate_alpha_boxplot <- function (data.obj, phylo.obj=NULL, rarefy=TRUE, depth=NULL, grp.name, strata=NULL, 
		measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'), gg.cmd=NULL, ann='', subject=NULL, p.size=2.5, l.size=0.5,
		hei = NULL, wid = NULL) {	
	# Rev: 2017_08_23
	if (is.null(phylo.obj)) {
		phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows=T), phy_tree(data.obj$tree), 
				tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))
	}
	# To be completed - jetter when strata is not null
	if (rarefy == TRUE) {
		if (is.null(depth)) {
			depth <- min(sample_sums(phylo.obj))
		} else {
			if (depth > min(sample_sums(phylo.obj))) {
				ind <- (sample_sums(phylo.obj) >= depth)
				cat(sum(!ind), " samples do not have sufficient number of reads!\n")
				
				sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ] 
				data.obj <- subset_data(data.obj, ind)
			}
		}

		phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed=12345)
		x <- estimate_richness(phylo.even, measures=measures)
	} else {
		x <- estimate_richness(phylo.obj, measures=measures)
	}
	
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	if (!is.null(subject)) {
		ID <- df[, subject]
	}
	
	if (is.null(hei) | is.null(wid))  {
		hei <- 5
		if (is.null(strata)) {
			wid <- 5
		} else {
			wid <- 5.5
		}
	}

	if (rarefy == T) {
		pdf(paste0('Alpha_diversity_boxplot_rarefied_', ann, '.pdf'), height=hei, width=wid)
	} else {
		pdf(paste0('Alpha_diversity_boxplot_unrarefied_', ann, '.pdf'), height=hei, width=wid)
	}
	if (is.null(subject)) {
		if (is.null(strata)) {
			for (measure in measures) {
				cat(measure, '\n')
				xx <- x[, measure]		
				df2 <- data.frame(Value=xx, Group=grp)
				dodge <- position_dodge(width=0.75)
				obj <- ggplot(df2, aes(x=Group, y=Value, col=Group)) +
						geom_boxplot(position=dodge,  outlier.colour = NA) + 
						geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) +
						labs(y=measure) +
						theme(legend.position="none")
				if (!is.null(gg.cmd)) {
					obj <- obj + eval(parse(text=gg.cmd))
				}
				
				print(obj)
			}	
		} else {
			for (measure in measures) {
				cat(measure, '\n')
				xx <- x[, measure]		
				grp2 <- df[, strata]
				df2 <- data.frame(Value=xx, Group=grp, Strata=grp2)
				
				dodge <- position_dodge(width=0.9)
				obj <- ggplot(df2, aes(x=Strata, y=Value, col=Group, fill=Group)) +
						geom_jitter(alpha=0.6, size=3.0,  position = position_jitterdodge(dodge.width=0.9)) +
						geom_boxplot(position=dodge,  outlier.colour = NA, fill='white', alpha=0.6) + 
						labs(y=measure, x=strata) 
				if (!is.null(gg.cmd)) {
					obj <- obj + eval(parse(text=gg.cmd))
				}
				print(obj)
			}
			
		}
		dev.off()
	} else {
		if (is.null(strata)) {
			for (measure in measures) {
				cat(measure, '\n')
				xx <- x[, measure]		
				df2 <- data.frame(Value=xx, Group=grp, subject=ID)
				dodge <- position_dodge(width=0.75)
			
				obj <- ggplot(df2, aes(x=Group, y=Value, shape=Group,  group=subject)) +
						geom_point(size=p.size) +
						geom_line(size=l.size) +
						labs(y=measure) +
						theme(legend.position="none")
				if (!is.null(gg.cmd)) {
					obj <- obj + eval(parse(text=gg.cmd))
				}
				
				print(obj)
			}	
		} else {
			warning('Stratum is disallowed when subject is specified!\n')
		}
		dev.off()
	}

}

# New: 2017_02_22
generate_alpha_scatterplot <- function (data.obj, phylo.obj=NULL, rarefy=TRUE, depth=NULL, grp.name, 
		strata=NULL, strata2=NULL, wrap.nrow=2, smooth.method='auto', hei0=NULL, wid0=NULL,
		measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'), gg.cmd=NULL, ann='') {	
	# To be completed - jetter when strata is not null
	
	# Rev: 2017_08_23
	if (is.null(phylo.obj)) {
		phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows=T), phy_tree(data.obj$tree), 
				tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))
	}
	
	if (rarefy == TRUE) {
		if (is.null(depth)) {
			depth <- min(sample_sums(phylo.obj))
		} else {
			if (depth > min(sample_sums(phylo.obj))) {
				ind <- (sample_sums(phylo.obj) >= depth)
				cat(sum(!ind), " samples do not have sufficient number of reads!\n")
				
				sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ] 
				data.obj <- subset_data(data.obj, ind)
			}
		}
		
		phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed=12345)
		x <- estimate_richness(phylo.even, measures=measures)
	} else {
		x <- estimate_richness(phylo.obj, measures=measures)
	}
	
	df <- data.obj$meta.dat
	grp <- df[, grp.name]

	if (is.null(hei0)) {
		hei <- 5
	} else {
		hei <- hei0
	}
	
	if (is.null(strata)) {
		if (is.null(wid0)) {
			wid <- 5
		} else {
			wid <- wid0
		}
		
	} else {
		if (is.null(wid0)) {
			wid <- 5.5
		} else {
			wid <- wid0
		}
	}
	
	if (rarefy == T) {
		pdf(paste0('Alpha_diversity_scatterplot_rarefied_', ann, '.pdf'), height=hei, width=wid)
	} else {
		pdf(paste0('Alpha_diversity_scatterplot_unrarefied_', ann, '.pdf'), height=hei, width=wid)
	}
	
	if (is.null(strata)) {
		for (measure in measures) {
			cat(measure, '\n')
			xx <- x[, measure]		
			df2 <- data.frame(Value=xx, Group=grp)
			obj <- ggplot(df2, aes(x=Group, y=Value)) +
					geom_point() +
					geom_smooth(method=smooth.method) +
					labs(y='Value', x=grp.name, title=measure) +
					theme(legend.position="none")
			if (!is.null(gg.cmd)) {
				obj <- obj + eval(parse(text=gg.cmd))
			}
			print(obj)
		}	
	} else {
		for (measure in measures) {
			cat(measure, '\n')
			
			xx <- x[, measure]		
			grp2 <- df[, strata]
			if (!is.null(strata2)) {
				grp3 <- df[, strata2]
			} else {
				grp3 <- grp2
			}
			
			df2 <- data.frame(Value=xx, Group=grp, Strata=grp2, Strata2=grp3)
			if (is.null(strata2)) {
				obj <- ggplot(df2, aes(x=Group, y=Value)) +
						geom_point() +
						geom_smooth(method=smooth.method) +
						labs(y='Value', x=grp.name, title=measure) + 
						facet_wrap(~ Strata, nrow=wrap.nrow)
			} else {
				obj <- ggplot(df2, aes(x=Group, y=Value)) +
						geom_point() +
						geom_smooth(method=smooth.method) +
						labs(y='Value', x=grp.name, title=measure) + 
						facet_grid(Strata2 ~ Strata)
			}		
			obj <- obj + theme(legend.title=element_blank())			
			if (!is.null(gg.cmd)) {
				obj <- obj + eval(parse(text=gg.cmd))
			}
			print(obj)
		}
		
	}
	dev.off()
	
}

perform_alpha_test_otu <- function (data.obj, formula, grp.name, ...) {
	# Assume two groups
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	dep <- colSums(data.obj$otu.tab)
	obs <- colSums(data.obj$otu.tab !=0)
	obj <- lm(as.formula(paste('log(obs)', formula, '+ log(dep)')), df, ...)
	sink('Alpha_diversity_test_results_observed_OTU.txt')
	date()
	cat('Linear model correcting for differential sequencing depth:\n')
	prmatrix(summary(obj)$coefficients)
	sink()
	
	pdf("Alpha_diversity_test_results_observed_OTU.pdf")
	plot(dep, obs, log='xy', xlab="Sequencing Depth", ylab="Observed OTU Number", 
			bg=c('blue', 'red')[grp], pch=21)
	grp.levels <- levels(grp)
	abline(lm(obs[grp==grp.levels[1]] ~ dep[grp==grp.levels[1]]), col='red')
	abline(lm(obs[grp==grp.levels[2]] ~ dep[grp==grp.levels[2]]), col='blue')
	legend('left', levels(grp), fill=c('blue', 'red'))
	dev.off()
}


# Rev: 2017_02_13 Individual label
# Rev: 2017_05_06 Add p value 
generate_ordination <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		grp.name, adj.name=NULL, emp.lev=NULL, indiv.lab=FALSE, indiv.lab.cex=0.5, is.lightcolor=TRUE,
		strata=NULL,  pc.separate, pca.method='cmd', ann=NULL, sub=NULL,
		clab=1.0, cex.pt=1.25, ellipse=T, cstar= 1, wid=5, hei=5, pdf=TRUE, is.pvalue=FALSE, ...) {
	# Implment strata
	# To be completed, add continuous case
	strata0 <- strata
	
	if (pdf) {
		if (is.null(ann)) {
			pdf(paste0('Beta_diversity_ordination_', pca.method, '_', grp.name, '.pdf'), width=wid, height=hei)
		} else {
			pdf(paste0('Beta_diversity_ordination_', pca.method, '_', ann, '.pdf'), width=wid, height=hei)
		}
	}

	df <- data.obj$meta.dat
	grp <- factor(df[, grp.name])
	
	if (!is.null(emp.lev)) {
		grp <- factor(grp, levels=c(setdiff(levels(grp), emp.lev), emp.lev))
	}
	
	if (!is.null(strata)) {
		strata <- factor(df[, strata])
	} else {
		strata <- factor(grp)
	}
	
	# To be revised
	darkcols <- hue_pal(l=40)(nlevels(grp))
	if (is.lightcolor) {
		lightcols <- hue_pal(c=45, l=80)(nlevels(grp))
	} else {
		lightcols <- darkcols
	}
	
	pchs <- rep(c(21, 22, 23, 24, 25), ceiling(nlevels(strata) / 5))[1:nlevels(strata)]
	
	for (dist.name in dist.names) {
		dist.temp <- dist.obj[[dist.name]]
		if (!is.null(adj.name)) {
			adj <- as.data.frame(df[, adj.name])
			obj <- cmdscale(as.dist(dist.temp), k=ncol(dist.temp)-1)
			dat2 <- apply(obj, 2, function(x) resid(lm(x ~ ., data=adj)))
			dist.temp <- dist(dat2)
		} 
		if (pca.method == 'cmd') {
			obj <- cmdscale(as.dist(dist.temp), k=2, eig=T)
			pve <- round(obj$eig[1:2]/sum(abs(obj$eig))*100, 1)
			y <- cbind(obj$points[, 1], obj$points[, 2])
			
			xlab <- paste0('PC1(', pve[1], '%)')
			ylab <- paste0('PC2(', pve[2], '%)')
		} 

		if (pca.method == 'nmds') {
				obj <- metaMDS(as.dist(dist.temp), k=2)
				y <- cbind(obj$points[, 1], obj$points[, 2])
				xlab <- 'NMDS1'
                ylab <- 'NMDS2'
		} 
		
		if (pca.method == 'pls') {
			require(mixOmics)
			# Test
			obj <- cmdscale(as.dist(dist.temp), k=ncol(dist.temp)-1)
			obj <- plsda(obj, grp, ncomp=2)
			y <- cbind(obj$variates$X[, 1], obj$variates$X[, 2])
			xlab <- 'PLS1'
			ylab <- 'PLS2'
		}
		
		if (is.pvalue == TRUE & nlevels(factor(grp)) > 1) {
			pvalue <- paste0('(P=', adonis(as.dist(dist.temp) ~ grp)$aov.tab[1, 6], ')')
		} else {
			pvalue <- ''
		}
		
		colnames(y) <- c("PC1", "PC2")		
		plot(y[, 1], y[, 2], type='n', xlim=range(y[, 1]) * 1.2, ylim=range(y[, 2])*1.2, 
				xlab=xlab, ylab=ylab)
#		points(y[, 1], y[, 2], bg='yellow', col=darkcols[grp], 
#				type='p', pch = pchs[grp], cex=cex.pt)
		col1 <- lightcols[1:nlevels(grp)]
		col2 <- darkcols[1:nlevels(grp)]
		col3 <- darkcols[grp]
		if (!is.null(emp.lev)) {
			
			col1 <- col2
			col1[which(levels(grp) != emp.lev)] <- rgb(0.5, 0.5, 0.5, 0.5)
			col1[which(levels(grp) == emp.lev)] <- rgb(1.0, 0, 0, 1.0)
			col2[which(levels(grp) != emp.lev)] <- rgb(0.5, 0.5, 0.5, 0.5)
			col2[which(levels(grp) == emp.lev)] <- rgb(1.0, 0, 0, 1.0)
            col3[grp !=  emp.lev] <- rgb(0.5, 0.5, 0.5, 0.5)
			col3[grp ==  emp.lev] <- rgb(1.0, 0, 0, 1.0)

		}
		if (ellipse == T) {
			s.class(y, 
					fac = grp,
					cstar = cstar,
					clab = 0,
					cpoint = 0,
					axesell = F,
					col = col1,
					grid = TRUE,
					add.plot=T,
					...
			)
		}
		
		points(y[, 1], y[, 2], bg=col3, col='black',
				type='p', pch = pchs[strata], cex=cex.pt)
		if (indiv.lab) {
			if (is.null(emp.lev)) {
				lab.temp <- rownames(dist.temp)
				text(y[, 1], y[, 2], labels=lab.temp, cex=indiv.lab.cex, col=rgb(1, 0, 0, 0.75))
			} else {
				lab.temp <- rownames(dist.temp)
				lab.temp[grp !=  emp.lev] <- ''
				text(y[, 1], y[, 2], labels=lab.temp, cex=indiv.lab.cex, col=rgb(1, 0, 0, 0.75))
			}
	
		}
		s.class(y, 
				fac = grp,
				cstar =0,
				cellipse = 0,
				clab = clab,
				cpoint = 0,
				axesell = F,
				col = col2,
				grid = TRUE,
				add.plot=T
		)
		if (!is.null(strata0)) {
			legend('topright', legend=(levels(strata)), pch=pchs[1:nlevels(strata)])	
		}

		title(main=paste(dist.name, "distance"), sub=paste0(sub, pvalue))
#		text(-0.25, 0.3, "PERMANOVA p=0.016")
	}
	
	if (pdf) {
		dev.off()
	}
}

# New: 2017_02_21  Separate PC plots for strata2
# Rev: 2017_05_04 Add p value
# Rev: 2017_06_01 Different pages for different distances
generate_ordination_separate <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		grp.name, adj.name=NULL, emp.lev=NULL, indiv.lab=FALSE, indiv.lab.cex=0.5, is.lightcolor=TRUE,
		strata=NULL,  separate=NULL,  layout=NULL, pc.separate, pca.method='cmd', ann=NULL, sub=NULL,
		clab=1.0, cex.pt=1.25, ellipse=T, cstar= 1, wid=5, hei=5, pdf=TRUE, is.pvalue = FALSE, ...) {
	# Implment strata
	# To be completed, add continuous case
	strata0 <- strata
	if (is.null(separate)){
		stop('The augment separate is required!\n')
	} else {
		data.obj$meta.dat[, separate] <- factor(data.obj$meta.dat[, separate])
		sep.levels <- levels(data.obj$meta.dat[, separate])
	}
	if (pdf) {

		if (is.null(ann)) {
			pdf(paste0('Beta_diversity_ordination_', pca.method, '_', grp.name, '_separate.pdf'), width=wid, height=hei)
		} else {
			pdf(paste0('Beta_diversity_ordination_', pca.method, '_', ann, '_separate.pdf'), width=wid, height=hei)
		}
	}
	

	df <- data.obj$meta.dat
	grp <- factor(df[, grp.name])
	
	if (!is.null(emp.lev)) {
		grp <- factor(grp, levels=c(setdiff(levels(grp), emp.lev), emp.lev))
	}
	
	if (!is.null(strata)) {
		strata <- factor(df[, strata])
	} else {
		strata <- factor(grp)
	}
	
	# To be revised
	darkcols <- hue_pal(l=40)(nlevels(grp))
	if (is.lightcolor) {
		lightcols <- hue_pal(c=45, l=80)(nlevels(grp))
	} else {
		lightcols <- darkcols
	}
	
	pchs <- rep(c(21, 22, 23, 24, 25), ceiling(nlevels(strata) / 5))[1:nlevels(strata)]
	
	for (dist.name in dist.names) {
		dist.temp <- dist.obj[[dist.name]]
		
		if (!is.null(layout)) {
			par(mfrow=layout)
		} else {
			par(mfrow=c(1, length(sep.levels)))
		}
		
		for (sep.level in sep.levels) {
			ind <- df[, separate] == sep.level
			df2 <- df[ind, ]
			grp2 <- grp[ind]
			strata2 <- strata[ind]
			dist.temp2 <- dist.temp[ind, ind]
			
			if (!is.null(adj.name)) {
				adj <- as.data.frame(df2[, adj.name])
				obj <- cmdscale(as.dist(dist.temp2), k=ncol(dist.temp2)-1)
				dat2 <- apply(obj, 2, function(x) resid(lm(x ~ ., data=adj)))
				dist.temp2 <- dist(dat2)
			} 
			if (pca.method == 'cmd') {
				obj <- cmdscale(as.dist(dist.temp2), k=2, eig=T)
				pve <- round(obj$eig[1:2]/sum(abs(obj$eig))*100, 1)
				y <- cbind(obj$points[, 1], obj$points[, 2])
				
				xlab <- paste0('PC1(', pve[1], '%)')
				ylab <- paste0('PC2(', pve[2], '%)')
			} 
			
			if (pca.method == 'nmds') {
				obj <- metaMDS(as.dist(dist.temp2), k=2)
				y <- cbind(obj$points[, 1], obj$points[, 2])
				xlab <- 'NMDS1'
				ylab <- 'NMDS2'
			} 
			
			if (pca.method == 'pls') {
				require(mixOmics)
				# Test
				obj <- cmdscale(as.dist(dist.temp2), k=ncol(dist.temp2)-1)
				obj <- plsda(obj, grp2, ncomp=2)
				y <- cbind(obj$variates$X[, 1], obj$variates$X[, 2])
				xlab <- 'PLS1'
				ylab <- 'PLS2'
			}
			
			if (is.pvalue == TRUE & nlevels(factor(grp2)) > 1) {
				pvalue <- paste0('(P=', adonis(as.dist(dist.temp2) ~ grp2)$aov.tab[1, 6], ')')
			} else {
				pvalue <- ''
			}
			colnames(y) <- c("PC1", "PC2")		
			plot(y[, 1], y[, 2], type='n', xlim=range(y[, 1]) * 1.2, ylim=range(y[, 2])*1.2, 
					xlab=xlab, ylab=ylab)
#		points(y[, 1], y[, 2], bg='yellow', col=darkcols[grp2], 
#				type='p', pch = pchs[grp2], cex=cex.pt)
			col1 <- lightcols[1:nlevels(grp2)]
			col2 <- darkcols[1:nlevels(grp2)]
			col3 <- darkcols[grp2]
			if (!is.null(emp.lev)) {
				
				col1 <- col2
				col1[which(levels(grp2) != emp.lev)] <- rgb(0.5, 0.5, 0.5, 0.5)
				col1[which(levels(grp2) == emp.lev)] <- rgb(1.0, 0, 0, 1.0)
				col2[which(levels(grp2) != emp.lev)] <- rgb(0.5, 0.5, 0.5, 0.5)
				col2[which(levels(grp2) == emp.lev)] <- rgb(1.0, 0, 0, 1.0)
				col3[grp2 !=  emp.lev] <- rgb(0.5, 0.5, 0.5, 0.5)
				col3[grp2 ==  emp.lev] <- rgb(1.0, 0, 0, 1.0)
				
			}
			if (ellipse == T) {
				s.class(y, 
						fac = grp2,
						cstar = cstar,
						clab = 0,
						cpoint = 0,
						axesell = F,
						col = col1,
						grid = TRUE,
						add.plot=T,
						...
				)
			}
			
			points(y[, 1], y[, 2], bg=col3, col='black',
					type='p', pch = pchs[strata2], cex=cex.pt)
			if (indiv.lab) {
				if (is.null(emp.lev)) {
					lab.temp <- rownames(dist.temp2)
					text(y[, 1], y[, 2], labels=lab.temp, cex=indiv.lab.cex, col=rgb(1, 0, 0, 0.75))
				} else {
					lab.temp <- rownames(dist.temp2)
					lab.temp[grp2 !=  emp.lev] <- ''
					text(y[, 1], y[, 2], labels=lab.temp, cex=indiv.lab.cex, col=rgb(1, 0, 0, 0.75))
				}
				
			}
			s.class(y, 
					fac = grp2,
					cstar =0,
					cellipse = 0,
					clab = clab,
					cpoint = 0,
					axesell = F,
					col = col2,
					grid = TRUE,
					add.plot=T
			)
			if (!is.null(strata0)) {
				legend('topright', legend=(levels(strata2)), pch=pchs[1:nlevels(strata2)])	
			}
			
			title(main=sep.level, sub=paste(dist.name, "distance", pvalue))

		}


#		text(-0.25, 0.3, "PERMANOVA p=0.016")
	}
	
	if (pdf) {
		dev.off()
	}
}


# Rev:2017_04_19 Add bootstrap standard error

generate_distance_barplot <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
		grp.name, strata=NULL, within=T, between=T, bt.no = 100, ann='') {
	strata.name <- strata
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	grp.levels <- levels(grp)
	grp.nlevels <- nlevels(grp)
	grp.btws <- outer(grp.levels, grp.levels, paste, sep="#")
	grp.btws <- grp.btws[lower.tri(grp.btws)]
	
	if (!is.null(strata.name)) {
		strata <- df[, strata.name]
	} else {
		strata <- factor(rep(1, nrow(df))) # pseudo strata
	}
	res.df <- NULL
	for (dist.name in dist.names) {
		for (stratum in levels(strata)) {
			ind <- strata %in% stratum
			dist.sub <- dist.obj[[dist.name]][ind, ind]
			df2 <- df[ind, , drop=FALSE]
			if (between) {
				for (grp.btw in grp.btws) {
					ind1 <- which(df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[1])
					ind2 <- which(df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[2])	
					temp <- as.vector(dist.sub[ind1, ind2])
					sem <- sd(sapply(1:bt.no, function (i) {
								mean(dist.sub[as.numeric(sample(paste(ind1), repl=TRUE)), 
												as.numeric(sample(paste(ind2), repl=TRUE))])
							}))
					res.df <- rbind(res.df, data.frame(DistanceMetric=dist.name, Strata=stratum, Compare='Between', 
									DistanceType=grp.btw, Distance=mean(temp), sd=sem))
				}	
			}

			if (within) {
				for (grp.wth in grp.levels) {
					ind1 <- which(df2[, grp.name] == grp.wth)
					temp <- dist.sub[ind1, ind1]			
					temp <- temp[lower.tri(temp)]
					sem <- sd(sapply(1:bt.no, function (i) {
										ind2 <- as.numeric(sample(paste(ind1), repl = TRUE))
										temp <- dist.sub[ind2, ind2]			
										temp <- temp[lower.tri(temp)]
										mean(temp)
									}))
					res.df <- rbind(res.df, data.frame(DistanceMetric=dist.name, Strata=stratum, Compare='Within',
									DistanceType=grp.wth, Distance=mean(temp), sd=sem))
				}	
			}
		}
		
	}
	if (between & within) {
		res.df$DistanceType <- factor(res.df$DistanceType, levels=c(grp.btws, grp.levels))
		levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'), paste('Within\n', grp.levels))
	} else {
		if (between) {
			res.df$DistanceType <- factor(res.df$DistanceType)
			levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'))
		} else {
			res.df$DistanceType <- factor(res.df$DistanceType)
			levels(res.df$DistanceType) <- c(paste('Within\n', grp.levels))
		}
	}
	
	if (is.null(strata.name)) {
		pdf(paste0("Beta_diversity_btw", between, "_wth", within, "no_strata_barplot_", ann, ".pdf"),
				width=5, height=5)
		
		limits <- aes(ymax = Distance + sd, ymin=Distance - sd)
		dodge <- position_dodge(width=0.9)
		for (dist.name in dist.names) {
			cat(dist.name, "...\n")
			temp <- res.df[res.df$DistanceMetric == dist.name, ]
			obj1 <- ggplot(temp, aes(x=DistanceType, y=Distance, fill=DistanceType)) + 
					geom_bar(position=dodge, stat="identity", width=0.75) + 
					geom_bar(position=dodge, stat="identity", width=0.75, colour="black", show_guide=FALSE, size=0.25) +
					geom_errorbar(limits, position=dodge, size=0.25, width=0.25) +
					labs(y=paste(dist.name, "Distance"), x='') +
					theme(legend.position="none") +
					theme(axis.text.x=element_text(angle=90, hjust=1))
					
			print(obj1)
			
		}
		dev.off()
	} else {
		pdf(paste0("Beta_diversity_btw", between, "_wth", within, "_strata", strata.name, "barplot_", ann, ".pdf"),
				width=2.5*(nlevels(strata)-1) + 5, height=5)
		
		limits <- aes(ymax = Distance + sd, ymin=Distance - sd)
		dodge <- position_dodge(width=0.9)
		for (dist.name in dist.names) {
			cat(dist.name, "...\n")
			temp <- res.df[res.df$DistanceMetric == dist.name, ]
			obj1 <- ggplot(temp, aes(x=Strata, y=Distance, fill=DistanceType)) + 
					geom_bar(position=dodge, stat="identity") + 
					geom_bar(position=dodge, stat="identity", colour="black", show_guide=FALSE, size=0.25) +
					geom_errorbar(limits, position=dodge, size=0.25, width=0.5) +
					labs(y=paste(dist.name, "Distance"), x=strata.name) +
					theme(axis.text.x=element_text(angle=90, hjust=1))
			print(obj1)
			
		}
		dev.off()
		
	}
	
}

generate_distance_boxplot <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
		grp.name, strata=NULL, within=F, between=T, ann='') {
	strata.name <- strata
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	grp.levels <- levels(grp)
	grp.nlevels <- nlevels(grp)
	grp.btws <- outer(grp.levels, grp.levels, paste, sep="#")
	grp.btws <- grp.btws[lower.tri(grp.btws)]
	
	if (!is.null(strata.name)) {
		strata <- df[, strata.name]
	} else {
		strata <- factor(rep(1, nrow(df))) # pseudo strata
	}
	res.df <- NULL
	for (dist.name in dist.names) {
		for (stratum in levels(strata)) {
			ind <- strata %in% stratum
			dist.sub <- dist.obj[[dist.name]][ind, ind]
			df2 <- df[ind, , drop=FALSE]
			if (between) {
				for (grp.btw in grp.btws) {
					ind1 <- df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[1]
					ind2 <- df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[2]	
					temp <- as.vector(dist.sub[ind1, ind2])
					n <- length(temp)
					res.df <- rbind(res.df, data.frame(DistanceMetric=rep(dist.name, n), Strata=rep(stratum, n), Compare=rep('Between', n),
									DistanceType=rep(grp.btw, n), Distance=temp))
				}	
			}
			
			if (within) {
				for (grp.wth in grp.levels) {
					ind1 <- df2[, grp.name] == grp.wth
					temp <- dist.sub[ind1, ind1]			
					temp <- temp[lower.tri(temp)]
					n <- length(temp)
					res.df <- rbind(res.df, data.frame(DistanceMetric=rep(dist.name, n), Strata=rep(stratum, n), Compare=rep('Within', n),
									DistanceType=rep(grp.wth, n), Distance=temp))
				}	
			}
		}
		
	}
	if (between & within) {
		res.df$DistanceType <- factor(res.df$DistanceType, levels=c(grp.btws, grp.levels))
		levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'), paste('Within\n', grp.levels))
	} else {
		if (between) {
			res.df$DistanceType <- factor(res.df$DistanceType)
			levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'))
		} else {
			res.df$DistanceType <- factor(res.df$DistanceType)
			levels(res.df$DistanceType) <- c(paste('Within\n', grp.levels))
		}
	}
	
	if (is.null(strata.name)) {
		pdf(paste0("Beta_diversity_btw", between, "_wth", within, "no_strata_boxplot_", ann, ".pdf"),
				width=5, height=5)
		
		for (dist.name in dist.names) {
			cat(dist.name, "...\n")
			temp <- res.df[res.df$DistanceMetric == dist.name, ]
			
			dodge <- position_dodge(width=0.95)		
			obj1 <- ggplot(temp, aes(x=DistanceType, y=Distance, col=DistanceType)) + 
					geom_boxplot(position=dodge, outlier.colour = NA) + 
#					geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) +
					labs(y=paste(dist.name, "Distance"), x='') +
					theme(legend.position="none")
			
			print(obj1)
			
		}
		dev.off()
	} else {
		pdf(paste0("Beta_diversity_btw", between, "_wth", within, "_strata", strata.name, "boxplot_", ann, ".pdf"),
				width=2.5*(nlevels(strata)-1) + 5, height=5)
		
		for (dist.name in dist.names) {
			cat(dist.name, "...\n")
			temp <- res.df[res.df$DistanceMetric == dist.name, ]
			dodge <- position_dodge(width=0.95)		
			obj1 <- ggplot(temp, aes(x=Strata, y=Distance, col=DistanceType)) + 
					geom_boxplot(position=dodge, outlier.colour = NA) + 
		#			geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) +
					labs(y=paste(dist.name, "Distance"), x=strata.name)
			print(obj1)
			
		}
		dev.off()
		
	}
	
}


# Rev:2016_11_25
generate_clustering <- function(data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), meta.info, cluster.method='average', 
		is.labRow=F, cex.lab=NULL, ann="", wid=10, hei=6, colFnsF=NULL, colFnsC=NULL) {
	if (is.null(colFnsC)) {
		colFnsC <- c(colorRampPalette(c('black', 'green')), colorRampPalette(c('black', 'blue')), colorRampPalette(c('black', 'red')))
	} 
	if (is.null(colFnsF)) {
		rainbow3 <- function(x) {
			rainbow(x + 2)[1:x]
		}
		jet3 <- function(x) {
			jet(x + 2)[1:x]
		}
		colFnsF <- c(rainbow3, jet3)
	}
	
	pdf(paste0("Beta_diversity_Hierachical_clustering_", ann, ".pdf"), width=wid, height=hei)
	for (dist.name in dist.names) {

		df <- data.obj$meta.dat
		dist.sub <- dist.obj[[dist.name]]
		dend <- hclust(as.dist(dist.sub), cluster.method)
		

		key.list <- list()
		mat <- NULL
		for (keyID in meta.info) {
			x <- df[, keyID]
			i <- 0
			j <- 0
			if (is.factor(x)) {
				key.list[[keyID]] <- list(breaks=levels(x), colors=(colFnsF[i+1][[1]])(nlevels(x)), base=NA, col.na=NA, right=F, include.lowest=F)	
				i <- (i + 1) %% length(colFnsF)
				mat <- cbind(mat, key.list[[keyID]]$colors[x])
			} else {
				key.list[[keyID]] <- makecmap(x, n=5, colFn=colFnsC[j+1][[1]])
				j <- (j + 1) %% length(colFnsC)
				mat <- cbind(mat, cmap(x, key.list[[keyID]]))
			}
		}
		colnames(mat) <- meta.info

		par(oma = c(1, 2, 1, 2))
		par(mar = c(5,4,4,3)+0.1)  # make space for color keys
		if (is.labRow == F) {
			labRow <- ""
		} else {
			labRow <- rownames(df)
		}
		if (is.null(cex.lab)) {
			cex.lab <- 25 / nrow(df) 
		}
		
		dendromat(dend, mat, labRow=labRow,
				ylab = 'Distance', main = paste(dist.name, "distance"), cex.lab=cex.lab)
		
		par(oma=c(0, 0, 0, 0))
		y.cord <- (1/length(meta.info)) * (0:(length(meta.info) - 1))
		k <- 1
		for (keyID in meta.info) {
			x <- df[, keyID]
			if (is.factor(x)) {
				vkey2(key.list[[keyID]], keyID, y=y.cord[k], stretch=1.2)
			} else {
				vkey(key.list[[keyID]], keyID, y=y.cord[k], stretch=1.2)
			}
			k <- k + 1
		}		
	}
	dev.off()
}


# Rev: 2016_09_10, Implement p value based omnibus test
PermanovaG2 <- function (formula, dat = NULL, ...) 
{
	save.seed <- get(".Random.seed", .GlobalEnv)
	lhs <- formula[[2]]
	lhs <- eval(lhs, dat, parent.frame())
	rhs <- as.character(formula)[3]
	p.perms <- list()
	p.obs <- list()
	for (i in 1:(dim(lhs)[3])) {
		assign(".Random.seed", save.seed, .GlobalEnv)
		Y <- as.dist(lhs[, , i])
		formula2 <- as.formula(paste("Y", "~", rhs))
		obj <- adonis(formula2, dat, ...)
		perm.mat <- obj$f.perms
		p.perms[[i]] <- 1 - (apply(perm.mat, 2, rank) - 1) / nrow(perm.mat)
		p.obs[[i]] <- obj$aov.tab[1:ncol(perm.mat), "Pr(>F)"]

	}
	
	omni.pv <- NULL
	indiv.pv <- NULL
	for (j in 1:ncol(perm.mat)) {
		p.perms.j <- sapply(p.perms, function (x) x[, j])
		p.obj.j <- sapply(p.obs, function (x) x[j])
		omni.pv <- c(omni.pv, mean(c(rowMins(p.perms.j ) <= min(p.obj.j), 1)))
		indiv.pv <- rbind(indiv.pv, p.obj.j)
	}
	colnames(indiv.pv) <- paste0('D', 1:ncol(indiv.pv), '.p.value')
	rownames(indiv.pv) <- 1:nrow(indiv.pv)
	
	aov.tab <- data.frame(indiv.pv, omni.p.value = omni.pv)
	rownames(aov.tab) <- rownames(obj$aov.tab)[1:ncol(perm.mat)]
	list(aov.tab = aov.tab)
}


perform_permanova_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		PermanovaG.dist=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
		formula=NULL,  grp.name=NULL, adj.name=NULL, pairwise=F, block.perm=F, strata=NULL, ann='', ...) {
	# PermanovaG not implemented for block permutation
	result <- list()
	
	df <- data.obj$meta.dat
	if (!is.null(strata)) {
		if (is.character(strata)) {
			strata <- df[, strata]
		}
	}
	if (!is.null(formula)) {
		ind <- apply(df[, gsub("^\\s+|\\s+$", "", strsplit(strsplit(formula, "~")[[1]][2], "\\+")[[1]]), drop=F], 1, function(x) sum(is.na(x))) == 0
		df <- df[ind, ]
		if (!is.null(strata)) {
			strata <- strata[ind]
		}
		sink(paste0('Beta_diversity_PERMANOVA_test_', ann, '.txt'))
		date()
		cat('\nPERMANOVA test: \n')
		permanova.obj <- list()
		for (dist.name in dist.names) {
			cat(dist.name, " distance: \n")
			dist.mat <- as.dist(dist.obj[[dist.name]][ind, ind])
			if (block.perm == F) {
				obj <- adonis(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
			} else {
				obj <- adonis2(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
			}
			
			prmatrix(obj$aov.tab)
			permanova.obj[[dist.name]] <- obj$aov.tab
			cat("\n")
		}
		result$permanova.obj <- permanova.obj
		permanovaG.obj <- NULL
		if (block.perm == F & !is.null(PermanovaG.dist)) {
			cat('\nPERMANOVA G test combining ', paste(PermanovaG.dist, collapse=','), '\n')
			response <- array(NA, c(sum(ind), sum(ind), length(PermanovaG.dist)), dimnames=list(NULL, NULL, PermanovaG.dist))
			for (dist.name in PermanovaG.dist) {
				response[, , dist.name] <- dist.obj[[dist.name]][ind, ind]
			}
			obj <- PermanovaG2(as.formula(paste("response", formula)), df,  strata=strata, ...)
			prmatrix(obj$aov.tab)
			permanovaG.obj <- obj$aov.tab
			cat("\n")
			result$permanovaG.obj <- permanovaG.obj
		}
		cat("\n")
		sink()
	} else {
		if (pairwise == F) {
			if (is.null(adj.name)) {
				formula <- paste('~', grp.name)
			} else {
				formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
			}
			
			ind <- apply(df[, gsub("^\\s+|\\s+$", "", strsplit(strsplit(formula, "~")[[1]][2], "\\+")[[1]]), drop=F], 1, function(x) sum(is.na(x))) == 0
			df <- df[ind, ]
			if (!is.null(strata)) {
				strata <- strata[ind]
			}
			sink(paste0('Beta_diversity_PERMANOVA_test_', ann, '.txt'))
			date()
			cat('\nPERMANOVA test: \n')
			permanova.obj <- list()
			for (dist.name in dist.names) {
				cat(dist.name, " distance: \n")
				dist.mat <- as.dist(dist.obj[[dist.name]][ind, ind])
				if (block.perm == F) {
					obj <- adonis(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
				} else {
					obj <- adonis2(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
				}
				prmatrix(obj$aov.tab)
				permanova.obj[[dist.name]] <- obj$aov.tab
				cat("\n")
			}
			result$permanova.obj <- permanova.obj
			permanovaG.obj <- NULL
			if (block.perm == F & !is.null(PermanovaG.dist)) {
				cat('\nPERMANOVA G test combining ', paste(PermanovaG.dist, collapse=','), '\n')
				response <- array(NA, c(sum(ind), sum(ind), length(PermanovaG.dist)), dimnames=list(NULL, NULL, PermanovaG.dist))
				for (dist.name in PermanovaG.dist) {
					response[, , dist.name] <- dist.obj[[dist.name]][ind, ind]
				}
				obj <- PermanovaG2(as.formula(paste("response", formula)), df,  strata=strata, ...)
				prmatrix(obj$aov.tab)
				permanovaG.obj <- obj$aov.tab
				cat("\n")
				result$permanovaG.obj <- permanovaG.obj
			}
			
			cat("\n")
			sink()
			
		} else {
			sink(paste0('Beta_diversity_PERMANOVA_test_', ann, '_pairwise.txt'))
			date()
			cat('\nPairwise PERMANOVA test: \n')
			grp <- factor(df[, grp.name])
			grp.levels <- levels(grp)
			grp.nlevels <- nlevels(grp)
			pmat.all <- NULL
			rmat.all <- NULL
			pmat.G <- matrix(NA, grp.nlevels, grp.nlevels)
			colnames(pmat.G) <- rownames(pmat.G) <- grp.levels
			for (dist.name in dist.names) {
				cat(dist.name, " distance: \n")
				pmat <- matrix(NA, grp.nlevels, grp.nlevels)
				colnames(pmat) <- rownames(pmat) <- grp.levels
				rmat <- matrix(NA, grp.nlevels, grp.nlevels)
				colnames(rmat) <- rownames(rmat) <- grp.levels
				for (i in 1:(grp.nlevels-1)) {
					grp.level1 <- grp.levels[i]
					for (j in (i+1):grp.nlevels) {
						
						grp.level2 <- grp.levels[j]
						cat(grp.level1, ' vs ', grp.level2, '\n')
						ind <- grp %in% c(grp.level1, grp.level2)
						df2 <- subset(df, ind)
						df2[, grp.name] <- factor(df2[, grp.name])
						dist.mat <- dist.obj[[dist.name]][ind, ind]
						strata2 <- strata[ind]
						
						if (is.null(adj.name)) {
							formula <- paste('~', grp.name)
						} else {
							formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
						}
						
						ind2 <- apply(df2[, gsub("^\\s+|\\s+$", "", strsplit(strsplit(formula, "~")[[1]][2], "\\+")[[1]]), drop=F], 1, function(x) sum(is.na(x))) == 0
						df2 <- df2[ind2, ]
						dist.mat2 <- as.dist(dist.mat[ind2, ind2])
						strata2 <- strata2[ind2]
						if (block.perm == F) {
							obj <- adonis(as.formula(paste("dist.mat2", formula)), df2,  strata=strata2, ...)
						} else {
							obj <- adonis2(as.formula(paste("dist.mat2", formula)), df2,  strata=strata2,  ...)
						}
						prmatrix(obj$aov.tab)
						cat("\n")
	
						if (block.perm == F) {
							pmat[i, j] <- pmat[j, i] <- obj$aov.tab[length(adj.name)+1, 6]
							rmat[i, j] <- rmat[j, i] <- obj$aov.tab[length(adj.name)+1, 5]
						} else {
							pmat[i, j] <- pmat[j, i] <- obj$aov.tab[1, 6]
							rmat[i, j] <- rmat[j, i] <- obj$aov.tab[1, 5]
						}

						# PERMANOVA G after last distance
						if (block.perm == F & !is.null(PermanovaG.dist)) {
							if (dist.name == dist.names[length(dist.names)]) {
								cat('\nPERMANOVA G test combining ', paste(PermanovaG.dist, collapse=','), '\n')
								response <- array(NA, c(sum(ind2), sum(ind2), length(PermanovaG.dist)), dimnames=list(NULL, NULL, PermanovaG.dist))
								for (dist.name in PermanovaG.dist) {
									response[, , dist.name] <- dist.mat[ind2, ind2]
								}
								obj <- PermanovaG2(as.formula(paste("response", formula)), df2,  strata=strata2, ...)
								prmatrix(obj$aov.tab)
								cat("\n")
								pmat.G[i, j] <- pmat.G[j, i] <- obj$aov.tab[length(adj.name)+1, 'omni.p.value']
							}
						}
					}
				}
				cat("\n")
				pmat.all <- rbind(pmat.all, c(dist.name, rep('', grp.nlevels-1)), formatC(pmat), rep("", grp.nlevels))
				rmat.all <- rbind(rmat.all, c(dist.name, rep('', grp.nlevels-1)), formatC(rmat), rep("", grp.nlevels))
			}
			cat("\n")
			sink()
			write.csv(pmat.all, paste0('Beta_diversity_PERMANOVA_test_', ann, '_pairwiseP.csv'))
			write.csv(rmat.all, paste0('Beta_diversity_PERMANOVA_test_', ann, '_pairwiseR.csv'))
			write.csv(pmat.G, paste0('Beta_diversity_PERMANOVA_G_test_', ann, '_pairwiseP.csv'))
		}
	}
	
	return(invisible(result))
}

randomForestTest <- function(x, y, perm.no=999, ...) {
	iris.rf <- randomForest(x=x, y=y, importance=FALSE)
	to1 <- mean(abs((as.numeric(y) - 1) - predict(iris.rf,  type='prob')[, 2])^2) #mean squarred error
	to2 <- mean(y !=  predict(iris.rf,   type='response'))  #mean prediction error
	tp <- sapply(1:perm.no, function(i) {
				if (i %% 10 == 0) cat('.')
				y.p <- sample(y)
				iris.rf <- randomForest(x=x, y=y.p, importance=FALSE)
				
				t1 <- mean(abs((as.numeric(y.p) - 1) - predict(iris.rf,  type='prob')[, 2])^2) 
				t2 <- mean(y.p !=  predict(iris.rf,   type='response'))
				c(t1, t2)
			}
	
	)
	pv1 <- (sum(tp[1, ] <= to1) + 1) / (perm.no + 1)
	pv2 <- (sum(tp[2, ] <= to2) + 1) / (perm.no + 1)
	c(pv1=pv1, pv2=pv2)
}

perform_rf_test <- function (data.obj, grp.name, taxa.level='Genus', perm.no=999, prev=0.1, minp=0.000, ann='',...) {
	if (taxa.level == 'Species') {
		if (taxa.level %in% names(data.obj$abund.list)) {
			ct <- data.obj$abund.list[[taxa.level]]
		} else {
			# Accomodate different version
			ct <- data.obj$otu.tab
			rownames(ct) <- paste0("OTU", rownames(ct), ":", data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
			data.obj$abund.list[['Species']] <- ct
		}
	} else {
		ct <- data.obj$abund.list[[taxa.level]]
	}
	
	prop <- t(t(ct) / colSums(ct))
	prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
	prop <- t(prop)
	
	grp <- data.obj$meta.dat[, grp.name]
	if (!is.factor(grp)) stop('Current random Forest test is designed to deal with binary factor data!\n')

	cat('RF test P values  ...\n')
	obj <- randomForestTest(prop, grp, perm.no, ...)
	sink(paste0('OverallAssoc_RF_test_', taxa.level, '_', ann, '.txt'))
	cat ('Prediction Probability:', obj[1], '\n')
	cat ('Binary Response:', obj[2], '\n')
	sink()
}

# Rev: 2016_12_02, MiKRAT for binary result, add out_type='D'
perform_mirkat_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
		grp.name=NULL, adj.name=NULL, pairwise=F,  ann='', ...) {
	
	# MiRKAT not implemented for correlated data
	df <- data.obj$meta.dat
	
	if (pairwise == F) {
		
		ind <- apply(df[, c(grp.name, adj.name), drop=F], 1, function(x) sum(is.na(x))) == 0
		df <- df[ind, ]
		
		grp <- df[, grp.name]
		
		if (is.character(grp)) {
			grp <- factor(grp)
		}
		# Rev: 2016_12_02
		if (is.factor(grp)) {
			if (nlevels(grp) > 2) {
				stop('Currently MiRKAT only supports binary outcome!')
			} else {
				grp <- as.numeric(grp) - 1
				out_type <- 'D'
			}
		} else {
			out_type <- 'C'
		}
		if (!is.null(adj.name)) {
			# No intercept
			adj <- model.matrix(~ ., data.frame(df[, adj.name]))
			# Remove collinear terms
			qadj <- qr(adj, tol = 1e-07)
			adj <- adj[, qadj$pivot, drop = FALSE]
			adj <- adj[, 1:qadj$rank, drop = FALSE]
			# Remove intercept
			adj <- adj[, colSums(adj==1) != nrow(adj)]
		} else {
			adj <- NULL
		}
		
		sink(paste0('Beta_diversity_MiRKAT_test_', ann, '.txt'))
		date()
		cat('\nMiRKAT test combining ', paste(dist.names, collapse=','), '\n')
		Ks <- list()
		for (dist.name in dist.names) {
			Ks[[dist.name]] <- D2K(dist.obj[[dist.name]][ind, ind])
		}
		# Rev: 2016_12_02
		obj <- MiRKAT(grp, X=adj, Ks, out_type=out_type)
		cat('Individual P value: ')
		prmatrix(t(obj$indivP))
		cat('\nOmnibus P value: ')
		cat(obj$omnibus_p)
		cat("\n")
		sink()
		
	} else {
		sink(paste0('Beta_diversity_MiRKAT_test_', ann, '_pairwise.txt'))
		date()
		cat('\nPairwise MiRKAT test: \n')
		grp <- factor(df[, grp.name])
		grp.levels <- levels(grp)
		grp.nlevels <- nlevels(grp)
		pmat.G <- matrix(NA, grp.nlevels, grp.nlevels)
		colnames(pmat.G) <- rownames(pmat.G) <- grp.levels
		parr <- array(NA, c(grp.nlevels, grp.nlevels, length(dist.names)), dimnames=list(grp.levels, grp.levels, dist.names))
		
		for (i in 1:(grp.nlevels-1)) {
			grp.level1 <- grp.levels[i]
			for (j in (i+1):grp.nlevels) {
				
				grp.level2 <- grp.levels[j]
				cat(grp.level1, ' vs ', grp.level2, '\n')
				ind <- grp %in% c(grp.level1, grp.level2)
				df2 <- subset(df, ind)
				df2[, grp.name] <- factor(df2[, grp.name])
				ind2 <- apply(df2[, c(grp.name, adj.name), drop=F], 1, function(x) sum(is.na(x))) == 0
				df2 <- df[ind2, ]
				
				grp2 <- as.numeric(df2[, grp.name]) - 1
				if (!is.null(adj.name)) {
					# No intercept
					adj <- model.matrix(~ ., data.frame(df2[, adj.name]))
					# Remove collinear terms
					qadj <- qr(adj, tol = 1e-07)
					adj <- adj[, qadj$pivot, drop = FALSE]
					adj <- adj[, 1:qadj$rank, drop = FALSE]
					# Remove intercept
					adj <- adj[, colSums(adj==1) != nrow(adj)]
				} else {
					adj <- NULL
				}
				
				Ks <- list()
				for (dist.name in dist.names) {
					Ks[[dist.name]] <- D2K(dist.obj[[dist.name]][rownames(df2), rownames(df2)])
				}
				# Rev: 2016_12_02
				obj <- MiRKAT(grp2, X=adj, Ks, out_type='D')
				pmat.G[i, j] <- pmat.G[j, i] <- obj$omnibus_p
				parr[i, j, ] <- parr[j, i, ] <- obj$indivP
				cat('Individual P value: ')
				prmatrix(t(obj$indivP))
				cat('\nOmnibus P value: ')
				cat(obj$omnibus_p)
				cat("\n")
				
			}
		}
		cat("\n")
		sink()
		
		pmat.all <- NULL
		for (dist.name in dist.names) {
			pmat <- parr[, , dist.name]
			pmat.all <- rbind(pmat.all, c(dist.name, rep('', grp.nlevels-1)), formatC(pmat), rep("", grp.nlevels))
		}
		
		write.csv(pmat.all, paste0('Beta_diversity_MiRKAT_test_', ann, '_pairwiseP.csv'))
		write.csv(pmat.G, paste0('Beta_diversity_MiRKAT_O_test_', ann, '_pairwiseP.csv'))
	}
}


perform_betadisper_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), grp.name) {
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	sink('Beta_diversity_BETADISPER_test.txt')
	date()
	cat('\nBetadisper test: \n')
	for (dist.name in dist.names) {
		cat(dist.name, " distance: \n")
		dist.mat <- as.dist(dist.obj[[dist.name]])
		obj <- betadisper(dist.mat, grp)
		prmatrix(anova(obj))
		cat("\n")
	}
	cat("\n")
	sink()
}

distance_compare_test <- function (dist.mat, ind1, ind2, ind3, ID2=NULL, ID3=NULL, alternative='greater', nperm=999) {
	
	ind23 <- c(ind2, ind3)
	if (!is.null(ID2) & !is.null(ID3)) {
		ID23 <- c(ID2, ID3)
		ID2u <- unique(ID2)
		ID3u <- unique(ID3)
		ID23u <- c(ID2u, ID3u)
		n2u <- length(ID2u)
		n3u <- length(ID3u)
	}
	n2 <- length(ind2)
	n3 <- length(ind3)
	
	dist12 <- dist.mat[ind1, ind2]
	dist13 <- dist.mat[ind1, ind3]
	
	stat.obs <- mean(dist12) - mean(dist13)
	stat.perm <- sapply(1:nperm, function(i) {
				if (is.null(ID2) & is.null(ID3)) {
					ind23.p <- sample(ind23)
					ind2.p <- ind23.p[1:n2]
					ind3.p <- ind23.p[(n2+1):(n2+n3)]
				} else {
					ID23u.p <- sample(ID23u)
					ID2.p <- ID23u.p[1:n2u]
					ID3.p <- ID23u.p[(n2u+1):(n2u+n3u)]
					ind2.p <- ind23[ID23 %in% ID2.p]
					ind3.p <- ind23[ID23 %in% ID3.p]
				}

				dist12.p <- dist.mat[ind1, ind2.p]
				dist13.p <- dist.mat[ind1, ind3.p]
				mean(dist12.p) - mean(dist13.p)
			})
	if (alternative == 'greater') {
		pv <- mean(c(stat.perm >= stat.obs, TRUE))
	} else {
		pv <- mean(c(stat.perm <= stat.obs, TRUE))
	}  
	pv
}

perform_distance_compare_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
		grp.name, level1, level2, level3, subject=NULL, alternative='greater', nperm=999, seed=123, ann='') {
	sink(paste0('Beta_diversity_test_', level1, '-', level2, '_', alternative, '_', level1, '-', level3, '_', ann, '.txt'))
	set.seed(seed)
	cat('Testing the distance ', level1, '-', level2, ' is ', alternative, ' than ', level1, '-', level3, '\n')
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	IDs <- df[, subject]
	grp.levels <- levels(grp)
	grp.nlevels <- nlevels(grp)
	ind1 <- which(grp == level1)
	ind2 <- which(grp == level2)
	ind3 <- which(grp == level3)
	if (is.null(subject)) {
		ID2 <- NULL
		ID3 <- NULL
	} else {
		ID2 <- as.character(IDs[grp == level2])
		ID3 <- as.character(IDs[grp == level3])
	}
	for (dist.name in dist.names) {
		cat(dist.name, ' Distance:p=') 
		pv <- distance_compare_test(dist.obj[[dist.name]], ind1, ind2, ind3, ID2, ID3, alternative, nperm)
		cat(pv, '\n')
	}
	sink()
	
}


perform_taxa_compare_test <- function (data.obj, grp.name, level1, level2, level3, alternative='greater', nperm=1000, seed=123,
		taxa.levels=c('Phylum', 'Family', 'Genus'), taxa.name='All', prev=0.1, minp=0.002, ann='') {
	
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	grp.levels <- levels(grp)
	grp.nlevels <- nlevels(grp)
	ind1 <- which(grp == level1)
	ind2 <- which(grp == level2)
	ind3 <- which(grp == level3)
	
	pv.list <- list()
	
	for (LOI in taxa.levels) {
		cat(LOI, "\n")
		prop <- data.obj$abund.list[[LOI]]
		prop <- t(t(prop) / colSums(prop))
		pv.vec <- NULL
		if (taxa.name == 'All') {
			prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
		} else {
			prop <- prop[taxa.name, , drop=FALSE]
		}
		
		for (taxon in rownames(prop)) {
			cat(".")
			dist.mat <- as.matrix(dist(rank(prop[taxon, ])))
			pv <- distance_compare_test(dist.mat, ind1, ind2, ind3, alternative, nperm)
			pv.vec <- c(pv.vec, pv)
		}
		names(pv.vec) <- rownames(prop)
		pv.list[[LOI]] <- pv.vec
		cat("\n")
	}
	
	for (LOI in taxa.levels) {
		pv.vec <- pv.list[[LOI]]

		qv.vec <- p.adjust(pv.vec, 'fdr')

		
		write.csv(data.frame("p value"=pv.vec, "q value"=qv.vec), paste0("Taxa_TrendAnalysis_", LOI, "_", ann, ".csv"))
	}
	
}

# New: 2017_11_03
# Now default is averaged distance for within or between-subject comparison
perform_site_correlation_test <- function (data.obj, dist.obj, site.name, sites, subject, 
		dist.names = c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), dist.ave = TRUE, nPerm = 999, ann = '') {

	pmat <- array(NA, c(length(dist.names), length(sites), length(sites)), dimnames=list(dist.names, sites, sites))
	for (site1 in sites) {
		for(site2 in sites) {
			if (site1 != site2 & is.na(pmat[1, site1, site2])) {
				cat('.')
				sampling.site <- c(site1, site2)
				samIDs <- data.obj$meta.dat[, site.name] %in% sampling.site 
				data.obj1 <- subset_data(data.obj, samIDs)
				dist.obj1 <- subset_dist(dist.obj, samIDs)
				
# Compute the observed distance within subject
				lIDs <- as.character(data.obj1$meta.dat[, site.name])
				pIDs <- as.character(data.obj1$meta.dat[, subject])
				uIDs <- as.character(unique(data.obj1$meta.dat[, subject]))
				
				
				pdf(paste('Permutation_test', site1, site2, 'Correlation_', ann, '.pdf', sep='_'), width=6, height=6)
				for (dist.name in dist.names) {
					cat(dist.name, '...\n')
					
					dist.mat <- dist.obj1[[dist.name]]					
# To self
					dist.o1 <- NULL
					for (uID in uIDs) {
						ind1 <- pIDs == uID & lIDs == site1
						ind2 <- pIDs == uID & lIDs == site2
						
						if (sum(ind1) != 0 & sum(ind2) != 0)
							if (dist.ave) {
								dist.o1 <- c(dist.o1, mean(dist.mat[ind1, ind2]))
							} else {
								dist.o1 <- c(dist.o1, dist.mat[ind1, ind2])
							}
							
					}
					
# To other people
					dist.o2 <- NULL
					for (uID in uIDs) {
						ind1 <- pIDs == uID & lIDs == site1
						ind2 <- pIDs != uID & lIDs == site2
						
						if (sum(ind1) != 0 & sum(ind2) != 0)
							if (dist.ave) {
								dist.o2 <- c(dist.o2, mean(dist.mat[ind1, ind2]))
							} else {
								dist.o2 <- c(dist.o2, dist.mat[ind1, ind2])
							}
							
					}
					
					if (!is.null(dist.o1) & !is.null(dist.o2)) {
						boxplot(list(ToSelf=dist.o1, ToOther=dist.o2), col='steelblue', ylab=paste(dist.name, 'distance'))
						stat.o <- mean(dist.o2) - mean(dist.o1)
						
						stat.p <- NULL
						for (i in 1:nPerm) {
							# Exchange the uterus between people :D
							pIDs.p <- pIDs
							pIDs.uterus <- pIDs[lIDs == site1]
							temp <- factor(pIDs.uterus)
							levels(temp) <- sample(levels(temp))
							pIDs.uterus.p <- as.character(temp)
							pIDs.p[lIDs == site1] <- pIDs.uterus.p
							
							# To increase the number of permutation
							pIDs.uterus <- pIDs[lIDs == site2]
							temp <- factor(pIDs.uterus)
							levels(temp) <- sample(levels(temp))
							pIDs.uterus.p <- as.character(temp)
							pIDs.p[lIDs == site2] <- pIDs.uterus.p
							
							# To self
							dist.o1 <- NULL
							for (uID in uIDs) {
								ind1 <- pIDs.p == uID & lIDs == site1
								ind2 <- pIDs.p == uID & lIDs == site2
								
								if (sum(ind1) != 0 & sum(ind2) != 0)
									if (dist.ave) {
										dist.o1 <- c(dist.o1, mean(dist.mat[ind1, ind2]))
									} else {
										dist.o1 <- c(dist.o1, dist.mat[ind1, ind2])
									}
							}
							
# To other people
							dist.o2 <- NULL
							for (uID in uIDs) {
								ind1 <- pIDs.p == uID & lIDs == site1
								ind2 <- pIDs.p != uID & lIDs == site2
								
								if (sum(ind1) != 0 & sum(ind2) != 0)
									if (dist.ave) {
										dist.o2 <- c(dist.o2, mean(dist.mat[ind1, ind2]))
									} else {
										dist.o2 <- c(dist.o2, dist.mat[ind1, ind2])
									}
							}
							
							# boxplot(list(ToSelf=dist.o1, ToOther=dist.o2))
							stat.p <- c(stat.p, mean(dist.o2) - mean(dist.o1))
						}
						
						pv <- mean(c(stat.p >= stat.o, TRUE))
						
						hist(stat.p, xlab=paste(dist.name, 'distance'), col='steelblue', xlim=c(min(stat.p) * 2, max(stat.p) * 2),
								main='Average distances observed under permutation', sub=paste('P < ', pv))
						
						abline(v=stat.o, col='red')
						pmat[dist.name, site1, site2] <- pmat[dist.name, site2, site1] <- pv
					}

				}
				dev.off()
			}
			
		}
	}
	cat('Finished!\n')
	
	pmat.csv <- NULL
	for (dist.name in dist.names) {
		temp <- as.matrix(pmat[dist.name, ,  ])
		mode(temp) <- 'character'
		pmat.csv <- rbind(pmat.csv, rbind(Distance=c(dist.name, rep('', ncol(temp) - 1)), temp, ''))
	}
	
	write.csv(pmat.csv, paste0('SiteCorrelationP_', ann, '.csv'))
	return(invisible(pmat))
}

# New: 2017_11_03
perform_site_taxa_correlation_test <- function (data.obj, site.name, sites, subject, 
		taxa.levels = c('Genus'), taxa.names = NULL, dist.func = function (x) {dist(sqrt(x))}, 
		prev=0.1, minp=0.002, medianp=NULL, ann = '', ...) {
	for (taxa.level in taxa.levels) {
		cat(taxa.level, '...\n')
		genus <- data.obj$abund.list[[taxa.level]]
		genus <- t(t(genus) / colSums(genus))
		if (is.null(taxa.names)) {
			ind <- rep(TRUE, nrow(genus))
			if (!is.null(prev)) ind <- ind & rowMeans(genus != 0) > prev
			if (!is.null(minp)) ind <- ind & rowMaxs(genus) > minp
			if (!is.null(medianp)) {
				nz.mean <- apply(genus, 1, function(x) median(x[x!=0]))
				ind <- ind & nz.mean > medianp
			}
			genus <- genus[ind, , drop = FALSE]
			gen.names <- rownames(genus)
		} else {
			if (sum(!(taxa.names %in% rownames(genus))) != 0) {
				stop('Could not find the taxa! Please check!\n')
			}
			genus <- genus[taxa.names, , drop = FALSE]
			gen.names <- taxa.names
		}

		# Convert into distance
	    dist.obj <- list()
		for (gen.name in gen.names) {
			dist.obj[[gen.name]] <- as.matrix(dist.func(genus[gen.name, ]))
		}
		
		perform_site_correlation_test(data.obj, dist.obj, site.name, sites, subject, 
				dist.names = gen.names,  ann = paste(taxa.level, ann, sep = '_'), ...)
		
	}
}

# Rev: 2017_02_19 Add hei and wid
generate_taxa_boxplot <- function (data.obj,  grp.name, strata=NULL, scale='P',  taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'), 
		taxa.name='All', pseudo.ct=0.5, rm.outlier=T, prev=0.1, minp=0.002, ann='All', subject=NULL, l.size=0.5, p.size=2.5, hei0=NULL, wid0=NULL) {
	# To be completed
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	if (!is.null(subject)) {
		ID <- df[, subject]
	}
	
	for (LOI in taxa.levels) {
		if (LOI == 'All') {
			if (taxa.name == 'All') stop('Taxa names are not required to be all when taxa level is also all!\n')
			headnames <- NULL
			prop <- NULL
			for (LOI2 in names(data.obj$abund.list)) {
				prop2 <- data.obj$abund.list[[LOI2]]
				taxa.name2 <- taxa.name[taxa.name %in% rownames(prop2)]
				if (length(taxa.name2) != 0) {
					if (scale == 'logP') {
						prop2 <- prop2 + pseudo.ct
					} 
					prop2 <- t(t(prop2) / colSums(prop2))			
					headnames2 <- sapply(strsplit(rownames(prop2), ";"), paste, collapse="\n")
					names(headnames2) <- rownames(prop2)
					prop <- rbind(prop, prop2[taxa.name2, , drop=FALSE])
					headnames <- c(headnames, headnames2)
				}
			}
			
		} else {
			cat(LOI, "\n")
			prop <- data.obj$abund.list[[LOI]]
			if (scale == 'logP') {
				prop <- prop + pseudo.ct
			} 
			prop <- t(t(prop) / colSums(prop))
			
			headnames <- sapply(strsplit(rownames(prop), ";"), paste, collapse="\n")
			names(headnames) <- rownames(prop)
			
			if (taxa.name == 'All') {
				prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
			} else {
				prop <- prop[taxa.name, , drop=FALSE]
			}
			
		}
		
		if (scale == 'logP') {
			prop <- log10(prop)
		}
		
		if (scale == 'sqrtP') {
			prop <- sqrt(prop)
		}
		
		if (scale == 'binary') {
			temp <- prop != 0
			prop[temp] <- 'Presence'
			prop[!temp] <- 'Absence'
		}
		if (is.null(hei0)) {
			hei <- 5
		} else {
			hei <- hei0
		}
		
		if (is.null(strata)) {
			if (is.null(wid0)) {
				wid <- 5
			} else {
				wid <- wid0
			}
			
		} else {
			if (is.null(wid0)) {
				wid <- 5.5
			} else {
				wid <- wid0
			}
		}
		
		if (scale == 'P') {
			ylab <- 'Proportion'
		} 
		if (scale == 'logP') {
			ylab <- 'log10(Proportion)'
		}
		if (scale == 'sqrtP') {
			ylab <- 'sqrt(Proportion)'
		}
		if (scale == 'binary') {
			ylab <- 'Count'
		} 
		
		pdf(paste("Taxa_Boxplot", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)
		
		if (is.null(strata)) {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				if (scale != 'binary') {
					if (scale == 'P') {
						if (rm.outlier == T) {
							ylims <- c(min(taxon.abund), quantile(taxon.abund, 0.95) * 1.25)
						} else {
							ylims <- range(taxon.abund)
						}
						
					} else {
						ylims <- range(taxon.abund)
					}
					
					if (is.null(subject)) {
						df2 <- data.frame(Value=taxon.abund, Group=grp)	
						dodge <- position_dodge(width=0.9)
						obj <- ggplot(df2, aes(x=Group, y=Value, col=Group)) +
								geom_boxplot(position=dodge, outlier.colour = NA) + 
								geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) +
								labs(y=ylab, title=headnames[taxon]) +
								ylim(ylims) +
								theme(legend.position="none")
						print(obj)
					} else {						
						df2 <- data.frame(Value=taxon.abund, Group=grp, subject=ID)	
						dodge <- position_dodge(width=0.9)
						obj <- ggplot(df2, aes(x=Group, y=Value, shape=Group, group=subject)) +
								geom_point(size=p.size) +
								geom_line(size=l.size) +
								labs(y=ylab, title=headnames[taxon]) +
							#	ylim(ylims) +
								theme(legend.position="none")
						print(obj)
						
					}

				} else {
					df2 <- data.frame(Value=taxon.abund, Group=grp)	
					
					obj <- ggplot(df2, aes(x=Group, fill=Value)) +
							geom_bar(width=.5) +
							labs(y=ylab, title=headnames[taxon]) +
							theme(legend.title=element_blank())
					
					print(obj)
				}
				
			}	
		} else {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				if (scale != 'binary') {
					if (scale == 'P') {
						if (rm.outlier == T) {
							ylims <- c(min(taxon.abund), quantile(taxon.abund, 0.95) * 1.25)
						} else {
							ylims <- range(taxon.abund)
						}
						
					} else {
						ylims <- range(taxon.abund)
					}
					grp2 <- df[, strata]
					df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2)
					dodge <- position_dodge(width=0.9)
					obj <- ggplot(df2, aes(x=Strata, y=Value, col=Group, fill=Group)) +
							geom_point(position=position_jitterdodge(dodge.width=0.9), size=3.0, alpha=0.6) +
							geom_boxplot(fill='white', alpha=0.6, position=dodge, outlier.colour = NA) + 
							labs(y=ylab, x=strata, title=headnames[taxon]) +
							ylim(ylims)
					print(obj)
				} else {
					grp2 <- df[, strata]
					df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2)
					obj <- ggplot(df2, aes(x=Group, fill=Value)) +
							geom_bar(width=.5) +
							labs(y=ylab, title=headnames[taxon]) +
							facet_wrap(~ Strata) + 
							theme(legend.title=element_blank())
					
					print(obj)
				}
			}
			
		}
		dev.off()
	}
	
}

# Rev: 2017_02_19 Add Pseudo.ct, hei0, wid0; strata2; smooth.method
generate_taxa_scatterplot <- function (data.obj,  grp.name, strata=NULL, strata2=NULL, scale='P',
		taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'), smooth.method='auto', wrap.nrow=2,
		pseudo.ct=0.5, taxa.name='All', rm.outlier=T, prev=0.1, minp=0.002, ann='All', hei0=NULL, wid0=NULL) {
	# To be completed
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	for (LOI in taxa.levels) {
		if (LOI == 'All') {
			if (taxa.name == 'All') stop('Taxa names are not required to be all when taxa level is also all!\n')
			headnames <- NULL
			prop <- NULL
			for (LOI2 in names(data.obj$abund.list)) {
				prop2 <- data.obj$abund.list[[LOI2]]
				taxa.name2 <- taxa.name[taxa.name %in% rownames(prop2)]
				if (length(taxa.name2) != 0) {
					if (scale == 'logP') {
						prop2 <- prop2 + pseudo.ct
					} 
					prop2 <- t(t(prop2) / colSums(prop2))			
					headnames2 <- sapply(strsplit(rownames(prop2), ";"), paste, collapse="\n")
					names(headnames2) <- rownames(prop2)
					prop <- rbind(prop, prop2[taxa.name2, , drop=FALSE])
					headnames <- c(headnames, headnames2)
				}
			}
			
		} else {
			cat(LOI, "\n")
			prop <- data.obj$abund.list[[LOI]]
			if (scale == 'logP') {
				prop <- prop + pseudo.ct
			} 
			prop <- t(t(prop) / colSums(prop))
			
			headnames <- sapply(strsplit(rownames(prop), ";"), paste, collapse="\n")
			names(headnames) <- rownames(prop)
			
			if (taxa.name == 'All') {
				prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
			} else {
				prop <- prop[taxa.name, , drop=FALSE]
			}
			
		}
		
		if (scale == 'logP') {
			prop <- log10(prop)
		}
		
		if (scale == 'sqrtP') {
			prop <- sqrt(prop)
		}
		
		if (scale == 'binary') {
			temp <- prop != 0
			prop[temp] <- 'Presence'
			prop[!temp] <- 'Absence'
		}
		
		if (is.null(hei0)) {
			hei <- 5
		} else {
			hei <- hei0
		}
		
		if (is.null(strata)) {
			if (is.null(wid0)) {
				wid <- 5
			} else {
				wid <- wid0
			}
			
		} else {
			if (is.null(wid0)) {
				wid <- 5.5
			} else {
				wid <- wid0
			}
		}
		
		
		if (scale == 'P') {
			ylab <- 'Proportion'
		} 
		if (scale == 'logP') {
			ylab <- 'log10(Proportion)'
		}
		if (scale == 'sqrtP') {
			ylab <- 'sqrt(Proportion)'
		}
		if (scale == 'binary') {
			ylab <- 'Count'
		} 
		
		pdf(paste("Taxa_Scatterplot", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)
		if (is.null(strata)) {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				if (scale != 'binary') {
					if (scale == 'P') {
						if (rm.outlier == T) {
							ylims <- c(min(taxon.abund), quantile(taxon.abund, 0.95) * 1.25)
						} else {
							ylims <- range(taxon.abund)
						}
						
					} else {
						ylims <- range(taxon.abund)
					}
					df2 <- data.frame(Value=taxon.abund, Group=grp)	
					dodge <- position_dodge(width=0.9)
					obj <- ggplot(df2, aes(x=Group, y=Value)) +
							geom_point() +
							geom_smooth(method=smooth.method) +
							labs(y=ylab, x=grp.name, title=headnames[taxon]) +
					#		ylim(ylims) +
							theme(legend.position="none")
					print(obj)
				} else {
					df2 <- data.frame(Value=taxon.abund, Group=grp)	
					
					obj <- ggplot(df2, aes(y=Group, x=Value)) +
							geom_boxplot() +
							geom_jitter() +
							labs(y=grp.name, x='Presence/Absence', title=headnames[taxon]) +
							theme(legend.title=element_blank())
					
					print(obj)
				}
				
			}	
		} else {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				if (scale != 'binary') {
					if (scale == 'P') {
						if (rm.outlier == T) {
							ylims <- c(min(taxon.abund), quantile(taxon.abund, 0.95) * 1.25)
						} else {
							ylims <- range(taxon.abund)
						}
						
					} else {
						ylims <- range(taxon.abund)
					}
					grp2 <- df[, strata]
					if (!is.null(strata2)) {
						grp3 <- df[, strata2]
					} else {
						grp3 <- grp2
					}
					
					# VERY bad names with only difference in capital letters
					df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2, Strata2=grp3)
					if (is.null(strata2)) {
						obj <- ggplot(df2, aes(x=Group, y=Value)) +
								geom_point() +
								geom_smooth(method=smooth.method) +
								labs(y=ylab, x=grp.name, title=headnames[taxon]) + 
								facet_wrap(~ Strata, nrow=wrap.nrow)
					} else {
						obj <- ggplot(df2, aes(x=Group, y=Value)) +
								geom_point() +
								geom_smooth(method=smooth.method) +
								labs(y=ylab, x=grp.name, title=headnames[taxon]) + 
								facet_grid(Strata2 ~ Strata)
					}		
					obj <- obj + theme(legend.title=element_blank())			
					print(obj)
					
				} else {				
					grp2 <- df[, strata]
					if (!is.null(strata2)) {
						grp3 <- df[, strata2]
					} else {
						grp3 <- grp2
					}
					df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2, Strata2=grp3)
					dodge <- position_dodge(width=0.9)

					if (is.null(strata2)) {
						obj <- ggplot(df2, aes(x=Strata, y=Group, col=Value, fill=Value)) +
								geom_point(position=position_jitterdodge(dodge.width=0.9), size=3.0, alpha=0.6) +
								geom_boxplot(fill='white', position=dodge, outlier.colour = NA) +
								labs(y=grp.name, x=strata, title=headnames[taxon]) 
					}  else {
						obj <- ggplot(df2, aes(x=Strata, y=Group, col=Value, fill=Value)) +
								geom_point(position=position_jitterdodge(dodge.width=0.9), size=3.0, alpha=0.6) +
								geom_boxplot(fill='white', position=dodge, outlier.colour = NA) +
								facet_wrap(~ Strata2, nrow=wrap.nrow) +
						        labs(y=grp.name, x=strata, title=headnames[taxon]) 
					}

					print(obj)
				}
			}
			
		}
		dev.off()
	}		
}

taxa_barplot_aggregate <- function (prop, df, grp.name, strata=NULL, scale='sqrt',  xsize=10, ylab='Proportion', error='se', cutoff=0.00001) {
	
	if (scale == 'log') {
		prop[prop <= cutoff] <- cutoff
	}
	
    grp <- factor(df[, grp.name])
	
	if (is.null(strata)) {
		df2 <- data.frame(Group=grp, t(prop))
		df2 <- melt(df2)
		colnames(df2) <- c('Group', 'Taxa', 'Value')
		
		# Could be revised
		temp1 <- aggregate(Value ~ Group + Taxa, df2, mean)

		if (error == 'se') {
			temp2 <- aggregate(Value ~ Group + Taxa, df2, function(x) sd(x) / sqrt(length(x)))
		}
		if (error == 'sd') {
			temp2 <- aggregate(Value ~ Group + Taxa, df2, function(x) sd(x))
		}
		
		
		df2 <- cbind(temp1, temp2[, 3])
		colnames(df2) <- c('Group', 'Taxa', 'Mean', 'SE')
		
		limits <- aes(ymax = Mean + SE, ymin = Mean - SE)
		dodge <- position_dodge(width=0.95)
		
		obj1 <- ggplot(df2, aes(x=Taxa, y=Mean, fill=Group)) + 
				geom_bar(position=dodge, stat="identity") + 
				geom_bar(position=dodge, stat="identity", colour="black", show_guide=FALSE, size=0.25) +
				geom_errorbar(limits, position=dodge, size=0.25, width=0.25) +
				labs(y=paste(ylab), x='') +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
				theme(legend.position="top", legend.title=element_blank())
		
		if (scale == 'sqrt') {
			# obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
			obj1 <- obj1 + scale_y_sqrt(
					breaks = trans_breaks("sqrt", function(x) x^2),
					labels = trans_format("sqrt", math_format(.x^2)))
		}
		
		# To be revised
		if (scale == 'log') {
			obj1 <- obj1 + scale_y_log10(
					breaks = trans_breaks("log10", function(x) 10^x),
					labels = trans_format("log10", math_format(10^.x)))
		}
		if (scale == 'boxcox') {
			obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
		}
		
	} else {
		grp2 <- factor(paste(strata, df[, strata]), levels=paste(strata, levels(df[, strata])))
		df2 <- data.frame(Group=grp, Strata=grp2, t(prop))
		df2 <- melt(df2)
		colnames(df2) <- c('Group', 'Strata', 'Taxa', 'Value')
		
		# Could be revised
		temp1 <- aggregate(Value ~ Group + Strata + Taxa, df2, mean)
		if (error == 'se') {
			temp2 <- aggregate(Value ~ Group + Strata + Taxa, df2, function(x) sd(x) / sqrt(length(x)))
		}
		if (error == 'sd') {
			temp2 <- aggregate(Value ~ Group + Strata + Taxa, df2, function(x) sd(x))
		}
		
		
		df2 <- cbind(temp1, temp2[, 4])
		colnames(df2) <- c('Group', 'Strata', 'Taxa', 'Mean', 'SE')
		
		limits <- aes(ymax = Mean + SE, ymin = Mean - SE)
		dodge <- position_dodge(width=0.95)
		
		obj1 <- ggplot(df2, aes(x=Taxa, y=Mean, fill=Group)) + 
				geom_bar(position=dodge, stat="identity") + 
				geom_bar(position=dodge, stat="identity", colour="black", show_guide=FALSE, size=0.25) +
				geom_errorbar(limits, position=dodge, size=0.25, width=0.25) +
				facet_wrap(~Strata, ncol=1) +
				labs(y=paste(ylab), x='') +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
				theme(legend.position="top", legend.title=element_blank())
		if (scale == 'sqrt') {
			# obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
			obj1 <- obj1 + scale_y_sqrt(
					breaks = trans_breaks("sqrt", function(x) x^2),
					labels = trans_format("sqrt", math_format(.x^2)))
		}
		if (scale == 'log') {
			obj1 <- obj1 + scale_y_log10(
					breaks = trans_breaks("log10", function(x) 10^x),
					labels = trans_format("log10", math_format(10^.x)))
		}
		if (scale == 'boxcox') {
			obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
		}
		
	}
	return(obj1)
	
}

taxa_boxplot_aggregate <- function (prop, df, grp.name, strata=NULL, scale='none', xsize=10, ylab='Proportion', cutoff=0.00001) {
	
	
	grp <- factor(df[, grp.name])
	
	if (scale == 'log') {
		prop[prop <= cutoff] <- cutoff
	}
	
	if (is.null(strata)) {
		df2 <- data.frame(Group=grp, t(prop))
		df2 <- melt(df2)
		colnames(df2) <- c('Group', 'Taxa', 'Value')
		

		dodge <- position_dodge(width=0.95)
		
		obj1 <- ggplot(df2, aes(x=Taxa, y=Value, fill=Group)) + 
				geom_boxplot(position=dodge) + 
				labs(y=paste(ylab), x='') +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
				theme(legend.position="top", legend.title=element_blank())
		if (scale == 'sqrt') {
			# obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
			obj1 <- obj1 + scale_y_sqrt(
					breaks = trans_breaks("sqrt", function(x) x^2),
					labels = trans_format("sqrt", math_format(.x^2)))
		}
		if (scale == 'log') {
			obj1 <- obj1 + scale_y_log10(
					breaks = trans_breaks("log10", function(x) 10^x),
					labels = trans_format("log10", math_format(10^.x)))
		}
		if (scale == 'boxcox') {
			obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
		}
		
	} else {
		grp2 <- factor(paste(strata, df[, strata]), levels=paste(strata, levels(df[, strata])))
		df2 <- data.frame(Group=grp, Strata=grp2, t(prop))
		df2 <- melt(df2)
		colnames(df2) <- c('Group', 'Strata', 'Taxa', 'Value')
		
		dodge <- position_dodge(width=0.95)
		
		obj1 <- ggplot(df2, aes(x=Taxa, y=Value, fill=Group)) + 
				geom_boxplot(position=dodge) + 
				facet_wrap(~Strata, ncol=1) +
				labs(y=paste(ylab), x='') +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
				theme(legend.position="top", legend.title=element_blank())
		if (scale == 'sqrt') {
			# obj1 <- obj1 + scale_y_sqrt(breaks=c(0, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0))
			obj1 <- obj1 + scale_y_sqrt(
					breaks = trans_breaks("sqrt", function(x) x^2),
					labels = trans_format("sqrt", math_format(.x^2)))
		}
		# To be revised
		if (scale == 'log') {
			obj1 <- obj1 + scale_y_log10(
					breaks = trans_breaks("log10", function(x) 10^x),
					labels = trans_format("log10", math_format(10^.x)))
		}
		if (scale == 'boxcox') {
			obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
		}
		
	}
	return(obj1)
	
}

generate_taxa_biplot <- function (data.obj, taxa, trans='sqrt', grp.name, ann='', ...) {
	
    grp <- data.obj$meta.dat[, grp.name]

	prop <- NULL
	for (LOI2 in names(data.obj$abund.list)) {
		ct <- data.obj$abund.list[[LOI2]]
		ct <- ct + 1
		prop0 <- t(t(ct) / colSums(ct))
		prop <- rbind(prop, prop0[intersect(rownames(prop0), taxa), , drop=FALSE])	
	}
	colnames(prop) <- colnames(prop0)
	if (nrow(prop) != length(taxa)) {
		warnings('Some taxa not found in abundance lists! Please check the names!\n')
	}
	
	if (trans == 'normal') 	prop <- t(apply(prop, 1, function(x) qqnorm(x, plot=F)$x))
	if (trans == 'log') prop <- log(prop)
	if (trans == 'sqrt') prop <- sqrt(prop)
	if (trans == 'rank') prop <- t(apply(prop, 1, rank))
	
	wine.pca <- prcomp(t(prop), scale. = TRUE)
	g <- ggbiplot(wine.pca, obs.scale = 1, var.scale = 1, 
			groups = grp, ellipse = TRUE, circle = FALSE, ...) 
	g <- g + scale_color_discrete(name = '') + theme_bw()
	g <- g + theme(legend.direction = 'horizontal', legend.position = 'top') 
	pdf(paste0("Taxa_Biplot_", ann, ".pdf"), height=6, width=6)
	print(g)
	dev.off()
}

# Rev: 2017_05_19: new formula for width
#  Ad presence and absence bar
generate_taxa_barplot_aggregate <- function (data.obj,  grp.name, strata=NULL, scale='sqrt', taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
		taxa.name='All',  prev=0.1, minp=0.002, ann='All', wids=NULL, heis=NULL, xsize=c(14, 12, 10, 9, 8)) {
	# To be completed
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	for (i in 1:length(taxa.levels)) {
		LOI <- taxa.levels[i]
		cat(LOI, "\n")
		prop <- data.obj$abund.list[[LOI]]
		prop <- t(t(prop) / colSums(prop))
				
		if (taxa.name == 'All') {
			prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
		} else {
			prop <- prop[taxa.name, , drop=FALSE]
		}
			
		# decreasing
		prop <- prop[rev(order(rowMeans(prop))), ]
		
		if (is.null(wids) | is.null(heis)) {
			wid <- 7 * ifelse(nrow(prop) / 30 < 1, 1, nrow(prop) / 30)
			wid <- sqrt(nlevels(grp) / 2) * wid
	        hei <- 7
		} else {
			wid <- wids[i]
			hei <- heis[i]
		}
			
		pdf(paste("Taxa_Barplot_Aggregate", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)		
		obj1 <- taxa_barplot_aggregate (prop, df, grp.name, strata, scale, xsize[i]) 
		print(obj1)
		dev.off()		
	}
}

# Rev: 2017_05_19: new formula for width
generate_taxa_boxplot_aggregate <- function (data.obj,  grp.name, strata=NULL, scale='sqrt', taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
		taxa.name='All',  prev=0.1, minp=0.002, ann='All', wids=NULL, heis=NULL, xsize=c(14, 12, 10, 9, 8)) {
	# To be completed
	df <- data.obj$meta.dat
	grp <- factor(df[, grp.name])
	
	for (i in 1:length(taxa.levels)) {
		LOI <- taxa.levels[i]
		cat(LOI, "\n")
		prop <- data.obj$abund.list[[LOI]]
		prop <- t(t(prop) / colSums(prop))
		
		if (taxa.name == 'All') {
			prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
		} else {
			prop <- prop[taxa.name, , drop=FALSE]
		}
		
		# decreasing
		prop <- prop[rev(order(rowMeans(prop))), ]
		
		if (is.null(wids) | is.null(heis)) {
			wid <- 7 * ifelse(nrow(prop) / 30 < 1, 1, nrow(prop) / 30)
			wid <- sqrt(nlevels(grp) / 2) * wid
		    
			hei <- 7
		} else {
			wid <- wids[i]
			hei <- heis[i]
		}
		
		pdf(paste("Taxa_Boxplot_Aggregate", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)		
		obj1 <- taxa_boxplot_aggregate (prop, df, grp.name, strata, scale, xsize[i]) 
		print(obj1)
		dev.off()		
	}
}

# Add combined barplot with error bar and presence/absence bar (currently presence/absence bar is in generate_taxa_boxplot
generate_taxa_barplot <- function (data.obj,  grp.name, strata=NULL, scale='P', taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
		taxa.name='All', rm.outlier=T, prev=0.1, minp=0.002, ann='All') {
	# To be completed
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	for (LOI in taxa.levels) {
		if (LOI == 'All') {
			if (taxa.name == 'All') stop('Taxa names are not required to be all when taxa level is also all!\n')
			headnames <- NULL
			prop <- NULL
			for (LOI2 in names(data.obj$abund.list)) {
				prop2 <- data.obj$abund.list[[LOI2]]
				taxa.name2 <- taxa.name[taxa.name %in% rownames(prop2)]
				if (length(taxa.name2) != 0) {
					if (scale == 'logP') {
						prop2 <- prop2 + 0.5
					} 
					prop2 <- t(t(prop2) / colSums(prop2))			
					headnames2 <- sapply(strsplit(rownames(prop2), ";"), paste, collapse="\n")
					names(headnames2) <- rownames(prop2)
					prop <- rbind(prop, prop2[taxa.name2, , drop=FALSE])
					headnames <- c(headnames, headnames2)
				}
			}
			
		} else {
			cat(LOI, "\n")
			prop <- data.obj$abund.list[[LOI]]
			if (scale == 'logP') {
				prop <- prop + 0.5
			} 
			prop <- t(t(prop) / colSums(prop))
			
			headnames <- sapply(strsplit(rownames(prop), ";"), paste, collapse="\n")
			names(headnames) <- rownames(prop)
			
			if (taxa.name == 'All') {
				prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
			} else {
				prop <- prop[taxa.name, , drop=FALSE]
			}
			
		}

		if (scale == 'logP') {
			prop <- log10(prop)
		}
		
		
		hei <- 5
		if (is.null(strata)) {
			wid <- 5
		} else {
			wid <- 4 * nlevels(df[, strata])
		}
		
		if (scale == 'P') {
			ylab <- 'Proportion'
		} else {
			ylab <- 'log10(Proportion)'
		}
		
		pdf(paste("Taxa_Barplot", LOI, scale, ann, ".pdf", sep="_"), height=hei, width=wid)
		if (is.null(strata)) {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				taxon.abund2 <- taxon.abund
				if (scale == 'P') {
					if (rm.outlier == T) {
						ylims <- c(0, quantile(taxon.abund, 0.95) * 1.25)
						taxon.abund[taxon.abund > quantile(taxon.abund, 0.95) * 1.25] <- quantile(taxon.abund, 0.95) * 1.25
					} else {
						ylims <- c(0, max(taxon.abund))
					}
					
				} else {
					ylims <- range(taxon.abund)
				}
				df2 <- data.frame(x=factor(1:length(taxon.abund)), Value=taxon.abund, Group=grp)			
				dodge <- position_dodge(width=0.99)
				obj <- ggplot(df2, aes(x=x, y=Value, fill=Group)) +
						geom_bar(position=dodge, stat='identity' ) + 
						facet_grid(. ~ Group, scales='free_x', space="free") +
						labs(y=ylab, title=headnames[taxon]) +
						ylim(ylims) +
						xlab('') +
						theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
						theme(legend.position="none")
				mean_df <- data.frame(Group=levels(grp), yint = tapply(taxon.abund2, grp, mean))
				obj <- obj + geom_hline(data = mean_df, aes(yintercept = yint))		
				mean_df <- data.frame(Group=levels(grp), yint = tapply(taxon.abund2, grp, median))
				obj <- obj + geom_hline(data = mean_df, aes(yintercept = yint), linetype=2)
				
				print(obj)
			}	
		} else {
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				taxon.abund2 <- taxon.abund
				if (scale == 'P') {
					if (rm.outlier == T) {
						ylims <- c(0, quantile(taxon.abund, 0.95) * 1.25)
						taxon.abund[taxon.abund > quantile(taxon.abund, 0.95) * 1.25] <- quantile(taxon.abund, 0.95) * 1.25
					} else {
						ylims <- c(0, max(taxon.abund))
					}
					
				} else {
					ylims <- range(taxon.abund)
				}
				
				grp2 <- df[, strata]
				obj.list <- list()
				for (level in levels(grp2)) {
					ind <- grp2 %in% level
					df2 <- data.frame(x=factor(1:length(taxon.abund[ind])), Value=taxon.abund[ind], Group=grp[ind])	
	
					dodge <- position_dodge(width=0.99)
					obj0 <- ggplot(df2, aes(x=x, y=Value, fill=Group)) +
							geom_bar(position=dodge, stat='identity' ) + 
							facet_grid(. ~ Group, scales='free_x', space="free") +
							labs(y=ylab, title=headnames[taxon]) +
							ylim(ylims) +
							xlab(paste(strata, level)) +
							theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
							theme(legend.position="none")
					mean_df <- data.frame(Group=levels(grp[ind]), yint = tapply(taxon.abund2[ind], grp[ind], mean))
					obj0 <- obj0 + geom_hline(data = mean_df, aes(yintercept = yint))		
					mean_df <- data.frame(Group=levels(grp[ind]), yint = tapply(taxon.abund2[ind], grp[ind], median))
					obj0 <- obj0 + geom_hline(data = mean_df, aes(yintercept = yint), linetype=2)
					obj.list[[level]] <- obj0		
				}
				multiplot(plotlist=obj.list, cols=nlevels(grp2))
			}	
			
		}
		dev.off()
	}		
}


twopart.test <- function(x1, x2, zero.p=0.2) {
	
	n1 <- length(x1)
	n2 <- length(x2)
	p1 <- mean(x1 != 0)
	p2 <- mean(x2 != 0)
	m1 <- sum(x1 != 0)
	m2 <- sum(x2 != 0)
	p12 <- (m1 + m2) / (n1 + n2)
	q12 <- 1 - p12
	
	if (q12 >= zero.p) {
		Z <- (abs(p1 - p2) - (1/(2*n1) + 1/(2*n2))) / sqrt(p12 * q12 * (1/n1 + 1/n2))
		x1 <- x1[x1!=0]
		x2 <- x2[x2!=0]
		R1 <- sum(rank(c(x1, x2))[1:length(x1)])
		ti <- as.vector(table(c(x1, x2)))
		W <- (abs(R1 - m1*(m1+m2+1)/2) - 1/2) / sqrt((m1*m2/12)*(m1+m2+1-sum(ti*(ti^2-1))/(m1+m2)/(m1+m2-1)))
		X2 <- Z^2 + W^2
		res <- list()
		res$stat <- X2
		res$p.value <- 1 - pchisq(X2, 2)
		res$Z <- Z
		res$W <- W
		res$test <- 'TwoPart'
	} else {
		res <- wilcox.test(x1, x2)
		res$test <- 'Wilcox'
	}
	res
}

getPermuteMatrix <- function (perm, N, strata = NULL) 
{
	if (length(perm) == 1) {
		perm <- how(nperm = perm)
	}
	if (!missing(strata) && !is.null(strata)) {
		if (inherits(perm, "how") && is.null(getBlocks(perm))) 
			setBlocks(perm) <- strata
	}
	if (inherits(perm, "how")) 
		perm <- shuffleSet(N, control = perm)
	if (is.null(attr(perm, "control"))) 
		attr(perm, "control") <- structure(list(within = list(type = "supplied matrix"), 
						nperm = nrow(perm)), class = "how")
	perm
}

perm_fdr_adj <- function (F0, Fp) {
	ord <- order(F0, decreasing = T)
	F0 <- F0[ord]
	perm.no <- ncol(Fp)
	Fp <- as.vector(Fp)
	Fp <- Fp[!is.na(Fp)]
	Fp <- sort(c(Fp, F0), decreasing = F)
	n <- length(Fp)
	m <- length(F0)
	FPN <- (n + 1) - match(F0, Fp) - 1:m
	p.adj.fdr <- FPN / perm.no / (1:m)
#		p.adj.fdr <- sapply(F0, function(x) sum(Fp >= 
#									x, na.rm=TRUE) / perm.no)/(1:length(F0))
	p.adj.fdr <- pmin(1, rev(cummin(rev(p.adj.fdr))))[order(ord)]
}

perm_fwer_adj <- function (F0, Fp) {
	ord <- order(F0, decreasing = T)
	F0 <- F0[ord]
	col.max <- colMaxs(Fp, na.rm=TRUE)
	p.adj.fwer <- sapply(F0, function(x) mean(col.max >= x))[order(ord)]
}

# Need to be further comprehensively tested
# a. Speed up, finished!
# b. Permutation method （response, covariate or residual permutation), residual will be default!
# c. Revise - add permutation stratified by subject, finished!
# Rev: 2017_02_02 permutation-based FDR control
# Rev: 2017_02_13 Add LMM-based permutation for block.perm=TRUE (type I error is controled, but power study hasn't been comprehensively studied)
# Rev: 2017_02_24 allow 'adj.name' to contain multiple covariates
# Still need to address NA's, currently simply remove NA's
permute_differential_analysis <- function (meta.dat, prop, grp.name, adj.name=NULL, strata=NULL, 
		block.perm=FALSE, sqrt.trans=TRUE, resid.perm=TRUE, perm.no=999) {
	# Square root transformation
	# User should take care of the normalization, transformation and addressing outliers
	if (sqrt.trans) {
		Y <- sqrt(prop)
	} else {
		Y <- prop
	}
	row.names <- rownames(Y)
	
	if (!is.null(strata)) {
		strata <- factor(strata)
	}
	
	# Prepare model matrix
	n <- ncol(prop)
	I <- diag(n)
	if (is.null(adj.name)) {
		M0 <- model.matrix(~ 1, meta.dat)
	} else {
		df0 <- meta.dat[, c(adj.name), drop=F]
		M0 <- model.matrix( ~., df0)
	}

#	P0 <- M0 %*% solve(t(M0) %*% M0) %*% t(M0)
	
	df1 <- meta.dat[, c(adj.name, grp.name), drop=F]
	M1 <- model.matrix( ~., df1)
	
	# QR decompostion
	qrX0 <- qr(M0, tol = 1e-07)
	Q0 <- qr.Q(qrX0)
	Q0 <- Q0[, 1:qrX0$rank, drop=FALSE]
	
	qrX1 <- qr(M1, tol = 1e-07)
	Q1 <- qr.Q(qrX1)
	Q1 <- Q1[, 1:qrX1$rank, drop=FALSE]
	
	# Got residual
	if (resid.perm) {
		if (!block.perm) {
			# Permute the residual
			if (is.null(adj.name)) {
				Y <- t(resid(lm(as.formula(paste('t(Y) ~ 1')), meta.dat)))
			} else {
				Y <- t(resid(lm(as.formula(paste('t(Y) ~ ', paste(adj.name, collapse='+'))), meta.dat)))
			}
			
		} else {
			
			if (is.null(strata)) {
				stop('Block permutation requires strata!\n')
			} else {
				Y.r <- matrix(NA, nrow(Y), nlevels(strata))
				Y.e <- Y
				cat('Fitting linear mixed effects model ...\n')
				for (j in 1:nrow(Y)) {
					# Linear mixed effects model
					yy <- Y[j, ]
					meta.dat$yy <- yy
					meta.dat$strata <- strata
					if (is.null(adj.name)) {
						obj <- lme(as.formula(paste('yy ~ 1')), random =~ 1 |  strata, data=meta.dat, method='ML')
						# The order is the same as the levels
						Y.r[j, ] <- random.effects(obj)[, 1]
						Y.e[j, ] <- resid(obj)
					} else {
						obj <- lme(as.formula(paste('yy ~ ', paste(adj.name, collapse='+'))), random =~ 1 |  strata, data=meta.dat, method='ML')
						Y.r[j, ] <- random.effects(obj)[, 1]
						Y.e[j, ] <- resid(obj)
					}
				}
				Y <- Y.r[, as.numeric(strata)] + Y.e
	#			Y <- Y - rowMeans(Y)
	#		    Y <- Y.e
			}
		}
	}
	

	
	TSS <- rowSums(Y^2)
	MSS1 <- rowSums((Y %*% Q1)^2)
	MSS0 <- rowSums((Y %*% Q0)^2)  # Not necessary, it's zero
	F0 <- (MSS1 - MSS0) /  (TSS - MSS1) 
	
#	P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
#	F0 <- diag(Y %*% (P1 - P0) %*% t(Y)) / diag(Y %*% (I - P1) %*% t(Y))
#	df3 <- df1
	
	if (block.perm == FALSE) {
		perm.ind <- vegan:::getPermuteMatrix(perm.no, n, strata = strata)
		perm.no <- nrow(perm.ind)
	}
	
	cat('Permutation test ....\n')
	Fp <- sapply(1:perm.no, function(i) {
				if (i %% 100 == 0) cat('.')
				if (block.perm == FALSE) {
					Yp <- Y[, perm.ind[i, ]]
				} else {
					# Double permutation
					strata.p <- factor(strata, levels=sample(levels(strata)))
				    Yp <- Y.r[, as.numeric(strata.p)] + Y.e[, sample(ncol(Y))]
#				    Yp <- Y.e[, sample(ncol(Y))]
#					Yp <- Yp - rowMeans(Yp)
				}
				
#				df3[, grp.name] <- sample(df1[, grp.name])
#				M1 <- model.matrix( ~., df3)
#				P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
				MSS1p <- rowSums((Yp %*% Q1)^2)
				MSS0p <- rowSums((Yp %*% Q0)^2)    
			    if (block.perm == FALSE) {
					TSSp <- TSS
				} else {
					TSSp <- rowSums(Yp^2)
				}
				(MSS1p - MSS0p) /  (TSSp - MSS1p) 
			})
	
	
	if (mean(is.na(F0)) >= 0.1) {
		warning('More than 10% observed F stats have NA! Please check! \n')
	}
	
	if (mean(is.na(Fp)) >= 0.1) {
		warning('More than 10% permuted F stats have NA! Please check! \n')
	}
	
	na.ind <- is.na(F0)
	F0 <- F0[!na.ind]
	Fp <- Fp[!na.ind, ]
	
	p.raw <- cbind(Fp >= F0, 1)
	p.raw <- rowMeans(p.raw)
#	p.raw[is.na(p.raw)] <- 1   
	
	p.adj.fdr <- perm_fdr_adj(F0, Fp)
	p.adj.fwer <- perm_fwer_adj(F0, Fp)
	
	# Pad back the NA values
	pad <- function (vec, ind) {
		vec0 <- numeric(length(ind))
		vec0[!ind] <- vec
		vec0[ind] <- NA
		vec0
	}
	
	F0 <- pad(F0, na.ind)
	p.raw <- pad(p.raw, na.ind)
	p.adj.fdr <- pad(p.adj.fdr, na.ind)
	p.adj.fwer <- pad(p.adj.fwer, na.ind)
	
	names(F0) <- names(p.raw) <- names(p.adj.fdr) <- names(p.adj.fwer) <- row.names
	return(list(F.stat=F0, p.raw=p.raw, p.adj.fdr=p.adj.fdr, p.adj.fwer=p.adj.fwer))
#	return(p.raw)
}
#


# New: 2017_08_17 Add a new variant of PERMANOVA with matrix decomposition
# Partially validated
permanova2 <- function (meta.dat, D, grp.name, adj.name=NULL, strata=NULL, 
		block.perm=FALSE, resid.perm=TRUE, perm.no=999, eig = c('All', 'Positive')) {
	
	eig <- match.arg(eig)
	# Square root transformation
	# User should take care of the normalization, transformation and addressing outliers
	
	D <- as.matrix(D)^2
	n <- nrow(D)
	G <- mean(D) + D - rowMeans(D) - matrix(rep(1, n), ncol=1) %*% colMeans(D)
	
	eig.obj <- eigen(G, symm=TRUE)
	
	if (eig == 'All') {
		Y <- eig.obj$vectors
		lambda <- eig.obj$values
	}
	
	if (eig == 'Positive') {
		lambda <- eig.obj$values
		ind <- lambda > 1e-6
		lambda <- lambda[ind]
		Y <- Y[, ind, drop = FALSE]
	}

    Y <- t(Y)
	
	if (!is.null(strata)) {
		strata <- factor(strata)
	}
	
	# Prepare model matrix
	n <- ncol(Y)
	I <- diag(n)
	if (is.null(adj.name)) {
		M0 <- model.matrix(~ 1, meta.dat)
	} else {
		df0 <- meta.dat[, c(adj.name), drop=F]
		M0 <- model.matrix( ~., df0)
	}
	
#	P0 <- M0 %*% solve(t(M0) %*% M0) %*% t(M0)
	
	df1 <- meta.dat[, c(adj.name, grp.name), drop=F]
	M1 <- model.matrix( ~., df1)
	
	# QR decompostion
	qrX0 <- qr(M0, tol = 1e-07)
	Q0 <- qr.Q(qrX0)
	Q0 <- Q0[, 1:qrX0$rank, drop=FALSE]
	
	qrX1 <- qr(M1, tol = 1e-07)
	Q1 <- qr.Q(qrX1)
	Q1 <- Q1[, 1:qrX1$rank, drop=FALSE]
	
	# Got residual
	if (resid.perm) {
		if (!block.perm) {
			# Permute the residual
			if (is.null(adj.name)) {
				Y <- t(resid(lm(as.formula(paste('t(Y) ~ 1')), meta.dat)))
			} else {
				Y <- t(resid(lm(as.formula(paste('t(Y) ~ ', paste(adj.name, collapse='+'))), meta.dat)))
			}
			
		} else {
			
			if (is.null(strata)) {
				stop('Block permutation requires strata!\n')
			} else {
				Y.r <- matrix(NA, nrow(Y), nlevels(strata))
				Y.e <- Y
#				cat('Fitting linear mixed effects model ...\n')
				for (j in 1:nrow(Y)) {
					# Linear mixed effects model
					yy <- Y[j, ]
					meta.dat$yy <- yy
					meta.dat$strata <- strata
					if (is.null(adj.name)) {
						obj <- lme(as.formula(paste('yy ~ 1')), random =~ 1 |  strata, data=meta.dat, method='ML')
						# The order is the same as the levels
						Y.r[j, ] <- random.effects(obj)[, 1]
						Y.e[j, ] <- resid(obj)
					} else {
						obj <- lme(as.formula(paste('yy ~ ', paste(adj.name, collapse='+'))), random =~ 1 |  strata, data=meta.dat, method='ML')
						Y.r[j, ] <- random.effects(obj)[, 1]
						Y.e[j, ] <- resid(obj)
					}
				}
				Y <- Y.r[, as.numeric(strata)] + Y.e
				#			Y <- Y - rowMeans(Y)
				#		    Y <- Y.e
			}
		}
	}
	
	
	
	TSS <- rowSums(Y^2)
	MSS1 <- rowSums((Y %*% Q1)^2)
	MSS0 <- rowSums((Y %*% Q0)^2)  # Not necessary, it's zero
	F0 <- sum(lambda * (MSS1 - MSS0)) /  sum(lambda * (TSS - MSS1))
	
#	P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
#	F0 <- diag(Y %*% (P1 - P0) %*% t(Y)) / diag(Y %*% (I - P1) %*% t(Y))
#	df3 <- df1
	
	if (block.perm == FALSE) {
		perm.ind <- vegan:::getPermuteMatrix(perm.no, n, strata = strata)
		perm.no <- nrow(perm.ind)
	}
	
#	cat('Permutation test ....\n')
	Fp <- sapply(1:perm.no, function(i) {
#				if (i %% 100 == 0) cat('.')
				if (block.perm == FALSE) {
					Yp <- Y[, perm.ind[i, ]]
				} else {
					# Double permutation
					strata.p <- factor(strata, levels=sample(levels(strata)))
					Yp <- Y.r[, as.numeric(strata.p)] + Y.e[, sample(ncol(Y))]
#				    Yp <- Y.e[, sample(ncol(Y))]
#					Yp <- Yp - rowMeans(Yp)
				}
				
#				df3[, grp.name] <- sample(df1[, grp.name])
#				M1 <- model.matrix( ~., df3)
#				P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
				MSS1p <- rowSums((Yp %*% Q1)^2)
				MSS0p <- rowSums((Yp %*% Q0)^2)    
				if (block.perm == FALSE) {
					TSSp <- TSS
				} else {
					TSSp <- rowSums(Yp^2)
				}
				sum(lambda * (MSS1p - MSS0p)) / sum(lambda * (TSSp - MSS1p)) 
			})
	
	p.value <- mean(c(Fp >= F0, 1))
  
	return(list(f0 = F0, f.perms = Fp, p.value = p.value))
#	return(p.raw)
}




# This function for nonparametric/permutaiton method
# Rev: 2017_02_16  Add normalization method; Add transformation; Remove rarefaction (only output warnings);
# Rev: 2017_10_30  Support filtering based on coefficient of variation
perform_differential_analysis <- function (data.obj, grp.name, adj.name=NULL, subject=NULL, 
		taxa.levels=c('Phylum', 'Order', 'Class', 'Family', 'Genus', 'Species'),
		method='perm', block.perm=FALSE, perm.no=999,
		norm='GMPR', norm.level='Species', intersect.no=4,
		transform='sqrt',
		prev=0.1, minp=0.002, medianp=NULL, cv=NULL,
		mt.method='fdr', cutoff=0.15, 
		ann='', seed=123, ...) {
	# To be completed
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	ind <- !is.na(grp)
	data.obj <- subset_data(data.obj, ind)
	grp <- grp[ind]
	df <- df[ind, ]
	
	if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
		data.obj$abund.list[['Species']] <- data.obj$otu.tab
		rownames(data.obj$abund.list[['Species']]) <- paste0("OTU", rownames(data.obj$otu.tab), ":", 
				data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
	}

	# Test for sequence-depth confounding
	dep <- colSums(data.obj$otu.tab)
	diff.seq.p <- summary(aov(dep ~ grp))[[1]][1, 'Pr(>F)']
	if (!is.na(diff.seq.p) & diff.seq.p <= 0.05) {
		warning(paste0(
		 '\nSignificant sequencing depth confounding with the variable of interest!\n',
		 'For nonparametric test/permutaiton test, there may be potentially many false postives (for those less prevalent taxa)!\n',
		 'Consider performing rarefaction first! (However, rarefaction will not completly solve the problem.)\n',
		 'May also try count-based models, which might have better false postive controls!\n'))
#		cat("For parametric test with sequence depth adjustment (DESeq2), please be cautious about the results!\n There may be potential residual sequence depth confounding!\n")
		# potential revising
		
#		# Causes prolbem for visualization (the names changed?)
#		otu.tab <- t(Rarefy(t(data.obj$otu.tab))$otu.tab.rff)
#		data.obj <- subset_data(data.obj, colnames(otu.tab))
#		data.obj$otu.tab <- otu.tab
#		otu.name <- data.obj$otu.name
#		# Re-summarize the taxa data
#		for (hierach in taxa.levels) {
#			if (hierach != 'Species') {
#				if (hierach != 'Phylum') {
#					tax.family <- paste(otu.name[, 'Phylum'], otu.name[, hierach], sep=";")
#					tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- paste0('Unclassified_', hierach)
#				} else {
#					tax.family <- otu.name[, 'Phylum']
#					tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- 'Unclassified_Phylum'
#				}
#				family <- aggregate(otu.tab, by=list(tax.family), FUN=sum)
#				rownames(family) <- family[, 1]
#				family <- as.matrix(family[, -1])
#				data.obj$abund.list[[hierach]] <- family
#			} else {
#				data.obj$abund.list[['Species']][, ] <- data.obj$otu.tab
#			}
#		}
	}
  
	# Calculate size.factor
	if (norm == 'Precalculated') {
		size.factor <- data.obj$size.factor
	}
	if (norm == 'GMPR') {
		if (norm.level %in% c('OTU', 'Species')) {
			tab <- data.obj$otu.tab
		} else {
			tab <- data.obj$abund.list[[norm.level]]
		}
		size.factor <- GMPR(tab, intersect.no)
	}	
	if (norm == 'TSS') {
		size.factor <- colSums(data.obj$otu.tab)
	}
	
	# Method-dependent processing
	if (is.null(method)) {
		if (nlevels(grp) == 2) {
			method <- 'wilcox'
		} else {
			method <- 'kruskal'
		}
	}
	if (method == 'wilcox.pair') {
		if (nlevels(grp) != 2) stop("Wilcox test requires two groups!\n")
		if (is.null(subject)) stop("Paired wilcox needs subject information!\n")
		subject <- factor(df[, subject])
		ind1 <- ind2 <- NULL
		for(sub in levels(subject)) {
			temp1 <- which(as.numeric(grp) == 1 & subject == sub)
			temp2 <- which(as.numeric(grp) == 2 & subject == sub)
			if (length(temp1) != 0 & length(temp2) != 0) {
				ind1 <- c(ind1, temp1[1])
				ind2 <- c(ind2, temp2[1])
			}
			
		}
		#	if (length(ind1) != length(ind2))  warning('Some subjects are not paired!\n')
	}
	
	if (method == 'perm.pair') {
		if (is.null(subject)) stop("Paired permutation test needs subject information!\n")
		subject <- factor(df[, subject])
	}
	
	if (method == 'perm') {
		if (!is.null(subject)) {
			subject <- factor(df[, subject])
		} 
	}
	
	if (is.factor(grp)) {
		
		pv.list <- qv.list <-  fc.list <- pc.list <- m.list <- nzm.list  <- prv.list <- list()
		res.final <- NULL
		for (LOI in taxa.levels) {
			cat(LOI, "\n")
			ct <- data.obj$abund.list[[LOI]]
			
			# Filtering
			prop0 <- t(t(ct) / colSums(ct))
			if (!is.null(prev)) {
				prop0 <- prop0[rowSums(prop0!=0) > prev * ncol(prop0), , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			if (!is.null(minp)) {
				prop0 <- prop0[rowMaxs(prop0) > minp, , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			
			if (!is.null(medianp)) {
				nz.mean <- apply(prop0, 1, function(x) median(x[x!=0]))
				prop0 <- prop0[nz.mean > medianp, , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			
			if (!is.null(cv)) {
				prop0 <- prop0[rowSds(prop0) / rowMeans(prop0) > cv, , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			
			# Normalization
	        prop <- t(t(ct) / size.factor)
			
			# Transformation - Others are possible/may be explored in the future
	        if (transform == 'sqrt') {
				prop <- sqrt(prop)
			}
			
#			if (method == 'DESeq2') {
#				if (is.null(adj.name)) {
#					formula <- as.formula(paste('~', grp.name))
#				} else {
#					formula <- as.formula(paste('~', paste(paste(adj.name, collapse='+'), '+', grp.name)))
#				}
#				dds <- DESeqDataSetFromMatrix(countData = ct + 1,
#						colData = df,
#						design = as.formula(formula))
#				dds <- DESeq(dds)
#				res <- results(dds)
#				pv.de2 <- res[, 'pvalue']
#				names(pv.de2) <- rownames(res)
#			}
			 		
			if (method == 'perm') {
				set.seed(seed)
				pv.de2 <- permute_differential_analysis(df, prop, grp.name, adj.name, strata=subject, block.perm=block.perm, sqrt.trans=FALSE, perm.no=perm.no)$p.raw
				names(pv.de2) <- rownames(prop)
			}
			
			# For legacy use
			if (method == 'perm.pair') {
				set.seed(seed)
				pv.de2 <- permute_differential_analysis(df, prop, grp.name, adj.name, strata=subject, block.perm=FALSE, sqrt.trans=FALSE, perm.no=perm.no)$p.raw
				names(pv.de2) <- rownames(prop)
			}
			
			pv.vec <- m.vec <- nzm.vec <- prv.vec <-  fc.vec <- pc.vec <-  NULL
			
			for (taxon in rownames(prop)) {
				taxon.abund <- prop[taxon, ]
				taxon.abund0 <- prop0[taxon, ]
				
				pv <- fc <- pc <- m <- nzm <- prv <- NULL
				if (method == 'wilcox') {
					if (nlevels(grp) != 2) stop("Wilcox test requires two groups!\n")
					pv <- wilcox.test(taxon.abund ~ grp)$p.value
					
				}
				if (method == 'wilcox.pair') {
					pv <- wilcox.test(taxon.abund[ind1], taxon.abund[ind2], paired=T)$p.value
				}
				if (method == 'twopart') {
					if (nlevels(grp) != 2) stop("Two part test requires two groups!\n")
					grp1 <- taxon.abund[as.numeric(grp)==1]
					grp2 <- taxon.abund[as.numeric(grp)==2]
					pv <- twopart.test(grp1, grp2)$p.value
				}	
				if (method == 'kruskal') {
					if (nlevels(grp) <= 2) warning("Kruskal-wallis test requires three or more groups!\n")
					pv <- kruskal.test(taxon.abund ~ grp)$p.value
				}
#				if (method == 'DESeq2') {
#					pv <- pv.de2[taxon]
#				}
				
				if (method == 'perm') {
					pv <- pv.de2[taxon]
				}
				
				if (method == 'perm.pair') {
					pv <- pv.de2[taxon]
				}
				m <- tapply(taxon.abund0, grp, function(x) mean(x))
				nzm <- tapply(taxon.abund0, grp, function(x) mean(x[x != 0]))
				prv <- tapply(taxon.abund0, grp, function(x) sum(x != 0))				
				
				# Rev: 2017_02_21 fc change baseline grp
				if (nlevels(grp) == 2) {
					grp.no <- table(grp)
					fc <- log2(m[2] / m[1])
					pc <- prv[2] / grp.no[2] / prv[1] * grp.no[1]
				} else {
					pc <- fc <- NA
				}
				
				pv.vec <- rbind(pv.vec, pv)
				fc.vec <- rbind(fc.vec, fc)
				m.vec <- rbind(m.vec, m)
				nzm.vec <- rbind(nzm.vec, nzm)
				pc.vec <- rbind(pc.vec, pc)
				prv.vec <- rbind(prv.vec, prv / table(grp))	
			}
			
			
			temp <- p.adjust(pv.vec[, 1], 'fdr')
			
			qv.vec <- matrix(temp, ncol=1)
			
			rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(pc.vec) <- rownames(m.vec) <- rownames(nzm.vec) <- rownames(prv.vec) <- rownames(prop)
			colnames(pv.vec) <- 'Pvalue'
			colnames(qv.vec) <- 'Qvalue'
			colnames(fc.vec) <- 'logFoldChange'
			colnames(pc.vec) <- 'PrevalChange'
			colnames(m.vec) <- paste(levels(grp), 'Mean')
			colnames(nzm.vec) <- paste(levels(grp), 'nzMean')
			colnames(prv.vec) <- paste(levels(grp), 'preval')
			
			pv.list[[LOI]] <- pv.vec
			qv.list[[LOI]] <- qv.vec
			fc.list[[LOI]] <- fc.vec
			pc.list[[LOI]] <- pc.vec
			m.list[[LOI]] <- m.vec
			nzm.list[[LOI]] <- nzm.vec
			prv.list[[LOI]] <- prv.vec
			
			res <- cbind(pv.vec, qv.vec, m.vec, nzm.vec, fc.vec, prv.vec, pc.vec)
			rownames(res) <- rownames(prop)
			write.csv(res, paste0("Taxa_DifferentialAbundanceAnalysis_", LOI, "_", ann, ".csv"))
			
			if (mt.method == 'fdr') {
				res.final <- rbind(res.final, res[res[, 'Qvalue'] <= cutoff, , drop=F])
			}
			if (mt.method == 'raw') {
				res.final <- rbind(res.final, res[res[, 'Pvalue'] <= cutoff, , drop=F])
			}
			
		}
		
		if (!is.null(res.final)) {
			colnames(res.final) <- colnames(res)
			write.csv(res.final, paste0("Taxa_DifferentialAbundanceAnalysis_AllLevels_", mt.method, '_', cutoff, "_", ann, ".csv"))
		}
		return(list(pv.list=pv.list, fc.list=fc.list, pc.list=pc.list, qv.list=qv.list, m.list=m.list))
	} else {
		if (is.null(method)) {
			method <- 'Spearman'
		}
		# Continuous case - currently only has DESeq2 
		pv.list <- qv.list <-  fc.list <-  m.list <- list()
		res.final <- NULL
		for (LOI in taxa.levels) {
			cat(LOI, "\n")
			ct <- data.obj$abund.list[[LOI]]
			
			# Filtering
			prop0 <- t(t(ct) / colSums(ct))
			if (!is.null(prev)) {
				prop0 <- prop0[rowSums(prop0!=0) > prev * ncol(prop0), , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			if (!is.null(minp)) {
				prop0 <- prop0[rowMaxs(prop0) > minp, , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			
			if (!is.null(medianp)) {
				nz.mean <- apply(prop0, 1, function(x) median(x[x!=0]))
				prop0 <- prop0[nz.mean > medianp, , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			
			if (!is.null(cv)) {
				prop0 <- prop0[rowSds(prop0) / rowMeans(prop0) > cv, , drop=FALSE]	
				ct <- ct[rownames(prop0), , drop=FALSE]
			}
			
			# Normalization
			prop <- t(t(ct) / size.factor)
			
			# Transformation - Others are possible/may be explored in the future
			if (transform == 'sqrt') {
				prop <- sqrt(prop)
			}
			
#			if (method == 'DESeq2') {
#				if (is.null(adj.name)) {
#					formula <- as.formula(paste('~', grp.name))
#				} else {
#					formula <- as.formula(paste('~', paste(paste(adj.name, collapse='+'), '+', grp.name)))
#				}
#				dds <- DESeqDataSetFromMatrix(countData = ct + 1,
#						colData = df,
#						design = as.formula(formula))
#				dds <- DESeq(dds)
#				res <- results(dds)
#				pv.de2 <- res[, 'pvalue']
##				fc.de2 <- res[, 'log2FoldChange']
#			# May need to revise
#				fc.de2 <- sapply(1:nrow(prop), function(i) cor(prop[i, ], df[, grp.name], method='spearman'))
#				names(pv.de2) <- rownames(res)
#			}
			
			if (method == 'perm') {
				set.seed(seed)
				pv.de2 <- permute_differential_analysis(df, prop, grp.name, adj.name, strata=subject, block.perm=block.perm, sqrt.trans=FALSE, perm.no=perm.no)$p.raw
				names(pv.de2) <- rownames(prop)
				# Place holder
				fc.de2 <- sapply(1:nrow(prop), function(i) cor(prop[i, ], df[, grp.name], method='spearman'))
			}
			
			if (method == 'Spearman') {
				if (!is.null(adj.name)) {
					stop("Spearman test can't adjust covariates!")
				}
				pv.de2 <- apply(prop, 1, function(x) {
							cor.test(x, grp, method='spearman')$p.value
						})
				names(pv.de2) <- rownames(prop)
				# Place holder
				fc.de2 <- sapply(1:nrow(prop), function(i) cor(prop[i, ], df[, grp.name], method='spearman'))
			}
			pv.vec <- matrix(pv.de2, ncol=1)	
			qv.vec <- matrix(p.adjust(pv.vec[, 1], 'fdr'), ncol=1)
			fc.vec <- matrix(fc.de2, ncol=1)
			m.vec <- matrix(rowMeans(prop0), ncol=1)
			
			rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(m.vec)  <- rownames(prop)
			colnames(pv.vec) <- 'Pvalue'
			colnames(qv.vec) <- 'Qvalue'
			colnames(fc.vec) <- 'SpearmanCorr'
			colnames(m.vec) <- 'Mean'
			
			
			pv.list[[LOI]] <- pv.vec
			qv.list[[LOI]] <- qv.vec
			fc.list[[LOI]] <- fc.vec
			m.list[[LOI]] <- m.vec
			
			
			res <- cbind(m.vec, fc.vec, pv.vec, qv.vec)
			rownames(res) <- rownames(prop)
			write.csv(res, paste0("Taxa_DifferentialAbundanceAnalysis_", LOI, "_", ann, ".csv"))
			
			if (mt.method == 'fdr') {
				res.final <- rbind(res.final, res[res[, 'Qvalue'] <= cutoff, , drop=F])
			}
			if (mt.method == 'raw') {
				res.final <- rbind(res.final, res[res[, 'Pvalue'] <= cutoff, , drop=F])
			}
		}
		
		if (!is.null(res.final)) {
			colnames(res.final) <- colnames(res)
			write.csv(res.final, paste0("Taxa_DifferentialAbundanceAnalysis_AllLevels_", mt.method, '_', cutoff, "_", ann, ".csv"))
		}
		return(list(pv.list=pv.list, fc.list=fc.list, qv.list=qv.list, m.list=m.list))
	}
	
}

#Owls <- transform(Owls, 
#		Nest=reorder(Nest,NegPerChick), 
#		logBroodSize=log(BroodSize), 
#		NCalls=SiblingNegotiation)
#
#fit_zinbinom <- glmmadmb(NCalls~(FoodTreatment+ArrivalTime)*SexParent+ 
#				offset(logBroodSize)+(1|Nest), 
#		data=Owls, 
#		zeroInflation=TRUE, 
#		family="nbinom")
confint.lme <- function (object, parm, level = 0.95, ...) {
	cf <- fixed.effects(object)
	pnames <- names(cf)
	if (missing(parm)) 
		parm <- pnames
	else if (is.numeric(parm)) 
		parm <- pnames[parm]
	a <- (1 - level)/2
	a <- c(a, 1 - a)
	pct <- stats:::format.perc(a, 3)
	fac <- qnorm(a)
	ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, 
					pct))
	ses <- sqrt(diag(vcov(object)))[parm]
	ci[] <- cf[parm] + ses %o% fac
	ci
}


# Rev: 2016_09_13 add glmmPQL-based overdipersed Poisson and binomial model
# glmmPQL default has overdispersion parameter. Regular Poisson and Binomial does not work.
perform_differential_analysis_para_single_RE <- function (taxon, ldep, grp.name, adj.name=NULL, subject, df, method='NB', LRT=FALSE) {
	# ldep: log depth (size factor)
	if (!is.null(adj.name) ) {
		if (sum(grepl(grp.name, c(adj.name)))) {
			stop('grp.name could not be part of adj.name, or there will be problem!\n')
		}
	}
	
	if (is.null(subject)) {
		stop('Random effects model require the specification of the subject parameter!\n')
	}	
	
	df$ldep <- ldep
	df$taxon <- taxon
	if (LRT) warning('Currently, only Wald test is implemented!\n')
	
	
	if (method %in% c('NB', 'ZINB1', 'ZINB0', 'B0')) {
		if (is.null(adj.name)) {
			grp.name.adj.name.subject <- grp.name.adj.name.subject <- paste(grp.name,  '+ (1|', subject, ')')
		} else {
			
			grp.name.adj.name.subject <- paste(grp.name, '+', adj.name, '+ (1|', subject, ')')
		}
		
		if (method == 'NB') {
			m1.nb <- glmmadmb(as.formula(paste('taxon ~', grp.name.adj.name.subject, '+ offset(ldep)')), data = df, 
					zeroInflation=FALSE, family='nbinom')
		}
		
		# 'B0' uses the glmmadmb, which will be compared to 'B' method of glmmPQL
		if (method == 'B0') {
			taxon2 <- as.numeric(taxon != 0)
			df$taxon2 <- taxon2
			m1.nb <- glmmadmb(as.formula(paste('taxon2 ~', grp.name.adj.name.subject, '+ ldep')), data = df, 
					zeroInflation=FALSE, family='binomial')
		}
		
		if (method == 'ZINB1') {
			m1.nb <- glmmadmb(as.formula(paste('taxon ~', grp.name.adj.name.subject, '+ offset(ldep)')), data = df, 
					zeroInflation=TRUE, family='nbinom')
		}
		
		if (method == 'ZINB0') {
			m1.nb <- glmmadmb(as.formula(paste('taxon ~', grp.name.adj.name.subject, '+ offset(ldep)')), data = df, 
					zeroInflation=TRUE, family='truncnbinom')
		}
		
		code <- list(m1.conv=m1.nb$convmsg)	
		pv.nb <- wald.test(b = coef(m1.nb), Sigma = vcov(m1.nb), Terms = grep(grp.name, names(coef(m1.nb))))$result$chi2['P']
		method <- paste(method, 'Wald')
		
		
		aic.nb <- -2 * m1.nb$loglik + 2 * m1.nb$npar	
		coef.nb <- coef(m1.nb)
		ci.nb <- confint.default(m1.nb)
	} 
	
	if (method %in% c('OP', 'B', 'QB')) {
		if (is.null(adj.name)) {
			grp.name.adj.name <- grp.name
		} else {
			grp.name.adj.name <- paste(grp.name, '+', adj.name)
		}
		
		subject <- paste('~ 1 |', subject)
		if (method == 'OP') {
			m1.nb <- glmmPQL(as.formula(paste('taxon ~', grp.name.adj.name, '+ offset(ldep)')), random = as.formula(subject), data = df, 
					verbose=F, family=quasipoisson())
		}
		
		# 'B0' uses the glmmadmb, which will be compared to 'B' method of glmmPQL
		if (method == 'B') {
			taxon2 <- as.numeric(taxon != 0)
			df$taxon2 <- taxon2
			m1.nb <- glmmPQL(as.formula(paste('taxon2 ~', grp.name.adj.name, '+ ldep')),  random = as.formula(subject), data = df, 
					verbose=F, family=binomial)
		}
		
		if (method == 'QB') {
			taxon2 <- as.numeric(taxon != 0)
			df$taxon2 <- taxon2
			m1.nb <- glmmPQL(as.formula(paste('taxon2 ~', grp.name.adj.name, '+ ldep')),  random = as.formula(subject), data = df, 
					verbose=F, family=quasibinomial())
		}
		
		pv.nb <- wald.test(b = fixed.effects(m1.nb), Sigma = vcov(m1.nb), Terms = grep(grp.name, names(coef(m1.nb))))$result$chi2['P']
		method <- paste(method, 'Wald')
		code <- NA   # placeholder
		aic.nb <- NA # NA placeholder
		coef.nb <- fixed.effects(m1.nb)
		ci.nb <- confint.lme(m1.nb)
		
	}
	
	fc.nb <- coef.nb[grep(grp.name, names(coef.nb))]
	obj <- ci.nb[grep(grp.name, rownames(ci.nb)), ]
	
	if (is.vector(obj)) {
		fc.lc.nb <- obj[1]
		fc.uc.nb <- obj[2]
	} else {
		fc.lc.nb <- obj[, 1]
		fc.uc.nb <- obj[, 2]
		names(fc.lc.nb) <- paste(names(fc.lc.nb), '2.5%')
		names(fc.uc.nb) <- paste(names(fc.uc.nb), '97.5%')
	}
	
	return(list(method=method, pv=pv.nb, lfc=fc.nb, lfc.lci=fc.lc.nb, lfc.uci=fc.uc.nb, aic=aic.nb, code=code))
}



perform_differential_analysis_para_single_FE <- function (taxon.abund, ldep, grp.name, adj.name=NULL, subject=NULL, df, method='NB', LRT=FALSE) {
	# ldep: log depth (size factor)
	if (!is.null(adj.name)) {
		if (sum(grepl(grp.name, c(adj.name)))) {
			stop('grp.name could not be part of adj.name or subject, or there will be problem!\n')
		}
	}
	
	if (!is.null(subject)) {
		warnings('Fixed effects model will ignore the subject variable! Please use randome effects model!\n')
	}
	if (LRT & method == 'OP') warning('Overdispersed Poisson does not support LRT! Wald test used!\n')
	
	if (is.null(adj.name)) {
		grp.name.adj.name <- grp.name
	} else {
		grp.name.adj.name <- paste(grp.name, '+', adj.name)
	}
	if (method == 'NB') {
		m1.nb <- glm.nb(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep)')), data = df)
		if (LRT) {
			m0.nb <- update(m1.nb, as.formula(paste('. ~ . -',  grp.name)))
			code <- list(m1.conv=m1.nb$converged, m1.bound=m1.nb$boundary, m0.conv=m0.nb$converged, m0.bound=m0.nb$boundary)
			pv.nb <- anova(m1.nb, m0.nb)['Pr(Chi)'][2, ]
			method <- paste(method, 'LRT')
		} else {
			code <- list(m1.conv=m1.nb$converged, m1.bound=m1.nb$boundary)	
			pv.nb <- wald.test(b = coef(m1.nb), Sigma = vcov(m1.nb), Terms = grep(grp.name, names(coef(m1.nb))))$result$chi2['P']
			method <- paste(method, 'Wald')
		}

		aic.nb <- summary(m1.nb)$aic
		
		coef.nb <- coef(m1.nb)
		fc.nb <- coef.nb[grep(grp.name, names(coef.nb))]
		
		ci.nb <- confint.default(m1.nb)
		obj <- ci.nb[grep(grp.name, rownames(ci.nb)), ]
		
		if (is.vector(obj)) {
			fc.lc.nb <- obj[1]
			fc.uc.nb <- obj[2]
		} else {
			fc.lc.nb <- obj[, 1]
			fc.uc.nb <- obj[, 2]
			names(fc.lc.nb) <- paste(names(fc.lc.nb), '2.5%')
			names(fc.uc.nb) <- paste(names(fc.uc.nb), '97.5%')
		}
		return(list(method=method, pv=pv.nb, lfc=fc.nb, lfc.lci=fc.lc.nb, lfc.uci=fc.uc.nb, aic=aic.nb, code=code))
	}
	if (method == 'B') {
		taxon.abund2 <- as.numeric(taxon.abund != 0)
		m1.b <- glm(as.formula(paste('taxon.abund2 ~', grp.name.adj.name, '+ ldep')), data = df, family=binomial)
		if (LRT) {
			m0.b <- update(m1.b, as.formula(paste('. ~ . -',  grp.name)))
			code <- list(m1.conv=m1.b$converged, m1.bound=m1.b$boundary, m0.conv=m0.b$converged, m0.bound=m0.b$boundary)
			pv.b <- pchisq(2 * (logLik(m1.b) - logLik(m0.b)), df = df.residual(m0.b) - df.residual(m1.b), lower.tail=FALSE)
			method <- paste(method, 'LRT')
		} else {
			code <- list(m1.conv=m1.b$converged, m1.bound=m1.b$boundary)	
			pv.b <- wald.test(b = coef(m1.b), Sigma = vcov(m1.b), Terms = grep(grp.name, names(coef(m1.b))))$result$chi2['P']
			method <- paste(method, 'Wald')
		}

		aic.b <- summary(m1.b)$aic
		coef.b <- coef(m1.b)
		fc.b <- coef.b[grep(grp.name, names(coef.b))]
		
		ci.b <- confint.default(m1.b)
		obj <- ci.b[grep(grp.name, rownames(ci.b)), ]
		
		if (is.vector(obj)) {
			fc.lc.b <- obj[1]
			fc.uc.b <- obj[2]
		} else {
			fc.lc.b <- obj[, 1]
			fc.uc.b <- obj[, 2]
			names(fc.lc.b) <- paste(names(fc.lc.b), '2.5%')
			names(fc.uc.b) <- paste(names(fc.uc.b), '97.5%')
		}
		return(list(method=method, pv=pv.b, lfc=fc.b, lfc.lci=fc.lc.b, lfc.uci=fc.uc.b, aic=aic.b, code=code))
	}
	
	# Rev: 2016_09_13 add 'QB', No likelihood ratio test
	if (method == 'QB') {
		taxon.abund2 <- as.numeric(taxon.abund != 0)
		m1.b <- glm(as.formula(paste('taxon.abund2 ~', grp.name.adj.name, '+ ldep')), data = df, family=quasibinomial)
		
		code <- list(m1.conv=m1.b$converged, m1.bound=m1.b$boundary)	
		pv.b <- wald.test(b = coef(m1.b), Sigma = vcov(m1.b), Terms = grep(grp.name, names(coef(m1.b))))$result$chi2['P']
		method <- paste(method, 'Wald')
		
		
		aic.b <- summary(m1.b)$aic
		coef.b <- coef(m1.b)
		fc.b <- coef.b[grep(grp.name, names(coef.b))]
		
		ci.b <- confint.default(m1.b)
		obj <- ci.b[grep(grp.name, rownames(ci.b)), ]
		
		if (is.vector(obj)) {
			fc.lc.b <- obj[1]
			fc.uc.b <- obj[2]
		} else {
			fc.lc.b <- obj[, 1]
			fc.uc.b <- obj[, 2]
			names(fc.lc.b) <- paste(names(fc.lc.b), '2.5%')
			names(fc.uc.b) <- paste(names(fc.uc.b), '97.5%')
		}
		return(list(method=method, pv=pv.b, lfc=fc.b, lfc.lci=fc.lc.b, lfc.uci=fc.uc.b, aic=aic.b, code=code))
	}
	
	if (method == 'OP') {
		# No LRT
		m1.op <- glm(as.formula(paste('taxon.abund ~', grp.name.adj.name)), offset=ldep, data = df, family=quasipoisson)
		code <- list(m1.conv=m1.op$converged, m1.bound=m1.op$boundary)

		# pv.op <- pchisq(2 * (logLik(m1.op) - logLik(m0.op)), df = df.residual(m0.op) - df.residual(m1.op), lower.tail=FALSE) # LRT not applicable
		coef.op <- coef(m1.op)		
		pv.op <- wald.test(b = coef.op, Sigma = vcov(m1.op), Terms = grep(grp.name, names(coef.op)))$result$chi2['P']
		method <- paste(method, 'Wald')
		fc.op <- coef.op[grep(grp.name, names(coef.op))]
		
		ci.op <- confint.default(m1.op)
		obj <- ci.op[grep(grp.name, rownames(ci.op)), ]
		
		if (is.vector(obj)) {
			fc.lc.op <- obj[1]
			fc.uc.op <- obj[2]
		} else {
			fc.lc.op <- obj[, 1]
			fc.uc.op <- obj[, 2]
			names(fc.lc.op) <- paste(names(fc.lc.op), '2.5%')
			names(fc.uc.op) <- paste(names(fc.uc.op), '97.5%')
		}
		return(list(method=method, pv=pv.op, lfc=fc.op, lfc.lci=fc.lc.op, lfc.uci=fc.uc.op, aic=NULL, code=code))
	}
	
	if (method == 'ZINB0') {
		m1.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep)')),
				data = df, dist = "negbin", EM = TRUE)
		if (LRT) {
			if (is.null(adj.name)) {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~ offset(ldep)')),
						data = df, dist = "negbin", EM = TRUE)
			} else {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', adj.name, '+ offset(ldep)')),
						data = df, dist = "negbin", EM = TRUE)
			}
			code <- list(m1.conv=m1.zinb$converged, m0.conv=m0.zinb$converged)
			# LRT
			pv.zinb <-  pchisq(2 * (logLik(m1.zinb) - logLik(m0.zinb)), df = df.residual(m0.zinb) - df.residual(m1.zinb), lower.tail=FALSE)
			method <- paste(method, 'LRT')
		} else {
			code <- list(m1.conv=m1.zinb$converged)	
			pv.zinb <- wald.test(b = coef(m1.zinb), Sigma = vcov(m1.zinb), Terms = grep(grp.name, names(coef(m1.zinb))))$result$chi2['P']
			method <- paste(method, 'Wald')
		}
		
		aic.zinb <- -2 * logLik(m1.zinb) + 2 * (m1.zinb$n - m1.zinb$df.residual)
		
		coef.zinb <- coef(m1.zinb)
		fc.zinb <- coef.zinb[grep(grp.name, names(coef.zinb))]
		
		ci.zinb <- confint.default(m1.zinb)
		obj <- ci.zinb[grep(grp.name, rownames(ci.zinb)), ]
		
		if (is.vector(obj)) {
			fc.lc.zinb <- obj[1]
			fc.uc.zinb <- obj[2]
		} else {
			fc.lc.zinb <- obj[, 1]
			fc.uc.zinb <- obj[, 2]
			names(fc.lc.zinb) <- paste(names(fc.lc.zinb), '2.5%')
			names(fc.uc.zinb) <- paste(names(fc.uc.zinb), '97.5%')
		}
		return(list(method=method, pv=pv.zinb, lfc=fc.zinb, lfc.lci=fc.lc.zinb, lfc.uci=fc.uc.zinb, aic=aic.zinb, code=code))
	}
	
	
	if (method == 'ZINB1') {
		m1.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep) | ldep')),
				data = df, dist = "negbin", EM = TRUE)
		if (LRT) {
			if (is.null(adj.name)) {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~ offset(ldep) | ldep')),
						data = df, dist = "negbin", EM = TRUE)
			} else {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', adj.name, '+ offset(ldep) | ldep')),
						data = df, dist = "negbin", EM = TRUE)
			}
			code <- list(m1.conv=m1.zinb$converged, m0.conv=m0.zinb$converged)
			# LRT
			pv.zinb <-  pchisq(2 * (logLik(m1.zinb) - logLik(m0.zinb)), df = df.residual(m0.zinb) - df.residual(m1.zinb), lower.tail=FALSE)
			method <- paste(method, 'LRT')
		} else {
			code <- list(m1.conv=m1.zinb$converged)	
			pv.zinb <- wald.test(b = coef(m1.zinb), Sigma = vcov(m1.zinb), Terms = grep(grp.name, names(coef(m1.zinb))))$result$chi2['P']
			method <- paste(method, 'Wald')
		}
		
		aic.zinb <- -2 * logLik(m1.zinb) + 2 * (m1.zinb$n - m1.zinb$df.residual)
		
		coef.zinb <- coef(m1.zinb)
		fc.zinb <- coef.zinb[grep(grp.name, names(coef.zinb))]
		
		ci.zinb <- confint.default(m1.zinb)
		obj <- ci.zinb[grep(grp.name, rownames(ci.zinb)), ]
		
		if (is.vector(obj)) {
			fc.lc.zinb <- obj[1]
			fc.uc.zinb <- obj[2]
		} else {
			fc.lc.zinb <- obj[, 1]
			fc.uc.zinb <- obj[, 2]
			names(fc.lc.zinb) <- paste(names(fc.lc.zinb), '2.5%')
			names(fc.uc.zinb) <- paste(names(fc.uc.zinb), '97.5%')
		}
		return(list(method=method, pv=pv.zinb, lfc=fc.zinb, lfc.lci=fc.lc.zinb, lfc.uci=fc.uc.zinb, aic=aic.zinb, code=code))
	}
	
	
	if (method == 'ZINB2') {
		m2.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', grp.name.adj.name, '+ offset(ldep) |', grp.name.adj.name, '+ ldep')),
				data = df, dist = "negbin", EM = TRUE)
		if (LRT) {
			if (is.null(adj.name)) {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~ offset(ldep) | ldep')),
						data = df, dist = "negbin", EM = TRUE)
			} else {
				m0.zinb <- zeroinfl(as.formula(paste('taxon.abund ~', adj.name, '+ offset(ldep) |', adj.name, ' + ldep')),
						data = df, dist = "negbin", EM = TRUE)
			}
			
			code <- list(m1.conv=m2.zinb$converged, m0.conv=m0.zinb$converged)
			# LRT
			pv2.zinb <-  pchisq(2 * (logLik(m2.zinb) - logLik(m0.zinb)), df = df.residual(m0.zinb) - df.residual(m2.zinb), lower.tail=FALSE)
			method <- paste(method, 'LRT')
		} else {
			code <- list(m2.conv=m2.zinb$converged)	
			pv2.zinb <- wald.test(b = coef(m2.zinb), Sigma = vcov(m2.zinb), Terms = grep(grp.name, names(coef(m2.zinb))))$result$chi2['P']
			method <- paste(method, 'Wald')
		}

		aic2.zinb <- -2 * logLik(m2.zinb) + 2 * (m2.zinb$n - m2.zinb$df.residual)
		
		coef.zinb <- coef(m2.zinb)
		fc2.zinb <- coef.zinb[grep(grp.name, names(coef.zinb))]
		
		ci.zinb <- confint.default(m2.zinb)
		obj <- ci.zinb[grep(grp.name, rownames(ci.zinb)), ]
		
		if (is.vector(obj)) {
			fc2.lc.zinb <- obj[1]
			fc2.uc.zinb <- obj[2]
		} else {
			fc2.lc.zinb <- obj[, 1]
			fc2.uc.zinb <- obj[, 2]
			names(fc2.lc.zinb) <- paste(names(fc2.lc.zinb), '2.5%')
			names(fc2.uc.zinb) <- paste(names(fc2.uc.zinb), '97.5%')
		}
		return(list(method=method, pv=pv2.zinb, lfc=fc2.zinb, lfc.lci=fc2.lc.zinb, lfc.uci=fc2.uc.zinb, aic=aic2.zinb, code=code))
	}

}

# Rev: 2016_09_14 Add Adaptive 0: switch OP and QB depending on the number of zeros. 
# Rev: 2017_10_30  Support filtering based on coefficient of variation
perform_differential_analysis_para <- function (data.obj,  grp.name, adj.name=NULL, subject=NULL, RE=FALSE, method='Adaptive0', zerop.cutoff=0.25, ZINB='ZINB1', LRT=FALSE, 
		taxa.levels=c('Phylum', 'Order', 'Class', 'Family', 'Genus', 'Species'), winsor=TRUE, winsor.qt=0.97, norm='GMPR', norm.level='Species', intersect.no=4,
		prev=0.1, minp=0.002, medianp=NULL, cv=NULL, mt.method='fdr', cutoff=0.15, ann='', ...) {
	# To be completed
	# subject holds the random effects formula
	if (!RE) {
		if (!(method %in% c('ZINB', 'B', 'QB', 'NB', 'OP', 'Adaptive0', 'Adaptive1', 'Adaptive2'))) stop('The speficied model is not supported!\n')
		perform_differential_analysis_para_single <- perform_differential_analysis_para_single_FE
		if (!is.null(subject)) warning('subject will not be used. Are you sure you want to run fixed effects model? ')
	} else {
		if (!(method %in% c('ZINB', 'B', 'B0', 'QB', 'NB', 'OP', 'Adaptive0', 'Adaptive1', 'Adaptive2'))) stop('The speficied model does not have random effects implementation!\n')
		if (ZINB != 'ZINB1') stop('Currently only ZINB1 is supported!\n')
		if (is.null(subject)) warning('subject is not supplied. Fixed effects model will be used instead!\n')
		perform_differential_analysis_para_single <- perform_differential_analysis_para_single_RE
	}

	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	ind <- !is.na(grp)
	data.obj <- subset_data(data.obj, ind)
	grp <- grp[ind]
	df <- df[ind, ]
	
	if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
		data.obj$abund.list[['Species']] <- data.obj$otu.tab
		rownames(data.obj$abund.list[['Species']]) <- paste0("OTU", rownames(data.obj$otu.tab), ":", 
				data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
	}
	
	dep <- colSums(data.obj$otu.tab)
	diff.seq.p <- summary(aov(dep ~ grp))[[1]][1, 'Pr(>F)']
	
	if (!is.na(diff.seq.p) & diff.seq.p <= 0.05) {
		cat("Signficant sequencing depth confounding!\n")
		cat("For parametric test with sequence depth adjustment, please be cautious about the results!\n")
		cat("There may be potential residual sequence depth confounding!\n")
	}
	
	pv.list <- qv.list <-  fc.list <- fc.lc.list <- fc.uc.list <- met.list <- list()
	res.final <- NULL
	
	if (norm == 'Precalculated') {
		dep <- data.obj$size.factor
	}
	
	if (norm == 'GMPR') {
		if (norm.level %in% c('OTU', 'Species')) {
			tab <- data.obj$otu.tab
		} else {
			tab <- data.obj$abund.list[[norm.level]]
		}
		dep <- GMPR(tab, intersect.no)
	}
	
	if (norm == 'TSS') {
		dep <- colSums(data.obj$otu.tab)
	}
	
	ldep <- log(dep)
	
	for (LOI in taxa.levels) {
		cat(LOI, "\n")

		taxon.ct <- data.obj$abund.list[[LOI]]
		

		if (winsor == TRUE) {
			# Addressing the outlier (97% percent) or at least one outlier
			
			taxon.ct.p <- t(t(taxon.ct) / dep)
			taxon.ct.p <- apply(taxon.ct.p, 1, function(x) {
						cutoff <- quantile(x, winsor.qt)
						x[x >= cutoff] <- cutoff
						x
					}
			)
			# column/row switch
			taxon.ct <- t(round(taxon.ct.p * dep))
			
		}
		
		prop <- t(t(taxon.ct) / colSums(taxon.ct))
#		if (!is.null(minp)) {
#			prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
#			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
#		}
#		
#		if (!is.null(medianp)) {
#			nz.mean <- apply(prop, 1, function(x) median(x[x!=0]))
#			prop <- prop[nz.mean > medianp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
#			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
#		}
		
		if (!is.null(prev)) {
			prop <- prop[rowSums(prop!=0) > prev * ncol(prop), , drop=FALSE]	
			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
		}
		if (!is.null(minp)) {
			prop <- prop[rowMaxs(prop) > minp, , drop=FALSE]	
			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
		}
		
		if (!is.null(medianp)) {
			nz.mean <- apply(prop, 1, function(x) median(x[x!=0]))
			prop <- prop[nz.mean > medianp, , drop=FALSE]	
			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
		}
		
		if (!is.null(cv)) {
			prop <- prop[rowSds(prop) / rowMeans(prop) > cv, , drop=FALSE]	
			taxon.ct <- taxon.ct[rownames(prop), , drop=FALSE]
		}
			
		pv.vec <-  fc.vec <- fc.lc.vec <- fc.uc.vec <- met.vec <- conv.vec <- NULL
		obj <- NULL
		for (taxon in rownames(taxon.ct)) {
			cat('.')
			taxon.abund <- taxon.ct[taxon, ]
			
			######## Logistic regression ###############
			if (method == 'B0') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='B0', LRT, ...)) 
			if (method == 'B') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='B', LRT, ...)) 
			if (method == 'QB') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='QB', LRT, ...)) 
			######## Overdispersed Poisson regression #########	
			if (method == 'OP') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='OP', LRT, ...)) 
			######## Negative binomial regression #########
			if (method == 'NB') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='NB', LRT, ...)) 
			######## Zeroinflated negbinomial regression 1 ########
			if (method == 'ZINB') error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method=ZINB, LRT, ...)) 
			
			# Adpative 0 selects OP and QB based on the zero proportion (Not optimal)
			if (method == 'Adaptive0') {
				temp <- mean(as.numeric(taxon.abund != 0))
				
				if (temp > zerop.cutoff) {
					error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='QB', LRT, ...)) 
				} else {
					error <- try(obj <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='OP', LRT, ...)) 
				}
			}
			
			# Adpative 1 selects NB and ZIB based on AIC
			if (method == 'Adaptive1') {
				error1 <- try(obj1 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='NB', LRT, ...)) 
				error2 <- try(obj2 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method=ZINB, LRT, ...)) 
				if (class(error1) != 'try-error' & class(error2) != 'try-error') {
					if (obj1$aic < obj2$aic) {
						obj <- obj1
					} else {
						obj <- obj2
					}
					error <- error1
				} else {
					# pv == 0 indicates some problems in fitting
					if (class(error1) != 'try-error' & obj1$pv != 0) {
						obj <- obj1
						error <- error1
					} else {
						if (class(error2) != 'try-error' & obj1$pv != 0) {
							obj <- obj2
							error <- error2
						} else {
							error <- error2
						}
					}
				}	
			}
			
			# Adaptive 2 starts with NB model, if it fails, it switches ZINB
			if (method == 'Adaptive2') {
				error1 <- try(obj1 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method='NB', LRT, ...)) 
				
				if (class(error1) == 'try-error' | obj1$pv == 0) {
					error2 <- try(obj2 <- perform_differential_analysis_para_single(taxon.abund, ldep, grp.name, adj.name, subject, df, method=ZINB, LRT, ...)) 
					if (class(error2) != 'try-error') {
						obj <- obj2
						error <- error2
					} else {
						error <- error2
					}
				} else {
					error <- error1
					obj <- obj1
				}
				
			}
			
			# Random Effects model 
			# ZINB, B, NB, Adpative1 is implemented based on glmmADMB
			
            # Set P value NA for those not makes sense
			if (class(error) == "try-error" | abs(obj$lfc) > 100) {
				obj$pv <- obj$lfc <- obj$lfc.lci <- obj$lfc.uci <- obj$method <- NA
			}
			
			
			pv.vec <- rbind(pv.vec, obj$pv)
			fc.vec <- rbind(fc.vec, obj$lfc)
			fc.lc.vec <- rbind(fc.lc.vec, obj$lfc.lci)
			fc.uc.vec <- rbind(fc.uc.vec, obj$lfc.uci)
			met.vec <- rbind(met.vec, obj$method)
			
		}
		cat('\n')

		qv.vec <- matrix(p.adjust(pv.vec[, 1], 'fdr'), ncol=1)
		
		rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(fc.uc.vec) <- rownames(fc.lc.vec) <- rownames(met.vec) <- rownames(prop)
		colnames(pv.vec) <- 'Pvalue'
		colnames(qv.vec) <- 'Qvalue'
        colnames(met.vec) <- 'Method'

		
		pv.list[[LOI]] <- pv.vec
		qv.list[[LOI]] <- qv.vec
		fc.list[[LOI]] <- fc.vec
		fc.lc.list[[LOI]] <- fc.lc.vec
		fc.uc.list[[LOI]] <- fc.uc.vec
		met.list[[LOI]] <- met.vec
 
		res <- cbind(pv.vec, qv.vec, fc.vec, fc.lc.vec, fc.uc.vec, met.vec)
		rownames(res) <- rownames(prop)
		write.csv(res, paste0("Taxa_DifferentialAbundanceAnalysis_", LOI, "_", ann, ".csv"))
		
		if (mt.method == 'fdr') {
			res.final <- rbind(res.final, res[as.numeric(res[, 'Qvalue']) <= cutoff, , drop=F])
		}
		if (mt.method == 'raw') {
			res.final <- rbind(res.final, res[ as.numeric(res[, 'Pvalue']) <= cutoff, , drop=F])
		}
	}
	
	if (!is.null(res.final)) {
		colnames(res.final) <- colnames(res)
		res.final <- res.final[rowSums(is.na(res.final)) == 0, , drop=F]
		write.csv(res.final, paste0("Taxa_DifferentialAbundanceAnalysis_AllLevels_", mt.method, '_', cutoff, "_", ann, ".csv"))
	}
	return(list(pv.list=pv.list, qv.list=qv.list, fc.list=fc.list, fc.uc.list=fc.uc.list, fc.lc.list=fc.lc.list, met.list=met.list))
}


plot_effect_size <- function (month,  value, pos.lab, neg.lab, ylab, hjust1=1.3, hjust2=-0.3, lab.size=3, xsize=10) {
	month <- factor(month, levels=month)
	dtm <- data.frame(month=month, value=value)
	dtm$colour <- factor(ifelse(dtm$value < 0, neg.lab, pos.lab), levels=c(pos.lab, neg.lab))
	dtm$hjust <- ifelse(dtm$value > 0, hjust1, hjust2)
	obj <- ggplot(dtm, aes(month, value, label = month, hjust = hjust)) + 
			geom_text(aes(y = 0, colour = colour), size=lab.size) + 
			geom_bar(stat = "identity", aes(fill = colour)) +
			theme(axis.text.x = element_text(size=xsize)) +
			ylim(c(-max(abs(value))*1.1, max(abs(value))*1.1)) +
			coord_flip() + 
			scale_x_discrete(breaks = NULL) +
			labs(x = "", y = ylab) +
			theme(legend.position="top", legend.title=element_blank())
}


plot_effect_size2 <- function (fold.dat.plot1, ylabel='log(Fold change)', is.ln=TRUE, ord=TRUE) {
	
	if (is.ln) {
		fold.dat.plot1[, c('Estimate', 'LCI', 'UCI')] <- fold.dat.plot1[, c('Estimate', 'LCI', 'UCI')] * log2(exp(1))
	}
	Alphas <- seq(1, 99, 2) / 100
#Multiplier <- qnorm(1 - Alphas / 2)
	Multiplier <- seq(1, 0.01, len=50)
	zzTransparency <<- 1/(length(Multiplier)/4)
	
#	fold.dat.plot1 <- data.frame(IV=rownames(fold.dat.plot), Estimate=fold.dat.plot[, 1], LCI=fold.dat.plot[, 2], UCI=fold.dat.plot[, 3])
	fold.dat.plot1 <- data.frame(cbind(fold.dat.plot1, Scalar=rep(Multiplier, each = nrow(fold.dat.plot1))))
#	fold.dat.plot1[,  -1] <- apply(fold.dat.plot1[, -1], 2, function(x){as.numeric(as.character(x))})
	fold.dat.plot1$Emphasis <- by(1 - seq(0, 1, length = length(Multiplier) + 1)[-(length(Multiplier) + 1)],
			as.character(round(Multiplier, 5)), mean)[as.character(round(fold.dat.plot1$Scalar, 5))]
	
	fold.dat.plot1$IV <- factor(fold.dat.plot1$IV, unique(fold.dat.plot1$IV))
	
	if (ord) {
		fold.dat.plot1 <- fold.dat.plot1[order(as.character(fold.dat.plot1$IV)), ]
	}
	
	OutputPlot <- ggplot2::qplot(data = fold.dat.plot1, x = IV, y = Estimate,
			ymin = Estimate - (Estimate -LCI)*Scalar, ymax = Estimate + (UCI - Estimate)*Scalar,
			ylab = NULL, xlab = NULL, alpha = I(zzTransparency), colour = I(gray(0)), geom = "blank")
	
	OutputPlot <- OutputPlot + geom_hline(yintercept = 0, lwd = I(7/12), colour = I(hsv(0/12, 7/12, 7/12)), alpha = I(5/12))
	OutputPlot <- OutputPlot + geom_linerange(data = fold.dat.plot1, aes(size = 1/Emphasis), alpha = I(zzTransparency), colour = I(gray(0)))
	OutputPlot <- OutputPlot + scale_size_continuous() + guides(size=FALSE)
#OutputPlot <- OutputPlot + facet_grid(~ ModelName)
	OutputPlot <- OutputPlot + coord_flip() + geom_point(aes(x = IV, y = Estimate), colour = I(gray(0))) + theme_bw() + ylab(ylabel)
}


# Rev: 2016_11_25 Uniqufy
# Rev: 2016_04_18 strata

visualize_differential_analysis <- function (data.obj, diff.obj,  grp.name=NULL, strata=NULL, test='Nonpara', mt.method='fdr', scale='sqrt', cutoff=0.15,
		taxa.levels=c('Phylum', 'Family', 'Genus'), ord=TRUE, eff.type='logP', indivplot=TRUE, colFnsC=NULL, colFnsF=NULL, subject=NULL,
		xsize=10, ann='', hei1=NULL, wid1=NULL, hei2=NULL, wid2=NULL) {
	
	# uniquefy names
	# For backward compatibility. Newer version will not need this and below. The old version has 'unclassified' which leads to duplicate names.
	# Newer version has 'Unclassified'. Case difference.
	
	# Check whether there is name duplication
	check.names <- NULL
	
	obj0 <- diff.obj[[1]]
	for (level in names(obj0)) {
		obj <- obj0[[level]]
		# rownames(obj) <- gsub('unclassified', paste0('Unclassified',substr(level, 1, 1)), rownames(obj))
		check.names <- c(check.names, rownames(obj))
	}
	
	if (sum(table(check.names) >= 2)) {
		data.obj <- uniquefy_taxa_names(data.obj)
		
		for (name1 in names(diff.obj)) {
			obj0 <- diff.obj[[name1]]
			for (level in names(obj0)) {
				obj <- obj0[[level]]
				# rownames(obj) <- gsub('unclassified', paste0('Unclassified',substr(level, 1, 1)), rownames(obj))
				rownames(obj) <- paste0(rownames(obj), substr(level, 1, 1))
				obj0[[level]] <- obj
			}
			diff.obj[[name1]] <- obj0
		}
		
	}

	#pv.list <- diff.obj$pv.list
	fc.list <- diff.obj$fc.list
	qv.list <- diff.obj$qv.list
	pv.list <- diff.obj$pv.list
	if (test == 'Para') {
		fc.lc.list <- diff.obj$fc.lc.list
		fc.uc.list <- diff.obj$fc.uc.list
	}
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	
	ind <- !is.na(grp)
	data.obj <- subset_data(data.obj, ind)
	grp <- grp[ind]
	df <- df[ind, ]
	
	prop <- NULL
	eff <- eff.lc <- eff.uc <- NULL
	taxa.names <- NULL
	if (is.null(taxa.levels)) {
		LOIs <- names(qv.list)
	} else {
		LOIs <- taxa.levels
		if (sum(!(taxa.levels %in% names(qv.list)))) {
			stop('Taxa levels are not contained in differential abundance analysis results!\n')
		}
	}
	for (LOI in LOIs) {
		pv.vec <- pv.list[[LOI]]
		fc.vec <- fc.list[[LOI]]
		#qv.vec <- qvalue(pv.vec[, 1])$qvalues
		qv.vec <- qv.list[[LOI]]
		
		if (test == 'Para') {
			fc.lc.vec <- fc.lc.list[[LOI]]
			fc.uc.vec <- fc.uc.list[[LOI]]
		}
		
		if (mt.method == 'fdr') {
			taxa.name <- rownames(qv.vec)[qv.vec <= cutoff]
			taxa.name <- taxa.name[!is.na(taxa.name)]
		}

		if (mt.method == 'raw') {
			taxa.name <- rownames(pv.vec)[pv.vec <= cutoff]
			taxa.name <- taxa.name[!is.na(taxa.name)]
		}
		
		if (length(taxa.name) != 0) {
			prop0 <- data.obj$abund.list[[LOI]]
			prop0 <- t(t(prop0) / colSums(prop0))
			prop0 <-  prop0[taxa.name, , drop=F]
			if (ord == TRUE) {
				prop0 <- prop0[rev(order(rowMeans(prop0))), , drop=F]
			}
			prop <- rbind(prop, prop0)
			# currently using fold change
		    if (test == 'Para') {
				eff <- rbind(eff, fc.vec[taxa.name, , drop=F])
				eff.lc <- rbind(eff.lc, fc.lc.vec[taxa.name, , drop=F])
				eff.uc <- rbind(eff.uc, fc.uc.vec[taxa.name, , drop=F])
			} else {
				if (eff.type == 'LFC') {
					eff <- c(eff, fc.vec[taxa.name, ])
				}
				if (eff.type == 'logP') {
					eff <- c(eff, sign(fc.vec[taxa.name, ]) * (-log10(pv.vec[taxa.name, ])))
				}
				
			}
			taxa.names <- c(taxa.names, taxa.name)
		}
	}
	
	if (length(taxa.names) == 0) {
		cat('No differential taxa! \n')
	} else {
		if (length(taxa.names) >= 2) {
			if (is.null(wid1) | is.null(hei1)) {
				wid1 <- 7 * ifelse(nrow(prop) / 30 < 1, 1, nrow(prop) / 30)
				hei1 <- 7
			}
	
			if (!is.null(grp.name)) {
				obj <- taxa_barplot_aggregate(prop, df, grp.name, strata, scale, xsize)
				pdf(paste("Taxa_DifferentialAbundance_AbundanceBarplot", scale, mt.method, cutoff, ann, ".pdf", sep="_"), height=hei1, width=wid1)
				print(obj)
				dev.off()
				
				obj1 <- taxa_boxplot_aggregate(prop, df, grp.name, strata, scale, xsize) 
				pdf(paste("Taxa_DifferentialAbundance_AbundanceBoxplot", scale, mt.method, cutoff, ann, ".pdf", sep="_"), height=hei1, width=wid1)
				print(obj1)
				dev.off()

			}

			# currently fold change
		    if (test == 'Para') {
				rownames(eff) <- rownames(eff.lc) <- rownames(eff.uc) <- taxa.names
				for (k in 1:ncol(eff)) {
					fold.dat.plot1 <- data.frame(Estimate=eff[, k], LCI=eff.lc[, k], UCI=eff.uc[, k], IV=taxa.names)
					obj <- plot_effect_size2(fold.dat.plot1)
					pdf(paste("Taxa_DifferentialAbundance_EffectBarplot", colnames(eff)[k], mt.method, cutoff, ann, ".pdf", sep="_"), height=hei2, width=wid2)
					print(obj)
					dev.off()
				}
			} else {
				if (!is.na(eff[1])) {
					names(eff) <- taxa.names
					eff <- eff[!is.na(eff) & is.finite(eff)]
					eff <- sort(eff)
					taxa.names2 <- names(eff)
					if (is.null(wid2) | is.null(hei2)) {
						hei2 <- 4 + length(taxa.names2) / 20 * 3 
						wid2 <- 6
					}
					if (eff.type == 'LFC') {
						obj <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='Log2 fold change')
						pdf(paste("Taxa_DifferentialAbundance_LFCBarplot", mt.method, cutoff, ann, ".pdf", sep="_"), height=hei2, width=wid2)
						print(obj)
						dev.off()	
					}

					if (eff.type == 'Spearman') {
						obj <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='Spearman correlation')
						pdf(paste("Taxa_DifferentialAbundance_LFCBarplot", mt.method, cutoff, ann, ".pdf", sep="_"), height=hei2, width=wid2)
						print(obj)
						dev.off()	
					}
					
					if (eff.type == 'logP') {
						obj <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='-log10(P)')
						pdf(paste("Taxa_DifferentialAbundance_logPBarplot", mt.method, cutoff, ann, ".pdf", sep="_"), height=hei2, width=wid2)
						print(obj)
						dev.off()	
					}					
				}
			}
			
			# create heatmp
#			taxa.names2 <- taxa.names[!grepl('unclassified', taxa.names, ignore.case=T)]
			if (!is.null(grp.name)) {
				
				generate_taxa_heatmap(data.obj, taxa.levels='All', sam.ord=order(grp), taxa=taxa.names, meta.info=c(grp.name, strata), Colv=F, dendrogram='row', 
						ann=paste0(mt.method, '_', cutoff, '_', ann, '_Unclustered'), colFnsC=colFnsC, colFnsF=colFnsF)
				generate_taxa_heatmap(data.obj, taxa.levels='All', taxa=taxa.names, meta.info=c(grp.name, strata),  ann=paste0(mt.method, '_', cutoff, '_', ann), colFnsC=colFnsC, colFnsF=colFnsF)
				generate_taxa_heatmap(data.obj, taxa.levels='All', taxa=taxa.names, meta.info=c(grp.name, strata), data.type='R', ann=paste0(mt.method, '_', cutoff, '_Rank_', ann), colFnsC=colFnsC, colFnsF=colFnsF)
				try(
						generate_taxa_biplot(data.obj, taxa=taxa.names, trans='sqrt', grp.name, ann=paste0(mt.method, '_', cutoff, '_', ann), varname.size = 1.5)	
				)	
			}
		}
		if (!is.null(grp.name)) {
			# Individual plots
	        if (indivplot == TRUE) {
				generate_taxa_boxplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, ann=paste0(mt.method, '_', cutoff, '_', ann))
				if (!is.null(subject)) {
					generate_taxa_boxplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, subject=subject, ann=paste0(mt.method, '_', cutoff, '_', ann, '_Paired'))
				}
				generate_taxa_boxplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, scale='binary', ann=paste0(mt.method, '_', cutoff, '_', ann))
				generate_taxa_barplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, ann=paste0(mt.method, '_', cutoff, '_', ann))
			}

		}
	}
}


# Rev: 2016_12_07 when the names of some differential taxa are not consistent
# Rev: 2017_02_21  Color by log fold change m.list ->> fc.list
# To-do-list:  size by effects, accomodate OTU level
create_lefse_format <- function(data.obj, diff.obj, grp.name, mt.method='fdr', cutoff=0.15, prev=0.1, minp=0.002, 
		lefse.dir="/data2/microbiome/jeff/tools/nsegata-lefse/", ann="") {
	
	# To be improved - currently no effect size
	df <- data.obj$meta.dat
	grp <- df[, grp.name]
	levels(grp) <- paste0(1:nlevels(grp), levels(grp))

	qv.list <- diff.obj$qv.list
	pv.list <- diff.obj$pv.list
	fc.list <- diff.obj$fc.list
	
	otu.name.12 <- data.obj$otu.name
	otu.tab.12 <- data.obj$otu.tab
	tax.family.a <- NULL
	for (i in 1:6) {
		tax.family <- apply(otu.name.12[, 1:i, drop=F], 1, paste, collapse=".")
		phlan.tab <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
		rownames(phlan.tab) <- phlan.tab[, 1]
		phlan.tab <- as.matrix(phlan.tab [, -1, drop=F])
		phlan.tab <- t(t(phlan.tab) / colSums(phlan.tab))	
		phlan.tab <- phlan.tab[rowSums(phlan.tab!=0) > prev*ncol(phlan.tab) & rowMaxs(phlan.tab) > minp, , drop=F]
		tax.family.a <- c(tax.family.a, rownames(phlan.tab))
		#tax.family.a <- c(tax.family.a, unique(tax.family))
	}
	alias.a <- sapply(strsplit(tax.family.a, "\\."), function(x) {
				if (length(x) >= 2) {
					if (length(x) == 2) {
						return(x[2])
					} else {
						return(paste0(x[2], ";", x[length(x)]))
					}
				} else {
					return(x)
				}
			})
	
	taxa.names <- NULL
	abundant.grp.names <- NULL
	
	for (LOI in setdiff(names(qv.list), 'Species')) {
		qv.vec <- qv.list[[LOI]]
		pv.vec <- pv.list[[LOI]]
		#qv.vec <- qvalue(pv.vec[, 1])$qvalues
		fc.vec <- fc.list[[LOI]]
		if (mt.method == 'fdr') {
			taxa.name <- rownames(qv.vec)[qv.vec <= cutoff]
		#	abundant.grp.name <- apply(m.vec, 1, function(x) levels(grp)[which.max(x)])[qv.vec <= cutoff]
			abundant.grp.name <- apply(fc.vec, 1, function(x) levels(grp)[as.numeric(x >= 0) + 1])[qv.vec <= cutoff]
		}

		if (mt.method == 'raw') {
			taxa.name <- rownames(pv.vec)[pv.vec <= cutoff]
		#	abundant.grp.name <- apply(m.vec, 1, function(x) levels(grp)[which.max(x)])[pv.vec <= cutoff]
		    abundant.grp.name <- apply(fc.vec, 1, function(x) levels(grp)[as.numeric(x >= 0) + 1])[pv.vec <= cutoff]
		}
		
		if (length(taxa.name) != 0) {
			taxa.names <- c(taxa.names, taxa.name)
			abundant.grp.names <- c(abundant.grp.names, abundant.grp.name)
		}
	}

	# remove 'unclassified'
    abundant.grp.names <- abundant.grp.names[!grepl('unclassified', taxa.names, ignore.case=T)]
	taxa.names <- taxa.names[!grepl('unclassified', taxa.names, ignore.case=T)]
	
	# Rev: 2016_12_07
#	taxa.names <- intersect(taxa.names, alias.a)

	ind <- match(taxa.names, alias.a)
	na.ind <- !is.na(ind)

	phlan <- cbind(tax.family.a, '1.5', '\t', '\t', '-')

	phlan[ind[na.ind], 2:5] <- cbind('3.5',  abundant.grp.names[na.ind], '3.0', '-')
	write.table(phlan,  paste0('Lefse.LDA.', ann, '.txt'), row.names=F, col.names=F, quote=F, sep='\t')
#	cmd1 <- paste0('python ', lefse.dir,  'plot_cladogram.py Lefse.LDA.txt Taxa_DifferentialAbundance_', mt.method, cutoff, '_Cladogram_', ann, '.pdf', ' --format pdf')
#	cmd2 <- paste0("python ", lefse.dir, 'plot_res.py Lefse.LDA.txt Taxa_DifferentialAbundance_', mt.method, cutoff, '_LDA_', ann, '.pdf --format pdf')
#	system(cmd1)
#	system(cmd2)
}

perform_lefse_analysis <- function (data.obj,  grp.name, sub.grp.name=NULL, prev=0.1, minp=0.002,
		lefse.dir="/data2/microbiome/jeff/tools/nsegata-lefse/", ann="") {
	otu.name.12 <- data.obj$otu.name
	otu.tab.12 <- data.obj$otu.tab
	meta.dat <- data.obj$meta.dat
	
	levels(meta.dat[, grp.name]) <- paste0(1:nlevels(meta.dat[, grp.name]), levels(meta.dat[, grp.name]))
	
	tax.family <- apply(otu.name.12, 1, paste, collapse="|")
	lefse.tab <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
	rownames(lefse.tab) <- lefse.tab[, 1]
	lefse.tab <- as.matrix(lefse.tab [, -1])

	lefse.tab <- t(t(lefse.tab) / colSums(lefse.tab))
	
	lefse.tab <- lefse.tab[rowMaxs(lefse.tab) > minp & rowSums(lefse.tab!=0) > prev*ncol(lefse.tab), , drop=FALSE]

	if (is.null(sub.grp.name)) {
		header <- rbind(class=as.character(meta.dat[colnames(lefse.tab), grp.name]), 
				id=colnames(lefse.tab))
	} else {
		header <- rbind(class=as.character(meta.dat[colnames(lefse.tab), grp.name]), 
				subclass=as.character(meta.dat[colnames(lefse.tab), sub.grp.name]), 
				id=colnames(lefse.tab))
	}

	
	lefse.tab <- rbind(header, lefse.tab)
	write.table(lefse.tab, "lefse.txt", sep="\t", col.names=F, quote=F)
	
#	system(paste0('mkdir LefSe_', ann))
#	output <- paste0('LefSe_', ann, '/')
#	if (is.null(sub.grp.name)) {
#		cmd1 <- paste0("python ", lefse.dir, "format_input.py lefse.txt ", output, "temp.in ", "-c 1 -u 2 -o 1000000")
#	} else {
#		cmd1 <- paste0("python ", lefse.dir, "format_input.py lefse.txt ", output, "temp.in ", "-c 1 -s 2 -u 3 -o 1000000")
#	}
#
#	cmd2 <- paste0("python ", lefse.dir, "run_lefse.py ", output, "temp.in ",  output, "lda.res")
#	cmd3 <- paste0("python ", lefse.dir, "plot_res.py ",  output, "lda.res ", output, "lda.pdf --format pdf")
#	cmd4 <- paste0("python ", lefse.dir, "plot_cladogram.py ",  output, "lda.res ", output, "cladogram.pdf --format pdf --labeled_stop_lev 6 --abrv_stop_lev  6")
#	
#	system(cmd1)
#	system(cmd2)
#	system(cmd3)
#	system(cmd4)
}

plot.Boruta2 <- function (x, colCode = c("green", "yellow", "red", "blue"), sort = TRUE,
		whichShadow = c(TRUE, TRUE, TRUE), col = NULL, xlab = "Attributes", ids=NULL,
		ylab = "Importance", ...)
{
	if (class(x) != "Boruta")
		stop("This function needs Boruta object as an argument.")
	lz <- lapply(1:ncol(x$ImpHistory), function(i) x$ImpHistory[is.finite(x$ImpHistory[, i]), i])
	names(lz) <- colnames(x$ImpHistory)
	numShadow <- sum(whichShadow)
	lz <- lz[c(rep(TRUE, length(x$finalDecision)), whichShadow)]
	col <- Boruta:::generateCol(x, colCode, col, numShadow)
	if (sort) {
		ii <- order(sapply(lz, median))
		lz <- lz[ii]
		col <- col[ii]
	}
	names(lz) <- gsub("^X",  "", names(lz))
	if (is.null(ids)) {
		len <- sum(x$finalDecision %in% c('Confirmed', 'Tentative'))
		ind <- (length(lz) - len) : length(lz)
	} else {
		ind <- match(ids, names(lz))
	}
	
	boxplot(lz[ind], xlab = xlab, ylab = ylab, col = col[ind], ...)
	invisible(x)
	names(lz[ind])
}



createROC <- function (pv.list, lab.list, pos.lab='1', file.name='ROC.pdf', width = 6, height = 6) {
	require(ROCR)
	n <- length(pv.list)
	aucs <- numeric(n)
	names(aucs) <- names(pv.list)

#	cols <- scales::hue_pal()(n)
	cols <- rep(c('red', 'blue', 'orange', 'cyan', 'purple'), ceiling(n/5))[1:n]
	ltys <- rep(c(1, 2), ceiling(n/2))[1:n]
	pdf(file.name, height=6, width=6)
	for (i in 1:n) {
		
  		    cat("*")
			pv.mat <- pv.list[[i]]
			lab.mat <- lab.list[[i]]
			
			pred <- prediction(pv.mat, lab.mat==pos.lab)
			perf <- performance(pred, "tpr", "fpr")
			aucs[i] <- mean(unlist(performance(pred, 'auc')@y.values))
			plot(perf, avg="threshold", spread.estimate = 'stddev', spread.scale = 2,  col=cols[i], lty=ltys[i], lwd=1.5,  add=ifelse(i==1, FALSE, TRUE),  main='ROC curve')
				
		}
	abline(0, 1, col='black')
	legend("right", legend=paste0(names(pv.list), "(AUC:", round(aucs, 3), ")"), col=cols, lty=ltys, lwd=2,  bty="n")
	dev.off()
}

#invlogit <- function(x) {
#	rbinom(length(x), 1, exp(x) / (1 + exp(x)))
#}
#pv.list <- list(x1=rnorm(100), x2=rnorm(100))
#lab.list <- list(x1=invlogit(pv.list[['x1']]), x2=invlogit(pv.list[['x2']]))
#createROC(pv.list, lab.list)

# Rev: 2017_02_17 supply aug.var to decision tree
# Rev: 2017_11_08 Change formula
predictionRF <- function (data.obj,  resp.name, formula=NULL, taxa.level='Species', binary=FALSE, prev=0.1, minp=0.002, B=50, seed=123, 
		boruta.level=c('Confirmed', 'Tentative'), ann='',...) {

	#sink(paste0("Taxa_RandomForest_", taxa.level, ".txt"))
	date()
	response <- data.obj$meta.dat[, resp.name]
	
	if (!is.null(formula)) {
		# adj.var <- unlist(strsplit(unlist(strsplit(formula, '\\s*~\\s*'))[2], '\\s*\\+\\s*'))
		# adj <- data.obj$meta.dat[, adj.var, drop=F]
		adj <- model.matrix(as.formula(formula), data.obj$meta.dat)
		adj.var <- colnames(adj)
	} 

	if (taxa.level == 'Species') {
		if (taxa.level %in% names(data.obj$abund.list)) {
			ct <- data.obj$abund.list[[taxa.level]]
		} else {
			# Accomodate different version
			ct <- data.obj$otu.tab
			rownames(ct) <- paste0("OTU", rownames(ct), ":", data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
			data.obj$abund.list[['Species']] <- ct
		}
	} else {
		ct <- data.obj$abund.list[[taxa.level]]
	}

	prop <- t(t(ct) / colSums(ct))
	prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
	prop <- t(prop)
	
	if (binary == TRUE) {
		prop <- (prop != 0)
	}
	
    original.names <- colnames(prop)
	set.seed(seed)
	
	if (is.null(formula)) {
		performance <- matrix(0, nrow=B, ncol=2)
		colnames(performance) <- c("RF_M","Guess")
		roc.list <- list(RF_M=NULL)
		lab.list <- roc.list
	} else {
		performance <- matrix(0, nrow=B, ncol=4)
		colnames(performance) <- c("RF_M", "RF_CF", "RF_M+CF", "Guess")
		roc.list <- list("RF_M"=NULL, "RF_CF"=NULL, "RF_M+CF"=NULL)
		lab.list <- roc.list
	}

	
	colnames(prop) <- gsub(";", "_", colnames(prop))
	colnames(prop) <- gsub(":", "_", colnames(prop))
	colnames(prop) <- gsub("-", "_", colnames(prop))
	colnames(prop) <- gsub("\\.", "_", colnames(prop))
	names(original.names) <- colnames(prop)
	if (!is.null(formula)) {
		names(adj.var) <- adj.var
		original.names <- c(original.names, adj.var)
	}
	
	
	if (!is.null(formula)) {
		padj <- cbind(prop, adj)
	} else {
		padj <- prop
	}
	
	I <- nrow(prop)
	cat('Begin to bootstrap ...\n')
	for(b in 1:B){
		if (b %% 50 == 0) cat (".")
		err <- try({
			bsample <- sample(1:I, I, replace=T)
			if (is.factor(response)) {
				rf1 <- randomForest(x=prop[bsample, ], y=response[bsample], xtest=prop[-bsample, ], ytest=response[-bsample], sampsize=table(response[bsample]), ...)
				performance[b, "RF_M"] <- mean(response[-bsample]!= rf1$test$predicted)
				performance[b, "Guess"]<- mean(response[-bsample] != levels(response)[which.max(tabulate(response[bsample]))])
				roc.list[['RF_M']] <- c(roc.list[['RF_M']], rf1$test$vote[, 1])
				lab.list[['RF_M']] <- c(lab.list[['RF_M']], response[-bsample])
				if (!is.null(formula)) {
					rf2 <- randomForest(x=adj[bsample, ], y=response[bsample], xtest=adj[-bsample, ], ytest=response[-bsample], sampsize=table(response[bsample]), ...)
					rf3 <- randomForest(x=padj[bsample, ], y=response[bsample], xtest=padj[-bsample, ], ytest=response[-bsample], sampsize=table(response[bsample]), ...)
					performance[b, "RF_CF"] <- mean(response[-bsample]!= rf2$test$predicted)
					performance[b, "RF_M+CF"] <- mean(response[-bsample]!= rf3$test$predicted)
					roc.list[['RF_CF']] <- c(roc.list[['RF_CF']], rf2$test$vote[, 1])
					roc.list[['RF_M+CF']] <- c(roc.list[['RF_M+CF']], rf3$test$vote[, 1])
					lab.list[['RF_CF']] <- c(lab.list[['RF_CF']], response[-bsample])
					lab.list[['RF_M+CF']] <- c(lab.list[['RF_M+CF']], response[-bsample])
				}
				
			} else {
				rf1 <- randomForest(x=prop[bsample, ], y=response[bsample], xtest=prop[-bsample, ], ytest=response[-bsample], ...)
				performance[b, "RF_M"] <- mean((response[-bsample] - rf1$test$predicted)^2)
				performance[b, "Guess"]<- mean((response[-bsample] - mean(response[bsample]))^2)
				
				if (!is.null(formula)) {
					
					rf2 <- randomForest(x=adj[bsample, ], y=response[bsample], xtest=adj[-bsample, ], ytest=response[-bsample],  ...)
					rf3 <- randomForest(x=padj[bsample, ], y=response[bsample], xtest=padj[-bsample, ], ytest=response[-bsample],  ...)
					performance[b, "RF_CF"] <- mean((response[-bsample] - rf2$test$predicted)^2)
					performance[b, "RF_M+CF"] <- mean((response[-bsample]- rf3$test$predicted)^2)
				}
			}
		})
		if (inherits(err, "try-error")) {
			next
		}

	}
	cat("\n")
	
	sink(paste0("Taxa_Random_forest_FriedmanTest_", taxa.level, '_', ann, ".txt"))
	if (is.null(formula)) {
		cat("Fridman.test p value: ", friedman.test(performance)$p.value, "\n")
	} else {
		cat("Fridman.test p value (M+CF vs CF): ", friedman.test(performance[, c('RF_M+CF', 'RF_CF')])$p.value, "\n")
	}
	sink()
	
	pdf(paste0("Taxa_Random_forest_misclassification_barplot_", taxa.level, '_', ann, ".pdf"), height=5, width=4)
	if (!is.null(formula)) {
		performance2 <- performance[, c('Guess', 'RF_CF', 'RF_M', 'RF_M+CF')]
	} else {
		performance2 <- performance
	}
	if (is.factor(response)) {
		boxplot(performance2, col="#4DAF4A", ylab='Classification error', las=2)
	} else {
		boxplot(performance2, col="#4DAF4A", ylab='PMSE', las=2)
	}	
	dev.off()
	
	# ROC curve - only compare to the first level of the factor 
    if (is.factor(response)) {
		lab.list <- lapply(lab.list, function(x) {
					y <- factor(x)
					levels(y) <- levels(response)
					y
				}
		)
		createROC(roc.list, lab.list, pos.lab=levels(response)[1], 
				file.name=paste0("Taxa_Random_forest_ROC_", taxa.level, '_', ann, ".pdf"))
	}
	
	if (is.factor(response)) {
		rf <- randomForest(x=padj, y=response, sampsize=table(response), importance=T, ...)
	} else {
		rf <- randomForest(x=padj, y=response,  importance=T, ...)
	}
	
	importance <- as.data.frame(rf$importance)
	
	if (is.factor(response)) {
		write.csv(importance[rev(order(importance$MeanDecreaseAccuracy)), ],
				paste0("Taxa_RandomForest_Ranking_MeanDecreaseAccuracy_", taxa.level,'_', ann, ".csv"))
		write.csv(importance[rev(order(importance$MeanDecreaseGini)), ], 
				paste0("Taxa_RandomForest_Ranking_MeanDecreaseGini_", taxa.level, '_', ann, ".csv"))
	} else {
		write.csv(importance[rev(order(importance[, '%IncMSE'])), ],
				paste0("Taxa_RandomForest_Ranking_IncMSE_", taxa.level, '_', ann, ".csv"))
		write.csv(importance[rev(order(importance$IncNodePurity)), ], 
				paste0("Taxa_RandomForest_Ranking_IncNodePurity_", taxa.level, '_', ann, ".csv"))
	}

	
# Boruta Feature Selection - All subset selection (can't solve the confounding problem)
	obj.Boruta <- Boruta(padj, response, doTrace = 2)	
	write.csv(obj.Boruta$finalDecision, paste0("Taxa_Random_forest_Boruta_Feature_Selection_", taxa.level, '_', ann, ".csv"))
	
	pdf(paste0("Taxa_Random_forest_Boruta_Feature_Selection_", taxa.level, '_', ann,  ".pdf"), height=6, width=10)
	par(mar=par('mar') + c(3, 0, 0, 0))
	plot(obj.Boruta, main = "Feature selection by Boruta", ylab="Importance z-score", lwd = 0.5, las = 3, xlab = "",
			cex=1 / (ncol(prop)/50), cex.axis=0.25*200/ncol(prop), yaxt='n')
	axis(2, cex.axis=1)
	dev.off()
	#sink()
	
	pdf(paste0("Taxa_Random_forest_Boruta_Feature_Selection_Significant_Only_", taxa.level, '_', ann, ".pdf"), height=6, width=6)
	par(mar=par('mar') + c(3, 0, 0, 0))
	otu.ids <- plot.Boruta2(obj.Boruta, main = "Feature selection by Boruta", lwd = 0.5, las = 3, 
			ylab="Importance z-score", xlab = "", cex.axis=1, yaxt='n')
	axis(2, cex.axis=1)
	dev.off()
	
	if ('Confirmed' %in% boruta.level) {
		taxa.names <- original.names[names(obj.Boruta$finalDecision)[obj.Boruta$finalDecision %in% c('Confirmed')]] 
		if (!is.null(formula)) {
			taxa.names <- setdiff(taxa.names, adj.var)
			aug.var <- adj.var
		} else {
			aug.var <- NULL
		}

		if (length(taxa.names) > 0) {
			if (length(taxa.names) > 1) {

				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12,
						margins=c(5, 15), ann=paste0('BorutaFeatures_P_Confirmed_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, data.type='R',
						margins=c(5, 15), ann=paste0('BorutaFeatures_R_Comfirmed_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, data.type='B',
						margins=c(5, 15), ann=paste0('BorutaFeatures_B_Comfirmed_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, Colv=F, dendrogram='row', sam.ord=order(response),
						margins=c(5, 15), ann=paste0('BorutaFeatures_P_Confirmed_Unclusterded_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, data.type='R', Colv=F, dendrogram='row', sam.ord=order(response),
						margins=c(5, 15), ann=paste0('BorutaFeatures_R_Comfirmed_Unclusterded_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,   meta.info=resp.name, width=12, data.type='B', Colv=F, dendrogram='row', sam.ord=order(response),
						margins=c(5, 15), ann=paste0('BorutaFeatures_B_Comfirmed_Unclusterded_', ann))
			}

			if (is.factor(response)) {
				generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
				generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Confirmed_', ann))
				generate_taxa_barplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
				if (length(taxa.names) > 1) {
					generate_taxa_barplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
					generate_taxa_boxplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
					build.decision.tree(data.obj,  resp.name=resp.name, aug.var=aug.var, taxa.level=taxa.level, binary=binary, taxa=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann)) 
				}
				gen <- data.obj$abund.list[[taxa.level]]
				gen <- t(t(gen) / colSums(gen))
				dat <- t(gen[taxa.names, ,drop=FALSE])
				colnames(dat) <- paste0('V', 1:ncol(dat))
				response2 <- factor(-as.numeric(response)+2)
				
				mylr <- function(formula, train, test){
					model <- randomForest(formula, data=train)
					x <- predict(model, newdata=test, type='prob')[, 2]
				}
				
				ACC <- Daim(response2 ~., model=mylr, data=as.data.frame(dat), labpos="1",
						control=Daim.control(method="boot", number=100), cutoff="0.632+")
				pdf(paste0('BorutaFeatures_Confirmed_ROC_', taxa.level, '_0.632+', ann, '.pdf'), height=6, width=6)
				plot(ACC, main='Boruta taxa', method='0.632+', lwd=1.5)	
				abline(0, 1, col='black')
				x <- ACC
				legend("bottomright", legend = paste("AUC:", 
								formatC(Daim::auc(x$roc$sens632p, x$roc$spec632p), 
										digits = max(3, getOption("digits") - 3))),
						inset = 0.01)
				dev.off()	
			} else {
				data.obj$meta.dat[, paste0(resp.name, '_b')] <- factor(data.obj$meta.dat[, resp.name] < median(data.obj$meta.dat[, resp.name]))
				levels(data.obj$meta.dat[, paste0(resp.name, '_b')]) <- c('High', 'Low')
				generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
				generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Confirmed_', ann))
				generate_taxa_barplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
				if (length(taxa.names) > 1) {
					generate_taxa_barplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
					generate_taxa_boxplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
				}
			}
		}
	}

	if ('Tentative' %in% boruta.level) {
		taxa.names <- original.names[names(obj.Boruta$finalDecision)[obj.Boruta$finalDecision %in% c('Confirmed',  'Tentative')]] 
		if (!is.null(formula)) {
			taxa.names <- setdiff(taxa.names, adj.var)
			aug.var <- adj.var
		} else {
			aug.var <- NULL
		}
		
		if (length(taxa.names) > 0) {
			
			if (length(taxa.names) > 1)  {
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12,
						margins=c(5, 15), ann=paste0('BorutaFeatures_P_Tentative_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, data.type='R',
						margins=c(5, 15), ann=paste0('BorutaFeatures_R_Tentative_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, data.type='B',
						margins=c(5, 15), ann=paste0('BorutaFeatures_B_Tentative_', ann))
				
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, Colv=F, dendrogram='row', sam.ord=order(response),
						margins=c(5, 15), ann=paste0('BorutaFeatures_P_Tentative_Unclustered_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, data.type='R', Colv=F, dendrogram='row', sam.ord=order(response),
						margins=c(5, 15), ann=paste0('BorutaFeatures_R_Tentative_Unclustered_', ann))
				generate_taxa_heatmap(data.obj, taxa.levels=taxa.level, taxa=taxa.names,  meta.info=resp.name, width=12, data.type='B', Colv=F, dendrogram='row', sam.ord=order(response),
						margins=c(5, 15), ann=paste0('BorutaFeatures_B_Tentative_Unclustered_', ann))
			}

			if (is.factor(response)) {		
				generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
				generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Tentative_', ann))
				generate_taxa_barplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
				if (length(taxa.names) > 1) {
					generate_taxa_barplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
					generate_taxa_boxplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
					build.decision.tree(data.obj,  resp.name=resp.name, aug.var=aug.var, taxa.level=taxa.level, binary=binary, taxa=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann)) 
				}
				
				gen <- data.obj$abund.list[[taxa.level]]
				gen <- t(t(gen) / colSums(gen))
				dat <- t(gen[taxa.names, ,drop=FALSE])
  			    colnames(dat) <- paste0('V', 1:ncol(dat))
				response2 <- factor(-as.numeric(response)+2)
				
				mylr <- function(formula, train, test){
					model <- randomForest(formula, data=train)
					x <- predict(model, newdata=test, type='prob')[, 2]
				}
				
				ACC <- Daim(response2 ~., model=mylr, data=as.data.frame(dat), labpos="1",
						control=Daim.control(method="boot", number=100), cutoff="0.632+")
				pdf(paste0('BorutaFeatures_Tentative_ROC_', taxa.level, '_0.632+', ann, '.pdf'), height=6, width=6)
				plot(ACC, main='Boruta taxa', method='0.632+', lwd=2)	
				abline(0, 1, col='black')
				x <- ACC
				legend("bottomright", legend = paste("AUC:", 
								formatC(Daim::auc(x$roc$sens632p, x$roc$spec632p), 
										digits = max(3, getOption("digits") - 3))),
						inset = 0.01)
				dev.off()	
				
			} else {
				data.obj$meta.dat[, paste0(resp.name, '_b')] <- factor(data.obj$meta.dat[, resp.name] < median(data.obj$meta.dat[, resp.name]))
				levels(data.obj$meta.dat[, paste0(resp.name, '_b')]) <- c('High', 'Low')
				generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
				generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Tentative_', ann))
				generate_taxa_barplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
				if (length(taxa.names) > 1) {
					generate_taxa_barplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
					generate_taxa_boxplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
				}
			}
		}
	}	
}

col.func <- colorRampPalette(c("blue", "cyan", "red", "yellow"))
#col.fun2 <- colorpanel(9, 'yellow', 'red')

# Rev: 2017_02_19 key alignment modified, rowsep, colsep bug
# Rev: 2017_12_19 row colors fixed
generate_taxa_heatmap <- function (data.obj, taxa.levels='Genus', taxa='All', meta.info, sam.ord=NULL, data.type='P',  prev=0.1, minp=0.002, 
		row.col.dat='Phyla', phy.no=4, sepwidth=0.01, colsep=NULL, rowsep=NULL, sepcolor='black',
		white='white', colFnsC=NULL, colFnsF=NULL, Rowv=T, Colv=T, dendrogram='both', margins=c(5, 15), in.grid=F,  is.labCol=T, cexCol=1, cexRow=NULL,
		omas=c(1, 1, 1, 8), width=12, height=6, ann='All', return.obj=FALSE, ...) {
	
	colsep0 <- colsep
	rowsep0 <- rowsep
	df <- data.obj$meta.dat

	# Determine the col/rowside color
	if (is.null(colFnsC)) {
		# colFnsC <- c(colorRampPalette(c('black', 'yellow', 'red'), colorRampPalette(c('black', 'green')), colorRampPalette(c('black', 'blue'))))
		# https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
		colFnsC <- c(colorRampPalette(colors = brewer.pal(9,"RdBu")), 
				colorRampPalette(colors = brewer.pal(9,"PRGn")), colorRampPalette(colors = brewer.pal(9,"RdYlBu")),  
		        colorRampPalette(colors = brewer.pal(9,"RdGy")), colorRampPalette(colors = brewer.pal(9,"PuOr")))
	} 
	if (is.null(colFnsF)) {
#		rainbow3 <- function(x) {
#			rainbow(x + 2)[1:x]
#		}
		jet3 <- function(x) {
			jet(x + 2)[1:x]
		}
#		colFnsF <- c(rainbow3, jet3)
#		colFnsF <- c(colorRampPalette(colors = brewer.pal(8, "Set1")), colorRampPalette(colors = brewer.pal(7, "Set2")), 
#				colorRampPalette(colors = brewer.pal(7,"Dark2")), colorRampPalette(colors = jet(8)))
		colFnsF <- function (x) {
			if (x <= 6) {
				return(brewer.pal(7, "Set2")[1:x])
			} else {
				return(colorRampPalette(colors = jet(8))(x))
			}
		}
		colFnsF <- c(colFnsF)
	}
	
	
	if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
		data.obj$abund.list[['Species']] <- data.obj$otu.tab
		rownames(data.obj$abund.list[['Species']]) <- paste0("OTU", rownames(data.obj$otu.tab), ":", 
				data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
	}
	for (LOI in taxa.levels) {
		cat(LOI, "\n")
		
		if (LOI == 'All') {
			if (taxa[1] == 'All') {
				stop("Please specify the taxa names that will be included in the heatmap!\n")
			} 
			prop <- NULL
			for (LOI2 in names(data.obj$abund.list)) {
				ct <- data.obj$abund.list[[LOI2]]
				prop0 <- t(t(ct) / colSums(ct))
				prop <- rbind(prop, prop0[intersect(rownames(prop0), taxa), , drop=FALSE])	
	
			}
			colnames(prop) <- colnames(prop0)
			if (nrow(prop) != length(taxa)) {
				warnings('Some taxa not found in abundance lists! Please check the names!\n')
			}
			
  		} else {
			ct <- data.obj$abund.list[[LOI]]
			prop <- t(t(ct) / colSums(ct))
			
			if (taxa[1] != 'All') {
				prop <- prop[taxa, , drop=FALSE]
			} else {
				prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
			}			
		}

		# Sort row and column
		if (is.null(sam.ord)) {
			prop <- prop[order(rownames(prop)), , drop=FALSE]
		} else {
			prop <- prop[order(rownames(prop)), sam.ord, drop=FALSE]
		}
		
		if (is.labCol) {
			labCol <- colnames(prop)
		}  else {
			labCol <- ''
		}
		# Deal with zeros
		if (data.type == 'B') {
			col.scheme <- c("lightyellow", "red")
			prop[, ] <- as.numeric(prop != 0)
			breaks <- c(-0.01, 0.01, 1.1)
		} 
		if (data.type == 'P'){
			#col.scheme <- c(white, colFunc(50))		
			col.scheme = c(white, brewer.pal(11, "Spectral"))
	        # col.scheme = c(white, colorRampPalette(c("blue", "yellow", "red"))(11))
			ind.temp <- prop != 0
			minp <- min(prop[prop!=0])/1.1
			prop[prop==0] <- minp
			prop <- log10(prop)
		#	breaks <- c(log10(minp)-0.01, seq(log10(minp)+0.01, 0, len=51))	
	        breaks <- c(log10(minp)-0.01, seq(log10(minp)+0.01, 0, len=12))
		}
		if (data.type == 'R'){
			col.scheme <- c(white, colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(ncol(prop)-1))
			#col.scheme <- c('white', colorRampPalette(c("green", "black", "red"))(ncol(prop)-1))
	
			prop <- t(apply(prop, 1, function(x) {
								temp <- rank(x[x!=0])
								s <- (ncol(prop) - 1) / (max(temp) - min(temp))
								temp <- 1 + (temp - min(temp)) * s
								x[x!=0] <- temp
								x
								}))
			breaks <- seq(0, ncol(prop), len=ncol(prop)+1)
		}
		
		if (nrow(prop) > 150) {
			labRow <- ''
			cexRow <- 1
		} else {
			labRow <- rownames(prop)
			if (is.null(cexRow)) {
				cexRow <- ifelse(0.5 * 60 / nrow(prop) > 1, 1, 0.5 * 60 / nrow(prop))
			}

		}
		
		if (is.null(colsep0)) {
			if (in.grid == T) {
				colsep <- c(0:ncol(prop))  
			} else {
				colsep <- c(0, ncol(prop))  
			}
		} else {
			colsep <- colsep0
		}
		
		if (is.null(rowsep0)) {
			if (in.grid == T) {
				rowsep <- c(0:nrow(prop))
			} else {
				rowsep <- c(0, nrow(prop))
			}
		} else {
			rowsep <- rowsep0
		}
		
		key.list <- list()
		colsidecol <- NULL
		i <- 0
		j <- 0
		for (keyID in meta.info) {
			if (is.null(sam.ord)) {
				x <- df[, keyID]
			} else {
				x <- df[sam.ord, keyID]
			}		

			if (is.factor(x) | is.character(x)) {
				x <- factor(x)
				key.list[[keyID]] <- list(breaks=levels(x), colors=(colFnsF[i+1][[1]])(nlevels(x)), base=NA, col.na=NA, right=F, include.lowest=F)	
				i <- (i + 1) %% length(colFnsF)
				colsidecol <- cbind(colsidecol, key.list[[keyID]]$colors[x])
			} else {
				key.list[[keyID]] <- makecmap(x, n=5, colFn=colFnsC[j+1][[1]])
				j <- (j + 1) %% length(colFnsC)
				colsidecol <- cbind(colsidecol, cmap(x, key.list[[keyID]]))
			}
		}
		colnames(colsidecol) <- meta.info
		
		# Rev: 2016_09_13
		colsidecol[is.na(colsidecol)] <- 'white'
		
		# add aunbdance key 
		if (data.type == 'B') {
			prop.cmap <- list(breaks=c("Absence", "Presence"), colors=c("lightyellow", "red"), base=NA, col.na=NA, right=F, include.lowest=F) 
			KeyName <- 'Abundance'
		} 
#		if (data.type == 'P'){
#			prop.cmap <- makecmap(prop[ind.temp], n = 5, colFn = colFunc)
#		}
#		if (data.type == 'R'){
#			prop.cmap <- makecmap(as.vector(prop), n = 5, colFn = colFunc)
#		}
		
		if (data.type == 'P'){
			KeyName <- 'log(Proportion)'	
		}
		if (data.type == 'R'){
			KeyName <- 'Rank'
		}
			
		# add row col key
        # Strong assumption of the structure of 'row' names
        if (row.col.dat == 'Phyla') {
			if (LOI %in% c('Class', 'Order', 'Family', 'Genus', 'Species')) {
				if (LOI == 'Species') {
					phy <- sapply(strsplit(rownames(prop), ":"), function(x) x[2])
					phy <- sapply(strsplit(phy, ";"), function(x) x[1])
				} else {
					phy <- sapply(strsplit(rownames(prop), ";"), function(x) x[1])
				}
				
				if (sum(is.na(phy)) == 0) {
					temp <- sort(table(phy), decr=T)
					if (length(temp) > phy.no) {
						rare.phy <- names(temp)[-(1:phy.no)]
						phy[phy %in% rare.phy] <- 'Other'
						phy <- factor(phy, levels=c(names(temp)[(1:phy.no)], 'Other'))
					} else {
						phy <- factor(phy)
					}
					rowsidecol<- rainbow(nlevels(phy))[phy]
					rowsidecol <- rbind(rowsidecol, rowsidecol)
					rownames(rowsidecol) <- c('', '')
					phy.cmap <- list(breaks=levels(phy), colors=rainbow(nlevels(phy)), base=NA, col.na=NA, right=F, include.lowest=F) 
				} else {
					rowsidecol <- NULL
				}

			} else {
				rowsidecol <- NULL
			}
		} else {
			rowsidecol <- NULL
		}
		
		pdf(paste0('Taxa_Heatmap_', LOI, '_', ann, '.pdf'), width=width, height=height)
		par(oma = omas)

#		if (data.type == 'R' | data.type == 'B') {
#			dist2 <- dist
#		}
#		if (data.type == 'P') {
#			# Better clustering of taxa
##			dist2 <- function(x) {    
##				as.dist((1-cor(t(x)))/2)
##			}
#	       dist2 <- dist
#		}
#
#		# Rev: 2016_09_13 handle zero sd cases
#        if (sum(rowSds(prop) == 0) != 0 | sum(colSds(prop) == 0) != 0)  {
#			dist2 <- dist
#			warning('Zero sd produced! Euclidean distance is used instead!\n')
#		}
     
		# Pearson correlation distance
		obj <- heatmap.3(prop, 
				Rowv=Rowv, 
				Colv=Colv, 
#				distfun = dist2,
				dendrogram=dendrogram,
				scale='none',
				col=col.scheme, 
				breaks=breaks, 
				symbreaks=F,
				trace='none',
				margins= margins, 
				colsep = colsep,  
				rowsep = rowsep,  
				sepcolor= sepcolor, 
				sepwidth=c(sepwidth, sepwidth),
				ColSideColors=colsidecol,
				RowSideColors=rowsidecol,
				cexRow=cexRow,
				labRow=labRow,
				labCol=labCol,
				cexCol=cexCol,
				key=(data.type != 'B'), density.info='none', symkey=F, KeyValueName=KeyName,	
				NumColSideColors= 0.5 *length(meta.info),
				NumRowSideColors= 0.5,
				...
		)

		par(cex=0.75)
		par(oma=c(0, 0, 1, 0)) 
		
		if (!is.null(rowsidecol) & row.col.dat == 'Phyla') {
			y.cord <- (1/(length(meta.info)+2)) * (0:((length(meta.info)+2) - 1))
			vkey2(phy.cmap, 'Phylum', x=0, y=-0.2, stretch=1.2)

		}

		y.cord <- (1/(length(meta.info))) * (0:(length(meta.info) - 1))
		k <- 1
		xs <- c(0.95, 1)
		for (keyID in meta.info) {
			x <- df[, keyID]
			if (is.factor(x) | is.character(x)) {
				x <- factor(x)
				vkey2(key.list[[keyID]], keyID, x=xs[(k-1) %% 2 + 1], y=y.cord[k], stretch=1.2)
			} else {
				vkey(key.list[[keyID]], keyID, x=xs[(k-1) %% 2 + 1], y=y.cord[k], stretch=1.2)
			}
			k <- k + 1
		}	
		
		if (data.type == 'B') {
			vkey2(prop.cmap, KeyName,  x=xs[(k-1) %% 2 + 1], y=1, stretch=1.2)		
		} 
#		if (data.type == 'P'){
#			vkey(prop.cmap, 'log10(Proportion)', y=1, stretch=1.2)
#		}	
#		if (data.type == 'R'){
#			vkey(prop.cmap, 'Rank', y=1, stretch=1.2)
#		}	
		dev.off()	
	}
	if (return.obj == TRUE) {
		return(obj)
	}
	
}

# New: 2017_03_02  General heatmap without taxa
# data.obj$data, data.obj$meta.dat
# Rev: 2017_11_20, Add 'LogC', taxa.as.row
generate_heatmap <- function (data.obj,  meta.info, sam.ord=NULL, data.type='P',
		row.col.dat=NULL,  sepwidth=0.01, colsep=NULL, rowsep=NULL, taxa.as.row = TRUE,
		colFnsC=NULL, colFnsF=NULL, Rowv=T, Colv=T, dendrogram='both', margins=c(5, 15),
		in.grid=F, sepcolor='black', is.labCol=T, cexCol=1, cexRow=NULL,
		omas=c(1, 1, 1, 8), width=12, height=6, ann='All', return.obj=FALSE, pdf=TRUE, ...) {
	
	colsep0 <- colsep
	rowsep0 <- rowsep
	df <- data.obj$meta.dat
	
	# Determine the col/rowside color
	# Determine the col/rowside color
	if (is.null(colFnsC)) {
		# colFnsC <- c(colorRampPalette(c('black', 'yellow', 'red'), colorRampPalette(c('black', 'green')), colorRampPalette(c('black', 'blue'))))
		# https://moderndata.plot.ly/create-colorful-graphs-in-r-with-rcolorbrewer-and-plotly/
		colFunsC <- c(colorRampPalette(colors = brewer.pal(9,"RdBu")), 
				colorRampPalette(colors = brewer.pal(9,"PRGn")), colorRampPalette(colors = brewer.pal(9,"RdYlBu")),  
				colorRampPalette(colors = brewer.pal(9,"RdGy")), colorRampPalette(colors = brewer.pal(9,"PuOr")))
	} 
	if (is.null(colFnsF)) {

		colFnsF <- function (x) {
			if (x <= 6) {
				return(brewer.pal(7, "Set2")[1:x])
			} else {
				return(colorRampPalette(colors = jet(8))(x))
			}
		}
		colFnsF <- c(colFnsF)
	}
	prop <- data.obj$data
	
	# Sort row and column
	if (!is.null(sam.ord)) {
		prop <- prop[, sam.ord, drop=FALSE]
	}
	
	if (is.labCol) {
		labCol <- colnames(prop)
	}  else {
		labCol <- ''
	}
	if (data.type == 'LogC') {
		col.scheme <- c('white', colorpanel(9, 'aliceblue', 'lightblue', 'blue'))
		breaks <- c(-0.1, 0.1, 2, 4, 6, 8, 10, 12, 14, 16, 18)
	}
	# Deal with zeros
	if (data.type == 'B') {
		col.scheme <- c("lightyellow", "red")
		breaks <- c(-0.01, 0.01, 1.1)
	} 
	if (data.type == 'P'){
		
		col.scheme = brewer.pal(11, "Spectral")
		breaks <- c(seq(min(prop), max(prop), len=12))
	}
	if (data.type == 'R'){
		col.scheme <- c(white, colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(ncol(prop)-1))
		breaks <- seq(0, ncol(prop), len=ncol(prop)+1)
	}
	
	if (nrow(prop) > 1500) {
		labRow <- ''
		cexRow <- 1
	} else {
		labRow <- rownames(prop)
		if (is.null(cexRow)) {
			cexRow <- ifelse(0.5 * 60 / nrow(prop) > 1, 1, 0.5 * 60 / nrow(prop))
		}
		
	}
	
	if (is.null(colsep0)) {
		if (in.grid == T) {
			colsep <- c(0:ncol(prop))  
		} else {
			colsep <- c(0, ncol(prop))  
		}
	} else {
		colsep <- colsep0
	}
	
	if (is.null(rowsep0)) {
		if (in.grid == T) {
			rowsep <- c(0:nrow(prop))
		} else {
			rowsep <- c(0, nrow(prop))
		}
	} else {
		rowsep <- rowsep0
	}
	
	key.list <- list()
	colsidecol <- NULL

	i <- 0
	j <- 0
	
	for (keyID in meta.info) {
		if (is.null(sam.ord)) {
			x <- df[, keyID]
		} else {
			x <- df[sam.ord, keyID]
		}		
		if (is.factor(x)) {
			key.list[[keyID]] <- list(breaks=levels(x), colors=(colFnsF[i+1][[1]])(nlevels(x)), base=NA, col.na=NA, right=F, include.lowest=F)	
			i <- (i + 1) %% length(colFnsF)
			colsidecol <- cbind(colsidecol, key.list[[keyID]]$colors[x])
		} else {
			key.list[[keyID]] <- makecmap(x, n=5, colFn=colFnsC[j+1][[1]])
			j <- (j + 1) %% length(colFnsC)
			colsidecol <- cbind(colsidecol, cmap(x, key.list[[keyID]]))
		}
	}
	colnames(colsidecol) <- meta.info
	
	# Rev: 2016_09_13
	colsidecol[is.na(colsidecol)] <- 'white'
	
	# add aunbdance key 
	if (data.type == 'B') {
		prop.cmap <- list(breaks=c("0", "1"), colors=c("lightyellow", "red"), base=NA, col.na=NA, right=F, include.lowest=F) 
		KeyName <- 'Value'
	} 
	
	if (data.type == 'P'){
		KeyName <- 'Value'	
	}
	if (data.type == 'R'){
		KeyName <- 'Rank'
	}
	
	if (data.type == 'LogC'){
		KeyName <- 'Log2(count+1)'
	}
	
	# add row col key
	if (!is.null(row.col.dat)) {
		
		phy <- df[, row.col.dat]
		phy <- factor(phy)
		
		rowsidecol<- rainbow(nlevels(phy))[phy]
		rowsidecol <- rbind(rowsidecol, rowsidecol)
		rownames(rowsidecol) <- c('', '')
		phy.cmap <- list(breaks=levels(phy), colors=rainbow(nlevels(phy)), base=NA, col.na=NA, right=F, include.lowest=F) 
	} else {
		rowsidecol <- NULL
	}
	if(pdf == TRUE) {
		pdf(paste0('Heatmap_', ann, '.pdf'), width=width, height=height)
	}
	
	par(oma = omas)
	
	
	# Pearson correlation distance
    if (taxa.as.row) {
		prop <- prop
	} else {
		prop <- t(prop)
	}
	obj <- heatmap.3(prop, 
			Rowv=Rowv, 
			Colv=Colv, 
#				distfun = dist2,
			dendrogram=dendrogram,
			scale='none',
			col=col.scheme, 
			breaks=breaks, 
			symbreaks=F,
			trace='none',
			margins= margins, 
			colsep = colsep,  
			rowsep = rowsep,  
			sepcolor= sepcolor, 
			sepwidth=c(sepwidth, sepwidth),
			ColSideColors=colsidecol,
			RowSideColors=rowsidecol,
			cexRow=cexRow,
			labRow=labRow,
			labCol=labCol,
			cexCol=cexCol,
			key=(data.type != 'B'), density.info='none', symkey=F, KeyValueName=KeyName,	
			NumColSideColors= 0.5 *length(meta.info),
			NumRowSideColors= 0.5,
			...
	)
	
	par(cex=0.75)
	par(oma=c(0, 0, 1, 0)) 
	
	if (!is.null(rowsidecol)) {
		y.cord <- (1/(length(meta.info)+2)) * (0:((length(meta.info)+2) - 1))
		vkey2(phy.cmap, rowsidecol, x=0, y=-0.2, stretch=1.2)
		
	}
	
	y.cord <- (1/(length(meta.info))) * (0:(length(meta.info) - 1))
	k <- 1
	xs <- c(0.95, 1)
	for (keyID in meta.info) {
		x <- df[, keyID]
		if (is.factor(x)) {
			vkey2(key.list[[keyID]], keyID, x=xs[(k-1) %% 2 + 1], y=y.cord[k], stretch=1.2)
		} else {
			vkey(key.list[[keyID]], keyID, x=xs[(k-1) %% 2 + 1], y=y.cord[k], stretch=1.2)
		}
		k <- k + 1
	}	
	
	if (data.type == 'B') {
		vkey2(prop.cmap, KeyName, x=xs[(k-1) %% 2 + 1], y=1, stretch=1.2)		
	} 
	
	if (pdf == TRUE) {
		dev.off()	
	}
	
	
	if (return.obj == TRUE) {
		return(obj)
	}
}



# Rev: 2016_06_13
# Still need to be revised: Multiple groups, the auto order should look similar for differnt groups;
# The legends do not work
# Rev: 2017_04_28 add col.map so the color scheme can be the same for different data sets
# Rev: 2017_11_09 add label for individual bar plot and change ordering 
generate_stacked_barplot <- function(data.obj, grp.name=NULL, taxa.levels=c('Phylum', 'Family', 'Genus'), agg.cutoff=0.005, 
		border=TRUE, order.auto=FALSE, order.method='abund', order.taxa=FALSE, cex.names=0.5, cex.top.lab = 0.75, separate=FALSE, cex.names2=1, indiv=TRUE, aggre=TRUE,
		hei1=6, wid1=9, hei2=6, wid2=9, margin=10, ann='', pdf=TRUE, col.map = NULL,...) {
	if (is.null(grp.name)) {
		grp <- rep(1, nrow(data.obj$meta.dat))
	} else {
		grp <- factor(data.obj$meta.dat[, grp.name])
	}
	
	name.list <- list()
	col.list <- list()
	
	# rearrange the order of taxa
	if (order.taxa) {
		for (taxa.level in taxa.levels) {
			abund0 <- data.obj$abund.list[[taxa.level]]
			abund0 <- t(t(abund0) / colSums(abund0))
			data.obj$abund.list[[taxa.level]] <- data.obj$abund.list[[taxa.level]][rev(order(rowMeans(abund0))), ]
		}
	}

	# create color mapping
	for (taxa.level in taxa.levels) {
		abund0 <- data.obj$abund.list[[taxa.level]]
		abund0 <- t(t(abund0) / colSums(abund0))
		
		abund1 <- abund0[rowMeans(abund0) >= agg.cutoff, , drop=F]
		abund2 <- abund0[rowMeans(abund0) < agg.cutoff, , drop=F]
		
		prop <- rbind(abund1, Other=colSums(abund2))
		colnames(prop) <- colnames(abund0)
		name.list[[taxa.level]] <- rownames(abund1)
		
		#rand.col <- c(sample(rainbow(nrow(prop)*2), nrow(prop)-1), 'gray')
		# Rev: 2017_04_28
		if (is.null(col.map[[taxa.level]])) {
			rand.col <- c(rep_len(brewer.pal(12, "Paired"), nrow(prop)-1), 'gray')
			names(rand.col) <- rownames(prop)
 			col.list[[taxa.level]] <- rand.col
		} else {
			if (sum(!(rownames(prop)%in% names(col.map[[taxa.level]]) ))) {
				stop('Color mapping not correct! Please check!\n')
			} else {
				col.list[[taxa.level]] <- col.map[[taxa.level]]
			}
			
		}

	}	
	
	if (indiv == TRUE) {
		if (pdf == TRUE)
			pdf(paste0("Taxa_Stacked_Barplot_Overall_Compo_", ann, ".pdf"), height=hei1, width=wid1)
		
		if (border == FALSE) {
			lty.o <- par("lty")
			par(lty = 0)
		}
		
		mar.o <- par(mar=par('mar') + c(0, margin, 0, 0))
		
		for (taxa.level in taxa.levels) {
			abund0 <- data.obj$abund.list[[taxa.level]]
			abund0 <- t(t(abund0) / colSums(abund0))
			
#			abund1 <- abund0[rowMeans(abund0) >= agg.cutoff, , drop=F]
#			abund2 <- abund0[rowMeans(abund0) < agg.cutoff, , drop=F]
			
			abund1 <- abund0[name.list[[taxa.level]], , drop=FALSE]
			abund2 <- 1 - colSums(abund1)
			
			prop <- rbind(abund1, Other=abund2)
			colnames(prop) <- colnames(abund0)
		
			#rand.col <- c(sample(rainbow(nrow(prop)*2), nrow(prop)-1), 'gray')
			#col.list[[taxa.level]] <- rand.col
			
			cex.legend = ifelse (nrow(prop) > 35, 35/nrow(prop)*0.75, 0.75)
			
			# Rev: 2016_09_13, better ordering within group based on single linkage
		    # Rev: 2017_11_10, ordering rule the same across group
			if (order.auto) {
				
				temp <- rev(order(rowMeans(prop)))
	
				prop <- do.call(cbind, tapply(1:length(grp), factor(grp), function(i){
									if (length(i) >= 2) {
										prop.sub <- prop[, i, drop=FALSE]
										if (order.method == 'single') {
											ord <- hclust(vegdist(t(prop.sub)), method='single')$order
										} else {
											prop.dummy <- prop.sub
											prop.dummy[prop.dummy < 0.05] <- 0
											if (length(temp) >=3) {
												ord <- (order(-prop.dummy[temp[1], ], -prop.dummy[temp[2], ], -prop.dummy[temp[3], ]))
											} else {
												ord <- (order(-prop.dummy[temp[1], ]))
											}
										}
										# temp <- rev(order(rowMeans(prop.sub)))
										prop.sub <- prop.sub[, ord, drop=FALSE]
										return(prop.sub)
									} else {
										return(prop[, i, drop=FALSE])
									}
								}))
			} else {
				prop <- prop[, order(grp), drop=FALSE]
			}
			
			grp2 <- sort(grp)
			if (border == FALSE) {
				if (is.null(grp.name)) {
					barplot(prop, col=col.list[[taxa.level]][rownames(prop)], ylab='Proportion', las=2, legend.text=rownames(prop), cex.names=cex.names, space=0,
							args.legend=list(x='left', bty='n',  cex=cex.legend, inset=c(-0.5, 0)), main=taxa.level, ...)
				}  else {
		
					coord <- barplot(prop, ylim = c(0, 1.05), col=col.list[[taxa.level]][rownames(prop)], ylab='Proportion', las=2, legend.text=rownames(prop), cex.names=cex.names, space=0,
							args.legend=list(x='left', bty='n',  cex=cex.legend, inset=c(-0.5, 0)), main=taxa.level, ...)
				
					coord.text <- tapply(coord, grp2, mean)
					coord.st <- tapply(coord, grp2, max)[-(length(unique(grp2)))]
					coord.end <- tapply(coord, grp2, min)[-1]
					coord.line <- (coord.st +  coord.end) / 2
	
					text.names <- names(coord.text)
	

					for (i in 1 : length(coord.line)) {
						abline(v = coord.line[i], lty = 1, lwd = 0.5)	
					}
					for (i in 1 : length(coord.text)) {
						if (i %% 2 == 0) {
							pos <- 3
						} else {
							pos <- 3
						}
						text(coord.text[i], 1, text.names[i], pos = pos, cex = cex.top.lab)
					}
					
				}
			} else {
				if (is.null(grp.name)) {
					barplot(prop, col=col.list[[taxa.level]][rownames(prop)], ylab='Proportion', las=2, legend.text=rownames(prop), cex.names=cex.names,
							args.legend=list(x='left', bty='n',  cex=cex.legend, inset=c(-0.5, 0)), main=taxa.level, ...)
				} else {
					coord <- barplot(prop, ylim = c(0, 1.05), col=col.list[[taxa.level]][rownames(prop)], ylab='Proportion', las=2, legend.text=rownames(prop), cex.names=cex.names,
							args.legend=list(x='left', bty='n',  cex=cex.legend, inset=c(-0.5, 0)), main=taxa.level, ...)
					
					coord.text <- tapply(coord, grp2, mean)
					coord.st <- tapply(coord, grp2, max)[-(length(unique(grp2)))]
					coord.end <- tapply(coord, grp2, min)[-1]
					coord.line <- (coord.st +  coord.end) / 2
					
					text.names <- names(coord.text)
					for (i in 1 : length(coord.line)) {
						abline(v = coord.line[i], lty = 1, lwd = 0.5)	
					}
					for (i in 1 : length(coord.text)) {
						if (i %% 2 == 0) {
							pos <- 3
						} else {
							pos <- 3
						}
						text(coord.text[i], 1, text.names[i], pos = pos, cex = cex.top.lab)
					}
				}
				
				
			}
			
		}
		if (border == FALSE) {
			par(lty = lty.o)
		}
		
		par(mar = mar.o)
		if (pdf == TRUE)
			dev.off()
		
	}

	# Generate averaged over stack barplot

	if (!is.null(grp.name) & aggre == TRUE) {
		if (separate != TRUE) {
			if (pdf == TRUE)
			pdf(paste0("Taxa_Stacked_Barplot_Grouped_Compo_", ann, "_combine.pdf"), height=hei2, width=wid2)
			mar.o <- par(mar=par('mar') + c(0, margin/length(taxa.levels), 0, 0))
			mfrow.o <- par(mfrow=c(1, length(taxa.levels)))
			cex.legend <- 10
			for (taxa.level in taxa.levels) {
				abund0 <- data.obj$abund.list[[taxa.level]]
				abund0 <- t(t(abund0) / colSums(abund0))
				
				# Time consuming for tapply
#				abund0 <- t(apply(abund0, 1, function(x) {
#									tapply(x, grp, mean)
#								}))
				
			    temp <- aggregate(t(abund0), by=list(grp), mean)
				abund0 <- as.matrix(temp[, -1])
				rownames(abund0) <- temp[, 1]
				abund0 <- t(abund0)
				
				abund1 <- abund0[name.list[[taxa.level]], , drop=FALSE]
				abund2 <- 1 - colSums(abund1)
				
				prop <- rbind(abund1, Other=abund2)
				colnames(prop) <- colnames(abund0)
				
				newsize <- ifelse (nrow(prop) > 35, 35/nrow(prop)*0.75, 0.75)
				cex.legend <- ifelse(cex.legend < newsize, cex.legend, newsize) 
				#prop <- prop[, order(grp)]
				barplot(prop, col=col.list[[taxa.level]][rownames(prop)], ylab='Proportion', las=2, cex.names=cex.names2, main=taxa.level, ...)
				#		legend.text=rownames(prop), args.legend=list(x='left', bty='n',  cex=cex.legend, inset=c(-2.2, 0)))
				
			}
			
            par(mar=c(0, 0, 0, 0))
			oma.o <- par(oma=c(0, 0, 0, 0))
			for (taxa.level in taxa.levels) {
				abund0 <- data.obj$abund.list[[taxa.level]]
				abund0 <- t(t(abund0) / colSums(abund0))
				
				abund0 <- t(apply(abund0, 1, function(x) {
									tapply(x, grp, mean)
								}))
				
				abund1 <- abund0[name.list[[taxa.level]], , drop=FALSE]
				abund2 <- 1 - colSums(abund1)
				
				prop <- rbind(abund1, Other=abund2)
				colnames(prop) <- colnames(abund0)
				
				#	cex.legend = ifelse (nrow(prop) > 35, 35/nrow(prop)*0.75, 0.75)
				plot(1, type="n", axes=FALSE, xlab="", ylab="")
				# Rev: 2016_09_10
				legend('left', legend=rev(rownames(prop)), bty='n', fill=rev(col.list[[taxa.level]][rownames(prop)]), cex=cex.legend)
				
			}
			par(mar=mar.o)
			par(oma=oma.o)
			par(mfrow=mfrow.o)
			if (pdf == TRUE)
			dev.off()
		} else {
			if (pdf == TRUE)
			pdf(paste0("Taxa_Stacked_Barplot_Grouped_Compo_", ann, "_Separate.pdf"), height=hei2, width=wid2)
		    
			mar.o <- par(mar=par('mar') + c(0, margin, 0, 0))
			for (taxa.level in taxa.levels) {
				abund0 <- data.obj$abund.list[[taxa.level]]
				abund0 <- t(t(abund0) / colSums(abund0))
				
				abund0 <- t(apply(abund0, 1, function(x) {
									tapply(x, grp, mean)
								}))
				
				abund1 <- abund0[name.list[[taxa.level]], , drop=FALSE]
				abund2 <- 1 - colSums(abund1)
				
				prop <- rbind(abund1, Other=abund2)
				colnames(prop) <- colnames(abund0)
				
				cex.legend <- ifelse (nrow(prop) > 35, 35/nrow(prop)*0.75, 0.75)
				barplot(prop, col=col.list[[taxa.level]][rownames(prop)], ylab='Proportion', las=2, cex.names=cex.names2, main=taxa.level,
						legend.text=rownames(prop), args.legend=list(x='left', bty='n',  cex=cex.legend, inset=c(-0.5, 0)), ...)
				
			}
			par(mar=mar.o)
			if (pdf == TRUE)
			dev.off()
			
		}

	}

}

build.decision.tree <- function(data.obj,  resp.name, taxa.level='Species', binary=FALSE, taxa, aug.var=NULL, ann='All') {
	ann <- paste(taxa.level, ann, sep="_")
	response <- data.obj$meta.dat[, resp.name]
	
	ct <- data.obj$abund.list[[taxa.level]]
	prop <- t(t(ct) / colSums(ct))
	prop <- prop[taxa, , drop=F] 	
	if (binary == TRUE) {
		prop <- (prop != 0)
	}
	
	dat <- as.data.frame(t(prop))
	# Rev: 2017_02_17 Add additional variables from meta dat
	if (!is.null(aug.var)) {
		dat <- cbind(dat, data.obj$meta.dat[, aug.var])
	}
	
	dat <- data.frame(dat, response)
	try(
		if (is.factor(response)) {
			fit <- rpart(response ~ ., method="class", data=dat)
			post(fit, file = paste0("Taxa_Unpruned_Classification_tree_", ann, ".ps"), title = "Unpruned Classification Tree")
			pfit<- prune(fit, cp= fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
			post(pfit, file = paste0("Taxa_Pruned_Classification_", ann, ".ps"), title = "Pruned Classification Tree")
			
		} else {
			fit <- rpart(response ~ ., method="anova", data=dat)
			post(fit, file = paste0("Taxa_Unpruned_Regression_tree_", ann, ".ps"), title = "Unpruned Regression Tree")
			pfit<- prune(fit, cp= fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
			post(pfit, file = paste0("Taxa_Pruned_Regression_tree_", ann, ".ps"), title = "Pruned RegressionTree")
		}
	)	
}

bootstrap.pwr <- function(pcs, formula, dat, ns=NULL, perm.no=199, iter.no=100) {
	if (is.null(ns)) {
		ns <- nrow(dat)*seq(1, 4, len=8)
	}
	pvs <- numeric(length(ns))
	names(pvs) <- paste(ns)
	for (n in ns) {
		cat('.')
		temp <- sapply(1:iter.no, function(i) {
			bt.ind <- sample(1:nrow(dat), n, repl=T)
			dat.bt <- dat[bt.ind, ]
			dist.bt <- dist(pcs[bt.ind, ])
			aov.tab <- adonis(as.formula(paste('dist.bt', formula)), dat=dat.bt, permutations=perm.no)$aov.tab
	        pv <- aov.tab[nrow(aov.tab)-2, ncol(aov.tab)]
		})
        pvs[paste(n)] <- mean(temp <= 0.05)
	}

	return(pvs)
}


perform_power_analysis <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), formula=NULL,  
		grp.name=NULL, adj.name=NULL, ann='', ...) {
	if (is.null(formula)) {
		if (is.null(adj.name)) {
			formula <- paste('~', grp.name)
		} else {
			formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
		}
	}	 
	df <- data.obj$meta.dat
    pvm <- NULL
	for (dist.name in dist.names) {
		cat('*')
		dist.mat <- dist.obj[[dist.name]]
		pcs <- cmdscale(dist.mat, k=nrow(dist.mat)-1)
		pvs <- bootstrap.pwr(pcs, formula, df, ...)
		pvm <- rbind(pvm, pvs)
	}
	colnames(pvm) <- names(pvs)
	rownames(pvm) <- dist.names
	
	pvdf <- melt(pvm)
	colnames(pvdf) <- c('Distance_type', 'Sample_size', 'Value')
	pdf(paste0("Power_curve_bootstrap_", ann, '.pdf')) 
	g.obj <- ggplot(pvdf, aes(x=Sample_size, y=Value)) +
				geom_point() +
				geom_line() +
				ylab('Power') +
				facet_wrap(~ Distance_type, ncol=2) +
				theme_bw()
	print(g.obj)
	dev.off()
	
	return(pvdf)
}

# More robust identification
k.nonoverlap.se <- function (tab, SE.factor=1) {
	max.v <- which.max(tab[, 'gap']) 
	
	for (i in (max.v-1):1) {
		if (tab[i, 'gap'] + SE.factor * tab[i, 'SE.sim'] < tab[i+1, 'gap'] - SE.factor * tab[i+1, 'SE.sim']) {
			break
		} 
	}
	
	return(i+1)
}


perform_cluster_analysis <- function (data.obj, dist.obj, dist.name=c('UniFrac'), k.best=NULL, method='pam', stat='gap', 
		grp.name=NULL, adj.name=NULL, subject=NULL, ann='', seed=1234) {
	
	df <- data.obj$meta.dat
	if (!is.null(grp.name)) {
		grp <- df[, grp.name]
	}
	if (!is.null(adj.name)) {
		adj <- df[, adj.name]
	}
	cat(dist.name, ' distance ...\n')
	mat <- dist.obj[[dist.name]]
	
	if (is.null(k.best)) {
		pc.obj <- cmdscale(as.dist(mat), k=ncol(mat) - 1, eig=T)
		eig <- pc.obj$eig
		eig <- eig[eig > 0]
		pvar <- eig / sum(eig)
		pc <- pc.obj$points
		
		cat('Assess cluster number by gap statistics ...\n')
		
		set.seed(seed)
		gs <-  gapstat_ord(pc, axes=1:ncol(pc), verbose=FALSE)
		plot_clusgap(gs)
		ggsave(paste0('cluster_assess_gap_statistic_', dist.name, ann, '.pdf') , width=6, height=6)
		
		
#		print(gs, method="Tibs2001SEmax")
#		print(gs, method="firstSEmax")
		
		tab <- gs$Tab
		write.csv(tab, paste0('gap_stat_', dist.name, '.csv'))
		
		k.best.gap <- k.nonoverlap.se(tab)
		
		cat('Gap statistic finds ', k.best.gap, ' clusters.\n')
		
		
		cat('Assess cluster number by average silhouette width ...\n')
		asw <- numeric(20)
		for (k in 2:20){
			asw[k] <- pam(as.dist(mat), k)$silinfo$avg.width
		}
		
		k.best.asw <- which.max(asw)
		cat("silhouette-optimal number of clusters:", k.best.asw, "\n")
		pdf(paste0('cluster_assess_asw_statistic_', dist.name, ann, '.pdf'), width=6, height=6)
		plot(1:20, asw, type= "h", main = "pam() clustering assessment",
				xlab= "k (# clusters)", ylab = "average silhouette width")
		axis(1, k.best.asw, paste("best",k.best.asw, sep="\n"), col = "red", col.axis = "red")
		dev.off()
		
		cat('ASW statistic finds ', k.best.asw, ' clusters.\n')
		
		if (stat == 'gap') {
			k.best <- k.best.gap
		} else {
			if (stat == 'asw') {
				k.best <- k.best.asw
			}
		}
		
	}
	
	if (k.best == 1) {
		cat('No robust cluster structure found!\n')
	} else {
		pam.obj <- pam(as.dist(mat), k=k.best)
		pam.class <- factor(pam.obj$clustering)
		write.csv(pam.class, paste0('Cluster.membership', dist.name, '.csv'))
		
		df$Cluster <- pam.class
		data.obj$meta.dat <- df
		if (is.null(grp.name)) {
			generate_ordination(data.obj, dist.obj, dist.name, grp.name='Cluster', ann=paste0('Cluster', dist.name, ann))
		} else {
			generate_ordination(data.obj, dist.obj, dist.name, grp.name='Cluster', strata=grp.name, ann=paste0('Cluster', dist.name, ann))
		}
		
		# Define cluster characteristics
	   diff.obj <- perform_differential_analysis(data.obj, grp.name='Cluster', 
			   taxa.levels=c('Genus'),
			   method='kruskal', mt.method='fdr', 
			   cutoff=0.01, prev=0.1, minp=0.002, ann='Cluster')
	   
	   visualize_differential_analysis(data.obj, diff.obj, grp.name='Cluster', taxa.levels=c('Genus'), indivplot=FALSE,
			   mt.method='fdr', cutoff=0.01, ann='Cluster')
		
		try(
				if (!is.null(grp.name)) {
					if (is.null(adj.name)) {
						form <- as.formula(paste('yy ~', grp.name))
					} else {
						form <- as.formula(paste('yy ~', paste(adj.name, collapse='+'), '+', grp.name))		
					}
					# Enrichment analysis - logistic regression - 1 vs other
					cat('Testing for association with the clusters - 1 vs other ...\n')
					sink(paste0('Cluster_association_test', dist.name, ann, '.txt'))
					if (is.null(subject)) {
						cat('Generalized linear model (logistic regression) is performed.\n')
						
						for (clus in levels(df$Cluster)[1:(nlevels(df$Cluster))]) {	
							cat('Test for enrichment in cluster', clus, '\n')
							y <- rep(0, length(df$Cluster)) 
							y[df$Cluster == clus] <- 1
							df$yy <- y
							prmatrix(summary(glm(form, data=df, family=binomial))$coefficients)	
						}
					} else {
						cat('Generalized linear mixed effects model (logistic regression) is performed.\n')
						for (clus in levels(df$Cluster)[1:(nlevels(df$Cluster))]) {
							cat('Test for enrichment in cluster', clus, '\n')
							y <- rep(0, length(df$Cluster)) 
							y[df$Cluster == clus] <- 1
							df$yy <- y
							prmatrix(summary(glmmPQL(form, data=df, random = as.formula(paste0('~ 1|', subject)), family=binomial, verbose=F))$tTable)
							
						}
					}
					sink()
				}
		)
		
	}
}


##########
B2M <- function(x) {
	x[x==0] <- min(x[x!=0])
	x[x==1] <- max(x[x!=1])
	log(x / (1-x))
}


fastDist <- function(X) {
	temp <- colSums(X^2)  
	D <- outer(temp, temp, "+") - 2 * t(X) %*% X
	diag(D) <- 0
	sqrt(D)
}

fastLM <- function(Y, M) {
	Y <- as.matrix(Y)
	XXI <- solve(t(M) %*% M)
	dof <- ncol(Y) - ncol(M)
	est <- XXI %*% t(M) %*% t(Y)
	resid <- t(Y) - M %*% est
	sigma <- sqrt(colSums(resid^2)/dof)
	Pvec <- 2*pt(-abs(t(est/(sqrt(diag(XXI))))/sigma), dof)
	return(Pvec)
}

matrix.paste0 <- function (m.list) {
	p <- length(m.list)
	
	m <- max(unlist(sapply(m.list, function (x) nrow(x))))
	n <- max(unlist(sapply(m.list, function (x) ncol(x))))
	res <- NULL
	for (i in 1:p) {
		mat <- m.list[[i]]
		if (!is.matrix(mat)) {
			mat <- matrix(mat, m, n)
		} else {
			if (!is.null(rownames(mat))) {
				row.names <- rownames(mat)
			}
			if (!is.null(colnames(mat))) {
				col.names <- colnames(mat)
			}
			
		}
		res <- paste0(res, mat)
	}
	return(matrix(res, m, n, dimnames=list(row.names, col.names)))
}

normalize.func <- function (x) {
	names <- paste(x)
	tab <- table(x)
	tab.name <- names(tab)
	tab.x <- qqnorm(as.numeric(tab.name), plot.it = F)$x
	names(tab.x) <- tab.name
	as.numeric(tab.x[names])
}


