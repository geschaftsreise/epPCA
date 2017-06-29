#' Easy Peasy PCA Plots
#'
#' This function generates PCA plots from a matrix.
#' @param mtx A matrix
#' @param lab A character label
#' @param ncp Number of components to examine. 0 (default) allows autodetection.
#' @param gradient A vector of colors as a 3 color gradient for displaying the component axis values from negative to positive
#' @param rowColors A vector of length 'r' of colors indicating group identity for each of 'r' rows of the matrix.
#' @param biplot A vector of 2 numbers indicating which components to use in the biplot (e.g. c(1,2)), default (0,0) will not generate a biplot figure.
#' @keywords PCA
#' @export
#' @examples
#' mtx <- t(as.matrix(iris[,1:4]))
#' rowColors <- rep("white",length(iris$Species));
#' rowColors[which(iris$Species=="setosa")] <- "red"
#' rowColors[which(iris$Species=="virginica")] <- "blue"
#' rowColors[which(iris$Species=="versicolor")] <- "light grey"
#' library(RColorBrewer); gradient <- rev(colorRampPalette(rev(brewer.pal(11,'BrBG')))(100))
#' epPCA(mtx,"irisData",ncp=3,gradient=gradient,rowColors=rowColors,biplot=c(1,2))

 
epPCA <- function(mtx,lab="epPCAplot",ncp=0,gradient="default",rowColors="dark grey",biplot=c(0,0)){
	require(RColorBrewer)
	require(grid)
	require(gridExtra)
	
	# dependency-specific code here	
	p <- prcomp(mtx)
	
	numColors <- 32
	cv.spect <- colorRampPalette(rev(brewer.pal(11,'BrBG')))(numColors)

	# prepare labels
	dateStamp <- format(Sys.time(),"%Y-%b-%d")
	
	# check inputs
	if (!is.matrix(mtx)) { 
		print("mtx is not a matrix object"); return();
	}
	if (ncp<1 | ncp>length(p$sdev) | ncp != round(ncp,0)) {
		ncp <- length(p$sdev)
	}
	if (!is.vector(gradient)) {
		gradient <- cv.spect
	}
	if (gradient[1]=="default") { gradient <- cv.spect }
	if (!is.vector(biplot)) {
		print("biplot variable must be a vector of length 2"); return()
	}
	if (length(biplot)!=2) {
		print("biplot variable must be a vector of length 2"); return()
	}
	if (any(biplot<1) | any(biplot>ncp)) {
		if (biplot[1]!=0 | biplot[2]!=0) {
			print("biplot variable must contain integer between 1 and ncp, referring to 2 valid component numbers"); return()
		}
	}
	
	######### generate main figure
	{
	# formatting variables
	mainFigFile <- paste(lab,"_mainPCAplot_",dateStamp,".pdf",sep="")
	evWidth   <- 1/nrow(p$rotation)
	evHeight  <- 1/ncp
	barGap    <- evHeight * 0.1
	barHeight <- evHeight * 0.8
	barWidth  <- p$sdev/max(p$sdev)
	gradValue <- t(ceiling(length(gradient)*(p$rotation[,1:ncp]/max(abs(p$rotation[,1:ncp]))+1)/2))
	gradValue[which(gradValue==0)] <- 1
	gradValue[which(gradValue>length(gradient))] <- length(gradient)

	pdf(mainFigFile,width=8,height=4,useDingbats=FALSE)
	vpMain <- viewport(x = 0.02, y = 0.02, w = .96, h = 0.96, just = c("left", "bottom"), name = "vpMain")
	pushViewport(vpMain)
	grid.text("principle components", x=0, y = 0.45,just="centre",gp=gpar(fontsize=10),rot=90)

	# vp0 contains the left label infoformation
	vp0 <- viewport(x = 0.05, y = 0.05, w = 0.045, h = .9, just = c("left", "bottom"), name = "vp0")
	pushViewport(vp0)
	
	vp0top <- viewport(x = 0, y = 0.875, w = 1, h = .1, just = c("left", "bottom"), name = "vp0top")
	pushViewport(vp0top)
	grid.text("group", x=0.5, y = 0.5,just="centre",gp=gpar(fontsize=10))
	upViewport()
	
	vp0bottom <- viewport(x = 0, y = 0.075, w = 1, h = .775, just = c("left", "bottom"), name = "vp0bottom")
	pushViewport(vp0bottom)
	for (i in 1:ncp) {
		grid.text(paste("PC",i,sep=""), x=0.5, y=1-i*evHeight+evHeight/2,just="centre",gp=gpar(fontsize=10))
	}
	upViewport()
	
	vp0foot <- viewport(x = 0, y = 0, w = 1, h = .05, just = c("left", "bottom"), name = "vp0foot")
	pushViewport(vp0foot)
	upViewport()
	upViewport()
	
	# vp1 contains the group ID plots as well as the raw component eigenvalues
	vp1 <- viewport(x = 0.105, y = 0.05, w = 0.685, h = .9, just = c("left", "bottom"), name = "vp1")
	pushViewport(vp1)

	# vp1top contains the group info
	vp1top <- viewport(x = 0, y = 0.875, w = 1, h = .1, just = c("left", "bottom"), name = "vp1top")
	pushViewport(vp1top)
	for (j in 1:nrow(p$rotation)) {
		grid.rect(gp=gpar(col=NA,fill=rowColors[j]),x=(j-1)*evWidth,width=evWidth,just="left")
	}
	upViewport()

	# vp1bottom contains the eigenvalue data
	vp1bottom <- viewport(x = 0, y = 0.075, w = 1, h = .775, just = c("left", "bottom"), name = "vp1bottom")
	pushViewport(vp1bottom)
	for (i in 1:ncp) {
		for (j in 1:nrow(p$rotation)) {
			grid.rect(gp=gpar(col=NA,fill=gradient[gradValue[i,j]]),x=(j-1)*evWidth,y=1-i*evHeight,height=evHeight,width=evWidth,just=c("left","bottom"))
		}
	}
	upViewport()
	
	# vp1foot contains the gradient scale
	vp1foot <- viewport(x = 0, y = 0, w = 1, h = .05, just = c("left", "bottom"), name = "vp1foot")
	pushViewport(vp1foot)
	binWidth <- .2/length(gradient)
	grid.text("-", x=0.08, y = 0.7,just="right",gp=gpar(fontsize=8))
	grid.text("+", x=0.32, y = 0.7,just="left",gp=gpar(fontsize=8))
	grid.text("eigenvalue", x=0.2, y = 0,just="centre",gp=gpar(fontsize=8))
	for (i in 1:length(gradient)) {
		grid.rect(gp=gpar(col=NA,fill=gradient[i]),x=.1+(i-1)*binWidth,y=0.4,height=.55,width=binWidth,just=c("left","bottom"))
	}
	upViewport()
	upViewport()

	vp2 <- viewport(x = 0.805, y = 0.05, w = 0.15, h = .9, just = c("left", "bottom"), name = "vp2")
	pushViewport(vp2)

	# vp2top contains the sdev axis
	vp2top <- viewport(x = 0, y = 0.875, w = 1, h = .1, just = c("left", "bottom"), name = "vp2top")
	pushViewport(vp2top)
	grid.lines(c(0,0),c(0,0.2))
	grid.lines(c(1,1),c(0,0.2))
	grid.lines(c(0,1),c(0.1,0.1))
	grid.text("0", x=0, y = 0.5,just="centre",gp=gpar(fontsize=8))
	grid.text(round(max(p$sdev),0), x=1, y = 0.5,just="centre",gp=gpar(fontsize=8))
	upViewport()

	# vp2bottom contains the sdev bars
	vp2bottom <- viewport(x = 0, y = 0.075, w = 1, h = .775, just = c("left", "bottom"), name = "vp2bottom")
	pushViewport(vp2bottom)
	for (i in 1:ncp) {
		grid.rect(gp=gpar(col=NA,fill="light grey"),x=0,y=1-((i-1)*evHeight+barGap)-barHeight/2,width=barWidth[i],height=barHeight,just="left")
	}
	upViewport()
	
	# vp2foot contains the sdev label
	vp2foot <- viewport(x = 0, y = 0, w = 1, h = .05, just = c("left", "bottom"), name = "vp2foot")
	pushViewport(vp2foot)
	grid.text("sdev", x=0.5, y = 0.5,just="centre",gp=gpar(fontsize=10))
	upViewport()
	
	popViewport(2)
	garbage <- dev.off()
	}
	
	######### generate biplot if requested
	if (biplot[1]!=0 & biplot[2]!=0) {
		biplotFigFile <- paste(lab,"_PC",biplot[1],"_PC",biplot[2],"_biplot_",dateStamp,".pdf",sep="")
		pcA <- p$rotation[,biplot[1]]
		pcB <- p$rotation[,biplot[2]]
		aMin <- round(min(pcA),3)
		aMax <- round(max(pcA),3)
		bMin <- round(min(pcB),3)
		bMax <- round(max(pcB),3)
		pdf(biplotFigFile,width=6,height=6,useDingbats=F)
		plot(pcA,pcB,pch=19,cex=1,
			xlim=c(aMin,aMax),
			ylim=c(bMin,bMax),
			xlab=paste("PC",biplot[1],sep=""),
			ylab=paste("PC",biplot[2],sep=""),
			col=rowColors,axes=F)
		axis(1,at=c(aMin,aMax))
		axis(2,at=c(bMin,bMax))
		garbage <- dev.off()
	}
}