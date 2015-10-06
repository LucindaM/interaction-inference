cal.logmeanlik <- function(x){
	b <- x*log10(exp(1))
	base <- 10^(b-as.integer(b))
	exponent <- as.integer(b)
	max.exp <- max(exponent)
	diff.exp <- exponent - max.exp
	sum.b <- sum(base*10^(diff.exp))
	mean.b <- sum.b/length(b)
	mean.loglik <- max.exp*log(10) + log(mean.b)
	mean.loglik
}

fpath <- "~/Documents/DenSerIntInf/filterstudy/par1/result_30000/"
fname <- "fs.delta.chi.rda"
fpath.case <- "~/Documents/DenSerIntInf/filterstudy/par1/data/"
fname.case <-"casedata.csv"

delta.real <- 365
imm.real <- 1/delta.real
chi.real <- 1

res.file <- paste(fpath,fname,sep='')
case.file <- paste(fpath.case,fname.case,sep='')

casedata <- read.csv(file=case.file, row=1)
casedata[1,] <- c(0, NA, NA, NA, NA, NA, NA, NA, NA)

load(res.file)
flt.stdy.delta.chi <- result

cols <- c('darkcyan','darkgoldenrod1','darkmagenta','chartreuse4')

delta <- sapply(flt.stdy.delta.chi, function(x)x$params['delta'])
imm <- 1/delta
chi <- sapply(flt.stdy.delta.chi, function(x)x$params['chi'])
loglik <- sapply(flt.stdy.delta.chi, function(x)x$pf$loglik)

delta.set <- sort(unique(delta), decreasing=F)
imm.set <- 1/delta.set
chi.set <- sort(unique(chi), decreasing=F)

chi.select <- which(chi.set >= chi.real-0.1 & chi.set <= chi.real+0.1)

result <- array(NA, c(5, length(delta.set), length(chi.set)))
result.mean <- array(NA, c(length(delta.set), length(chi.set)))

for(jd in 1:length(delta.set)){
	for(jc in 1:length(chi.set)){
		subset <- which(delta == delta.set[jd] & chi == chi.set[jc])
		result[,jd,jc] <- loglik[subset]
		result.mean[jd,jc] <- cal.logmeanlik(loglik[subset])
	}
}

dum <- matrix(NA,3,1)
dum[1,1] <- 1
dum[2,1] <- 2
dum[3,1] <- 3

graph.name <- "profile_par1.pdf"

pdf(width=6, height=10, file=graph.name)

layout(dum,width=5.5,height=c(2.5, 3.25,3.25),TRUE)
par(mar=c(1,4,4,2))
case.range <- range(c(casedata$y1, casedata$y2, casedata$y3, casedata$y4), na.rm=T)/1000
plot(range(casedata$time), range(c(casedata$y1, casedata$y2, casedata$y3, casedata$y4), na.rm=T)/1000, type='n', xlab='',ylab='', cex.axis=1.25, axes=F)
lines(casedata$time, casedata$y1/1000, lty=1, lwd=1.5, col=cols[1])
lines(casedata$time, casedata$y2/1000, lty=1, lwd=1.5, col=cols[2])
lines(casedata$time, casedata$y3/1000, lty=1, lwd=1.5, col=cols[3])
lines(casedata$time, casedata$y4/1000, lty=1, lwd=1.5, col=cols[4])

title(ylab=list('sim. cases (000s)',cex=1.25))
mtext(text='time (years)', side=3, cex=0.95, padj=-2.5)
axis(3, seq(0,40,by=10), cex.axis=1.25)
axis(2, seq(0, round(max(casedata$y1, casedata$y2, casedata$y3, casedata$y4, na.rm=T)/1000), length=3), cex.axis=1.25)
box()

delta.expand <- seq(0.1,365,by=0.001)
imm.expand <- rev(1/delta.expand)
prof1.loess <- loess(rev(apply(result.mean, 1, max)) ~ rev(imm.set), span=0.7)
prof1.pred <- predict(prof1.loess, rev(imm.expand))
ci.95 <- qchisq(df=1, p=0.95)/2
baseline.imm <- max(prof1.pred)-ci.95
l1 <- imm.expand[min(which(prof1.pred > baseline.imm))]
l2 <- imm.expand[max(which(prof1.pred > baseline.imm))]
par(mar=c(4,4,1,2))

plot(rev(imm.set), rev(apply(result.mean, 1, max) -baseline.imm), pch=19, col='grey', cex=2, xlab='',ylab='', cex.axis=1.25, ylim=c(-500,5), xlim=c(0,2))
lines(rev(imm.expand), prof1.pred-baseline.imm, lwd=2)

abline(h=0, lty=2)
abline(v=imm.real, lty=2, col='red')
title(xlab=list(expression(paste(1/delta, ', duration of temporary imm. (years)', sep='')), cex=1.25), ylab=list(expression(paste(Delta, ' loglik', sep='')),cex=1.25))

chi.expand <- seq(0,2,by=0.001)
prof2.loess <- loess(apply(result.mean, 2, max) ~ chi.set, span=0.5)
prof2.pred <- predict(prof2.loess, chi.expand)

ci.95 <- qchisq(df=1, p=0.95)/2
baseline.chi <- max(prof2.pred)-ci.95
l1 <- chi.expand[min(which(prof2.pred > baseline.chi))]
l2 <- chi.expand[max(which(prof2.pred > baseline.chi))]


par(mar=c(4,4,1,2))

plot(chi.set[], apply(result.mean, 2, max)[]-baseline.chi, pch=19, col='grey', cex=2, xlab='',ylab='', cex.axis=1.25, ylim=c(-500,5), xlim=c(0.5,1.5))
lines(chi.expand, prof2.pred-baseline.chi, lwd=2)
abline(h=0, lty=2)
abline(v=chi.real, lty=2, col='red')

title(xlab=list(expression(paste(chi,', enhancement', sep='')), cex=1.25), ylab=list(expression(paste(Delta, ' loglik', sep='')),cex=1.25))


op <- par(fig=c(0.63,0.96,.08,.32), new = TRUE)
plot(chi.set, apply(result.mean, 2, max)-baseline.chi, pch=19, col='grey', cex=1.5, xlab='',ylab='', axes=F)
ll.range <- round(range(apply(result.mean, 2, max)-baseline.chi))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "azure")
points(chi.set, apply(result.mean, 2, max)-baseline.chi, pch=19, col='grey', cex=0.95)
lines(chi.expand, prof2.pred-baseline.chi, lwd=1)
axis(1, seq(0, 2, by=0.2),cex.axis=0.85, tcl=-0.25, padj= -1.5)
axis(4, seq(ll.range[1], ll.range[2], length=5),cex.axis=0.85, tcl=-0.25, hadj=0.5, padj=-1)

abline(h=0, lty=2)
abline(v=chi.real, lty=2, col='red')
box(lwd=1)
par(op)

op <- par(fig=c(0.63,0.96,.4,.64), new = TRUE)
plot(rev(imm.set), rev(apply(result.mean, 1, max)-baseline.imm), pch=19, col='grey', cex=1.5, xlab='',ylab='', axes=F)
ll.range <- round(range(apply(result.mean, 1, max)-baseline.imm))
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "azure")
points(rev(imm.set), rev(apply(result.mean, 1, max)-baseline.imm), pch=19, col='grey', cex=0.95)
lines(rev(imm.expand), prof1.pred-baseline.imm, lwd=1)

axis(1, seq(0, 10, by=2),cex.axis=0.85, tcl=-0.25, padj= -1.5)
axis(4, seq(ll.range[1], ll.range[2], length=5),cex.axis=0.85, tcl=-0.25, hadj=0.5, padj=-1)

abline(h=0, lty=2)
abline(v=imm.real, lty=2, col='red')
box(lwd=1)

dev.off()


