#R script

draw <- function(filename,main){
    a1 <- read.table(filename)
    plot(a1[,1],a1[,2],type='h',lwd=6,lend="square",
    xlim=c(0,40),ylim=c(0,3000000),
    xlab="group size",ylab="group number",main="")
}

png("./figure7.png",res=200,width=1800,height=1800)
par(mfrow=c(2,2),mar=c(4,4,1,1))
draw("WTreplicate1Forward.cov","ec1, replicate1")
draw("WTreplicate1.cov","ec2,replicate1")
draw("WTreplicate2Forward.cov","ec1, replicate2")
draw("WTreplicate2.cov","ec2,replicate2")



