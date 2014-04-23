
pathToTheOmsProject="~/Desktop/IS/oms3.prj.NewAge_EGU.SNOW/"

DegreeDay=read.table(paste(pathToTheOmsProject,"output/Jho_Wrigth_SWE_LAST_DAILY_DD.csv",sep=""),sep=",",skip=7)
CazorziDallaFontana=read.table(paste(pathToTheOmsProject,"output/Jho_Wrigth_SWE_LAST_DAILY_CD.csv",sep=""),sep=",",skip=7)
Measured=read.table(paste(pathToTheOmsProject,"data/JoeWrigth/Daily/MEASURED_2000_2014.csv",sep=""))


plot(Measured[,1],pch=16,col="gray",cex=0.8,ylim=c(0,1600),xlab="Time [Days]",ylab="SWE [mm]")
lines(DegreeDay[,3],col="blue")
lines(CazorziDallaFontana[,3],col="red")
legend(1,1600, c("Measured","Cazorzi-Dalla Fontana", "Degree-Day"), col=c("gray","red","blue"),lty = c(-1, 1, 1), pch = c(16, NA, NA),merge = TRUE)
