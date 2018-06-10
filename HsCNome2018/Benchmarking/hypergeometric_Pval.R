cat("enter the entire filepath incl name (col labels MUST BE (CASE-SENSITIVE) a, b, k and x):\n")
infile<-readLines(con = file("stdin"),n=1,warn=TRUE)
cat(infile,"\n")

data<-read.table(infile,header=T)
print(data)

cat("\nHypergeometric test will be formulated as follows:\n\nphyper(x-1,a,b,k,x,lower.tail=F)\nWhere:\nx=overlap\na:members of the type of category you wish to detect in the population\nb:entire population or the background size\nk:sample size\n\n")

x<-data$x
a<-data$a
b<-data$b
k<-data$k

print(phyper(x-1,a,b,k,lower.tail=F))
