library(stats)
args = commandArgs(trailingOnly=TRUE)
print(args[[1]])

df = read.csv(args[[1]],sep=",",header=TRUE)
total = 50280

res <- phyper(q = df$k, m = df$n, n = total - df$n, k=df$N,lower.tail=F)
df$res <- res
write.csv(df,"tmp/phyper.csv")
