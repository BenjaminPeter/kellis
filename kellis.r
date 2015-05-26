data = t(read.table("s.txt"))[,4:6]
to.switch <- data[,2] > data[,3]
data[to.switch,2:3] <-data[to.switch,3:2] 
ps <- rowSums(data)


llmultinomial <- function(x,p){                                      
    p <- p/sum(p)                                                    
    ll <- lfactorial(sum(x)) - sum(lfactorial(x))  + sum(x * log(p)) 
    return(ll)                                                       
}                                                                    


ll.full <- sapply(1:76, function(x)llmultinomial(data[x,], data[x,]/ps[x]))
p.merged <- rowSums(data[,1:2])/2
p.alt <- cbind(p.merged, p.merged, data[,3])

ll.2 <- sapply(1:76, function(x)llmultinomial(data[x,],
                                              p.alt[x,]/ps[x]))

deccel <- data[,1] > data[,2]
p.deccel <- data
p.deccel[deccel,1:2] <- p.alt[deccel,1:2]

ll.3 <- sapply(1:76, function(x)llmultinomial(data[x,],
                                              p.deccel[x,]/ps[x]))

ll.full <- sapply(1:76, function(x)llmultinomial(data[x,], data[x,]/ps[x]))
p.merged <- rowSums(data[,2:3])/2
p.alt <- cbind(data[,1], p.merged, p.merged)

ll.4 <- sapply(1:76, function(x)llmultinomial(data[x,],
                                              p.alt[x,]/ps[x]))



aic.full <- 6-2*ll.full
aic.2 <- 4-2*ll.2      

bic.full <- log(rowSums(data)) * 3 - 2*ll.full
bic.2 <- log(rowSums(data)) * 2 - 2*ll.2
bic.3 <- log(rowSums(data)) * 3 - 2*ll.2
bic.4 <- log(rowSums(data)) * 3 - 2*ll.4

pdf("kellis1.pdf",width=8,height=4)
par(mfcol=c(1,2))             
plot(data[,1:2], xlab="S_k", ylab="S_y")              
abline(0,1)                   
plot(data[,1], data[,3], xlab="S_k", ylab="S_z")              
abline(0,1)                   
dev.off()

pdf("kellis2.pdf",width=8,height=4)
par(mfcol=c(1,2))             
plot(bic.3, bic.2, xlab="BIC H_0", ylab="BIC H_1")              
abline(0,1)                   
hist(2*(ll.3 - ll.2), xlab="2*(L(H_0)-L(H_1))")              
#abline(0,1)                   
dev.off()
