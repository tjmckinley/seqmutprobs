#functions for calculating normalised K-L entropy
entropy.fn<-function(z1,z2,...)
{
	s1<-sum(z1)
	s2<-sum(z2)
	log(s1-1)+log(s1-z1[4]+2)+log(s1-z1[4]+1)-log(s2-1)-log(s2-z2[4]+2)-log(s2-z2[4]+1)+lgamma(s1+1)-sum(lgamma(z1+1))-lgamma(s2+1)+sum(lgamma(z2+1))+sum((z1-z2)*digamma(z1+1))-(s1-z1[4]-s2+z2[4])*(digamma(s1-z1[4]+3)-digamma(s1-z1[4]+1))-(s1-s2)*digamma(s1+2)
}

#function for calculating Shannon entropy
shannon.fn <- function(z1, ...)
{
	z1 <- z1 / sum(z1)
	z1[z1 == 0] <- 1
	z1 <- z1 * log(z1)
	z1 <- -sum(z1)
	z1
}

