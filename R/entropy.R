#functions for calculating K-L entropy
entropy.fn<-function(z1,z2)
{
#	s1<-sum(z1)
#	s2<-sum(z2)
#	log(s1-1)+log(s1-z1[4]+2)+log(s1-z1[4]+1)-log(s2-1)-log(s2-z2[4]+2)-log(s2-z2[4]+1)+lgamma(s1+1)-sum(lgamma(z1+1))-lgamma(s2+1)+sum(lgamma(z2+1))+sum((z1-z2)*digamma(z1+1))-(s1-z1[4]-s2+z2[4])*(digamma(s1-z1[4]+3)-digamma(s1-z1[4]+1))-(s1-s2)*digamma(s1+2)
	s1 <- sum(z1)
	s2 <- sum(z2)
	lgamma(s2 + 1) - sum(lgamma(z2 + 1)) + log(s2 + 1) + log(s2 - z2[4] + 2) + log(s2 - z2[4] + 1) - lgamma(s1 + 1) + sum(lgamma(z1 + 1)) - log(s1 + 1) - log(s1 - z1[4] + 2) - log(s1 - z1[4] + 1) + sum((z2 - z1) * digamma(z2 + 1)) - (s2 - z2[4] - s1 + z1[4]) * (digamma(s2 - z2[4] + 3) - digamma(s2 - z2[4] + 1)) - (s2 - s1) * digamma(s2 + 2)
}

#function for calculating Shannon entropy
shannon.fn <- function(z1)
{
	z1 <- z1 / sum(z1)
	z1[z1 == 0] <- 1
	z1 <- z1 * log(z1)
	z1 <- -sum(z1)
	z1
}

#function for calculating K-L entropy relative to the prior
klprior.fn <- function(z)
{
	s <- sum(z)
	lgamma(s + 1) - sum(lgamma(z + 1)) + log(s + 3) + log(s + 2) + log(s + 1) - log(6) + sum(z * digamma(z + 1)) - s * digamma(s + 4)
}

#klprior.fn <- function(z)
#{
#	s <- sum(z)
#	lgamma(s + 1) - sum(lgamma(z + 1)) + log(s + 1) + log(s - z[4] + 2) + log(s - z[4] + 1) + sum((z * digamma(z + 1))) - (s - z[4]) * (digamma(s - z[4] + 3) - digamma(s - z[4] + 1)) - s * digamma(s + 2)
#}

#klprior.fn <- function(z1)
#{
#	s <- sum(z1)
#	pd <- log(2) - log((s + 1) * (s - z1[4] + 2) * (s - z1[4] + 1))
#	a <- lfactorial(s) - sum(lfactorial(z1))
#	pd + (5 / 2) * (s + z1[4]) - a
#}




