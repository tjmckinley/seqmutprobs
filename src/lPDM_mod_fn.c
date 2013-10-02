#include "functions.h"
#include <R.h>
#include <Rmath.h>

//function to calculate l[P'(D|M)] for a given distribution of bases
double lPDM_mod_fn(int *z, int ind, double pstar)
{	
	/*'z' is a pointer to a vector of length 5, where the
		first 4 elements correspond to each base (with the consensus
		as the fourth element z[3]). The final element is 
		S-z[3]=sum(z[0:2])
	'ind' denotes which model (from 0:9) is to be calculated
	'pstar' is the overall mutation rate*/
	
	double lPDM = 0.0;
/*	switch(ind)*/
/*	{*/
/*		//Null p1=p2=p3=p3*/
/*		case 0 :lPDM=z[4]*log(pstar/3.0)+z[3]*log(1.0-pstar);*/
/*			break;*/
/*		//Alt that one free pi is different: e.g. p1!=p2=p3 etc. but mutation rate constrained to p**/
/*		case 1 :lPDM=-(z[1]+z[2])*log(2.0)+lfactorial(z[0])+lfactorial(z[1]+z[2])-lfactorial(z[4]+1)+z[4]*log(pstar)+z[3]*log(1.0-pstar);*/
/*			break;*/
/*		case 2 :lPDM=-(z[0]+z[2])*log(2.0)+lfactorial(z[1])+lfactorial(z[0]+z[2])-lfactorial(z[4]+1)+z[4]*log(pstar)+z[3]*log(1.0-pstar);*/
/*			break;*/
/*		case 3 :lPDM=-(z[0]+z[1])*log(2.0)+lfactorial(z[2])+lfactorial(z[0]+z[1])-lfactorial(z[4]+1)+z[4]*log(pstar)+z[3]*log(1.0-pstar);*/
/*			break;*/
/*		//Alt that all pis are different but mutation rate constrained to p**/
/*		case 4 :lPDM=log(2.0)+lfactorial(z[0])+lfactorial(z[1])+lfactorial(z[2])-lfactorial(z[4]+2)+z[4]*log(pstar)+z[3]*log(1.0-pstar);*/
/*			break;*/
/*		//Alt that pis are uniform but not constrained to sum to p**/
/*		case 5 :lPDM=-z[4]*log(3.0)+lfactorial(z[4])+lfactorial(z[3])-lfactorial(z[3]+z[4]+1);*/
/*			break;*/
/*		//Alt that one free pi is different: e.g. p1!=p2=p3 etc.*/
/*		case 6 :lPDM=-(z[1]+z[2])*log(2.0)+lfactorial(z[1]+z[2])+lfactorial(z[0])+lfactorial(z[3])-log(z[4]+1)-lfactorial(z[4]+z[3]+1);*/
/*			break;*/
/*		case 7 :lPDM=-(z[0]+z[2])*log(2.0)+lfactorial(z[0]+z[2])+lfactorial(z[1])+lfactorial(z[3])-log(z[4]+1)-lfactorial(z[4]+z[3]+1);*/
/*			break;*/
/*		case 8 :lPDM=-(z[0]+z[1])*log(2.0)+lfactorial(z[0]+z[1])+lfactorial(z[2])+lfactorial(z[3])-log(z[4]+1)-lfactorial(z[4]+z[3]+1);*/
/*			break;*/
/*		//Alt that all pis are different*/
/*		case 9 :lPDM=log(2.0)+lfactorial(z[0])+lfactorial(z[1])+lfactorial(z[2])+lfactorial(z[3])-log(z[4]+2)-log(z[4]+1)-lfactorial(z[4]+z[3]+1);*/
/*			break;*/
/*	}*/
	int i;
	double z1[5];
	for(i=0;i<5;i++) z1[i]=(double) z[i];
	switch(ind)
	{
		//Null p1=p2=p3=p/3 where p<=p*
		case 0 :lPDM=-z1[4]*log(3.0)-log(pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,1,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		//Alt that one free pi is different: e.g. p1!=p2=p3 etc. but mutation rate constrained to be <= p*
		case 1 :lPDM=-(z[1]+z[2])*log(2.0)+lfactorial(z[0])+lfactorial(z[1]+z[2])-lfactorial(z[4]+1)-log(pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,1,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		case 2 :lPDM=-(z[0]+z[2])*log(2.0)+lfactorial(z[1])+lfactorial(z[0]+z[2])-lfactorial(z[4]+1)-log(pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,1,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		case 3 :lPDM=-(z[0]+z[1])*log(2.0)+lfactorial(z[2])+lfactorial(z[0]+z[1])-lfactorial(z[4]+1)-log(pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,1,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		//Alt that all pis are different but mutation rate constrained to be <= p*
		case 4 :lPDM=log(2.0)+lfactorial(z[0])+lfactorial(z[1])+lfactorial(z[2])-lfactorial(z[4]+2)-log(pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,1,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		//Alt p1=p2=p3=p/3 where p>p*
		case 5 :lPDM=-z1[4]*log(3.0)-log(1.0-pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,0,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		//Alt that one free pi is different: e.g. p1!=p2=p3 etc. but mutation rate constrained to be > p*
		case 6 :lPDM=-(z[1]+z[2])*log(2.0)+lfactorial(z[0])+lfactorial(z[1]+z[2])-lfactorial(z[4]+1)-log(1.0-pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,0,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		case 7 :lPDM=-(z[0]+z[2])*log(2.0)+lfactorial(z[1])+lfactorial(z[0]+z[2])-lfactorial(z[4]+1)-log(1.0-pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,0,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		case 8 :lPDM=-(z[0]+z[1])*log(2.0)+lfactorial(z[2])+lfactorial(z[0]+z[1])-lfactorial(z[4]+1)-log(1.0-pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,0,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
		//Alt that all pis are different but mutation rate constrained to be > p*
		case 9 :lPDM=log(2.0)+lfactorial(z[0])+lfactorial(z[1])+lfactorial(z[2])-lfactorial(z[4]+2)-log(1.0-pstar)+pbeta(pstar,z1[4]+1,z1[3]+1,0,1)+lbeta(z1[4]+1,z1[3]+1);
			break;
	}
	return lPDM;
}

