#Modified from Balbuena et al. 2013 - retrieved from https://www.uv.es/cophylpaco/index.html

#PACo simulations Type I error as per Balbuena et al. **** - comparison with HCT [1] and Parafit [2]#
#[1] Hommola et al. (2009) Mol. Biol. Evol. 26: 1457-1468.
#[2] Legendre et al. (2002) Syst. Biol. 51: 217-234.
#load packages: ape, vegan#
library(ape)
library(vegan)
#
#function adjusting the sizes of host and parasite patristic distances matrices & performing Procrustes fit#
PACo <- function (H.dist, P.dist, HP.bin)
{ 
HP.bin <- which(HP.bin > 0, arr.in=TRUE)
H.PCo <- pcoa(H.dist, correction="cailliez")$vectors #Performs PCo of Host distances
P.PCo <- pcoa(P.dist, correction="cailliez")$vectors #Performs PCo of Parasite distances
H.PCo <-H.PCo[HP.bin[,1],] #adjust Host PCo vectors 
P.PCo <-P.PCo[HP.bin[,2],]  ##adjust Parasite PCo vectors
m2 <- procrustes(H.PCo, P.PCo)$ss
}
#The following function implements HCT. Code written by Kerstin Hommola
#obtained at http://www1.maths.leeds.ac.uk/~kerstin/#  
# 
#############################
### test for cospeciation ###
#############################

sim_cosp=function(d_h, d_p, rel, nsample){

#d_h: host distance matrix
#d_p: parasite distance matrix
#rel: matrix with host-parasite interaction

#check for correct arguments
if(is.data.frame(d_h))
	d_h=as.matrix(d_h)
if(is.data.frame(d_p))
	d_p=as.matrix(d_p)
if(is.matrix(d_h)==FALSE)
	stop(paste(sQuote("d_h"),"is not a matrix"))
if(is.matrix(d_p)==FALSE)
	stop(paste(sQuote("d_p"),"is not a matrix"))
if(dim(rel)[1]!=dim(d_p)[1]|| dim(rel)[2]!=dim(d_h) [2])
	stop(paste(sQuote("rel"),"has not correct dimensions"))
if(sum(rel!=0)+sum(rel!=1)!=dim(rel)[1]*dim(rel)[2])
	stop(paste("Entries in", sQuote("rel"), "not correct. Please check input."))
if(is.numeric(nsample)==FALSE)
	stop(paste(sQuote("nsample"),"is not numeric"))

m=sum(rel) #number of interaction links

#rearrange relation matrix 
#first column contains host (labeled by number)
#second column contains parasite infecting host in column 1
rel_1=matrix(nrow=m,ncol=2)
s=1;
while(s<=m){
	for(i in 1:dim(rel)[1]){
		for(j in 1:dim(rel)[2]){
			if(rel[i,j]==1){
				rel_1[s,1]=j
				rel_1[s,2]=i
				s=s+1
			}
		}
	}
}

#get distance of hosts and parasites for each pair of edges
#get distance of hosts in rel_1[i,1] and rel_1[j,1] + parasites respectively

getdist=function(v,matrix,row){
	#v...vector containing label of relation matrix (either host or parasite)
	#matrix... distance matrix
	#row... label of rows of distance matrix (either host or parasite)
	#equals column labels cause distance matrix is symmetric

	vec=vector(length=choose(m,2))     #vector contain distances
	o.r=order(row)
			
	a=1
	for(i in 1:(m-1)){
		k=o.r[v[i]]
		for(j in (i+1):m){
			t=o.r[v[j]]
			vec[a]=matrix[k,t]
			a=a+1
		}
	}
	return(vec)
}

x=getdist(rel_1[,1],d_h,as.integer(rownames(d_h)))	#distance between two host
y=getdist(rel_1[,2],d_p,as.integer(rownames(d_p)))	#distance between two parasites

#Pearson correlation coefficient for real data
r_real=cor(x,y)   			

############	
#simulation#
############

m.h=dim(d_h)[1]   #number of hosts
m.p=dim(d_p)[1]	#number of parasites
r_s=vector(length=nsample)	#vector with simulated correlation coefficients r*

for(i in 1:nsample){
	simlabel.h=sample(1:m.h)      #random permutation of numbers 1 to m.h
	simlabel.p=sample(1:m.p)	#random permutation of numbers 1 to m.p
	
	x_s=getdist(rel_1[,1],d_h,simlabel.h)	#distance between two host
	y_s=getdist(rel_1[,2],d_p,simlabel.p)	#distance between two parasites

	#correlation coefficient for simulated data
	r_s[i]=cor(x_s,y_s)
}
p=sum(r_s>=r_real)/nsample
#plot(r_s,xlab="number of permutations", ylab="r*")  #plot r* obtained after each randomization
return(p)
}
#HHHHHHHHHHHHHHHHHHHHHHHH** PACo SIMUL STARTS HERE **HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
ptm <- proc.time()
N.runs =  9 #######             Change number of simulations as required           ######
N.perm = 9 #######               Change number of permutation for tests as required ######
NH = 10 #######   Change number of hosts (NH), parasites (NP) and H-P associtations (N.assoc) as required  ######
NP = 10
N.assoc = 10
#
seed1 <- .Random.seed[trunc(runif(1,1,626))]
seed2 <- .Random.seed[trunc(runif(1,1,626))]
#
for (i in c(1:N.runs))
{
set.seed(i+seed1)
TreeH <- rtree (n= NH, rooted=TRUE, tip.label=NULL, br = runif)			
host.D <- cophenetic (TreeH) #Host patristic distance matrix
TreeP <- rtree (n= NP, rooted=TRUE, tip.label=NULL, br = runif)
para.D <- cophenetic (TreeP)  #Parasite patristic distance matrix
if (N.assoc <= NH | N.assoc <= NP) #control statement to avoid all parasites beig associated to 
{	flag <- TRUE					#a single host and vice versa
	while (flag == TRUE)	{ 
	HP <- matrix ((sample (c(rep(1, each=N.assoc), rep(0, each=((NH*NP)-N.assoc))))), NH)
		if(any(rowSums(HP) == N.assoc)) flag <- TRUE 
		else if(any(colSums(HP) == N.assoc)) flag <- TRUE 
		else flag <- FALSE
							}					
} 
else HP <- matrix((sample (c(rep(1, each=N.assoc), rep(0, each=((NH*NP)-N.assoc))))), NH)
############ PARAFIT SIMUL ######################
#
P.parafit <- parafit(host.D, para.D, HP, nperm = N.perm, test.links = FALSE,
        seed = NULL, correction = "cailliez", silent = TRUE)$p.global
#
############ HOMMOLA ET AL. SIMUL ###############
#	

P.hommola <-  suppressWarnings(sim_cosp (host.D, para.D, t(HP), N.perm)) #uses the transposed HP matrix
													  # i.e., P in rows and H in columns
############ PACo SIMUL ###################	
m2.obs <- suppressWarnings(PACo (host.D, para.D, HP)) #computes sum of sq residuals
j = 0	
	for (n in c(1:N.perm))
	{set.seed(n+seed2)
	if (N.assoc <= NH | N.assoc <= NP) 	#control statement to avoid all parasites beig associated to 
		{	flag2 <- TRUE 					#a single host 
			while (flag2 == TRUE)	{ 
		HP.perm <- t(apply(HP,1,sample))
		if(any(colSums(HP.perm) == N.assoc)) flag2 <- TRUE else flag2 <- FALSE
									}  
		} else { HP.perm <- t(apply(HP,1,sample))} #permutes each row independently
		m2.perm <- suppressWarnings (PACo (host.D, para.D, HP.perm))
		if (m2.perm <= m2.obs)
		{ j = j + 1} 
	}
j <- j/N.perm
#
#######################
P.matrix101010 <- c(P.parafit, P.hommola, j)
#write (P.matrix101010, file = "TypeI_sim_example.txt", sep ="\t", append =TRUE)
}                    #writes to disk p-values of Parafit, HCT and PACo
proc.time() - ptm
#end.of.program

######Add PACo_function 
PACo_N.perm <- function(host.D, para.D, HP, N.perm) {
  
m2.obs <- suppressWarnings(PACo (host.D, para.D, HP)) #computes sum of sq residuals
j = 0	
NH = length(colnames(host.D))
NP = length(colnames(para.D))
N.assoc = sum(HP)
for (n in c(1:N.perm))
{set.seed(n+seed2)
  if (N.assoc <= NH | N.assoc <= NP) 	#control statement to avoid all parasites beig associated to 
  {	flag2 <- TRUE 					#a single host 
  while (flag2 == TRUE)	{ 
    HP.perm <- t(apply(HP,1,sample))
    if(any(colSums(HP.perm) == N.assoc)) flag2 <- TRUE else flag2 <- FALSE
  }  
  } else { HP.perm <- t(apply(HP,1,sample))} #permutes each row independently
  m2.perm <- suppressWarnings (PACo (host.D, para.D, HP.perm))
  if (m2.perm <= m2.obs)
  { j = j + 1} 
}
j <- j/N.perm
return(j)
}

PACo.P = PACo_N.perm(host.D, para.D, HP, N.perm)
