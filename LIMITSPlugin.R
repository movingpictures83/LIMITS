#' @title LIMITS
#'
#' @description LIMITS is an algorithm developed by Fisher & Mehta to estimate the interaction matrix from time series data
#' assuming a Ricker model.
#'
#' @param x time series with taxa as rows and time points as columns
#' @param indices.perturb indices of species considered as environmental parameters
#' @param bagging.iter the number of iterations used for bagging
#' @param verbose print which taxon LIMITS is processing
#' @return a list with the estimated interaction matrix as Aest and the estimation error as error
#' @references Fisher & Mehta (2014). Identifying Keystone Species in the Human Gut Microbiome from Metagenomic Timeseries using Sparse Linear Regression. PLoS One 9:e102451
#' @examples
#' \dontrun{
#' N=20
#' A=generateA(N,c=0.1)
#' ts=ricker(N=N,A=A)
#' Aest=limits(ts,verbose=TRUE)$Aest
#' par(mfrow=c(2,1))
#' plotA(A,header="original")
#' plotA(Aest,header="estimated")
#' par(mfrow=c(1,1))
#' }
#' @export

#' @title Generate an interaction matrix
#'
#' @description Generate an interaction matrix, either randomly from a uniform distribution or
#' using Klemm-Eguiluz algorithm to generate a modular and scale-free interaction matrix.
#'
#' @param N number of species
#' @param type random (sample a uniform distribution), klemm (generate a Klemm-Eguiluz matrix) or empty (zero everywhere, except for diagonal which is set to d)
#' @param pep desired positive edge percentage (only for klemm)
#' @param d diagonal values (should be negative)
#' @param min.strength random: minimal off-diagonal interaction strength (only for random)
#' @param max.strength random: maximal off-diagonal interaction strength, klemm: maximal absolute off-diagonal interaction strength
#' @param c desired connectance (interaction probability)
#' @param ignore.c do not adjust connectance
#' @param negedge.symm set symmetric negative interactions (only for klemm)
#' @param clique.size modularity parameter (only for klemm)
#' @param groups vector of group memberships for each species, assign NA if species does not belong to any group (only for random)
#' @param intra.group.strength interaction strength between members of the same group
#' @param inter.group.strength interaction strength between members of different groups (if not defined, will be assigned randomly)
#' @return the interaction matrix
#' @examples
#' klemm=generateA(N=10,type="klemm",c=0.5)
#' groups=c(rep(NA,5),rep(1,10),rep(2,5),rep(3,10),rep(4,10))
#' Agroup=generateA(N=40,groups=groups,c=0.5,intra.group.strength=0.1,inter.group.strength=-0.5, d=-1)
#' @references Klemm & Eguiluz, Growing Scale-Free Networks with Small World Behavior \url{http://arxiv.org/pdf/cond-mat/0107607v1.pdf}
#' @export

#' @title Modify the interaction matrix
#'
#' @description The interaction matrix can be manipulated in the following ways:
#' adjustc: Adjust the connectance to the given value
#' schur: apply a Schur decomposition to remove positive eigenvalues that lead to explosions in simulations with Ricker or gLV
#' negpercent: convert positive into negative interactions, such that the interaction matrix reaches the specified percentage of negative edges/arcs (diagonal is not specially treated)
#' tweak: multiply a randomly chosen positive entry in the interaction matrix with -1. This is useful when searching for an interaction matrix that does not lead to explosions in Ricker and gLV simulations.
#' enforceneg: multiply off-diagonal negative entries in interaction matrix with the given factor
#' removeorphans: remove all taxa that do not interact with other taxa and thus represent orphan nodes in the network representation of the interaction matrix
#' mergeposlinks: remove redundant taxon names by summing replicate positive arcs and keeping the sum as entries in the matrix (negative arcs are ignored)
#' mergeneglinks: same as mergeposlinks, but only negative entries are kept
#' mergelinks: All entries are kept. Negative entries are subtracted, positive entries are added.
#' The modes mergeposlinks, mergeneglinks and mergelinks expect a first column with the taxon names. They will return a matrix with row and column names representing unique taxa.
#'
#' @details An entry in the interaction matrix represents an arc in a directed network. An interaction involves two
#' entries in the interaction matrix, which represent the influence of species A on species B and vice versa.
#' Mode negpercent: By default, the negative arc percentage is adjusted, i.e. the percentage of
#' negative entries in the interaction matrix. If symmetric is true, both the forward and the reverse arc of an interaction
#' will be set to a negative value (corresponding to a competition edge), else only one arc will be set to a negative value
#' (corresponding to an edge representing parasitism, predation or amensalism).
#'
#' @param A the interaction matrix
#' @param mode modification mode, values: adjustc, schur, negpercent, tweak, enforceneg, removeorphans, mergeposlinks, mergeneglinks, mergelinks
#' @param strength interaction strength, binary (0/1) or uniform (sampled from uniform distribution from minstrength to 1)
#' @param factor multiplication factor for enforceneg mode
#' @param minstrength minimum interaction strength for uniform mode (maximum is 1)
#' @param c the target connectance (only for mode adjustc)
#' @param perc negative edge percentage (only for mode negpercent)
#' @param symmetric only introduce symmetric negative interactions (only for mode negpercent)
#' @return the modified interaction matrix
#' @export

#' @title Connectance
#' @description Compute connectance of an interaction matrix.
#' @details The connectance is defined as \eqn{c=E/(N*(N-1))}, where \eqn{E} is the number
#' of realized arcs (the number of non-zero entries in the interaction matrix)
#' and \eqn{N*(N-1)} the number of possible arcs.
#' The diagonal (self-arcs) is excluded.
#'
#' @param A an interaction matrix
#' @return the connectance
#' @examples
#'   A <- cbind(c(-1,0,1),c(0,-1,0),c(-1,0,-1))
#'   x <- getConnectance(A)
#' @export


ricker<-function(N, A, K=rep(0.1,N), y=runif(N), sigma=0.05, K.trend=NA, tend=100, death.t=NA, tskip=0, explosion.bound=10^8, perturb=NULL){
  if(length(y) != N){
    stop("y needs to have N entries.")
  }
  if(nrow(A)!=N || ncol(A)!=N){
    stop("A needs to have N rows and N columns.")
  }
  if(length(K)!=N){
    stop("K needs to have N entries.")
  }
  if(length(K.trend)>1 && length(K.trend)!=N){
    stop("K.trend needs to have N entries.")
  }
  if(tskip>=tend){
    stop("There are as many or more time points to skip than time points specified.")
  }
  out=matrix(nrow=N, ncol=tend-tskip)
  out[,1]=y
  perturbCounter=1
  durationCounter=1
  K.copy=K
  perturbationOn=FALSE
  # simulate difference equation
  for(t in 2:tend){
    if(sigma > 0){
      b=rlnorm(N,meanlog=0,sdlog=sigma)
    }else{
      b=rep(1,N)
    }
    if(!is.null(perturb)){
      if(perturb$times[1]==1){
        stop("Please do not specify a perturbation at the first time point.")
      }
      applied=applyPerturbation(perturb,t=t,perturbCounter=perturbCounter,durationCounter=durationCounter,perturbationOn=perturbationOn,ori.growthrates=K.copy,abundances=y)
      y=applied$abundances
      K=applied$growthrates
      durationCounter=applied$durationCounter
      perturbCounter=applied$perturbCounter
      perturbationOn=applied$perturbationOn
      #print(perturbCounter)
      #print(perturbationOn)
    }
    if(length(K.trend)>1){
      # calculate percentages to be added
      K.onepercent=K.copy/100
      K.percent=K.trend*100*K.onepercent
      K=K+K.percent
      # set negative carrying capacities to zero
      negative.K.indices=which(K<0)
      K[negative.K.indices]=0
    }
    y=b*y*exp(A%*%(y-K))
    if(max(y) > explosion.bound){
      # report which species explodes
      print("Explosion!")
      res=c(-1,which(y==max(y)))
      return(res)
    }
    if(!is.na(death.t) && death.t>0){
      y[y<death.t]=0 # kill species below the threshold
    }
    else if(length(y[y<0]) > 0){
      stop("Species below 0!")
    }
    if(t > tskip){
      out[,t-tskip]=y
    }
  }
  return(out)
}

getConnectance <- function(A){

	 N <- nrow(A)

	 # exclude diagonal from observed and possible interactions
	 c <- (length(A[A!=0])-N)/(ncol(A)*ncol(A)-N)

	 return(c)
 }

modifyA<-function(A, mode="adjustc", strength="binary", factor=2, minstrength=0.1, c=0.2, perc=50, symmetric=FALSE){
  edgeNumAdded = 0
  print(paste("Initial edge number", length(A[A!=0])))
  c_obs = getConnectance(A)
  print(paste("Initial connectance", c_obs))
  if(mode == "adjustc"){
    if(c_obs < c){
      while(c_obs < c){
        # randomly select source node of edge
        xpos=sample(c(1:ncol(A)))[1]
        # randomly select target node of edge
        ypos=sample(c(1:ncol(A)))[1]
        # avoid diagonal
        if(xpos != ypos){
          # count as added if there was no edge yet
          if(A[xpos,ypos]==0){
            edgeNumAdded = edgeNumAdded+1
          }
          # add edge
          A[xpos,ypos]=getStrength(strength=strength,pos=TRUE, minstrength=minstrength)
          c_obs=getConnectance(A=A)
        }
      }
      print(paste("Number of edges added", edgeNumAdded))
    }else if(c_obs > c){
      edgeNumRemoved = 0
      while(c_obs > c){
        xpos=sample(c(1:ncol(A)))[1]
        ypos=sample(c(1:ncol(A)))[1]
        # avoid diagonal
        if(xpos != ypos){
          # count as removed if there was an edge before
          if(A[xpos,ypos]!=0){
            edgeNumRemoved = edgeNumRemoved+1
          }
          # remove edge
          A[xpos,ypos]=0
          c_obs = getConnectance(A)
        }
      }
      print(paste("Number of edges removed", edgeNumRemoved))
    }
    print(paste("Final connectance", c_obs))
  }else if(mode=="mergeposlinks" || mode=="mergeneglinks" || mode=="mergelinks"){
    taxonnames=as.character(A[,1])
    entries=unique(taxonnames)
    # remove taxon name column
    A=A[,2:ncol(A)]
    mergedlinks=matrix(0,nrow=length(entries),ncol=length(entries))
    rownames(mergedlinks)=entries
    colnames(mergedlinks)=entries
    for(i in 1 : nrow(A)){
      for(j in 1 : ncol(A)){
        if(!is.na(A[i,j]) && A[i,j]!=0){
          merge=FALSE
          xIndex=which(entries==taxonnames[i])
          yIndex=which(entries==taxonnames[j])
          #print(paste("x index:",xIndex))
          #print(paste("y index:",yIndex))
          if(mode=="mergeposlinks" && A[i,j]>0){
            merge=TRUE
          }else if(mode=="mergeneglinks" && A[i,j]<0){
            merge=TRUE
          }else if(mode=="mergelinks"){
            merge=TRUE
          }
          if(merge==TRUE){
            if(mode=="mergelinks"){
              if(A[i,j]<0){
                mergedlinks[xIndex,yIndex]=mergedlinks[xIndex,yIndex]-1
              }else{
                mergedlinks[xIndex,yIndex]=mergedlinks[xIndex,yIndex]+1
              }
            }else{
              mergedlinks[xIndex,yIndex]=mergedlinks[xIndex,yIndex]+1
            }
          }
        } # interaction is not zero
      } # column loop
    } # row loop
    A=mergedlinks
  }else if(mode=="removeorphans"){
    # since A can be asymmetric, only those species can be removed for which rows and columns are simultaneously zero (except for diagonal)
    toKeep=c()
    diagvals=diag(A)
    diag(A)=0
    for(i in 1:nrow(A)){
      rowsum=sum(abs(A[i,]))
      colsum=sum(abs(A[,i]))
      if(rowsum != 0 || colsum!=0){
        toKeep=append(toKeep,i)
      }
    }
    A=A[toKeep,toKeep]
    diag(A)=diagvals[toKeep]
  }else if(mode == "negpercent"){
    # arc number: number of non-zero entries in the interaction matrix
    num.edge=length(A[A!=0])
    num.perc=(num.edge/100)*perc
    # subtract the negative edges that are already there
    num.neg.edge=length(A[A<0])
    print(paste("Number of negative edges already present:",num.neg.edge))
    if(num.neg.edge>num.perc){
      warning("The matrix has more negative edges than are required to reach the desired negative edge percentage!")
    }else if(num.neg.edge==num.perc){
      print("The matrix has already the desired negative edge percentage.")
    }else{
      # those negative edges already present do not need to be added
      num.perc=num.perc-num.neg.edge
      # symmetric interactions: we will count negative edges, not arcs
      if(symmetric == TRUE){
        num.perc=round(num.perc/2)
      }else{
        num.perc=round(num.perc)
      }
      print(paste("Converting",num.perc,"edges into negative edges",sep=" "))
      indices=which(A>0,arr.ind=T)
      # randomly select indices
      rand=sample(c(1:nrow(indices)))
      xyposAlreadySeen = c()
      counter = 0
      # loop over number of negative edges to introduce
      for(i in 1:num.perc){
        xpos=indices[rand[i],1]
        ypos=indices[rand[i],2]
        xypos=paste(xpos,"_",ypos, sep="")
        yxpos=paste(ypos,"_",xpos,sep="")
        # if we find an index pair that was already used, we have to look for another index pair,
        # since using the same index pair means to use the same arc or the same arc in reverse direction
        if(symmetric == TRUE && is.element(xypos,xyposAlreadySeen) == TRUE){
          xpos = indices[rand[nrow(indices)-counter],1]
          ypos = indices[rand[nrow(indices)-counter],2]
          counter = counter + 1
          if((num.perc + counter) > nrow(indices)){
            stop("More negative edges requested than can be set!")
          }
        }
        xyposAlreadySeen = c(xypos, yxpos, xyposAlreadySeen)
        # print for tests
        # print(paste("xpos",xpos,"ypos",ypos,"value:",A[xpos,ypos],sep=" "))
        negval=getStrength(strength=strength,pos=FALSE,minstrength=minstrength)
        A[xpos,ypos]=negval
        if(symmetric == TRUE){
          A[ypos,xpos]=negval
        }
        #print(paste("xpos",xpos,"ypos",ypos,"value:",A[xpos,ypos],sep=" "))
        #print(paste("reverse value:",A[ypos,xpos],sep=" "))
      }
    }
  }else if(mode == "tweak"){
    # check that positive arcs are present
    if(length(A[A>0]) > 0){
      # row and column indices of positive arcs
      indices.pos = which(A>0,arr.ind=TRUE)
      # randomly select a positive arc
      x=sample(1:nrow(indices.pos),1)
      # convert positive arc into negative one, keeping the same interaction strength
      A[indices.pos[x,1],indices.pos[x,2]]=A[indices.pos[x,1],indices.pos[x,2]]*(-1)
    }else{
      warning("Cannot tweak. No positive arc in the given matrix.")
    }
  }else if(mode == "enforceneg"){
    diag=diag(A)
    indices.neg = which(A<0,arr.ind=TRUE)
    # multiply negative entries by given factor
    A[indices.neg]=A[indices.neg]*factor
    # keep original diagonal
    diag(A)=diag
  }else if(mode == "schur"){
    # remove positive real parts of eigenvalues if any (using schur decomposition)
    sd<-dim(A)

    if(max(Re(eigen(A)$values))){
      # division by max.A helps removing positive eigenvalues
      max=max(A)
      A=A/max

      diagt<-diag(sd[2])+0i

      # Computes the generalized eigenvalues and Schur form of a pair of matrices.
      # R=imaginary part identical to 0 with a tolerance of 100*machine_precision as determined by Lapack
      schur.A<-geigen::gqz(A,diagt,"R")
      # generalized inverse of a matrix
      T<-schur.A$S%*%MASS::ginv(schur.A$T)
      rediag<-Re(diag(T))
      imdiag<-Im(diag(T))

      indicesP=rediag>0
      listind=1:sd[2]

      for(k in listind[indicesP]){
        T[k,k]<- complex(real=-Re(T[k,k]),imaginary=Im(T[k,k]))
      }

      A <- schur.A$Q %*% T %*% MASS::ginv(schur.A$Q)
      A<-Re(A)
      A=A*max
    }

  }else{
    stop(paste("Mode",mode,"not known."))
  }
  c=getConnectance(A)
  print(paste("Final connectance:",c))
  return(A)
}

################## helper functions ################

# Get the interaction strength.
getStrength<-function(strength="binary", minstrength=0.1, pos=TRUE){
  value = NA
  if(strength=="binary"){
    value = 1
  }else if(strength == "uniform"){
    value = runif(1,min=minstrength,max=1)
  }
  if(!pos){
    value = -1*value
  }
  return(value)
}

# Check whether the interaction matrix is fully connected (entirely filled with 1)
isFullyconnected<-function(A){
  if(length(A[A!=0])==nrow(A)*nrow(A)){
    return(TRUE)
  }else{
    return(FALSE)
  }
}

generateA<-function(N=100, type="random",pep=50, d=-0.5, min.strength=-0.5, max.strength=0.5, c=0.02, ignore.c=FALSE, negedge.symm=FALSE, clique.size=5, groups=c(), intra.group.strength=0.5, inter.group.strength=NA){
  A=matrix(0,nrow=N,ncol=N)      # init species interaction matrix
  if(type=="random"){
    if(length(groups)>0){
      if(length(groups)!=nrow(A)){
        stop("Please define a group membership for each species.")
      }
    }
    for (i in 1:N){
      for(j in 1:N){
        if(i==j){
          A[i,j]=d
        }else{
          if(length(groups)==0){
            A[i,j] = runif(1,min=min.strength,max=max.strength)
          }else{
            group1=groups[i]
            group2=groups[j]
            if(!is.na(group1) && !is.na(group2) && group1==group2){
              A[i,j] = intra.group.strength
            }else{
              # assign interaction strength between groups randomly
              if(is.na(inter.group.strength)){
                A[i,j] = runif(1,min=min.strength,max=max.strength)
              }else{
                # assign selected interaction strength between groups
                A[i,j] = inter.group.strength
              }
            }
          }
        }
      }
    }
  }else if(type=="empty"){
    diag(A)=d
  }else if(type=="klemm"){
    g<-klemm.game(N,verb=FALSE, clique.size)
    A=get.adjacency(g)
    A=as.matrix(A)
    diag(A)=d
  }

  if(ignore.c==FALSE){
    print(paste("Adjusting connectance to",c))
    A=modifyA(A,c=c, mode="adjustc")
  }

  # for Klemm-Eguiluz: inroduce negative edges and set random interaction strengths
  if(type=="klemm"){
    if(pep < 100){
      A=modifyA(A=A, perc=(100-pep), symmetric=negedge.symm, mode="negpercent")
    }
    # excluding edges on the diagonal
    print(paste("Final arc number (excluding self-arcs)", length(A[A!=0])-N ))
    # excluding negative edges on the diagonal
    print(paste("Final negative arc number (excluding self-arcs)", length(A[A<0])-N ))

    # check PEP and number of asymmetric negative interactions
    # assuming diagonal values are negative
    pep = getPep(A)
    print(paste("PEP:",pep))

    # convert binary interaction strengths (-1/1) into continuous ones using uniform distribution
    # zero would remove the edge, so the minimum strength is small, but non-zero
    min.klemm.strength=0.00001
    for(i in 1:nrow(A)){
      for(j in 1:nrow(A)){
        # skip diagonal
        if(i != j){
          A[i,j]=A[i,j]*runif(1,min=min.klemm.strength,max=max.strength)
        }
      }
    }
  }
  return(A)
}

limits<-function(x, indices.perturb = NA, bagging.iter=2, verbose=TRUE){
  x=t(x)
  print(paste("Time series has",ncol(x),"taxa"))
  if(verbose==TRUE){
    print("Processing first taxon.")
  }
  results<-limitscolumnwise(x, 1, indices.perturb, r=bagging.iter)
  Aest<-t(results[[1]])
  error<-results[[2]]
  # loop taxa
  for(i in 2:ncol(x)){
    if(verbose==TRUE){
      print(paste("Processing taxon",i))
    }
    results<-limitscolumnwise(x, i, indices.perturb, r=bagging.iter)
    Aest.temp<-results[[1]]
    error.temp<-results[[2]]
    Aest<-cbind(Aest,t(Aest.temp))
    error<-append(error,error.temp)
  }

  res=list(t(Aest),error)
  names(res)=c("Aest","error")

  return(res)
}

#### R is a time series with the columns corresponding to the species N, and
#### lines to different time points Ntp. If the first column of your time series
#### i.e. dim(R)=Ntp N
#### labels time, remove it before applying limits (ie, do R<-R[:,2:end] ).
#### i is the column of the interaction matrix we try to reconstruct.

#### output
#### Beval : the output of the function is the estimation of the "ith" line of the interaction matrix
#### error : mean of the errors made on evaluation of "y" using the estimated B matrix normalized by the variance of y ("one-time-step evalutaion").
#### errorF : idem as error but without the normalization by the variance.

# by Sophie de Buyl, translated from a Mathematica code provided by Charles Fisher and Pankaj Mehta.
limitscolumnwise <- function(R,i,indices.perturb, r=100){

  listnumbkeysp<-c(); #list of number species kept. NOT MANDATORY

  # choices to be made:
  thresh <- .5 # orignially put to 5 (diminish to increase precision)

  # manipulating R to put the data in the appropriated form to apply limits:

  R[R <= 0] <- 2.22507e-308  # we will have to take the log, so let's remove zero's.
  sd<-dim(R)
  N<-sd[2] #number of species
  Ntp<-sd[1]-1 #number of time points

  # comput medians column-wise
  colMedians=apply(R,2,median)
  # the median should not be substracted for the perturbations
  if(all(is.na(indices.perturb)) == FALSE){
    colMedians[indices.perturb]=rep(0,length(indices.perturb))
  }
  # formulation with dependency on matrixStats
  #data<- R-t(kronecker(matrix(1,1,sd[1]),t(t(colMedians(R)))));#first N column of data matrix needed for limits
  data<- R-t(kronecker(matrix(1,1,sd[1]),colMedians))
  data<-data[1:(sd[1]-1),]
  data<-cbind(data,(log(R[2:sd[1],i])-log(R[1:(sd[1]-1),i])))
  if(sum(is.na(log(R)))>0){
    print(which(is.na(log(R)),arr.ind = T))
  }
  # variable initiation
  res<-array(0,c(N,1)) # array for storing results
  errorlist<-c()
  listspecies<-seq(1,N) # to construct the choices of excluded/includes species

  for(k in 1:r){ # we do r times the same thing to get r estimations the interaction matrix B

    #initialize covariates to be included/excluded from regression model
    if (i!= 1 & i!=N){
      c1<-1:(i-1)
      c2<-(i+1):N
      excluded <-  c(c1,c2)
      included <- i
    }

    if(i==1){
      excluded<-2:N
      included<-i
    }
    if(i==N){
      excluded<-1:(N-1)
      included<-N
    }

    #randomly partition data into training and test sets

    trainset <- as.matrix(sample(1:Ntp,round(Ntp/2)))
    testset <- as.matrix(setdiff(1:Ntp,trainset))

    data1 <- data[trainset,]
    data2 <- data[testset,]

    test <- included
    results<- restrictedleastsquare(data1,data2,test) #perform the estimation
    errorEt <-results[[1]]
    B1t<-results[[2]]

    # loop that adds covariates(=species) as long as prediction error decreases by thresh
    xxx=1
    yyy=1

    while(xxx == 1 && yyy != N){

      yyy=yyy+1

      #loop to create list of performances for all regression models including an additional covariate

      #initial the loop
      test <- c(included,excluded[1])
      errorlist<-restrictedleastsquare(data1,data2,test)[[1]]
      #the loop
      for(kk in 2:(length(excluded)-1)){
        test = c(included,excluded[kk])
        errorlist<-c(errorlist,restrictedleastsquare(data1,data2,test)[[1]])
      }

      #sort the list so that prev[[1]] has the lowest prediction error
      ind<-which(errorlist == min(errorlist), arr.ind = TRUE)
      if(dim(as.matrix(ind))[1]>1){
        ind<-ind[[1]]
      }
      test <- c(included,excluded[ind])  #we choose the test set with the smallest error
      resulttemp <- restrictedleastsquare(data1, data2,test) #we re-run limits for that test set since we didn't store results
      errorEtemp <-resulttemp[[1]]
      B1temp<-resulttemp[[2]]
      # redefine included excluded covariates to account for including new covariate
      included <- test
      excluded <- setdiff(union(listspecies,test),intersect(listspecies,test))

      # if improvement over the current best model by thresh, include new species, otherwise exit loop
      if(is.na(errorEt)==FALSE && errorEt!=0 && (100.0*(errorEtemp - errorEt)/errorEt< -thresh)){ #we keep adding species
        errorEt <- errorEtemp# we update the error to be compared with.
        B1t <- B1temp #useless
        errorEtempt <- errorEtemp #useless
        testt <- test
      } else{ # we stop adding species
        xxx<-0
        listnumbkeysp<-c(listnumbkeysp,yyy)
      }
    } # end of while loop

    # store final regression coefficients in res
    B <- matrix(0,N,1)
    B[test]<-t(B1temp)

    res<-cbind(res, B)
    errorlist<-c(errorlist,errorEtemp)

  } #end or r loop

  # Bagging step: output median of res

  res<-res[,-1]

  # compute medians row-wise
  rowMedians=apply(res,1,median)
  Beval<-t(as.matrix(rowMedians))
  error.median<-median(errorlist)
  #listnumbkeysp

  return(list(Beval,error.median))
}

##############################################################
# restrictedleastsquare
# function need for limits.r
#
# by Sophie de Buyl, translated from a Mathematica code provided by Charles Fisher and Pankaj Mehta.
################################################################
restrictedleastsquare <- function(inbag,outbag,test){

  # This subroutine performs the linear regression and computes the error on the test set
  #"inbag" is a matrix containing the training set. One row of the matrix looks like {x_1(t) - u_1, ..., x_N(t) -
  # u_N, ln x_i(t+1) - ln x_i(t) },
  #where u_i is the median abundance of species i. ;
  #"outbag" is a matrix containing the test set.
  #It is in the same form as inbag.;
  #"test" is a vector of integers specifying which covariates are included in the regression. I.e. if species 1,2,3 are included then test = {1,2,3}.

  #perform linear regression on training set
  temp<-dim(inbag)
  lastcol<-temp[2]
  X1 <- as.matrix(inbag[,test])
  y1 <- as.matrix(inbag[,lastcol])
  B1 <- as.matrix(MASS::ginv(X1,tol=2e-308)%*%y1)

  # calculate prediction error on test set
  X2 <- outbag[test]#as.matrix(outbag[,test])
  y2 <- outbag[lastcol]#as.matrix(outbag[,lastcol])
  errorE <- mean((y2-X2%*%B1)*(y2-X2%*%B1))/var(y2)
  result <- list(errorE,B1)
  return(result)
}

#N=40
#A=generateA(N,c=0.1)
#ts=ricker(N=N,A=A)
#write.table(ts, "timeseries.csv", sep=',', row.names=TRUE, col.names=TRUE);
#Aest=limits(ts,verbose=TRUE)$Aest
#print(Aest)


input <- function(inputfile) {
   pc <<- read.csv(inputfile, header = TRUE);
}

run <- function() {
cn <<- colnames(pc);
cn <<- cn[2:length(cn)];
pc <<- pc[,-1];
pc <<- apply(pc, 1, as.numeric);
result <<- limits(pc);
}

output <- function(outputfile) {
   write.table(result$Aest, file=outputfile, sep=",", append=FALSE, row.names=unlist(cn), col.names=unlist(cn), na="")
}

#pc <<- read.csv("Early.abund.norm.csv", header = TRUE);
#cn <<- colnames(pc);
#cn <<- cn[2:length(cn)];
#pc <<- pc[,-1];
#pc <<- apply(pc, 1, as.numeric);
#result <<- limits(pc);
#write.table(result, file="output.csv", sep=",", append=FALSE, row.names=unlist(cn),  na="")
#print(result);


#  pc <<- t(pc);
#  correlations <<- rcorr(pc[,], type=c("pearson"));
#  pc <<- as.matrix(correlations$r);
#  pc[is.na(pc)] <<- 0;
#  empty <- c("");
#  pc[which(correlations$P>p_value)] <<- 0;
#}

#output <- function(outputfile) {
#   write.table(pc, file=outputfile, sep=",", append=FALSE, row.names=unlist(cn), col.names=unlist(cn), na="");
#}
  
