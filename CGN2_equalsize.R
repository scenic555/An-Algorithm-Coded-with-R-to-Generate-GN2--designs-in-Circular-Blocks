
################################################################################
# CGN2_equalsize: Circular Generalized neighbor design for block of equal 
# size(K)
################################################################################

# Algorithm from paper:

# Hurria Ali, Jamshaid ul Hassan, Khadija Noreen, Muhammad Sajid Rashid,
# Mahmood Ul Hassan and Rashid Ahmed (2022). Efficient Balanced and Strongly 
# Balanced Neighbor Designs through Circular GN2 -designs Generated with R 


# Coded by Ali et al., 01-03-2022 to 05-07-2022
# Version 1.4.0  (2022-07-05)
################################################################################




################################################################
# Division of adjusted A in i-1 groups of size k and one
# of size K-2 to get the set(s) of shifts
################################################################
grouping1<-function(A,k,v,i){
  bs<-c()
  z=0;f=1
  A1=A
  while(f<=(i-1)){
    
    for(y in 1:5000){
      comp<-sample(1:length(A1),k)
      com<-A1[comp]
      cs<-sum(com)
      if(cs%%v==0){
        bs<-rbind(bs,com)
        A1<-A1[-comp]
        z<-z+1
        f=f+1
      }
      if(z==(i-1)) break
    }
    if(z<(i-1)) {bs<-c();z=0;f=1;A1=A}  
    
  }
  
 
  bs1<-t(apply(bs,1,sort))
  bs1<-cbind(bs1,rowSums(bs),rowSums(bs)/(v))
  rownames(bs1)<-paste("G",1:(i-1), sep="")
  colnames(bs1)<-c(paste(1:k, sep=""),"sum" ,"sum/(v-1)")
  
  bs2<-t(apply(as.matrix(A1),1,sort))
  bs2<-cbind(bs2,rowSums(bs2),rowSums(bs2)/(v))
  rownames(bs2)<-paste("G",i, sep="")
  colnames(bs2)<-c(paste(1:(k-2), sep=""),"sum" ,"sum/(v-1)")
  
  
  fs1<-t(apply(bs,1,sort))
  fs1<-delmin(fs1)
  rownames(fs1)<-paste("S",1:(i-1), sep="")
  colnames(fs1)<-rep("",(k-1))
  
  fs2<-t(apply(as.matrix(A1),1,sort))
  rownames(fs2)<-paste("S",i, sep="")
  colnames(fs2)<-rep("",(k-2))
  
  
  list(B1=list(fs1,fs2),B2=list(bs1,bs2))
  }


#######################################################################
# Obtain set(s) of shifts by deleting smallest value of group
#######################################################################

delmin<-function(z){
  fs<-c()
  n<-nrow(z)
  c<-ncol(z)-1
  for(i in 1:n){
    z1<-z[i,]
    z2<-z1[z1!=min(z1)]
    fs<-rbind(fs,z2)
  }
  return(fs)
}


################################################################################
# Selection of adjusted A and the set(s) of shifts to obtain Circular 
# Generalized neighbor design for block of equal size.
################################################################################
#   K: Block sizes
#   i: Number of set of shifts for K


CGN2_equalsize<-function(k,i){
  
if(k<=3) stop("k= Block size: Block size must be greater than 4")
if(i<=1) stop("i= Must be a greater than 1")

setClass( "stat_test", representation("list"))
  
setMethod("show", "stat_test", function(object) {

row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")

cat("Following are required sets of shifts to obtain the 
minimal CGN2 for", "v=" ,object[[3]][1], "and","k=",object[[3]][2], "\n")
    

row <- paste(rep("=", 51), collapse = "")
      cat(row, "\n")
      print(object$S[[1]])
      cat("\n")
      print(object$S[[2]])
})


v=2*i*k-2; m=(v-2)/2

if(m%%8==0){
  j=m/8
  if(j<1) {return("Conditions are not satisfied for CGN2")}
  A=c(1:m,(2*m-j+1));A=A[-(j)]
   A1<-grouping1(A,k,(v-1),i)
   A2<-c(v,k);names(A2)<-c("V","K")
   x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}

if(m%%8==1){
  j=(m-1)/8
  if(j<1) {return("Conditions are not satisfied for CGN2")}
  A=c(1:(m-1),(m+1),(2*m-3*j));A=A[-(3*j+1)]
  A1<-grouping1(A,k,v=(v-1),i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
 }

if(m%%8==2){
  j=(m-2)/8
  if(j<1) {return("Conditions are not satisfied for CGN2")}
  A=c(1:(m-1),(m+1),(2*m-5*j-1));A=A[-(5*j+2)]
  A1<-grouping1(A,k,v=(v-1),i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}

if(m%%8==3){
  j=(m-3)/8
  if(j<0) {return("Conditions are not satisfied for CGN2")}
  A=c(1:m,(m+j+1));A=A[-(m-j)]
  A1<-grouping1(A,k,v=(v-1),i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}


if(m%%8==4){
  j=(m-4)/8
  if(j<0) {return("Conditions are not satisfied for CGN2")}
  A=c(1:(m-1),(m+1),(2*m-j));A=A[-(j+1)]
  A1<-grouping1(A,k,v=(v-1),i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}

if(m%%8==5){
  j=(m-5)/8
  if(j<0) {return("Conditions are not satisfied for CGN2")}
  A=c(1:m,(2*m-3*j-1));A=A[-(3*j+2)]
  A1<-grouping1(A,k,v=(v-1),i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}

if(m%%8==6){
  j=(m-6)/8
  if(j<0) {return("Conditions are not satisfied for CGN2")}
  A=c(1:m,(2*m-5*j-3));A=A[-(5*j+4)]
  A1<-grouping1(A,k,v=(v-1),i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}

if(m%%8==7){
  j=(m-7)/8
  if(j<1) {return("Conditions are not satisfied for CGN2")}
  A=c(1:(m-1),(m+1),(m+j+1));A=A[-(m-j)]
  A1<-grouping1(A,k,v=(v-1),i)
  A2<-c(v,k);names(A2)<-c("V","K")
  x<-list(S=A1$B1,G=A1$B2,R=A2,A=A)
}


new("stat_test", x)

}

##################################################################
# Generation of design using sets of cyclical shifts
###################################################################
# H is an output object from CGN2_equalsize
# The output is called using the design_CWBND to generate design

design<-function(H){
  
  setClass( "CWBND_design", representation("list"))
  setMethod("show", "CWBND_design", function(object) {
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    cat("Following is minimal CGN2 for", "v=" ,object$R[1], "and","k=",object$R[2], "\n")
    row <- paste(rep("=", 51), collapse = "")
    cat(row, "\n")
    for(i in 1:length(ss)){
      W<-ss[[i]]
      nr<-dim(W)[1]
      for(j in 1:nr){
        print(object$Design[[i]][[j]])
        cat("\n\n")
      }}
  })  
  
  v<-H$R[1]
  k<-H$R[2]
  ss<-H$S  
  treat<-(1:v)-1
  fn<-(1:v)
  G<-list()
  
  
  for(j in 1:length(ss)){ 
    W<-ss[[j]]
    nr<-dim(W)[1]
    nc<-dim(W)[2]
    D<-list()
    
    for(i in 1:nr){
      dd<-c()
      d1<-matrix(treat,(nc+1),v,byrow = T)
      ss1<-cumsum(c(0,W[i,]))
      dd2<-d1+ss1
      dd<-rbind(dd,dd2)
      rr<-dd[which(dd>=v)]%%v
      dd[which(dd>=v)]<-rr
      colnames(dd)<-paste("B",fn, sep="")
      rownames(dd)<-rep("",(nc+1))
      fn<-fn+v
      D[[i]]<-dd
    }
    G[[j]]<-D
    
  }
  
  x<-list(Design=G,R=H$R)
  new("CWBND_design", x)
}



################################################################################
# Examples: Using CGN2_equalsize function to obtain the set(s) of shifts
# for construction of Circular Generalized neighbor design for equal block  
# sizes (k)
################################################################################


# example#1
(H<-CGN2_equalsize(k=5,i=3))
(D<-design(H))


# example #2
(H<-CGN2_equalsize(k=7,i=3))
design(H)



# example #3
(H<-CGN2_equalsize(k=20,i=5))
design(H)


