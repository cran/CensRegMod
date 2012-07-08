#### EM Algorithm  ####

em.cens <- function(cc,x,y,nu,type,diagnostic,typediag){

  # cc: vector of censoring indicators: 0 if non-censored, 1 if censored
  # y: vector of responde variable
  # x: matrix or vector of covariables
  # nu: initial value for the degree of freedon in case of Student-t model
  # type: Normal or T distribution to be used
  # diagnostic: TRUE ou FALSE - should compute diagnostic measures? 
  # typediag: in case diagnostic=TRUE, what kind of diagnostic?
  	# typediag=1 => Case delection measure: Generalized Cook Distance (GD) 
  	# typediag=2 => Local influence: case-weight perturbation
  	# typediag=3 => Local influence: scale perturbation
	

	## Verify error at parameters specification

	if(ncol(as.matrix(y)) > 1) stop("Only univariate linear regression supported!")
	if( length(y) != nrow(as.matrix(x)) ) stop("X variable does not have the same number of lines than y")
	if( (length(x) == 0) | (length(y) == 0) ) stop("All parameters must be provided.")
	if( (type != "T") && (type != "Normal") ) stop("Distribution family not supported. Check documentation!") 
	#if( (length(nu) != 0) && (type != "T") ) stop("For the Normal distribution the Nu parameters must not be provided.") 
	if( type == "T" ){
  		if(length(nu) > 1) stop("Nu parameters must have only one parameter")
  		if(length(nu) == 0) stop("Nu parameters must be provided.")
  		if(nu <= 0) stop("Nu parameters must be positive.")
  	} 
	if( (diagnostic=="TRUE") && (length(typediag)==0) ) stop("typediag parameter must be provided.")
	if( (diagnostic=="TRUE") && ( (typediag!=1) && (typediag!=2) && (typediag!=3) )) stop("typediag must be 1, 2 or 3")



	x<-cbind(1,x)
 	p <- ncol(x)
	n <- nrow(x)
	reg <- lm(y ~ x[,2:p])

	#Parameters
	ERRO<-1e-6
	TOLERANCIA<-1e-6
	MAX_NU<-150
	MIN_NU <- 2.1
	
	#Initialize beta and sigma2 with least squares estimatives
	
	beta<- as.vector(coefficients(reg),mode="numeric")
	sigma2 <- sum((y-x%*%beta)^2)/(n-p)
	
	if(type=="T"){
 
		teta_velho <- matrix(c(beta,sigma2,nu),ncol=1)
 		cont <- 0
	
		repeat{
        		mu<-x%*%beta
  		  	u1<-y*(nu+1)/(nu+(y-mu)^2/(sigma2)) 
        		u2<-y^2*(nu+1)/(nu+(y-mu)^2/(sigma2))
        		u0<-(nu+1)/(nu+(y-mu)^2/(sigma2))  
               
      			aux1<-CalMom(mu,sigma2,nu,y,type)
      			aux2t<-CalMomMIt(mu,sigma2,nu,y)
      
      			aux2t0<-u0^2
      			aux2t0[cc==1]<- aux2t$Ey0[cc==1]
      
      			aux2t1<-u0^2*y
      			aux2t1[cc==1]<- aux2t$Ey1[cc==1]
              
     				aux2t2<-u0^2*y^2
      			aux2t2[cc==1]<-aux2t$Ey2[cc==1]
      
      			vari<- (aux2t2-2*mu*aux2t1+mu^2*aux2t0)-((aux2t1-mu*aux2t0)^2)
      
      
      			u1[cc==1]<- aux1$Ey1[cc==1]
      			u2[cc==1]<- aux1$Ey2[cc==1]
      			u0[cc==1]<- aux1$Ey0[cc==1]
      
  				suma1<-matrix(0,p,p)
      			suma2<-matrix(0,p,1)
      			suma3<-matrix(0,p,p)

  			for(i in 1:n){
  				suma1<-suma1+u0[i]*((x[i,])%*%t(x[i,]))
      				suma2<-suma2+(x[i,])*u1[i]
     				suma3 <- suma3 + (((nu+1)/(nu+3))*(x[i,]%*%t(x[i,]))/sigma2) - (vari[i]*x[i,]%*%t(x[i,])/(sigma2^2))
  			}

  			beta<-solve(suma1)%*%suma2
  			sigma2<-sum(u2-2*u1*mu+mu^2*u0)/n
  		  		
  	 		auxf<-(y-x%*%beta)/sqrt(sigma2)
 		   	f<-function(nu){sum(log(dt(auxf[cc==0],df=nu)/sqrt(sigma2)))+sum(log(pt(-auxf[cc==1],df=nu)))}
       			nu<-optimize(f, c(MIN_NU,MAX_NU), tol = TOLERANCIA, maximum = TRUE)$maximum 

    			teta_novo<-matrix(c(beta,sigma2,nu),ncol=1)

			#In case of convergence, stop the iterations

			if (sqrt(sum((teta_velho-teta_novo)^2)) < ERRO ){break}
		
			#If there´s no convergence, update the parameters estimatives
	
   			teta_velho <- teta_novo
			
		}
  
 	 	logver<-f(nu)
	}
	
	
	if(type=="Normal"){
   
		teta_velho <- matrix(c(beta,sigma2),ncol=1)
 		cont <- 0
	
		repeat{
        		mu<-x%*%beta
 		  
        		u1<-y 
        		u2<-y^2  
       		u0<-rep(1,n)       
      		aux1<-CalMom(mu=mu,sigma2=sigma2,a=y,type="Normal")
      		vari <- (aux1$Ey2 - ((aux1$Ey1)^2))
      		vari[cc==0]<-0
      
      		u1[cc==1]<- aux1$Ey1[cc==1]
      		u2[cc==1]<- aux1$Ey2[cc==1]
      		u0[cc==1]<- aux1$Ey0[cc==1]

      	  	suma1<-matrix(0,p,p)
      		suma2<-matrix(0,p,1)
      		suma3<-matrix(0,p,p)

  			for(i in 1:n){
  		      		suma1<-suma1+u0[i]*((x[i,])%*%t(x[i,]))
            			suma2<-suma2+(x[i,])*u1[i]
            			suma3<-suma3+(x[i,])%*%t(x[i,])/sigma2-((x[i,])%*%t(x[i,])/sigma2^2)*vari[i]
  			}

  			beta<-solve(suma1)%*%suma2
  			sigma2<-sum(u2-2*u1*mu+mu^2*u0)/n

  			teta_novo<-matrix(c(beta,sigma2),ncol=1)
  		
  		 	auxpdf<-dnorm(y,x%*%beta,sqrt(sigma2))
 		  	auxcdf<-pnorm(-(y-x%*%beta)/sqrt(sigma2))
 		   	logver<-sum(log(auxpdf[cc==0]))+sum(log(auxcdf[cc==1]))

			#In case of convergence, stop the iterations

			if (sqrt(sum((teta_velho-teta_novo)^2)) < ERRO ){break}
		
			#If there´s no convergence, update the parameters estimatives
	
			teta_velho <- teta_novo
		
				
		} 
	}
	
	if(diagnostic=="TRUE"){
		
		Matriz<-HessianaQ(cc=cc,x=x,y=y,beta=beta,sigma2=sigma2,u0=u0,u1=u1,u2=u2)
 		MatrizQ<-Matriz$MatrizQ
 		MatrizQbeta<-Matriz$MatrizQbeta
 		MatrizQsigma2<-Matriz$MatrizQsigma2

		if(typediag==1){

			CDm<-CD(cc=cc,x=x,y=y,beta=beta,sigma2=sigma2,u0=u0,u1=u1,u2=u2,MatrizQ=MatrizQ,MatrizQbeta=MatrizQbeta,MatrizQsigma2=MatrizQsigma2)
			GD<-CDm$GD
			GDbeta<-CDm$GDbeta
			GDsigma2<-CDm$GDsigma2

			medida<-data.frame(GD,GDbeta,GDsigma2)

		}

		if(typediag==2){

			If1<-If(cc=cc,x=x,y=y,beta=beta,sigma2=sigma2,u0=u0,u1=u1,u2=u2,MatrizQ=MatrizQ,Perturb=1)
			
			medida<-diag(If1)/sum(diag(If1))
		}

		if(typediag==3){

			If2<-If(cc=cc,x=x,y=y,beta=beta,sigma2=sigma2,u0=u0,u1=u1,u2=u2,MatrizQ=MatrizQ,Perturb=2)
			medida<-diag(If2)/sum(diag(If2))
		}


		if(type=="T"){
			return(list(beta=teta_novo[1:p],sigma2=teta_novo[p+1],nu=round(teta_novo[p+2]),logver=logver,dp=sqrt(diag(solve(suma3))),measure=medida))
		}

		if(type=="Normal"){
			return(list(beta=teta_novo[1:p],sigma2=teta_novo[p+1],logver=logver,dp=sqrt(diag(solve(suma3))),measure=medida))
		}





	}

	else{

		if(type=="T"){
			return(list(beta=teta_novo[1:p],sigma2=teta_novo[p+1],nu=round(teta_novo[p+2]),logver=logver,dp=sqrt(diag(solve(suma3)))  ))
		}

		if(type=="Normal"){
			return(list(beta=teta_novo[1:p],sigma2=teta_novo[p+1],logver=logver,dp=sqrt(diag(solve(suma3)))  ))
		}
	}


}





