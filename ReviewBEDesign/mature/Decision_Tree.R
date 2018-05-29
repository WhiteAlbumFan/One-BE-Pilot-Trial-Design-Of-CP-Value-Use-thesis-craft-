# TODO: Add comment
# 
# Author: Ted under credit of Young Lee
########################################################################################################

for(prepare_const in 1:1){
	
	GMR=0.95;	
	delta_true=0;	delta_nomial=delta_true
	#+0.05*rnorm(1);
	
	Period_dif=7e-1;				#	futility-para_1
	alpha=0.10;						#	default CI-para
	theta=0.80;						#	tolerance para
	U=-log(theta);	L=log(theta);	#	tolerance border
	
	sigma_square_true=0.25;			#	major para
	gamma=1.1;						#	rate of est sigma_s to real sigma_s
	sigma_nomial=gamma*sigma_square_true;
	
	beta=0.80;						#	power
	lambda=-0.1;					#	assist para
	pilot_power=0.5;	
	pilot_size=0.3;					#	pilot-magnitude
	pilot_dtr=0.5;
	pivotal_dtr=0.5;
};rm(prepare_const);

for(easy_func in c(1)){
	CV_W<-function(s_square_local=0){sqrt(exp(s_square_local)-1);}
	
	CV_W_Inverse<-function(cv_local=0.5){log(cv_local^2+1);}
	
	
	
};rm(easy_func);

########################################################################################################

power_tree<-function(choice="N", sigma=0.5, mean=1, S1=5, S2=5, lamb=0, alpha_power=alpha){
	lamb=min(lamb, 0);
	if(choice=="N"){
		return( pnorm( (U-mean-lamb)/(sqrt(sigma*(1/S1+1/S2) )) - qnorm(1-alpha_power/2)) -
						pnorm( (L-mean+lamb)/(sigma*sqrt(1/S1+1/S2))+qnorm(1-alpha_power/2)) );}
	else{
		if(choice=="T"){
			return( pnorm( (U-mean-lamb)/(sqrt(sigma*(1/S1+1/S2) )) - qt(p=1-alpha_power/2, df=S1+S2-2) )  -
							pnorm( (L-mean+lamb)/(sigma*sqrt(1/S1+1/S2))+qt(p=1-alpha_power/2, df=S1+S2-2)));}
		else{return(NA);}
	}
}

#Empherical_power_size_loc


########################################################################################################

Fuglsang_Simulation_core<-function(Order=c(1000, "Yes", "Yes", "No", "Yes"),
                                   num=1, Max=300, Para=1, app=1){
  
  N1_dtr=c(num, pilot_dtr); N2_dtr=pivotal_dtr;
	iter_inside=max(as.integer(Order[1]), 5000);
	pilot_permission=Order[2];	pilot_veto=Order[3];
	GMR_obs=Order[4];			CV_obs=Order[5];
	
	true_counter=0;
	Acc<-matrix(rep(0,10), nrow=5, ncol=2,
			dimnames=list(c("Pilot_True","Pilot_False", "Routine", 
			                "design", "rate"),
			              c("Pivotal_Positive","Pivotal_Negative")));
	
	TPNF<-matrix(rep(0,4), nrow=2, ncol=2,
			dimnames=list(c("True","False"), c("Positive","Negative")));
	
	pilot_indicator=0;		pivotal_indicator=1;		
	N=N1_dtr[1]*6;			N1=as.integer(N*N1_dtr[2]);	N2=N-N1;
	sigma_obs=sigma_nomial;	delta_est=log(GMR);			N_consume=0;
	
	Variance_inside<-function(l1, l2){
		return(( var(l1)*(length(l1)-1)+var(l2)*(length(l2)-1))/(length(l1)+length(l2)-2));
	}
	N_old=N;
	
	for(index in 1:iter_inside){
		N=N1_dtr[1]*6;			N1=as.integer(N*N1_dtr[2]);	N2=N-N1;
		N_consume=N_consume+N;	
		pilot_indicator=3;			pivotal_indicator=1;
		seq1=c( rnorm(n=N1, mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
		seq2=c( rnorm(n=N2, mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
		
		delta_hat=mean(seq1)+mean(seq2);	sigma_s_hat=Variance_inside(seq1, seq2)
		
		d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=1-alpha/2, df=N-2);
		d_U=2*delta_hat-d_L;
		
		if( (d_L > L & d_U<U) & (pilot_permission=="Yes") ){
			pilot_indicator=1;}
		else{
			if( (d_L > U | d_U<L) & (pilot_veto=="Yes") ){
			  pilot_indicator=2;
			}
			else{pilot_indicator=3;}
		}
		
		if(CV_obs=="Yes"){sigma_obs=2*sigma_s_hat;}
		if(GMR_obs=="Yes"){delta_est=delta_hat;}
			
		for(prepare_again in c(1)){
			N=12;	N1=as.integer(N*N2_dtr);	N2=N-N1;
			while(N<=Max & power_tree(choice = "T", sigma=sigma_obs, mean=delta_est, S1=N1, S2=N2)<=beta){
				N=N+2;	N1=as.integer(N*N2_dtr);	N2=N-N1;
			}
		};rm(prepare_again);
		
		if(pilot_indicator==3){N_consume=N_consume+N;}
		N1=as.integer(N*N2_dtr);			N2=N-N1;
		seq1=c( rnorm(n=N1, mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
		seq2=c( rnorm(n=N2, mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
		
		delta_hat=mean(seq1)+mean(seq2);	sigma_s_hat=Variance_inside(seq1, seq2)

		d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=1-alpha/2, df=N-2);
		d_U=2*delta_hat-d_L;
		if(  ( (d_L>L & d_U<U) & pilot_indicator==3) | (pilot_indicator==1) ){
			true_counter=true_counter+1;
			pivotal_indicator=1;
		}
		else{pivotal_indicator=2;}
		
		Acc[pilot_indicator, pivotal_indicator]=Acc[pilot_indicator, pivotal_indicator]+1;
	}
	
	
	Fact=1;		if(abs(delta_true)>=U){Fact=2;}
	TPNF[Fact, 1]=true_counter;TPNF[Fact, 2]=iter_inside - true_counter;
	Acc[5,1]=N_consume/iter_inside; Acc[5,2]=app;
	Acc[4,1]=N_old;                 Acc[4,2]=sigma_square_true;
	
	name="D:/Young_Lee_Current/ZXH_PD/MDZZ/Fuglsang_delta";
	name=paste(name, Para, sep="");
	name=paste(name, ".csv", sep="");
	write.table(Acc, file=name, append=TRUE, sep=",", col.name=FALSE, row.names=TRUE);

	print("N_Consume");
	print(N_consume/iter_inside);
	return(1);
}

#Fuglsang_Simulation_core(N1_dtr = c(2, 0.5))

Fuglsang_Simulation_Console<-function(method=1, nums=1, iter=10000, sample_sup=300, Paras=1){
	iter_inside=as.integer(iter);
	sample_sup=max(sample_sup, 100);
	Ord=c(iter_inside, "Yes", "Yes", "No", "Yes");
	
	if(method==1){Ord=c(iter_inside, "No", "No", "No", "Yes");}
	if(method==2){Ord=c(iter_inside, "Yes", "No", "No", "Yes");}
	if(method==3){Ord=c(iter_inside, "No", "No", "Yes", "Yes");}
	if(method==4){Ord=c(iter_inside, "Yes", "No", "Yes", "Yes");}
	if(method==5){Ord=c(iter_inside, "No", "Yes", "No", "Yes");}
	if(method==6){Ord=c(iter_inside, "Yes", "Yes", "No", "Yes");}
	if(method==7){Ord=c(iter_inside, "No", "Yes", "Yes", "Yes");}
	if(method==8){Ord=c(iter_inside, "Yes", "Yes", "Yes", "Yes");}
	
	return(Fuglsang_Simulation_core(Order = Ord, num=nums, Max=sample_sup, Para = Paras, app=method));
}

Fuglsang_Simulation_Console(method=8, nums=1, Paras=2)
#     Paras means that delta's select; and nums mean that sample size  


for(paras in c(10,5,0)){
  delta_true=paras*U/10;
  for(den in 1:5){
    sigma_square_true=den*0.03;
    sigma_nomial=gamma*sigma_square_true;
    for(nnn in c(6,8) ){
      for(nnns in 1:4){
        Fuglsang_Simulation_Console(Paras=paras, method=nnn, nums = nnns, sample_sup = 120);
      }
    }
  }
}
pilot_size=0.5





delta_true=0
delta_true=2*U
sigma_square_true=0.2;sigma_nomial=gamma*sigma_square_true;


########################################################################################################

Potvin_Simulation_Console<-function(method="C", iter=10000, S1_dtr=0.5, S2_dtr=0.5, number=2, Paras=1){
	iter_inside=as.integer(iter);
	if(method=="A"){
		return(Potvin_method_A(iter=iter_inside, N1_dtr = S1_dtr, N2_dtr = S2_dtr, K=number, Para=Paras));}
	if(method=="B"){
		return(Potvin_method_B(iter=iter_inside, N1_dtr = S1_dtr, N2_dtr = S2_dtr, K=number, Para=Paras));}
	if(method=="C"){
		return(Potvin_method_C(iter=iter_inside, N1_dtr = S1_dtr, N2_dtr = S2_dtr, K=number, Para=Paras));}
	if(method=="D"){
		return(Potvin_method_D(iter=iter_inside, N1_dtr = S1_dtr, N2_dtr = S2_dtr, K=number, Para=Paras));}
	else{
		print("method alg error!");
		return(NA);}
}
#   only need to select the value of number and paras

Potvin_method_C<-function(iter=5000, N1_dtr=0.3, N2_dtr=0.3, K=2, Para=1){
	iter_inside=as.integer(iter);	N_consume=0;	true_counter=0;
	pilot_indicator=3;			pivotal_indicator=1;
	
	Acc<-matrix(rep(0,10), nrow=5, ncol=2,
			dimnames=list(c("Pilot_True", "Pilot_False", "Routine",
			                "design", "rate"),
			              c("Pivotal_Positive", "Pivotal_Negative")));
	TPNF<-matrix(rep(0,4), nrow=2, ncol=2,
			dimnames=list(c("True","False"), c("Positive","Negative")));
	
	sigma_obs=sigma_nomial;		delta_est=log(GMR);	alpha_Wang=0.0294;
	
	Variance_inside<-function(l1, l2){
		return( ( var(l1)*(length(l1)-1)+var(l2)*(length(l2)-1))/(length(l1)+length(l2)-2)  );
	};pilot_ind=0;
	for(index in 1:iter_inside){
		N=6*K;					N1=as.integer(N*N1_dtr);	N2=N-N1;
		N_consume=N_consume+N;	pilot_indicator=3;			pivotal_indicator=1;
		
		seq1=c( rnorm(n=N1, mean=(delta_true+Period_dif)/2, sd=sqrt(sigma_square_true/2) ) );
		seq2=c( rnorm(n=N2, mean=(delta_true-Period_dif)/2, sd=sqrt(sigma_square_true/2) ) );
		
		delta_hat=mean(seq1)+mean(seq2);	sigma_s_hat=Variance_inside(seq1, seq2);
		d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=1-alpha/2.0, df=N-2);
		d_U=delta_hat*2 - d_L;
		
		if(power_tree(choice="T", sigma=sigma_s_hat, mean=delta_est, S1=N1, S2=N2)>=beta){
			pilot_ind=pilot_ind+1;
			if(d_L>L & d_U<U){
				true_counter=true_counter+1;	pilot_indicator=1;
			}
			else{pilot_indicator=2;};
		}
		
		else{
			d_L=delta_hat-sqrt(sigma_s_hat*(1/length(seq1)+1/length(seq2)) )*qt(p=1-alpha_Wang/2, df=N+N_old-2);
			d_U=delta_hat*2-d_L;
			if(d_L>L & d_U<U){
				true_counter=true_counter+1;
				pilot_indicator=1;
			}
			else{
				pilot_indicator=3;
			}
		}
		
		for(prepare_again in c(1)){
			sigma_obs=sigma_s_hat;		N_old=N;	N=12;
			N1=as.integer(N*N2_dtr);	N2=N-N1;
			while(N<=100 & power_tree(choice = "T", sigma=sigma_obs, mean=delta_est, S1=N1, S2=N2, alpha_power = alpha_Wang)<=beta){
				N=N+2;	N1=as.integer(N*N2_dtr);	N2=N-N1;
			};
		};rm(prepare_again);
		if(pilot_indicator==3){N_consume=N_consume+N;}
		#N=as.integer(2*sigma_s_hat*(qnorm(0.5+beta/2)+qnorm(1-alpha_local/2))^2*1.0/U^2);
		N1=as.integer(N*N2_dtr);	N2=N-N1;
		
		seq1=c(seq1, rnorm(n=N1,mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
		seq2=c(seq2, rnorm(n=N2,mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
		
		delta_hat=mean(seq1)+mean(seq2);	sigma_s_hat=Variance_inside(seq1, seq2);
		
		d_L=delta_hat-sqrt(sigma_s_hat*(1/length(seq1)+1/length(seq2)) )*qt(p=1-alpha_Wang/2, df=N+N_old-2);
		d_U=delta_hat*2-d_L;
		if(d_L>L*1.0 & d_U<U ){
			pivotal_indicator=1;
			if(pilot_indicator==3){true_counter=true_counter+1;}
		}
		else{pivotal_indicator=2;}
		Acc[pilot_indicator, pivotal_indicator]=Acc[pilot_indicator, pivotal_indicator]+1;
	}
	Acc[4,1]=K*6;Acc[4,2]=sigma_square_true;
	Acc[5,1]=N_consume/iter_inside;Acc[5,2]=3
	
	Fact=1;		if(abs(delta_true)>=U){Fact=2;}
	TPNF[Fact, 1]=true_counter;	TPNF[Fact, 2]=iter_inside - true_counter;
	
	name="D:/Young_Lee_Current/ZXH_PD/MDZZ_Re/Potvin_delta";
	name=paste(name, Para, sep="");
	name=paste(name, ".csv", sep="");
	write.table(Acc, file=name, append=TRUE, sep=",", col.name=FALSE, row.names=TRUE);
	
	print("N_Consume");
	#print(N_consume/iter_inside);
	return(N_consume/iter_inside);
}



for(paras in c(10,5,0)){
  delta_true=paras*U/10;
  for(den in 1:10){
    sigma_square_true=den*0.02;
    sigma_nomial=gamma*sigma_square_true;
    for(nnn in (1:6) ){
      Potvin_Simulation_Console(Paras=paras, number=nnn);
    }
  }
}




















Potvin_method_A<-function(iter=5000, N1_dtr=0.5, N2_dtr=0.5, K=2, Para=2){
  iter_inside=as.integer(iter);	N_consume=0;		true_counter=0;
  pilot_indicator=3;				pivotal_indicator=1;
  
  Acc<-matrix(rep(0,8), nrow=4, ncol=2,
              dimnames=list(c("Pilot_True", "Pilot_False", "Routine"), c("Pivotal_Positive", "Pivotal_Negative")));
  TPNF<-matrix(rep(0,4), nrow=2, ncol=2,
               dimnames=list(c("True","False"), c("Positive","Negative")));
  
  sigma_obs=sigma_nomial;			delta_est=log(GMR);	alpha_Wang=0.05;
  
  Variance_inside<-function(l1, l2){
    return(( var(l1)*(length(l1)-1)+var(l2)*(length(l2)-1))/(length(l1)+length(l2)-2));
  };
  pilot_ind=0;
  for(index in 1:iter_inside){
    N=12*K;					N1=as.integer(N*N1_dtr);	N2=N-N1;
    N_consume=N_consume+N;	pilot_indicator=3;			pivotal_indicator=1;
    
    seq1=c( rnorm(n=N1, mean=(delta_true+Period_dif)/2, sd=sqrt(sigma_square_true/2) ) );
    seq2=c( rnorm(n=N2, mean=(delta_true-Period_dif)/2, sd=sqrt(sigma_square_true/2) ) );
    
    delta_hat=mean(seq1)+mean(seq2);	sigma_s_hat=Variance_inside(seq1, seq2);
    
    if(power_tree(choice="T", sigma=sigma_s_hat, mean=delta_est, S1=N1, S2=N2)>=beta){
      d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=1-alpha_Wang/2.0, df=N-2);
      d_U=delta_hat*2-d_L;
      pilot_ind=pilot_ind+1;
      if(d_L>L & d_U<U){
        true_counter=true_counter+1;
        pilot_indicator=1;
      }
      else{pilot_indicator=2;};
    }
    
    for(prepare_again in c(1)){
      sigma_obs=sigma_s_hat;		N_old=N;	N=12;
      N1=as.integer(N*N2_dtr);	N2=N-N1;
      while(N<=100 & power_tree(choice = "T", sigma=sigma_obs, mean=delta_est, S1=N1, S2=N2)<=beta){
        N=N+2;	N1=as.integer(N*N2_dtr);	N2=N-N1;
      };
    };rm(prepare_again);
    if(pilot_indicator==3){N_consume=N_consume+N;}
    #N=as.integer(2*sigma_s_hat*(qnorm(0.5+beta/2)+qnorm(1-alpha_local/2))^2*1.0/U^2);
    N1=as.integer(N*N2_dtr);	N2=N-N1;
    
    seq1=c(seq1, rnorm(n=N1,mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
    seq2=c(seq2, rnorm(n=N2,mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
    
    delta_hat=mean(seq1)+mean(seq2);	sigma_s_hat=Variance_inside(seq1, seq2);
    
    d_L=delta_hat-sqrt(sigma_s_hat*(1/length(seq1)+1/length(seq2)) )*qt(p=1-alpha_Wang/2, df=N+N_old-2);
    d_U=delta_hat*2-d_L;
    if(d_L>L*1.0 & d_U<U ){
      pivotal_indicator=1;
      if(pilot_indicator==3){true_counter=true_counter+1;}
    }
    else{pivotal_indicator=2;}
    Acc[pilot_indicator, pivotal_indicator]=Acc[pilot_indicator, pivotal_indicator]+1;
  }
  
  Fact=1;		if(abs(delta_true)>=U){Fact=2;}
  TPNF[Fact, 1]=true_counter;	TPNF[Fact, 2]=iter_inside - true_counter;
  Acc[4,1]=N, Acc[4,2]=sigma_square_true
  
  name="D:/Young_Lee_Current/ZXH_PD/MDZZ_Re/TF_delta";
  name=paste(name, Para, sep="");
  name=paste(name, ".csv", sep="");
  write.table(Acc, file=name, append=TRUE, sep=",", col.name=FALSE, row.names=TRUE);
  
  print(N_consume/iter_inside);
  print("pilot_ind");print(pilot_ind);
  return(TPNF);
}

Potvin_method_B<-function(iter=5000, N1_dtr=0.3, N2_dtr=0.3, K=2){
  iter_inside=as.integer(iter);	N_consume=0;	true_counter=0;
  pilot_indicator=3;			pivotal_indicator=1;
  
  Acc<-matrix(rep(0,6), nrow=3, ncol=2,
              dimnames=list(c("Pilot_True", "Pilot_False", "Routine"), c("Pivotal_Positive", "Pivotal_Negative")));
  TPNF<-matrix(rep(0,4), nrow=2, ncol=2,
               dimnames=list(c("True","False"), c("Positive","Negative")));
  
  sigma_obs=sigma_nomial;		delta_est=log(GMR);	alpha_Wang=0.0294;
  
  Variance_inside<-function(l1, l2){
    return( ( var(l1)*(length(l1)-1)+var(l2)*(length(l2)-1))/(length(l1)+length(l2)-2)  );
  };pilot_ind=0;
  for(index in 1:iter_inside){
    N=12*K;					N1=as.integer(N*N1_dtr);	N2=N-N1;
    N_consume=N_consume+N;	pilot_indicator=3;			pivotal_indicator=1;
    
    seq1=c( rnorm(n=N1, mean=(delta_true+Period_dif)/2, sd=sqrt(sigma_square_true/2) ) );
    seq2=c( rnorm(n=N2, mean=(delta_true-Period_dif)/2, sd=sqrt(sigma_square_true/2) ) );
    
    delta_hat=mean(seq1)+mean(seq2);	sigma_s_hat=Variance_inside(seq1, seq2);
    
    d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=1-alpha_Wang/2.0, df=N-2);
    d_U=delta_hat*2-d_L;
    
    if(d_L>L & d_U<U){
      pilot_ind=pilot_ind+1;	pilot_indicator=1;
      true_counter=true_counter+1;
    }
    else{
      if(power_tree(choice="T", sigma=sigma_s_hat, mean=delta_est, S1=N1, S2=N2, alpha_power = alpha_Wang)>=beta){
        pilot_indicator=2;
      }
      else{pilot_indicator=3;};
    }
    
    for(prepare_again in c(1)){
      sigma_obs=sigma_s_hat;		N_old=N;	N=12;
      N1=as.integer(N*N2_dtr);	N2=N-N1;
      while(N<=100 & power_tree(choice = "T", sigma=sigma_obs, mean=delta_est, S1=N1, S2=N2, alpha_power = alpha_Wang)<=beta){
        N=N+2;	N1=as.integer(N*N2_dtr);	N2=N-N1;
      };
    };rm(prepare_again);
    if(pilot_indicator==3){N_consume=N_consume+N;}
    #N=as.integer(2*sigma_s_hat*(qnorm(0.5+beta/2)+qnorm(1-alpha_local/2))^2*1.0/U^2);
    N1=as.integer(N*N2_dtr);	N2=N-N1;
    
    seq1=c(seq1, rnorm(n=N1,mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
    seq2=c(seq2, rnorm(n=N2,mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
    
    delta_hat=mean(seq1)+mean(seq2);	sigma_s_hat=Variance_inside(seq1, seq2);
    
    d_L=delta_hat-sqrt(sigma_s_hat*(1/length(seq1)+1/length(seq2)) )*qt(p=1-alpha_Wang/2, df=N+N_old-2);
    d_U=delta_hat*2-d_L;
    if(d_L>L*1.0 & d_U<U ){
      pivotal_indicator=1;
      if(pilot_indicator==3){true_counter=true_counter+1;}
    }
    else{pivotal_indicator=2;}
    Acc[pilot_indicator, pivotal_indicator]=Acc[pilot_indicator, pivotal_indicator]+1;
  }
  
  Fact=1;		if(abs(delta_true)>=U){Fact=2;}
  TPNF[Fact, 1]=true_counter;	TPNF[Fact, 2]=iter_inside - true_counter;
  
  print(Acc);	print("N_Consume");
  print(N_consume/iter_inside);
  print("pilot_ind");print(pilot_ind);
  return(TPNF);
}

Potvin_method_D<-function(iter=5000, N1_dtr=0.3, N2_dtr=0.3, K=2){
	iter_inside=as.integer(iter);	N_consume=0;	true_counter=0;
	pilot_indicator=3;			pivotal_indicator=1;
	
	Acc<-matrix(rep(0,6), nrow=3, ncol=2,
			dimnames=list(c("Pilot_True", "Pilot_False", "Routine"), c("Pivotal_Positive", "Pivotal_Negative")));
	TPNF<-matrix(rep(0,4), nrow=2, ncol=2,
			dimnames=list(c("True","False"), c("Positive","Negative")));
	
	sigma_obs=sigma_nomial;		delta_est=log(GMR);	alpha_Wang=0.028;
	
	Variance_inside<-function(l1, l2){
		return( ( var(l1)*(length(l1)-1)+var(l2)*(length(l2)-1))/(length(l1)+length(l2)-2)  );
	};pilot_ind=0;
	for(index in 1:iter_inside){
		N=12*K;					N1=as.integer(N*N1_dtr);	N2=N-N1;
		N_consume=N_consume+N;	pilot_indicator=3;			pivotal_indicator=1;
		
		seq1=c( rnorm(n=N1, mean=(delta_true+Period_dif)/2, sd=sqrt(sigma_square_true/2) ) );
		seq2=c( rnorm(n=N2, mean=(delta_true-Period_dif)/2, sd=sqrt(sigma_square_true/2) ) );
		
		delta_hat=mean(seq1)+mean(seq2);	sigma_s_hat=Variance_inside(seq1, seq2);
		d_L=delta_hat-sqrt(sigma_s_hat*(1.0/N1+1.0/N2))*qt(p=1-alpha/2.0, df=N-2);
		d_U=delta_hat*2 - d_L;
		
		if(power_tree(choice="T", sigma=sigma_s_hat, mean=delta_est, S1=N1, S2=N2)>=beta){
			pilot_ind=pilot_ind+1;
			if(d_L>L & d_U<U){
				true_counter=true_counter+1;	pilot_indicator=1;
			}
			else{pilot_indicator=2;};
		}
		
		else{
			d_L=delta_hat-sqrt(sigma_s_hat*(1/length(seq1)+1/length(seq2)) )*qt(p=1-alpha_Wang/2, df=N+N_old-2);
			d_U=delta_hat*2-d_L;
			if(d_L>L & d_U<U){
				true_counter=true_counter+1;
				pilot_indicator=1;}
			else{pilot_indicator=3;};
		}
		
		for(prepare_again in c(1)){
			sigma_obs=sigma_s_hat;		N_old=N;	N=12;
			N1=as.integer(N*N2_dtr);	N2=N-N1;
			while(N<=100 & power_tree(choice = "T", sigma=sigma_obs, mean=delta_est, S1=N1, S2=N2, alpha_power = alpha_Wang)<=beta){
				N=N+2;	N1=as.integer(N*N2_dtr);	N2=N-N1;
			};
		};rm(prepare_again);
		if(pilot_indicator==3){N_consume=N_consume+N;}
		#N=as.integer(2*sigma_s_hat*(qnorm(0.5+beta/2)+qnorm(1-alpha_local/2))^2*1.0/U^2);
		N1=as.integer(N*N2_dtr);	N2=N-N1;
		
		seq1=c(seq1, rnorm(n=N1,mean=(delta_true+Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
		seq2=c(seq2, rnorm(n=N2,mean=(delta_true-Period_dif)/2.0, sd=sqrt(sigma_square_true/2.0)));
		
		delta_hat=mean(seq1)+mean(seq2);	sigma_s_hat=Variance_inside(seq1, seq2);
		
		d_L=delta_hat-sqrt(sigma_s_hat*(1/length(seq1)+1/length(seq2)) )*qt(p=1-alpha_Wang/2, df=N+N_old-2);
		d_U=delta_hat*2-d_L;
		if(d_L>L*1.0 & d_U<U ){
			pivotal_indicator=1;
			if(pilot_indicator==3){true_counter=true_counter+1;}
		}
		else{pivotal_indicator=2;}
		Acc[pilot_indicator, pivotal_indicator]=Acc[pilot_indicator, pivotal_indicator]+1;
	}
	
	Fact=1;		if(abs(delta_true)>=U){Fact=2;}
	TPNF[Fact, 1]=true_counter;	TPNF[Fact, 2]=iter_inside - true_counter;
	
	print(Acc);	print("N_Consume");
	print(N_consume/iter_inside);
	print("pilot_ind");print(pilot_ind);
	return(TPNF);
}


delta_true=log(GMR)
delta_true=U
Potvin_Simulation_Console(method="D", iter=5000, number=4)









death_list=as.character(ls() );death_list_1=death_list[1]
for(name in death_list[2:length(death_list)]){
	death_list_1<-paste(death_list_1, name, sep=",")
}
death_list_1;rm(death_list, name, death_list_1)




