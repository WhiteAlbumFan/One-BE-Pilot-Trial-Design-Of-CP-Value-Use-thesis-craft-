# TODO: Add comment
# 
# Author: Ted under credit of Young Lee
#getwd()
#"D:/R-3.4.4"
file_path_1="D:/JAVA_Palace/LAKILL/result/"
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
	CV_W<-function(s_square_local=0){sqrt(exp(s_square_local)-1);};
	CV_W_Inverse<-function(cv_local=0.5){log(cv_local^2+1);};
};rm(easy_func);

###############################################################################
N1=60;	N11=N1/2;	N12=N11
seq1=rnorm(n=N11, mean=Period_dif+delta_true, sd=sqrt(sigma_square_true));
seq2=rnorm(n=N12, mean=Period_dif-delta_true, sd=sqrt(sigma_square_true));

for(running in c(1)){
	Rec=c();
	Rec1=c();
	for(fl in 4:N11){
		seq_1=seq1[1:fl];	seq_2=seq2[1:fl];
		sigma_s_hat=(sum((seq_1-mean(seq_1))^2)+sum((seq_2-mean(seq_2))^2) )/(length(seq_1)+length(seq_2)-2);
		Rec1=c(Rec1, sigma_s_hat/sigma_nomial -1);
		CP_=( (sigma_s_hat/sigma_nomial) - 1 )*(2*fl-2)/sqrt(2*(N1-2) );
		Rec=c(CP_, Rec);
	};rm(fl);
	print(optimize(Old, interval=c(-5,5), maximum=TRUE) );
	print( (mean(Rec1)+1)*gamma );
}

gamma=1.0;sigma_nomial=gamma*sigma_square_true;

#	0.8 - 1.237297


Old<-function(k){
	waits=1;
	for(fl in 2:length(Rec)){
		waits=waits*dnorm( x=(Rec[fl]-Rec[fl-1])*sqrt(2*N1-2)-k, mean=0, sd=1 );
		#print((Rec[fl]-Rec[fl-1])*sqrt(2*N1-2));
	}
	return(waits);
}








