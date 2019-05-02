function fval=FitFreq(D,Freq,alpha,lamda)
% a=alpha(1);
% b=alpha(2);

beta1=1/alpha;
 beta0=0;%-log(a)/b;

F_cap=exp(beta0+beta1*log(D));

scale=norm(Freq(D>0&Freq>0))/norm(F_cap(D>0&Freq>0));
F_cap=F_cap*scale;
F_cap(D==0)=1;
Freq(D==0)=1;

fval=sum(sum(abs(F_cap(Freq>0)-Freq(Freq>0))));

%  fval=sum(sum(abs(F_cap(Freq>0)-Freq(Freq>0))./(Freq(Freq>0))));
 
% F_cap_0=exp(beta0+beta1*log(D(Freq==0)));
%  fval=fval+sum(F_cap_0(:))*lamda;
 
%  [~,fval]=matCorr(F_cap,Freq);
%  fval=-fval;
end