function [YY,consensusIndex]=linear_SDP(D_eq,rho)
%D is a sparse matrix, given the upper bound distance constrainst
%this version will also maximize the local distance
num_points=size(D_eq,2);

 denseRatio=sum(sum(D_eq>0))/(num_points*num_points)
 scale=sqrt(max(D_eq(:)));
% flag_mat=zeros(num_points);

sizeD=size(D_eq);

 K = sdpvar(num_points,num_points);
 Constraints=[sum(sum(K))==0,K>=0];



e=0;
weightVec=[];
if ~isempty(D_eq)
[a,b,c] = find(D_eq);
 flagmat=zeros(num_points);
 odd=[];
 for i=1:length(a)
    flagmat(a(i),b(i))=1; 
    if  flagmat(b(i),a(i))==0
        odd=[odd i];
    end
 end
 
 a=a(odd);
 b=b(odd);
 c=c(odd);
e1=sdpvar(length(odd),1);
e2=sdpvar(length(odd),1);
  Iaa=sub2ind(sizeD,a,a);
  Iba=sub2ind(sizeD,b,a);
  Iab=sub2ind(sizeD,a,b);
  Ibb=sub2ind(sizeD,b,b);
 
  e=sum(K(Iaa)-K(Iab)-K(Iba)+K(Ibb));
  
  Constraints=[Constraints,K(Iaa)-K(Iab)-K(Iba)+K(Ibb)+e1-e2==(c/scale).^2];
  d=(c/scale);
  
  weightVec=(1./d);

  [B IX]=sort(weightVec,'descend');
  cutoff=B(size(D_eq,1));
  weightVec(weightVec>cutoff)=cutoff;

  Constraints=[Constraints,e1>=0,e2>=0];
  
a=0.58;

weightVec2=weightVec;
weightVec1=weightVec;

end


lamda=rho*4*denseRatio;
regT=-lamda*(num_points*trace(K)-e);
DIAGNOSTIC=solvesdp(Constraints,regT+weightVec1'*e1+weightVec2'*e2,sdpsettings('solver','sdpt3','csdp.maxiter',1000,'debug',1,'verbose',1,'dualize',1)); 
regterm=abs(double(regT))
errterm=double(weightVec1'*e1+weightVec2'*e2)
errterm/regterm
Corr_Mat=double(K);
 satisfiedRatio=sum(min(d./sqrt(double(K(Iaa)-K(Iab)-K(Iba)+K(Ibb))),sqrt(double(K(Iaa)-K(Iab)-K(Iba)+K(Ibb)))./d))/length(e1);
sprintf('satisfiedRatio: %f (%d, %d)',satisfiedRatio,length(e1),min(d))
%from Gram matrix to 2D coordinate
    [eig_vec,~]=eig(Corr_Mat);
    eig_val=eig(Corr_Mat); 
    [B,IX]=sort(eig_val,'descend');
    YY=[sqrt(eig_val(IX(1)))*eig_vec(:,IX(1)),  sqrt(eig_val(IX(2)))*eig_vec(:,IX(2)), sqrt(eig_val(IX(3)))*eig_vec(:,IX(3))];
    YY=real(YY)*scale;
  
   consensusIndex=sum( ( (eig_val(IX(1:3)))))/sum( ( (eig_val(eig_val>0))))*satisfiedRatio;
    sprintf('Consensus Ratio: %f',consensusIndex)
end

