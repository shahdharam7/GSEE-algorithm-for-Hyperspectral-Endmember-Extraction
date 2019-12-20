function [cov_ED]=GSEE(Y,n)
% GSEE : 
% Finds the location of the pure endmember from the data 
%
% Reference: 
% Shah D, Zaveri T. A Novel Geo-Stat Endmember Extraction Algorithm. In TENCON 2019-2019 IEEE Region 10 Conference (TENCON) 2019 Oct 17 (pp. 2685-2689). IEEE.
%
% Inputs: Y  -- A 2D array of data. 
%			    The first dimension is of length contains spectral data of length num.
%               The second dimension is number of pixels in the image i.e. P=height*width.    
%         n  -- The number of pure endmembers (hyperspectral subspace dimension)
%
% Output: cov_ED -- Location of n number of pure endmembers.
       
%%% Step-A : Pre-processing %%%
Y1=Y;
[Bands,num] = size(Y);
for i=1:1:Bands
    R=Y1(i,:);
    meanr=mean(R);
    stdr=std(R);
    R1=(R-meanr)/stdr;
    Y(i,:)=R1;
end

%%% Step-B : Band selection using statistical feature %%%
R = (Y*Y')/num;
u = mean(Y,2);
K = R-u*u';
K1 = K.*eye(Bands);
K = K - K1;
imagesc(K);
[q1,q2]=max(K);
[q3,q4]=max(q1);
n_y1=q2(q4);
n_x1=q4;
band_x1=Y(n_x1,:)';
band_y1=Y(n_y1,:)';

%%% Step-C : Convex set geometry %%%
p11=[band_x1,band_y1];
for i=0:0.001:1
    conv_points = boundary(p11,i);
    [m1,n1]=size(conv_points);
    if(m1>n)
        break
    end
end

%%% Step-D : Removal of extra points %%%
cov_ED=conv_points;
z=(size(cov_ED)-1);
for i=1:z(1)
    [t1]=p11(cov_ED(i),:);
    [t2]=p11(cov_ED(i+1),:);
    dist(i,1)=sqrt(sum((t1-t2).^2));
end
e=[p11(cov_ED),p11(cov_ED,2)];
e(end,:)=[];
cov_ED(end,:)=[];
for k=1:1:(m1-n-1)
    [d1,d2]=min(dist);
    e(d2,:)=[];
    dist(d2,:)=[];
    cov_ED(d2,:)=[];
end

end