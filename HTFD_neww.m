function [Inew,orient]= HTFD_neww(s,alpha1,alpha2,N)

times=1;
if mod(length(s),2)==0

    I = quadtfd(s, length(s)-1, 1, 'wvd',length(s)*times);
else
    I = quadtfd(s, length(s), 1, 'wvd',length(s)*times);
end

    %I = quadtfd(s, length(s)-1, 1, 'mb',0.25,length(s));
    %length
%I=wvd(s,length(s)-1);
alpha=0.6;
%alpha=0.9;

if length(s)>200
   % alpha=0.25;
end
%alpha=1;
%alpha=0.25;
%I=sign(I).*abs(I).^0.25;

t=-floor(N/2):floor(N/2);
h1=1./((cosh(t).^2).^alpha);
h1=h1/sum(h1);
t=-floor(N/2):floor(N/2);
h2=1./((cosh(t).^2).^alpha);
h2=h2/sum(h2);
B2=h1'*h2;
%mesh(abs(fftshift(fft2(B2))))
I=filter2(B2,I);
M=180;
R=1;
R=3;

[X,Y]=meshgrid(-1:2/N:1,-1:2/N:1);
BB = cell(M/R);
for kk=0:M/R-1
    angle=pi*kk*(R)/M;
    
    X1=X*cos(angle)-Y*sin(angle);
    Y1=X*sin(angle)+Y*cos(angle);
    
    
    A=exp((-1/2)*(((alpha1*X1).^2)+(alpha2*Y1).^2));
   % BB1(kk+1,:,:)=exp((-1/2)*(((alpha1*X1).^2)+(2*alpha2*Y1).^2));
   A=A.*(1-alpha2*alpha2*Y1.^2);
% A=A.*alpha2*Y1;
    A=A/sum(sum(abs(A)));
%     BB(kk+1,:,:)=A;
BB{kk+1} = A;
end
%  A1(:,:)=BB(10,:,:);
%  mesh(A1)

[M1,N1]=size(I);
PQ = paddedsize(size(I));
F = fft2(double(I), PQ(1), PQ(2));

FIabs = fft2(double((abs(I)).^1), PQ(1), PQ(2));
 II=zeros([PQ(1)-length(A) PQ(2)-length(A) M/R]);
% 
%size(II)
 II3=II;
% II=zeros(M1,N1,M/R);
% II3=II;
jjjj=0;
for ii=0:M/R-1
%     B(:,:)=BB(ii+1,:,:);
    if or(ii*R<75,ii*R>105)
    B = BB{ii+1};
    %B1(:,:)=BB1(ii+1,:,:);
    H = fft2(double(B), PQ(1), PQ(2));
    H1 = fft2(double(B), PQ(1), PQ(2));
    
    F_fH = H.*F;
    ffi = real(ifft2(F_fH));
    II3(:,:,jjjj+1)=(ffi(round(length(B)/2):end-round(length(B)/2), round(length(B)/2):end-round(length(B)/2)));
    %F_fH = H.*FIabs;
    
    F_fH = H1.*FIabs;
    
    ffi =( abs(ifft2(F_fH))).^2;
    
    II(:,:,jjjj+1)=(ffi(round(length(B)/2):end-round(length(B)/2), round(length(B)/2):end-round(length(B)/2)));
    end
    jjjj=jjjj+1;
end
% Inew=zeros(size(ffi(1:end/2, 1:end/2)));
% [M1,N1]=size(ffi(1:end/2, 1:end/2));

Inew=I;
lag=0;
 [b,a]=max(II,[],3);


for m=1:M1+lag
    for n=1:N1+lag
        %a=find(max(II(m,n,:))==II(m,n,:));
        
        %aa=find(min(II(m,n,:))==II(m,n,:));
        %orient(m,n)=aa(1);  % Direction of each pixel
        Inew(m,n)=II3(m,n,a(m,n));
        orient(m,n)=(a(m,n)-1)*(R);
    end
end

Inew=filter2(B2,Inew);
%Inew=sign(Inew).*abs(Inew).^4;

Inew(Inew<0)=0;
end