function mpsnr_tr= MSI_denoising(ori, kappa, atom, rhos)

if ~exist('rhos', 'var')
    rhos=0.30; % percentage of random noise
end
if ~exist('kappa', 'var')
    kappa = 0.7; % 0.2 0.4 0.7 1
end
if ~exist('atom', 'var')
    atom = 500; % 500 250 100
end

path=strcat('./data/', ori);
T=load(path);
T=T.f1;
sizT=size(T);
n1=sizT(1);
n2=sizT(2);
n3=sizT(3);
Xn=reshape(T,[n1*n2,n3]);
PSNRtpcpsf=zeros(n3,1);
clear f0 path
ind=find(rand(n1*n2*n3,1)<rhos);
Xn(ind)=rand(length(ind),1);
TXn=reshape(Xn,[n1 n2 n3]);
clear ind Xn

TW = zeros(sizT);
for i=1:1:n3
    TW(:,:,i) =  medfilt2(T(:,:,i));
end

lambda=1/sqrt(n3*max(n1,n2));
PC = TW;
D0 = ones(size(PC,1), atom);
D0 = repmat(D0 / sqrt(size(PC, 3)), [1 1 size(PC, 3)]);    
D = ktsvd(PC, D0, 0.1, 10);
D0Y = ones(size(PC, 2),atom);
D0Y = repmat(D0Y / sqrt(size(PC, 3)),[1 1 size(PC, 3)]);    
DY = ktsvd(tran(PC), D0Y, 0.1, 10);
tX=D;
tY=DY;
X=fft(tX,[],3);
Y=fft(tY,[],3);
tXz = zeros(size(X));
tYz = zeros(size(Y));
for i = 1:size(X,3)
    a = orth(X(:,:,i));
    [m1,m2] = size(a);
    tXz(1:m1,1:m2,i) = a;
    b = orth(Y(:,:,i));
    [m3,m4] = size(b);
    tYz(1:m3,1:m4,i) = b;
end
tXnew=ifft(tXz,[],3);
tYnew=ifft(tYz,[],3);
clear m1 m2 m3 m4 a b i tXz tYz tX tY D0 D0Y Xs Side_origion

[L_tpcpsf_tr,]=tpcpsf(TXn,TW,tXnew,tYnew,lambda,kappa);
for i=1:n3
    PSNRtpcpsf(i)=PSNR_RGB(L_tpcpsf_tr(:,:,i),T(:,:,i));
end
mpsnr_tr=mean(PSNRtpcpsf);
end