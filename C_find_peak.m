clear
clc

load('data_heart_a_03.mat')
N=length(P0);

s=1;
for i=2:N-1
    if P0(i)>P0(i-1) && P0(i)>P0(i+1)
        Peak_t(1,s)=i;
        s=s+1;
    end
end

M=s-1;
for i=1:M
    t1=Peak_t(1,i)-1;
    while P0(t1)<P0(t1+1) && t1~=1 
        t1=t1-1;
    end
    if P0(t1)<P0(Peak_t(1,i))/2
            Peak_t(2,i)=1;
    end
end

for i=1:M
    t1=Peak_t(1,i)+1;
    while P0(t1)<P0(t1-1) && t1~=N 
        t1=t1+1;
    end
    if P0(t1)<P0(Peak_t(1,i))/2
        Peak_t(3,i)=1;
    end
end

s=1;
for i=1:M
    if Peak_t(2,i)*Peak_t(3,i)==1
        Peak(s)=Peak_t(1,i);
        s=s+1;
    end
end

semilogy(P0);hold on;semilogy(Peak_t(1,:),ones(1,M),'r*');hold on
semilogy(Peak(:),ones(1,length(Peak)),'bo');hold on
%axis([400 600 0 1000])

for i=1:length(Peak)
    Peak_V(i)=V2R(Peak(i));
end

for i=1:length(Peak)
    t1=Peak(i)-1;
    t2=Peak(i)+1;
    while P0(t1)>P0(Peak(i))/2
        t1=t1-1;
    end
    while P0(t2)>P0(Peak(i))/2
        t2=t2+1;
    end
    Peak_L(i)=V2R(t2)-V2R(t1);
end
save([pwd,'/Peak_V.mat'],'Peak_V')
save([pwd,'/Peak_L.mat'],'Peak_L')