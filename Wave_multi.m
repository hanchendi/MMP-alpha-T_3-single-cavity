clear
clc

load('Peak_V.mat')
N=length(Peak_V);

for i=1:N
	V_in=Peak_V(i);
     [Psi0,error,X,Y,zX,zY,Pd,Ux,Uy]=MMP(V_in);
     save([pwd,'/Wave/alpha_050_',num2str(i),'.mat'],'Pd','Ux','Uy','zX','zY')

    figure();
    pcolor(zX, zY,  (Pd).^(1/2)); 
    hold on;
    plot(X, Y, 'k--')
    shading flat; axis equal; axis tight
    title(['VR=',num2str(V_in)])
   
    saveas(gca,num2str(i),'bmp');
    close()
    
    delta=(zX(1,2)-zX(1,1))*(zY(2,1)-zY(1,1));
    in = inpolygon(zX,zY,X,Y);
    C(i)=sum(sum(abs(zX).* (Pd).*in))/(sum(sum(abs(zX).*in))*sum(sum((Pd).*in))*delta);
    miuB(i)=abs(sum(sum((Ux.*zY-Uy.*zX).*in))/sum(sum(abs(Pd).*in)));
    
    label(1,i)=Psi0;
    label(2,i)=error;
    disp(i)
end


save([pwd,'/miuB.mat'],'miuB')
save([pwd,'/C.mat'],'C')
save([pwd,'/label.mat'],'label')