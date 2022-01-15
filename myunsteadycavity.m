%Solve the two-dimensional unsteady Navier-Stokes equations for an elliptic
%flow in a square cavity with top wall sliding at a given velocity
clc;
clear;
Lx=1;%The width of the domain of interest is 1m.
Ly=1;%The height of the domain of interest is 1m.
Re=50;%Reynolds number
Nx=50;%Number of nodes in the horizontal direction
Ny=50;%Number of nodes in the vertical direction
DeltaX=Lx/(Nx-2);%Spatial increment in the horizontal 
DeltaY=Ly/(Ny-2);%Spatial increment in the vertical direction
DeltaT=0.002;
%Time increment that yields to the stability limit
TotalT=100;%Total simulation time
u(:,:,1)=zeros(Ny,Nx);
v(:,:,1)=zeros(Ny,Nx); 

%Initialization of the velocity field
p(:,:,1)=zeros(Ny,Nx);%Initialization of the pressure field
u_hat=zeros(Ny,Nx);
v_hat=zeros(Ny,Nx);
u_hat_mid=zeros(Ny,Nx-1);
v_hat_mid=zeros(Ny-1,Nx);
CE=zeros(Ny,Nx);%Matrix to store continuity errors.
for n=1:round(TotalT/DeltaT)
    for i=2:(Nx-1)
    u(Ny,i,1)=2*sin(n*DeltaT*0.1)-u(Ny-1,i,1);%Initialization of ghost nodes value
    end
    %Solve for u_hat and v_hat at each node (i,j) inside the domain.
    for i=2:(Nx-1)
        for j=2:(Ny-1)
            C_x=u(j,i,n)/(2*DeltaX)*(u(j,i+1,n)-u(j,i-1,n))+v(j,i,n)/(2*DeltaY)*(u(j+1,i,n)-u(j-1,i,n));%Horizontal convective term
            C_y=u(j,i,n)/(2*DeltaX)*(v(j,i+1,n)-v(j,i-1,n))+v(j,i,n)/(2*DeltaY)*(v(j+1,i,n)-v(j-1,i,n));%Vertical convective term
            D_x=1/(Re*DeltaX^2)*(u(j,i+1,n)-2*u(j,i,n)+u(j,i-1,n))+1/(Re*DeltaY^2)*(u(j+1,i,n)-2*u(j,i,n)+u(j-1,i,n));%Horizontal diffusive term
            D_y=1/(Re*DeltaX^2)*(v(j,i+1,n)-2*v(j,i,n)+v(j,i-1,n))+1/(Re*DeltaY^2)*(v(j+1,i,n)-2*v(j,i,n)+v(j-1,i,n));%Vetical diffusive term
            u_hat(j,i)=u(j,i,n)+DeltaT*(-1*C_x+D_x);
            v_hat(j,i)=v(j,i,n)+DeltaT*(-1*C_y+D_y);
        end
    end
    
    %Evaluate u_hat and v_hat at midpoints.
    for i=2:(Nx-2)
        for j=1:(Ny)
            u_hat_mid(j,i)=(u_hat(j,i+1)+u_hat(j,i))/2;
        end
    end
    for i=1:(Nx)
        for j=2:(Ny-2)
            v_hat_mid(j,i)=(v_hat(j+1,i)+v_hat(j,i))/2;
        end
    end
    for i=2:(Nx-1)
        for j=2:(Ny-1)
            CE(j,i)=(1/DeltaT)*((u_hat_mid(j,i)-u_hat_mid(j,i-1))/DeltaX+(v_hat_mid(j,i)-v_hat_mid(j-1,i))/DeltaY);
            %Continuty error
        end
    end
    
    %Using Successive Over-relaxation Scheme to solve the Discrete
    %Pressure-Poisson Equations.
    omega=1.4;
    p(:,:,n+1)=p(:,:,n);
    
    for iii=1:200
        for i=2:(Nx-1)
            p(1,i,n+1)=p(2,i,n+1);
            p(Ny,i,n+1)=p(Ny-1,i,n+1);
        end
        for j=2:(Ny-1)
            p(j,1,n+1)=p(j,2,n+1);
            p(j,Nx,n+1)=p(j,Nx-1,n+1);
        end
        for i=2:(Nx-1)
            for j=2:(Ny-1)
                R=((p(j,i+1,n+1)+p(j,i-1,n+1))/(DeltaX^2)+(p(j+1,i,n+1)+p(j-1,i,n+1))/(DeltaY^2)-CE(j,i))/(2/(DeltaX^2)+2/(DeltaY^2));%R=P(j,i)
                p(j,i,n+1)=omega*R+(1-omega)*p(j,i,n+1);
            end
        end
        
    end
    
    for i=2:(Nx-1)
        p(1,i,n+1)=p(2,i,n+1);
        p(Ny,i,n+1)=p(Ny-1,i,n+1);
    end
    for j=2:(Ny-1)
        p(j,1,n+1)=p(j,2,n+1);
        p(j,Nx,n+1)=p(j,Nx-1,n+1);
    end
    
    %Correction of grid node velocities
    for i=2:(Nx-1)
        for j=2:(Ny-1)
            u(j,i,n+1)=u_hat(j,i)-DeltaT*(p(j,i+1,n+1)-p(j,i-1,n+1))/(2*DeltaX);
            v(j,i,n+1)=v_hat(j,i)-DeltaT*(p(j+1,i,n+1)-p(j-1,i,n+1))/(2*DeltaY);
        end
    end
    
    %Update ghost node values.
    for i=2:(Nx-1)
        u(Ny,i,n+1)=2*sin(n*DeltaT*0.1)-u(Ny-1,i,n+1);
        u(1,i,n+1)=-u(2,i,n+1);
        v(Ny,i,n+1)=-v(Ny-1,i,n+1);
        v(1,i,n+1)=-v(2,i,n+1);
    end
    for j=2:(Ny-1)
        u(j,1,n+1)=-u(j,2,n+1);
        u(j,Nx,n+1)=-u(j,Nx-1,n+1);
        v(j,1,n+1)=-v(j,2,n+1);
        v(j,Nx,n+1)=-v(j,Nx-1,n+1);
    end
    
    if  n==15000 %solution at time 30 s
        break;
    end
 
end

%Streamline plot
startx=[];
starty1=[];
starty2=[];
starty3=[];
FigHandle = figure('Position', [100, 40, 800, 600]);
Jump=1;%Coefficient of pathline skipping
startx=DeltaX/2:Jump*DeltaX:Lx-DeltaX/2;
Jump=3;
for i=DeltaY/2:Jump*DeltaY:Ly-DeltaY/2
    starty1(round(1-0.5/Jump+i/DeltaY/Jump),:)=i*ones(size(startx));
end
for i=DeltaY/2+Jump/3*DeltaY:Jump*DeltaY:Ly-DeltaY/2
    starty2(round(2/3-0.5/Jump+i/DeltaY/Jump),:)=i*ones(size(startx));
end
for i=DeltaY/2+Jump*2/3*DeltaY:Jump*DeltaY:Ly-DeltaY/2
    starty3(round(1/3-0.5/Jump+i/DeltaY/Jump),:)=i*ones(size(startx));
end

for i=DeltaY/2:Jump*DeltaY:Ly-DeltaY/2
    h=streamline([DeltaX/2:DeltaX:Lx-DeltaX/2],[DeltaY/2:DeltaY:Ly-DeltaY/2],u(2:Ny-1,2:Nx-1,n),v(2:Ny-1,2:Nx-1,n),startx,starty1(round(1-0.5/Jump+i/DeltaY/Jump),:),[0.1,500]);set(h,'Color','red');
end
for i=DeltaY/2+Jump/3*DeltaY:Jump*DeltaY:Ly-DeltaY/2
    h=streamline([DeltaX/2:DeltaX:Lx-DeltaX/2],[DeltaY/2:DeltaY:Ly-DeltaY/2],u(2:Ny-1,2:Nx-1,n),v(2:Ny-1,2:Nx-1,n),startx,starty2(round(2/3-0.5/Jump+i/DeltaY/Jump),:),[0.1,500]);set(h,'Color','blue');
end
for i=DeltaY/2+Jump*2/3*DeltaY:Jump*DeltaY:Ly-DeltaY/2
    h=streamline([DeltaX/2:DeltaX:Lx-DeltaX/2],[DeltaY/2:DeltaY:Ly-DeltaY/2],u(2:Ny-1,2:Nx-1,n),v(2:Ny-1,2:Nx-1,n),startx,starty3(round(1/3-0.5/Jump+i/DeltaY/Jump),:),[0.1,500]);set(h,'Color','green');
end

axis([0 Lx 0 Ly]);
axis square;
xlabel('X');
ylabel('Y');

%Contour of u-velocity 
FigHandle = figure('Position', [100, 40, 800, 600]);
contourf([DeltaX/2:DeltaX:Lx-DeltaX/2],[DeltaY/2:DeltaY:Ly-DeltaY/2],u(2:Ny-1,2:Nx-1,n));
axis([0 Lx 0 Ly]);
axis square;
h=colorbar;
title(h,'U');
xlabel('X');
ylabel('Y');

%Contour of v-velocity 
FigHandle = figure('Position', [100, 40, 800, 600]);
contourf([DeltaX/2:DeltaX:Lx-DeltaX/2],[DeltaY/2:DeltaY:Ly-DeltaY/2],v(2:Ny-1,2:Nx-1,n));
axis([0 Lx 0 Ly]);
axis square;
h=colorbar;
title(h,'V');
xlabel('X');
ylabel('Y');

%u-velocity profile on vertical center line

FigHandle = figure('Position', [100, 40, 800, 600]);
plot(DeltaX/2:DeltaX:Lx-DeltaX/2,u(round((1+Ny)/2),2:Nx-1,n+1));
xlabel('X');
ylabel('velocity u at vertical centerline');

%v-velocity profile on horizontal center line
FigHandle = figure('Position', [100, 40, 800, 600]);
plot(DeltaX/2:DeltaX:Lx-DeltaX/2,v(round((1+Ny)/2),2:Nx-1,n+1));
xlabel('X');
ylabel('velocity v at vertical centerline');

