function [ ] = three_body()

G=1; %gravitational constant

%masses of bodies in [kg]
M1=1;
M2=1;
M3=1;

%initial coordinates in [m]
x1=1.2;
y1=1.2;
x2=0.5;
y2=0.7;
x3=-0.4;
y3=0.4;

%initial velocities (x and y components)
v1x=0;
v1y=0;
v2x=0;
v2y=0;
v3x=0;
v3y=0;

%matrices for storing trajectory data
r1_matrix=zeros(2,1);
r2_matrix=zeros(2,1);
r3_matrix=zeros(2,1);

%coordinates of the system's center of mass
center_of_mass=zeros(2,1);

dt=0.0001;     %delta t used in iteration [seconds]

for i=1:30000
    %differencies of position vectors (e.g. r1_2 means r1-r2)
    r1_2=sqrt((x1-x2)^2+(y1-y2)^2);
    r1_3=sqrt((x1-x3)^2+(y1-y3)^2);
    r2_3=sqrt((x2-x3)^2+(y2-y3)^2);
    
    dv1x=-G*((M2)/r1_2^3)*(x1-x2)*dt-G*((M3)/r1_3^3)*(x1-x3)*dt;
    dv1y=-G*((M2)/r1_2^3)*(y1-y2)*dt-G*((M3)/r1_3^3)*(y1-y3)*dt;
    
    dv2x=-G*((M1)/r1_2^3)*(x2-x1)*dt-G*((M3)/r2_3^3)*(x2-x3)*dt;
    dv2y=-G*((M1)/r1_2^3)*(y2-y1)*dt-G*((M3)/r2_3^3)*(y2-y3)*dt;
    
    dv3x=-G*((M1)/r1_3^3)*(x3-x1)*dt-G*((M2)/r2_3^3)*(x3-x2)*dt;
    dv3y=-G*((M1)/r1_3^3)*(y3-y1)*dt-G*((M2)/r2_3^3)*(y3-y2)*dt;
    
    v1x=v1x+dv1x;
    v1y=v1y+dv1y;
    
    v2x=v2x+dv2x;
    v2y=v2y+dv2y;
    
    v3x=v3x+dv3x;
    v3y=v3y+dv3y;
    
    x1=x1+v1x*dt;
    y1=y1+v1y*dt;
    r1_matrix(1, i)=x1;
    r1_matrix(2, i)=y1;
    
    x2=x2+v2x*dt;
    y2=y2+v2y*dt;
    r2_matrix(1, i)=x2;
    r2_matrix(2, i)=y2;
    
    x3=x3+v3x*dt;
    y3=y3+v3y*dt;
    r3_matrix(1, i)=x3;
    r3_matrix(2, i)=y3;
    
    center_of_mass(1, i)=(M1*x1+M2*x2+M3*x3)/(M1+M2+M3);
    center_of_mass(2, i)=(M1*y1+M2*y2+M3*y3)/(M1+M2+M3);
    
end

p1=plot(r1_matrix(1,:),r1_matrix(2,:),'-', 'Linewidth', 1.5);
hold on
p2=plot (r2_matrix(1,:),r2_matrix(2,:),'-', 'Linewidth', 1.5);
hold on
p3=plot (r3_matrix(1,:),r3_matrix(2,:),'-', 'Linewidth', 1.5);
hold on
p4=plot (center_of_mass(1,:),center_of_mass(2,:),'o');
set(p4, 'MarkerFaceColor', get(p4,'Color'));
legend([p1 p2 p3 p4],'M1','M2', 'M3', 'Center of mass');
title('Trajectories');

end
