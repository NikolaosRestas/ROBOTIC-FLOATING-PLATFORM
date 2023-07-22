%%%% CLEAR ALL PROGRAM %%%%%%
clear;
clc;
%%% ALL PROGRAM DATA AND MATRIXES FROM PAPER FOR THE CALCULATIONS %%%
m=425000.000;
Izz = (25.73*(10^7));
lab=45.00;
lac=45.00;
lbc=35.00;
dae=27.64;
dge=0.00;
ruc=2.20;
rlc=3.50;
huc=6.50;
hlc=3.00;
ca=0.80;
cd=0.80;
pw=1024.00;
dad=(3/2)*dae;
L=dad;
dbf=lbc/2;
p=1.22;
h = huc-(1.0/(ruc*ruc))*(m/(3*pi*pw)-rlc*rlc*hlc);
ma = -ca*pi*p*(ruc*ruc*huc+ rlc*rlc*hlc); 
k = ((-2*dad)+(4*dad*dae)-(3*dae*dae)-((5/4)*(lbc*lbc))+(3*lbc*dbf)-(3*dbf*dbf));
m33 = (Izz+(ma*k));
D=cd*pw*(ruc*(huc-h)+rlc*hlc);
%%%%%%%% MATRIXES %%%%%%%%
J = [1 0 1 0 1 0;
    0 -1 0 -1 0 -1;
    -(dbf-(lbc/2)) -dae -dbf (dad-dae) (lbc-dbf) (dad-dae)];
    
Sag=[dae; (dbf-lbc/2); 0];
Sbg=[-(dad-dae); dbf; 0];
Scg=[-(dad-dae); dbf-lbc; 0];

M = [(m-(3*ma)) 0 (3*(dbf-(lbc/2))*ma);
0 (m-(3*ma)) (((2*dad)-(3*dae))*ma);
(3*(dbf-(lbc/2))*ma) (((2*dad)-(3*dae))*ma) m33];

x(1)=0;y(1)=100;theta(1)=0;
u(1)=0;v(1)=0;r(1)=0;
gonia_anemou(1)=30;
gonia(1)=0;
v_anemou=10;


dt=0.5; 
t=3600; %1000 %3600

xdes=500; 
ydes=0;

kx=300; 
ky=300;      
ktheta=20000;

index=1;
sum_avg_jet_force=0;

%%%%%%%OPTIMIZATION%%%%%%
lb = -200000*ones(6,1);
ub = 200000*ones(6,1);

%objective_function = @(fx) abs((atan2(fx(1,1),fx(2,1))) + (atan2(fx(3,1),fx(4,1))) + (atan2(fx(5,1),fx(6,1))));
%objective_function = @(fx) ((abs(atan2(fx(1,1),fx(2,1)))) + (abs(atan2(fx(3,1),fx(4,1)))) + (abs(atan2(fx(5,1),fx(6,1)))));
%objective_function = @(fx) norm(fx,2);
%objective_function = @(fx) sqrt(fx(1,1)*fx(1,1)+fx(2,1)*fx(2,1))+sqrt(fx(3,1)*fx(3,1)+fx(4,1)*fx(4,1))+sqrt(fx(5,1)*fx(5,1)+fx(6,1)*fx(6,1));

f0 = [0; 0; 0; 0; 0; 0];

for i=1:dt:t
  time(index)=i;
  
  gonia(index)=(gonia_anemou(index)-theta(index)); 
  theta_des(index)=gonia_anemou(index); 
  
  error_x(index)=(xdes-x(index));
  error_y(index)=(ydes-y(index));
  error_theta(index)=(theta_des(index)-theta(index));
  
  b_qc=[kx*error_x(index);
        ky*error_y(index);
        ktheta*error_theta(index)];

  %g = @(fy) J*fy-b_qc;
  %[bf, obj, info, iter, nf, lambda] = sqp (f0, objective_function, g, [], lb, ub,1000);
  %f = bf;
  
  f = ((pinv(J)*(b_qc)));
  
  fa(index)=sqrt(((f(1,1)*f(1,1))+(f(2,1)*f(2,1)))); %THE DISIRED JET THRUSTS
  fb(index)=sqrt(((f(3,1)*f(3,1))+(f(4,1)*f(4,1))));
  fc(index)=sqrt(((f(5,1)*f(5,1))+(f(6,1)*f(6,1))));
  
  theta_a(index)=atan2(f(1,1),f(2,1));
  theta_b(index)=atan2(f(3,1),f(4,1));
  theta_c(index)=atan2(f(5,1),f(6,1));
  
  jet_forces_avg(index)=((fa(index)+fb(index)+fc(index))/3); %CALCULATE AVERAGE JET THRUST
  sum_avg_jet_force=sum_avg_jet_force+jet_forces_avg(index);
  

  %%%%%%%%%%%% ANTISTASH %%%%%%%%%%%%  
  VA=[(u(index)+r(index)*(lbc/2-dbf));
      (v(index)+r(index)*dae);
      0];
      
  VB=[(u(index)-r(index)*dbf);
      (v(index)-r(index)*(dad-dae));
      0];
      
  VC=[(u(index)+r(index)*(lbc-dbf));
      (v(index)-r(index)*(dad-dae));
      0];
   
  fha=D*[abs(VA(1,1))*(-VA(1,1));
         abs(VA(2,1))*(-VA(2,1));
         abs(VA(3,1))*(-VA(3,1))];
       
  fhb=D*[abs(VB(1,1))*(-VB(1,1));
         abs(VB(2,1))*(-VB(2,1));
         abs(VB(3,1))*(-VB(3,1))]; 
        
  fhc=D*[abs(VC(1,1))*(-VC(1,1));
         abs(VC(2,1))*(-VC(2,1));
         abs(VC(3,1))*(-VC(3,1))];
          
  k1=cross(Sag,fha);
  k2=cross(Sbg,fhb);
  k3=cross(Scg,fhc);
  
  qha=[fha(1,1);
      fha(2,1);
      k1(3,1)];
      
  qhb=[fhb(1,1);
      fhb(2,1);
      k2(3,1)];
      
  qhc=[fhc(1,1);
      fhc(2,1);
      k3(3,1)];
      
  sum=(qha+qhb+qhc);
  %%%%%%% WAVES %%%%%%%%%
  f_water1(index)=25000*cosd(3*index*dt); %MAKING THE WAVE FORCES 
  fa_water_x(index)=f_water1(index)*sind(gonia_anemou(index));
  fa_water_y(index)=f_water1(index)*cosd(gonia_anemou(index));
  fb_water_x(index)=f_water1(index)*sind(gonia_anemou(index));
  fb_water_y(index)=f_water1(index)*cosd(gonia_anemou(index));
  fc_water_x(index)=f_water1(index)*sind(gonia_anemou(index));
  fc_water_y(index)=f_water1(index)*cosd(gonia_anemou(index));
  
  f_water=[fa_water_x(index);
           fa_water_y(index);
           fb_water_x(index);
           fb_water_y(index);
           fc_water_x(index);
           fc_water_y(index)];
  
  fa_water(index)=sqrt(((f_water(1,1)*f_water(1,1))+(f_water(2,1)*f_water(2,1))));%WAVE CYLINDER FORCES
  fb_water(index)=sqrt(((f_water(3,1)*f_water(3,1))+(f_water(4,1)*f_water(4,1))));
  fc_water(index)=sqrt(((f_water(5,1)*f_water(5,1))+(f_water(6,1)*f_water(6,1))));
  
  b_water=(J*f_water);
  %%%%%%%%%% WIND %%%%%%%%%%
   R=[cosd(theta(index)) -sind(theta(index)) 0;
     sind(theta(index)) cosd(theta(index)) 0;
     0 0 1]; 
  
  cx(index)= ((-0.9)*cosd(gonia(index)));
  cy(index)= ((-0.85)*sind(gonia(index)));
  cn(index)= ((-0.125)*sind(2*gonia(index)));
  
  at = 35.00*2.50;
  al = 45.00*2.50;
  
  x_wind(index)=(0.5*cx(index)*p*(v_anemou*v_anemou)*at);
  y_wind(index)=(0.5*cy(index)*p*(v_anemou*v_anemou)*al);
  n_wind(index)=(0.5*cn(index)*p*(v_anemou*v_anemou)*al*L); 
  WIND=[x_wind(index);y_wind(index);n_wind(index)];
  
  b_wind=((inv(R)*WIND));
  
  %%%%%%%%%% DYNAMIC %%%%%%%%%%%
  q_dist = b_wind + b_water;
  
  q_c    = J*f;
  
  vdot= inv(M)*(sum + inv(R)*q_c + q_dist); %CALCULATING THE PLATFORMS ACCELARATION IN {B}
  u(index+1) = (u(index)+vdot(1,1)*dt);
  v(index+1) = (v(index)+vdot(2,1)*dt); %VELOCITIES IN {B}
  r(index+1) = (r(index)+vdot(3,1)*dt);
  b_velocities=[u(index);v(index);r(index)]; 
  
  
  xdot = (R*b_velocities); %VELOCITIES IN {I}
  
  
  x(index+1)=(x(index)+(xdot(1,1)*dt)); % X,Y,THETA IN {I}
  y(index+1)=(y(index)+(xdot(2,1)*dt));
  theta(index+1)=(theta(index)+(xdot(3,1)*dt));
  
  
  if(mod(index,100)==0)  %MAKING THE WIND
    gonia_anemou(index+1)=mod(10*cosd(3*index*dt)+30,360); 
  else
    gonia_anemou(index+1)=gonia_anemou(index);
  endif
  
  index=index+1;
endfor
printf("sum_avg_jet_force = %d\n ",sum_avg_jet_force);
%%%%%%% PLOT DATA %%%%%%
j=1:1:index-1;
plot(x,y);title("X/Y");xlabel("X");ylabel("Y");xlim([0 600]); ylim([-100 500]);
figure
plot(time(j),(x_wind(j)));title("X WIND");xlabel("TIME");ylabel("X WIND POWER"); 
figure
plot(time(j),(y_wind(j)));title("Y WIND");xlabel("TIME");ylabel("Y WIND POWER"); 
figure
plot(time(j),(n_wind(j)));title("N WIND");xlabel("TIME");ylabel("N WIND POWER");
%figure
%plot(time(j),x(j));title("X");xlabel("TIME");ylabel("X"); ylim([0 600]);
%figure
%plot(time(j),y(j));title("Y");xlabel("TIME");ylabel("Y"); ylim([-100 500]);
%figure
%plot(time(j),theta(j));title("THETA");xlabel("TIME");ylabel("THETA");
figure
plot(time,fa,time,fb,time,fc);title("JET FORCES");xlabel("TIME");ylabel("JET FORCE");
%figure
%plot(time,jet_forces_avg);title("JET AVERAGE FORCES");xlabel("TIME");ylabel("JET AVERAGE FORCE");


%%% MAKING THE MATRIX PINAKAS FOR THE ANIMATION WITH THE REAL DATA %%%%%%
pinakas(1,:) = x(1:index); 
pinakas(2,:) = y(1:index);
pinakas(3,:)= theta(1:index);
%%%%%%% SAVE THE DATA %%%%%%%%%
save result1.mat pinakas;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%