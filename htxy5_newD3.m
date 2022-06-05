%for validation
num_classes=5;
v=10; 
plot_profile=0; save_files=1; depth=10; salinity_a=500*ones(v,1); T=298.15*ones(v,1); rhow=1021*ones(v,1); %D=0e-4; 
% resus_ratio=0.001;

save_files=0;
%salinity_a=(600-linspace(0,1,v)*100)';
%ones(v,1)*500; %salinity_a(5,1)=650;
%T=((273.15+25)-linspace(0,1,v)*21)';
ccc=40; Frac_arr=[]; sites_covered=[]; time_elapsed=[]; N_arr=[]; Nagg_arr=[]; Nhet_arr=[]; Nij_aa=[]; N_t=zeros(v,1,100);  N_np=zeros(v,1,100); %num_classes=10; %figure;
%save_files=1;
%T=ones(v,1)*298.15; T(end-3:end)=273.15+4;
delta_h=ones(v,1); output=0;
avgconc=0;
while ccc<=50
    N_temp=[]; Ndom_arr=[]; Ndom_att_arr=[];
h=1;

N=zeros(num_classes,v); turnoffdom=0;
alpha_ss=zeros(v,1); alpha=zeros(v,1);
Nss=0.9999+0.0002*rand(num_classes,v); Nss=ones(num_classes,v); N_het=zeros(num_classes,v); N_sed=zeros(num_classes,1);
Nijss=zeros(num_classes,v);
Nnp_ss=zeros(num_classes,v);
Ndom=zeros(2,v); %1 is free, 2 is attached to NP
Rdom=10e-9;%+ccc*1e-9; 
rho_dom=1770; 
depth=1000;
%depth(1:2,:)=10; depth(h-1:h,:)=10; depth(:,1:2)=10; depth(:,v-1:v)=10;
dx=depth/(v); delta_h(:,:)=dx; delta_h(v,1)=dx*2;
size_p=100*1e-9;
arr_sizes=2.^(linspace(0,num_classes-1,num_classes));
arr_sizes=arr_sizes*size_p;
arr_sizes=arr_sizes';
%arr_sizes_ss=[1 1 2 2 400 400 800 800 1000 1000];  
arr_sizes_ss=[1 1 2 2 4 4 8 8 16 16];
%arr_sizes_ss=ones(1,num_classes);%[1 1 2 2 4 4 8 8 16 16];
arr_sizes_ss=[1 2 4 8 16]; 
arr_sizes_ss=500*1e-9*arr_sizes_ss;
%arr_sizes*5;2
arr_sizes=arr_sizes*1;
%max_size=1000; 
numberx=10^(log10(1000)/50);
%arr_sizes=arr_sizes*(1e-2*numberx^ccc); 
%sites=4*pi*(arr_sizes_ss).^2/(pi*(arr_sizes(1))^2);
sites=round(4*pi*(arr_sizes(1)/2)^2/(pi*(Rdom)^2)); sites(sites==0)=1; sites=0;
inital_area_ratio=1;

%inital_c=1.82e10*(100*1e-9)^3/arr_sizes(1)^3;
D_mat=zeros(5,v);

j=1;
b=depth/2; c=depth/10; a=1/(c*sqrt(2*pi));
%dx=depth/v;
%b=0; c=1; a=1/(c*sqrt(2*pi));
mu=8.9e-4;
g=9.81;
%T=273.15+25; 
mu=-9.80215e-9*T.^3+9.16232E-6*T.^2-2.86927E-3*T+3.01685E-1;
%salinity_a=500;
s=0; s=salinity_a/1000*58.44;
df=1.7;%+ccc*0.03;
mu=mu.*(1.9609e-3*s+9.9276e-1); %sites=zeros(5,1); 
rho_p=1050;
rho_ps=2600; %Cnp=0.0001;
inital_c=Cnp*1e-6/(rho_p*4/3*pi*(arr_sizes(1)/2)^3)*1; %ug/L
inital_c_add=Cnp_add*1e-6/(rho_p*4/3*pi*(arr_sizes(1)/2)^3)*1; %ug/L
input=inital_c_add(1)/(60*60*24)*1/dx*1e-3*numberx^ccc;
rho_w=rhow;
%rho_w=rho_w*(0.7223*s+998.29);
%rho_w=rho_w*(0.7223*s+4.497e-5*T.^3-4.518E-2*T.^2+1.469E1*T-5.576E2);
%rhow=rho_w;
u=(2*(arr_sizes(1)/2)^2*(rho_p-rhow)*9.81/(9*mu));
while j<=v
   %if j*dx<=depth/2+
     % N(1,j)=(inital_c*depth)*a*exp(-(j*dx-b)^2/(2*c^2)); 
   %end
    %  N(1,j)=-inital_c*2/depth*j*dx+2*inital_c;
     %N(1,j)=inital_c*2/depth*j*dx;
    % N(1,j)=inital_c*1e-6;
     %  N(1,j)=-inital_c*depth*u/(exp(-u*depth)-1)*exp(-u*dx*j);
     N(1,j)=inital_c;
     j=j+1;
end
%N(:,:)=0;
avgconc=sum(sum(N))*dx/depth; Ndom_conc=Cdom; 
Ndom=ones(v,1)*Ndom_conc*0.001/((4/3*pi()*(Rdom).^3)*rho_ps)*dx*1;
%Ndom=[Ndom zeros(v,1)]; 

Ndom=Ndom'; Ndom_att=(zeros(v,1))';
Ndom_inital=Ndom;
%Ndom=(Nss_conc*0.001/(sum(4/3*pi()*(arr_sizes_ss/2).^3)*rho_ps))*500^3/500^3;
j=1;
%N(1,1)=inital_c*depth/dx; 
N00=N; %initial conditions
%N(:,1,:)=0; N(:,:,1)=0; N(:,h,:)=0; N(:,:,v)=0;
%N(1,5,15)=1.82E13;
%Nss_conc=1; %mg/L
Nss(:,:)=Nss(:,:)*(Nss_conc*0.001/(sum(4/3*pi()*(arr_sizes_ss/2).^3)*rho_ps))*500^3/500^3;%*(1e-3*1.5^ccc); %Nss(:,1)=Nss(:,1)/2; %mg/L
%Nss(:,1,:)=0; Nss(:,:,1)=0; Nss(:,h,:)=0; Nss(:,:,v)=0;
% n0=n;
t=0;
K=zeros(num_classes,num_classes,v);
Kss=zeros(num_classes,num_classes,v);
hydro_size=[];
v_dndt=zeros(numel(arr_sizes),1);
kb=1.380648813131313131313131e-23;

t_end=5*12*30*24*60*60;
%t_end=60*60*1*1; %t_end=1e-5; 
%t_end=60*60*200*1; %t_end=43.6344; 
%t_end=60*2;
 
Hamaker_A_a=5.91E-21; % PS
Hamaker_Ass_a=(sqrt(3.7e-20)+sqrt(3.7e-20))^2; % Kao
densityNP_a=1050; %
densityNC_a=2600; %
densityNOM_a=1.33;
OM_size_a=0e-9;
OM_packing_a=phi;  
Temperature_a=298.15;
phi_pa=OM_packing_a; L=OM_size_a; Mw=500000;
N_sum=(2.^(linspace(0,num_classes-1,num_classes))*(N(1:end,:)));

sites_occupied=mean((Ndom_inital-Ndom)./(N_sum+1e-10));
%phi=sites_occupied/sites;


rho_pa=densityNOM_a*1000;
rho_L=rho_pa*phi_pa+rhow*(1-phi_pa);
V_L=4/3*pi*((arr_sizes(1)/2+L)^3-(arr_sizes(1)/2)^3);
rho_pL=(rho_p*4/3*pi*(arr_sizes(1)/2)^3+rho_L*V_L)/(4/3*pi*((arr_sizes(1)/2+L)^3));
arr_sizes=arr_sizes*(size_p+L*2)/size_p; %adjusted size to take into account DOM layer
rho_p=rho_pL; %adjusted density based on packing of organic ligands.

alpha=pb_xDLVO(salinity_a(1), arr_sizes(1)/2+2*Rdom*phi, T(1), 0.485, 0, 0*2*Rdom, rho_dom,Mw,Hamaker_A_a)*1;
alpha_ss=pb_het_DLVO(salinity_a(1),T(1),0.485, 0, 0*2*Rdom, rho_dom, Mw,(sqrt(Hamaker_A_a)+sqrt(3.7e-20))^2,Hamaker_Ass_a,arr_sizes(1)/2+2*Rdom*phi,arr_sizes_ss(1)/2)*1;

if isnan(alpha) alpha=0; end
if isnan(alpha_ss) alpha_ss=0; end
%alpha=0;
%


%alpha_ss=0.0;
t=0;
i=1;
j=1;
rhop=zeros(num_classes,v);
rhops=zeros(num_classes,1);
a1=arr_sizes(1)/2;
dt=10;
%D=kb*T/(6*pi*mu); 
%D=1*(0.01)^2; 
%D=D+D_out; %D=0;
while i<=num_classes
    ai=arr_sizes(i)/2;
    rhop(i,1:end)=((ai/a1)^df*4/3*pi*(a1+2*Rdom*phi)^3*mean(rho_p)+mean(rho_w)*(4/3*pi*ai^3-(ai/a1)^df*4/3*pi*(a1+2*Rdom*phi)^3))/(4/3*pi*(ai+2^i*Rdom*phi)^3);
    rhops(i)=rho_ps;%((ai/a1)^df*4/3*pi*(a1*5)^3*rho_ps+rho_w*(4/3*pi*(ai*5)^3-(ai/a1)^df*4/3*pi*(a1*5)^3))/(4/3*pi*(ai*5)^3);
    i=i+1;
end

for k=1:v
i=1;
while i<=num_classes
    j=1;
    while j<=num_classes
    %assuming no shear rate:
    ai=arr_sizes(i)/2+2^(i)*Rdom*phi;
    aj=arr_sizes(j)/2+2^(j)*Rdom*phi;
    K(i,j,k)=((2*kb*T(k))/(3*mu(k))*(ai+aj)^2/(ai*aj)+...
        2*pi*g/(9*mu(k))*(rhop(j,k)-rho_w(k))*(ai+aj)^3*(ai-aj));
    j=j+1;
    end
    i=i+1;
end
end    

for k=1:v
i=1;
j=1;
while i<=num_classes
    j=1;
    while j<=num_classes
        ai=arr_sizes_ss(i)/2;
        aj=arr_sizes(j)/2+2^(j)*Rdom*phi;
        vsj=2*aj^2*(rhop(j,k)-rho_w(k))*g/(9*mu(k));
        vsi=2*ai^2*(rho_ps-rho_w(k))*g/(9*mu(k));   
        Kss(i,j,k)=((2*kb*T(k))/(3*mu(k))*(ai+aj)^2/(ai*aj)+...
        pi*(aj+ai)^2*abs(vsj-vsi));
        j=j+1; 

    end
    i=i+1;
end
end
% dom+np
for k=1:v
i=1;
Kdom=zeros(num_classes,1);
sites_occupied=(Ndom_inital-inital_c*sites)/(Ndom_inital+1e-10);
while i<=num_classes
 
   % while j<=num_classes
        hydro_radius=arr_sizes(i)/2*((arr_sizes(i)/arr_sizes(1))/0.64)^(1/df);
        %ai=arr_sizes(i)/2;
        aj=Rdom;
        vsi=2*hydro_radius^2*(rhop(i)-rhow(k))*g/(9*mu(k));
        vsj=2*Rdom^2*(rho_dom-rhow(k))*g/(9*mu(k));   
        Kdom(i,1)=((2*kb*T(k))/(3*mu(k))*(hydro_radius+aj)^2/(hydro_radius*aj)+...
        pi*(aj+hydro_radius)^2*abs(vsj-vsi));
        %j=j+1; 

   % end
    i=i+1;
end
end
%
sites_occupied=mean((Ndom_inital-real(Ndom))./(real(N_sum)+1e-10));
%phi=sites_occupied/sites;
% if Ndom_inital-inital_c*sites>0
% t_initial=log((Ndom_inital-inital_c*sites)/Ndom_inital)/(-1*Kdom(1)*inital_c); %t_end=t_initial;
% else
% t_initial=0;    
% end
%alpha=pb_xDLVO(salinity_a, arr_sizes(1)/2, T, 0.485, phi, 2*Rdom, rho_dom,Mw,Hamaker_A_a)*1;
%alpha_ss=pb_het_DLVO(salinity_a,T,0.485, phi, 2*Rdom, rho_dom, Mw,(sqrt(Hamaker_A_a)+sqrt(3.7e-20))^2,Hamaker_Ass_a,arr_sizes(1)/2,arr_sizes_ss(1)/2)*1;

num_steps=120;
dt=t_end/num_steps; %t_end=t_initial*1.01;
ii=1;
%
skip_initial=1;
tt=0; tti=0;
exp_v=t_end^(1/num_steps);
while ii<=num_steps
%     if phi<0.9999 && t_initial>0 %&& skip_initial==1
%        dt=t_initial/5; turnoffdom=0; num_steps=num_steps+1;
%     else
%        dt=(t_end-t_initial)/(num_steps); turnoffdom=1*(t_initial>0);
%         %if ii==1
%         %t_initial=t_end/10000;   skip_initial=0; num_steps=num_steps+1;
%         %dt=t_initial;
%         %end
%     end


Nss(:,:)=Nss(:,:)+resus_ratio*ones(num_classes,v).*(Nss_conc*0.001/(sum(4/3*pi()*(arr_sizes_ss/2).^3)*rho_ps))*500^3/500^3;

dt=((exp_v)^ii-(exp_v)^(ii-1)); %
dt=t_end/num_steps;
if phi>0.999
    turnoffdom=1;
end

% if ii>30
%    input=0; 
% end
time_elapsed=[time_elapsed; tt];    
tti=tt;    
tt=tt+dt;


N(N<0)=0;
N_total=[N(:); Nss(:); N_het(:); N_sed(:); Nijss(:); Nnp_ss(:); Ndom(:); Ndom_att(:)];
N_temp=[N_temp;(N(:))'];
N0=N_total;
%[Cnp ccc ii log10(N(1,1))]
%N0=N;
%input=inital_c*(4.03221e-9);
%output=1*input;
m=1; I=0; VZ=0; K_in=0; K_out=0; Vfish=0 ;%inital_c=0;

%if ii<=5
N_sum=(2.^(linspace(0,num_classes-1,num_classes))*(N(1:end,:)));%+sum(Nnp_ss);
N_lost=inital_c-(N_sum); sites=0;
%sites_occupied=((Ndom_inital-Ndom))./((N_sum+1e-10)); 
%sites_occupied=Ndom_att./((N_sum+1e-10));
%phi=sites_occupied/sites; sites_covered=[sites_covered; mean(phi)];
%theta_e=(min((Ndom./N_sum),sites))/sites;
%Ce=Ndom_inital-N_sum.*min((Ndom./N_sum),sites);
%theta=(sites_occupied)/sites;
%beta=(Ndom_inital-Ce)./theta_e;

Mw=33395*Rdom^2-92817*Rdom+73501; 
%phi=0;
for j=1:v
alpha(j)=pb_xDLVO2(salinity_a(j), arr_sizes(1)/2+2*Rdom*phi, T(j), 0.485, phi, 2*Rdom, rho_dom,Mw,Hamaker_A_a,4/3*pi*Rdom^3)*1;
alpha_ss(j)=1*pb_het_DLVO(salinity_a(j),T(j),0.485, phi, 2*Rdom, rho_dom, Mw,(sqrt(Hamaker_A_a)+sqrt(3.7e-20))^2,Hamaker_Ass_a,arr_sizes(1)/2+2*Rdom*phi,arr_sizes_ss(1)/2)*1;
end
if phi>1
   phi=1; 
end

%end
num_np=Nnp_ss./(abs(Nijss)+1e-10);


[t,y]=ode15s(@(t,N_total)htxy5_ode_input_old2(N_total,K,Kss,Kdom,dt,alpha,alpha_ss,arr_sizes,arr_sizes_ss,rhow,rhop,rho_ps,mu,m,I,VZ,K_in,K_out,Vfish,inital_c,dx,input,output,delta_h,v,D,depth,D_mat,sites,num_classes,turnoffdom,phi,Rdom),[tti tt]',N0);
%[t,y]=ode15s(@(t,N_total)htxy5_ode4(N_total,K(:,:,1),Kss(:,:,1),Kdom(:,:,1),dt,alpha(1),alpha_ss(1),arr_sizes,arr_sizes_ss,rhow(1),rhop(:,1),rho_ps,mu(1),m,I,VZ,K_in,K_out,Vfish,inital_c,dx,input,output,delta_h,v,D,depth,D_mat,sites,num_classes,turnoffdom,phi,Rdom),[tti tt]',N0);



a1=arr_sizes(1)/2+2*Rdom;
while i<=num_classes
    ai=arr_sizes(i)/2;
    rhop(i,1)=((ai/a1)^df*4/3*pi*(a1+2*Rdom*phi)^3*rho_p+rho_w*(4/3*pi*ai^3-(ai/a1)^df*4/3*pi*(a1+2*Rdom)^3)+i*sites_occupied*4/3*pi*Rdom^3*rho_dom)/(4/3*pi*(ai+2^i*Rdom*phi)^3);
    i=i+1;
end

for k=1:v
i=1;
while i<=num_classes
    j=1;
    while j<=num_classes
    %assuming no shear rate:
    ai=arr_sizes(i)/2+2^(i)*Rdom*phi;
    aj=arr_sizes(j)/2+2^(j)*Rdom*phi;
    K(i,j,k)=((2*kb*T(k))/(3*mu(k))*(ai+aj)^2/(ai*aj)+...
        2*pi*g/(9*mu(k))*(rhop(j,k)-rho_w(k))*(ai+aj)^3*(ai-aj));
    j=j+1;
    end
    i=i+1;
end
end  
    
for k=1:v
i=1;
j=1;
while i<=num_classes
    j=1;
    while j<=num_classes
        ai=arr_sizes_ss(i)/2;
        aj=arr_sizes(j)/2+2^(j)*Rdom*phi;
        vsj=2*aj^2*(rhop(j,k)-rho_w(k))*g/(9*mu(k));
        vsi=2*ai^2*(rho_ps-rho_w(k))*g/(9*mu(k));   
        Kss(i,j,k)=((2*kb*T(k))/(3*mu(k))*(ai+aj)^2/(ai*aj)+...
        pi*(aj+ai)^2*abs(vsj-vsi));
        j=j+1; 

    end
    i=i+1;
end
end

for k=1:v
i=1;
Kdom=zeros(num_classes,1);
sites_occupied=(Ndom_inital-inital_c*sites)/(Ndom_inital+1e-10);
while i<=num_classes
 
   % while j<=num_classes
        hydro_radius=arr_sizes(i)/2*((arr_sizes(i)/arr_sizes(1))/0.64)^(1/df);
        %ai=arr_sizes(i)/2;
        aj=Rdom;
        vsi=2*hydro_radius^2*(rhop(i)-rhow(k))*g/(9*mu(k));
        vsj=2*Rdom^2*(rho_dom-rhow(k))*g/(9*mu(k));   
        Kdom(i,1)=((2*kb*T(k))/(3*mu(k))*(hydro_radius+aj)^2/(hydro_radius*aj)+...
        pi*(aj+hydro_radius)^2*abs(vsj-vsi));
        %j=j+1; 

   % end
    i=i+1;
end
end
%
clf;
%hold on;
% plot(Nijss(1,:));
% plot(Nijss(2,:));
% plot(Nijss(3,:));
% plot(Nijss(4,:));
% plot(Nijss(5,:));
%plot((2.^(linspace(0,num_classes-1,num_classes)))*N+sum(Nnp_ss));
plot(2.^(linspace(1,num_classes-1,num_classes-1))*(N(2:end,:))./((2.^(linspace(0,num_classes-1,num_classes)))*N+sum(Nnp_ss)));
pause(0.0001); hold off;

[y_iter,last_d]=size(y);
N=reshape(y(y_iter,1:num_classes*v),num_classes,v);
Nss=reshape(y(y_iter,num_classes*v+1:num_classes*2*v),num_classes,v);
N_het=reshape(y(y_iter,2*num_classes*v+1:3*num_classes*v),num_classes,v);
N_sed=reshape(y(y_iter,3*num_classes*v+1:3*num_classes*v+num_classes),num_classes,1);
Nijss=reshape(y(y_iter,3*num_classes*v+num_classes+1:4*num_classes*v+num_classes),num_classes,v);
Nnp_ss=reshape(y(y_iter,4*num_classes*v+num_classes+1:5*num_classes*v+num_classes),num_classes,v);
Ndom=reshape(y(y_iter,5*num_classes*v+num_classes+1:5*num_classes*v+num_classes+v),1,v);
Ndom_att=reshape(y(y_iter,5*num_classes*v+num_classes+v+1:5*num_classes*v+num_classes+2*v),1,v);
Nss(Nss<=0)=0;% Nijss(Nijss<=0)=1e-20;  Nnp_ss(Nnp_ss<=0)=1e-20; 
%Ndom_att(Ndom_att<0)=0;
avgconc=avgconc+input*dt*dx/depth;
%
N_t(:,ccc+1,ii)=(2.^(linspace(0,num_classes-1,num_classes)))*N+sum(Nnp_ss);
N_np(:,ccc+1,ii)=arr_sizes(1)*(2.^(linspace(0,num_classes-1,num_classes))/0.64).^(1/df)*N./(sum(N)+1e-10);
Ndom_arr=[Ndom_arr; sum(Ndom)]; Ndom_att_arr=[Ndom_att_arr; sum(Ndom_att)];
%(2.^(linspace(0,num_classes-1,num_classes)))*N(1,:)+N(2,:)*2+N(3,:)*4+N(4,:)*8+N(5,:)*16+sum(Nnp_ss);
%

ii=ii+1;
end
N_arr=[N_arr; N(1,:)]; %Nagg_arr=[Nagg_arr; N(2,:)*2+N(3,:)*4+N(4,:)*8+N(5,:)*16];
Nagg_arr=[Nagg_arr; 2.^(linspace(1,num_classes-1,num_classes-1))*(N(2:end,:))];
Nhet_arr=[Nhet_arr; sum(Nnp_ss)];

%Nij_aa=[Nij_aa; (squeeze(sum(Nijss(1,:))+sum(Nijss(2,:))*2+sum(Nijss(3,:))*4+sum(Nijss(4,:,:))*8+sum(Nijss(5,:,:)*16)))'];
F_single=sum(N(1,:))*1/v; F_2=sum(N(2,:))/v; F_4=sum(N(3,:))/v; F_8=sum(N(4,:))/v; F_16=sum(N(5,:))/v;
F_homo=sum((2.^(linspace(1,num_classes-1,num_classes-1)))*N(2:end,:))*1/v;
F_het=sum(sum(Nnp_ss))*1/v;
F_sed=0;%(sum(N_sed(1,:))+sum(N_sed(2,:))*2+sum(N_sed(3,:))*4+sum(N_sed(4,:))*8+sum(N_sed(5,:))*16)*1/v;
F_total=F_sed+F_het+F_homo+F_single;

Fracs=[F_homo/avgconc; F_het/avgconc; F_single/avgconc;1-F_total/avgconc; std(N(1,:)/F_single); std(N(2,:)/F_2);...
    std(N(3,:)/F_4); std(N(4,:)/F_8); std(N(5,:)/F_16)];
Frac_arr=[Frac_arr Fracs];
[resus_idx ccc (Fracs(1:4))']

%hold on;
% for time1=1:ii-1
% plot(N_t(:,ccc+1,time1));
% %plot(N_t(:,ccc+1,time1));
% ylim([min(min(min(N_t(:,ccc+1,time1))),1) max(max(N_t(:,ccc+1,time1)))]); xlim([1 v]); set(gca, 'YScale', 'linear'); 
% pause(0.2);
% %error_mb=F_total/avgconc
% %dt/dx^2
% end

ccc=ccc+1;
%

%

end

if save_files==1
%properties
properties=[v; num_classes; depth; max_size; Cnp; D; t_end; salinity_a; rho_p; T; Nss_conc];
prop_name=strcat('Input_RR=',num2str(resus_ratio),'_Phi=',num2str(phi),'_D=',num2str(D),'m2s-1_Cnp=',num2str(Cnp),'_Nss_conc=',num2str(Nss_conc),'_num_classes=',num2str(num_classes),'_T=',num2str(T(1)),'_salinity_a=',num2str(salinity_a(1)),'_depth=',num2str(depth));
mkdir(prop_name);
%writetable(N_t,strcat(prop_name,'/NP_concentration.xlsx'),N_t,'delimiter','\t');
dlmwrite(strcat(prop_name,'/NP_concentration.xlsx'),N_t,'delimiter','\t');
%writetable(N_np,strcat(prop_name,'/NP_size.xlsx'),N_t,'delimiter','\t');
dlmwrite(strcat(prop_name,'/NP_size.xlsx'),N_np,'delimiter','\t');

dlmwrite(strcat(prop_name,'/N_arr.xlsx'),N_arr,'delimiter','\t');
dlmwrite(strcat(prop_name,'/Nagg_arr.xlsx'),Nagg_arr,'delimiter','\t');
dlmwrite(strcat(prop_name,'/Nhet_arr.xlsx'),Nhet_arr,'delimiter','\t');
dlmwrite(strcat(prop_name,'/Fracs.xlsx'),Frac_arr(1:4,:),'delimiter','\t');
end

if plot_profile==1
    
   for time1=1:ii-1
    conc_np=Cnp.*1e-6./(rho_p.*4/3*pi.*((1e-2*(10^(log10(max_size)/50)).^linspace(1,50,51)*100*1e-9)./2).^3);
    %prop_name=strcat('Cnp=',num2str(Cnp),'_Nss_conc=',num2str(Nss_conc),'_num_classes=',num2str(num_classes),'_T=',num2str(T),'_salinity_a=',num2str(T),'_D=',num2str(D),'_depth=',num2str(depth));
    Nt=zeros(10,11,3);
    N_t_plot=N_t(:,:,time1)./(repmat(conc_np,10,1)*4);
   % Nt(:,:,1)=1-0.5*squeeze(N_t_plot(:,end-50:end))./conc_np; Nt(:,:,2)=1-0.5*squeeze(N_t_plot(:,end-50:end))./conc_np; Nt(:,:,3)=1-0.5*squeeze(N_t_plot(:,end-50:end))./conc_np;
   % top_mat(:,iii,1)=(squeeze(N_t_plot(1,end-50:end))./conc_np)';
    N_t_plot = imresize(N_t_plot, 8);
    imshow(N_t_plot);
   end

end

