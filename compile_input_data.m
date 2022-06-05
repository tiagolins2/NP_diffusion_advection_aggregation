v=10; num_classes=5; depth=1000; max_size=1000; Cnp=1; t_end=60*60*24*30; salinity_a=500; rho_p=1050; T=298.15; Nss_conc=1; D=0; 
save_files=1;
Cnp=1; Nss_conc=1; Cdom=0; Cnp_add=1;
% for D=[0 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2]
% %htxy5_newD
%     htxy5_8_2
% end
D=1E-5; Rdom=1e-8;
numberphi=10^(log10(1000)/10);
numberx=10^(log10(1000)/50);
%arr_sizes=arr_sizes*(1e-2*numberx^ccc);
resus_ratio_vec=0.001*numberphi.^linspace(0,10,11);
%phi_idx=9:51
%resus_ratio_vec=0.001+linspace()
phi=0;
inital_c_add=Cnp_add*1e-6/(rho_p*4/3*pi*((100e-9)/2)^3)*1;
top_mat=zeros(51,10,3); bottom_mat=zeros(51,10,3);  mid_mat=zeros(51,10,3); iii=1;
input=inital_c_add(1)/(60*60*24)*1/dx;
iii=1;
Fracs_v=zeros(4,51,51); single_vec=zeros(51,51,1); homo_vec=zeros(51,51,1); het_vec=zeros(51,51,1);
for resus_idx=1:10
    resus_ratio=resus_ratio_vec(resus_idx);
    prop_name=strcat('Input_RR=',num2str(resus_ratio),'_Phi=',num2str(phi),'_D=',num2str(D),'m2s-1_Cnp=',num2str(Cnp),'_Nss_conc=',num2str(Nss_conc),'_num_classes=',num2str(num_classes),'_T=',num2str(T(1)),'_salinity_a=',num2str(salinity_a(1)),'_depth=',num2str(depth));
    %for phi=0.1*1.1482.^linspace(0,50,51)
    conc_np=Cnp.*1e-6./(rho_p.*4/3*pi.*((100*1e-9)/2).^3);
    avgconc=conc_np+input*(tt)*dx/depth*1e-3*numberx.^linspace(0,50,51);
%pause(0.5); 
    %prop_name=strcat('Phi=',num2str(phi),'_D=',num2str(D),'m2s-1_Cnp=',num2str(Cnp),'_Nss_conc=',num2str(Nss_conc),'_num_classes=',num2str(num_classes),'_T=',num2str(T),'_salinity_a=',num2str(salinity_a),'_depth=',num2str(depth));
    Nt=zeros(10,51,3);
    N_t_plot=dlmread(strcat(prop_name,'\NP_concentration.xlsx'));
    Nt(:,:,1)=squeeze(N_t_plot(:,end-50:end))./avgconc; Nt(:,:,2)=1-squeeze(N_t_plot(:,end-50:end))./avgconc; Nt(:,:,3)=1-squeeze(N_t_plot(:,end-50:end))./avgconc;
    Nt1=zeros(10,51,3); 
    N_arr_plot=dlmread(strcat(prop_name,'\N_arr.xlsx'));  
    N_agg_arr_plot=dlmread(strcat(prop_name,'\Nagg_arr.xlsx')); 
    N_het_arr_plot=dlmread(strcat(prop_name,'\Nhet_arr.xlsx'));
    Ntotal=N_arr_plot+N_agg_arr_plot+N_het_arr_plot;
    N_mat(:,:,iii)=N_arr_plot'; top_mat(:,iii,1)=N_arr_plot(:,1)./Ntotal(:,1); bottom_mat(:,iii,1)=N_arr_plot(:,10)./Ntotal(:,10); mid_mat(:,iii,1)=mean(N_arr_plot(:,5:6),2)./mean(Ntotal(:,5:6),2);
    Nagg_mat(:,:,iii)=N_agg_arr_plot'; top_mat(:,iii,2)=N_agg_arr_plot(:,1)./Ntotal(:,1); bottom_mat(:,iii,2)=N_agg_arr_plot(:,10)./Ntotal(:,10); mid_mat(:,iii,2)=mean(N_agg_arr_plot(:,5:6),2)./mean(Ntotal(:,5:6),2);
    Nhet_mat(:,:,iii)=N_het_arr_plot'; top_mat(:,iii,3)=N_het_arr_plot(:,1)./Ntotal(:,1); bottom_mat(:,iii,3)=N_het_arr_plot(:,10)./Ntotal(:,10); mid_mat(:,iii,3)=mean(N_het_arr_plot(:,5:6),2)./mean(Ntotal(:,5:6),2);
    Nt1(:,:,1)=squeeze(N_arr_plot'./Ntotal');
    Nt1(:,:,2)=squeeze(N_agg_arr_plot'./Ntotal');
    Nt1(:,:,3)=squeeze(N_het_arr_plot'./Ntotal');
    Nt1 = imresize(Nt1, 8);
    
    
    
    %Nt(:,:,3)=squeeze(N_t_plot(:,end-50:end))./conc_np;
    %top_mat(:,iii,1)=(squeeze(N_t_plot(1,end-50:end))./conc_np)';
    Nt = imresize(Nt, 8);
imshow(Nt1); axis on;
hold on; axis on;
    [C,h] = contour((squeeze(Nt(:,:,1))),[0.01 linspace(0.05,0.95,19) 0.99]); 
    %[C,h] = contour((squeeze(Nt(:,:,1)))); 
    colormap(cmap);
    clabel(C,h); set(gcf,'color','w');
    
    numberx=10^(log10(max_size)/50);
    sizes=(numberx.^linspace(0,50,51)); 
    xticks(8*[1 1+log(0.1/(1e-2))/log(numberx) 1+log(1/(1e-2))/log(numberx) 1+log(10/(1e-2))/log(numberx)])
    set(gca, 'xTickLabel', {'1' '10' '100' '1000'}); 
    yticks([1 8*5 8*10]);
    set(gca, 'yTickLabel', {'0' '500' '1000'}); pause(0.1);
    xlabel('Inflow (ng/L)'); ylabel('Depth (m)');
    
    img = getframe(gcf);
   % writeVideo(vidWriter,img);
   title(strcat('Resuspension rate=',num2str(resus_ratio)));
    set(gcf,'color','w');
    
    if iii==1 || iii==3 || iii==6 || iii==10
       saveas(gca,strcat('RR=',num2str(resus_ratio),'_profile.png'));  
    end
    Fracs_v(:,:,iii)=dlmread(strcat(prop_name,'\Fracs.xlsx'));
    single_vec(:,iii)=Fracs_v(3,:,iii);
    homo_vec(:,iii)=Fracs_v(1,:,iii);
    het_vec(:,iii)=Fracs_v(2,:,iii);
    %homo_vec=zeros(51,51,1); het_vec=zeros(51,51,1);
iii=iii+1;
end
Cnp=1;


img1=zeros(51,51,3);
img1(:,:,1)=single_vec;
img1(:,:,2)=homo_vec;
img1(:,:,3)=het_vec;
img1 = imresize(img1, 8);
imshow(img1); axis on;
%close(vidWriter);

%iii=1;
%figure;
% N_mat=zeros(10,51,51); Nagg_mat=zeros(10,51,51); Nhet_mat=zeros(10,51,51); iii=1; top_mat=zeros(51,51,3); bottom_mat=zeros(51,51,3);  mid_mat=zeros(51,51,3);
% for Nss_conc=0.1*1.1482.^linspace(0,50,51)
%     conc_np=Cnp.*1e-6./(rho_p.*4/3*pi.*((1e-2*(10^(log10(max_size)/50)).^linspace(0,50,51)*100*1e-9)./2).^3);
%     %pause(0.5); 
%     prop_name=strcat('Cnp=',num2str(Cnp),'_Nss_conc=',num2str(Nss_conc),'_num_classes=',num2str(num_classes),'_T=',num2str(T),'_salinity_a=',num2str(T),'_D=',num2str(D),'_depth=',num2str(depth));
%     Nt=zeros(10,51,3);
%     N_arr_plot=dlmread(strcat(prop_name,'\N_arr.xlsx'));  N_mat(:,:,iii)=N_arr_plot'; top_mat(:,iii,1)=N_arr_plot(:,1); bottom_mat(:,iii,1)=N_arr_plot(:,10); mid_mat(:,iii,1)=mean(N_arr_plot(:,5:6),2);
%     N_agg_arr_plot=dlmread(strcat(prop_name,'\Nagg_arr.xlsx')); Nagg_mat(:,:,iii)=N_agg_arr_plot'; top_mat(:,iii,2)=N_agg_arr_plot(:,1); bottom_mat(:,iii,2)=N_agg_arr_plot(:,10); mid_mat(:,iii,2)=mean(N_agg_arr_plot(:,5:6),2);
%     N_het_arr_plot=dlmread(strcat(prop_name,'\Nhet_arr.xlsx')); Nhet_mat(:,:,iii)=N_het_arr_plot'; top_mat(:,iii,3)=N_het_arr_plot(:,1); bottom_mat(:,iii,3)=N_het_arr_plot(:,10); mid_mat(:,iii,3)=mean(N_het_arr_plot(:,1),2);
%     Nt(:,:,1)=squeeze(N_arr_plot');
%     Nt(:,:,2)=squeeze(N_agg_arr_plot');
%     Nt(:,:,3)=squeeze(N_het_arr_plot');
%     Nt = imresize(Nt, 8);
% %imshow(Nt); 
% iii=iii+1;
% end


%for top/middle/bottom plot only
% hold on; axis on;
% xticks(8*[1 1+16.6618 1+33.3236 1+50]);
% set(gca, 'XTickLabel', {'0.1' '1' '10' '100'});
% numberx=10^(log10(max_size)/50);
% yticks(8*[1 1+log(0.1/(1e-2))/log(numberx) 1+log(1/(1e-2))/log(numberx) 1+log(10/(1e-2))/log(numberx)])
% set(gca, 'YTickLabel', {'1' '10' '100' '1000'});
% xlabel('Concentration'); ylabel('Primary particle size');
% set(gcf,'color','w');

%for top/middle/bottom plot only
top=imresize(top_mat, [408 200]); img=figure; imshow(top);
hold on; axis on;
numberx=10^(log10(max_size)/50);
xticks(20*[1 10]);
set(gca, 'XTickLabel', {'0.001' '1'});
yticks(8*[1 1+log(0.1/(1e-2))/log(numberx) 1+log(1/(1e-2))/log(numberx) 1+log(10/(1e-2))/log(numberx)]);
set(gca, 'YTickLabel', {'1' '10' '100' '1000'});
ylabel('Inflow (ng/L)'); xlabel('Resuspension rate');
set(gcf,'color','w');
saveas(img,'RR_top.emf'); set(gca,'FontSize',16)
% 
% 
mid=imresize(mid_mat, [408 200]); img=figure; imshow(mid);
hold on; axis on;
xticks(8*[1 1+16.6618 1+33.3236 1+50]);
set(gca, 'XTickLabel', {'0.1' '1' '10' '100'});
numberx=10^(log10(max_size)/50);
yticks(8*[1 1+log(0.1/(1e-2))/log(numberx) 1+log(1/(1e-2))/log(numberx) 1+log(10/(1e-2))/log(numberx)])
xticks(20*[1 10]);
set(gca, 'XTickLabel', {'0.001' '1'});
set(gca, 'YTickLabel', {'1' '10' '100' '1000'});
ylabel('Inflow (ng/L)'); xlabel('Resuspension rate');
set(gcf,'color','w');
saveas(img,'RR_mid.emf'); set(gca,'FontSize',16)
% 
% 
bottom=imresize(bottom_mat, [408 200]); img=figure; imshow(bottom);
hold on; axis on;
numberx=10^(log10(max_size)/50);
xticks(8*[1 1+log(0.01/0.001)/log(numberx) 1+log(0.1/0.001)/log(numberx) 1+log(1/0.001)/log(numberx)]);
set(gca, 'XTickLabel', {'0.001' '0.01' '0.1' '1'});
yticks(8*[1 1+log(0.1/(1e-2))/log(numberx) 1+log(1/(1e-2))/log(numberx) 1+log(10/(1e-2))/log(numberx)]);
set(gca, 'YTickLabel', {'1' '10' '100' '1000'});
yticks(8*[1 1+log(0.1/(1e-2))/log(numberx) 1+log(1/(1e-2))/log(numberx) 1+log(10/(1e-2))/log(numberx)])
xticks(20*[1 10]);
set(gca, 'XTickLabel', {'0.001' '1'});
set(gca, 'YTickLabel', {'1' '10' '100' '1000'});
ylabel('Inflow (ng/L)'); xlabel('Resuspension rate');
set(gcf,'color','w');
saveas(img,'RR_bottom.emf'); set(gca,'FontSize',16)
% 
% 
% 
