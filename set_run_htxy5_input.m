v=10; num_classes=10; depth=10; max_size=1000; Cnp=1; t_end=60*60*24*30; salinity_a=500; rho_p=1050; T=298.15; Nss_conc=1; D=0; 
save_files=1;
Cnp=1; Nss_conc=1; Cdom=0; Cnp_add=1;
% for D=[0 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2]
% %htxy5_newD
%     htxy5_8_2
% end
D=1E-5; Rdom=1e-8;
numberphi=10^(log10(1000)/10);
%arr_sizes=arr_sizes*(1e-2*numberx^ccc);
resus_ratio_vec=0.001*numberphi.^linspace(0,10,11);
%phi_idx=9:51
%resus_ratio_vec=0.001+linspace()
phi=0;
for resus_idx=1:10
    resus_ratio=resus_ratio_vec(resus_idx);
    htxy5_newD3;
end