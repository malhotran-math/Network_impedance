clc
clearvars;
%for p= 10.^[-4:1:0]
    
    %% ParamËtres de simulation
    nombre_de_cas = 100;
    Fichier = 0;
    
    Nombre_de_producteurs  = 3;
    Nombre_de_consomateurs = 3;
    Bus = Nombre_de_producteurs + Nombre_de_consomateurs + 3;
    
    %% Initialisations des puissances avec impÈdances connues
    r = [0      0.17    0.39    0       0.0119 0.0085   0       0.032   0.01];
    x = [0.0576 0.092   0.17    0.0586  0.1008 0.072    0.0625  0.161   0.085];
    b = [ 0     0.1580  0.3580  0       0.2090 0.1490   0       0.3060  0.1760];
    
    r_ = zeros(nombre_de_cas,length(r));
    x_ = zeros(nombre_de_cas,length(x));
    b_ = zeros(nombre_de_cas,length(b));
    
    N = 0; %ParamËtre de la loi normale des paramÈtres du rÈseau (ImpÈdeance)
    
    mpc = loadcase(case9);
    mpc.branch(:,3) = r';
    mpc.branch(:,4) = x';
    mpc.branch(:,5) = b';
    
    results = runopf(mpc);
    
    %On Stocke les valeurs des puissances dans le cas o˘ l'on connait
    %prÈcisemment les valeurs des impÈdances
    results.gen(:,2); %Producteurs en W    [P]
    results.gen(:,3); %Producteurs en  Var [Q]
    ProducteursW_init = results.gen(:,2);
    ProducteursQ_init = results.gen(:,3);
    
    results.bus(:,3); %Consomateur en W    [P]
    results.bus(:,4); %Consomateurs en var [Q]
    ConsomateursW_init = results.bus(:,3);
    ConsomateursQ_init = results.bus(:,4);
    %% DÈclaration tableau calculs
    Puissances_gen        = zeros(Nombre_de_producteurs,nombre_de_cas);
    Puissances_gen_test   = zeros(Nombre_de_producteurs,nombre_de_cas);
    Puissances_bus        = zeros(Nombre_de_producteurs,nombre_de_cas);
    
    Pertes_ligne_W_inc = zeros(Nombre_de_producteurs,nombre_de_cas);
    Pertes             = zeros(Nombre_de_producteurs,nombre_de_cas);
    
    
    Puissance_prod_W = mpc.gen(:,2);
    Puissance_prod_Q = mpc.gen(:,3);
    
    %% Programme OPF -> Puissance -> PF -> Pertes (avec la puissance initiale)
    %%Variation lois logarythmique
for N = 0:0.1:1
    for i = 1:nombre_de_cas
        for j = 1:length(r)
            r_(i,j)=r(j)*(1+normrnd(0,N));
            x_(i,j)=x(j)*(1+normrnd(0,N));
            b_(i,j)=b(j)*(1+normrnd(0,N));
        end
        %%% OPF
        mpc = loadcase(case9);
        mpc.branch(:,3) = r_(i,:)';
        mpc.branch(:,4) = x_(i,:)';
        mpc.branch(:,5) = b_(i,:)';
        
        results = runopf(mpc);
        Puissances_gen_test(:,i)= results.gen(:,2);
        
        
        %%% PF
        mpcpf = loadcase(case9);
        mpcpf.branch(:,3) = r';
        mpcpf.branch(:,4) = x';
        mpcpf.branch(:,5) = b';
        
        mpcpf.gen(:,2)  =results.gen(:,2);
        mpcpf.gen(:,3)  =results.gen(:,3);
        
        mpcpf.bus(:,3)  =results.bus(:,3);
        mpcpf.bus(:,4)  =results.bus(:,4);
        
        results_pf = runpf(mpcpf);
        
        %%Calculs des pertes
        Puissances_gen(:,i) = results_pf.gen(:,2);
        ind                 = (results_pf.bus(:,3)~=0); %Enlever les zÈro
        Puissances_bus(:,i) = -results_pf.bus(ind,3);
        %-- Pertes lignes Puissance rÈseau inconu retransposÈ sur le rÈseau
        %connu (OPF puis PF)
        Pertes_ligne_W_inc(:,i) = Puissances_gen(:,i);   %- Puissance_prod_W;
    end
    %--Pertes lignes rÈseau connu (OPF)
    %Pertes_ligne_W          = ProducteursW_init     - Puissance_prod_W;
    
    %%--Pertes entre un rÈseau avec impÈdance connus et impÈdances mal
    %%rÈglÈe
    for m = 1:nombre_de_cas
        Pertes(:,m) = Pertes_ligne_W_inc(:,m) - ProducteursW_init;
        %plot(Pertes(:,m)','x')
    end
    
    str = strcat('ZZ_',num2str(Fichier),'_Pertes_.mat');
    save(str,'Pertes',num2str(Fichier));
    
    str = strcat('ZZ_',num2str(Fichier),'_Pertes_ligne_W_inc.mat');
    save(str,'Pertes_ligne_W_inc',num2str(Fichier));
    
    str = strcat('ZZ_',num2str(Fichier),'_ProducteursW_init.mat');
    save(str,'ProducteursW_init',num2str(Fichier));
    
    Fichier = Fichier+1;
end
%      figure()
%      boxplot(Pertes')
%      ylabel('Pertes (MW)')
%      xlabel('Producteurs')
%     
    
