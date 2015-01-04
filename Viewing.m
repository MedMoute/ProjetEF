Uex=textread('u_exact');
Uglo=textread('u_global');
Errloc=abs(Uex-Uglo);
msh=load_gmsh4('fichierTest/testpart.msh',[]);
X=msh.POS(:,1);
Y=msh.POS(:,2);

S=size(X);
s=S(1);
cdata_uex=horzcat((Uex+1)/(max(Uex)+abs(min(Uex))),...
    zeros(s,1)...
    ,(-Uex+1)/(max(Uex)+abs(min(Uex))));

cdata_uglo=horzcat((Uglo+1)/(max(Uex)+abs(min(Uex))),...
    zeros(s,1),...
    (-Uglo+1)/(max(Uglo)+abs(min(Uglo))));

cdata_err=horzcat((Errloc)/(max(Errloc)),...
    zeros(s,1),...
    zeros(s,1));

color=horzcat((((0:1:s-1))/(s-1))',...
    zeros(s,1),...
    zeros(s,1));
clear s S;
subplot(2,2,1);
scatter3(X,Y,Uex,'.','cdata',cdata_uex);
view([15 40]);
subplot(2,2,2);
scatter3(X,Y,Uglo,'.','cdata',cdata_uglo);
view([15 40]);
subplot(2,2,[3 4]);
scatter3(X,Y,Errloc,'.','cdata',cdata_err);
view([0 90]);
colormap(color);
caxis([0,max(Errloc)]);
colorbar;
