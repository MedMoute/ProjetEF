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

subplot(2,2,1);
scatter3(X,Y,Uex,'.','cdata',cdata_uex);
title('Solution exacte');
view([15 40]);
subplot(2,2,2);
scatter3(X,Y,Uglo,'.','cdata',cdata_uglo);
title('Solution calculée par la Méthode de Jacobi');
view([15 40]);
subplot(2,2,[3 4]);
scatter3(X,Y,Errloc,'o','filled','cdata',cdata_err);
title('Erreur absolue de la Méthode de Jacobi');
xlabel('X');
ylabel('Y');
view([0 90]);
colormap(color);
caxis([0,max(Errloc)]);
colorbar;

nrelm=msh.nbElm;
nnod=s;

Edof=zeros(nrelm,4);
Ex=zeros(nrelm,3); Ey=Ex; Ez=Ex;
Edof(:,1)=sort(1:nrelm);
Elm(:,1)=sort(1:nrelm);

for i=1058:nrelm 
    Ey(i,:)=[msh.POS(msh.ELE_NODES(i,2),2) msh.POS(msh.ELE_NODES(i,3),2) msh.POS(msh.ELE_NODES(i,4),2)];
    Ez(i,:)=[msh.POS(msh.ELE_NODES(i,2),3) msh.POS(msh.ELE_NODES(i,3),3) msh.POS(msh.ELE_NODES(i,4),3)];
    Ex(i,:)=[msh.POS(msh.ELE_NODES(i,2),1) msh.POS(msh.ELE_NODES(i,3),1) msh.POS(msh.ELE_NODES(i,4),1)];    
end


figure(2);
hold on;
for i=1058:nrelm;
    plot(Ex(i,[1:end 1]),Ey(i,[1:end 1]));
end
axis equal

clear s S;