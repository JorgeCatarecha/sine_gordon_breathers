nombre = 'Autovalores_Clásico.avi';
%v = VideoWriter(nombre);
%open(v);
n=41;
cs = linspace(0.79,0.83,n);
shg
for i=1:n
    c = cs(i);
    texto = ['.\Soluciones\N41C',num2str(c),'.csv'];
    u =table2array(readtable(texto));
    [~,D]=estab(41,c,u,2*pi/0.8,0);
    txt = ['.\Autovalores\n41C',num2str(c),'.csv'];
    writematrix(diag(D),txt)
%frame = getframe(gcf);
%writeVideo(v,frame);
end
%close(v);