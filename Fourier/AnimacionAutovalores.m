% nombre = 'Autovalores_Fourier0.5.avi';
% v = VideoWriter(nombre);
% open(v);
n=485;
cs = linspace(0.16,5,n);
shg
for i=1:n
    c = cs(i);
    texto = ['.\Solucion_Fourier\m750s0.5km20\750C',num2str(c),'.csv'];
    u =table2array(readtable(texto));
    [~,D]=estabfracfourier2(750,c,u,2*pi/0.8,1,0.5);
    txt = ['.\Autovalores\s0.5m750\m750C',num2str(c),'.csv'];
    writematrix(diag(D),txt)
% frame = getframe(gcf);
% writeVideo(v,frame);
end
% close(v);
