nombre = 'AnimaAuto.avi';
v = VideoWriter(nombre);
open(v);
t = linspace(0,2*T,200);
cs = linspace(0,3,61);
shg
for i=1:61
    c = cs(i)
    u = table2array(readtable(['.\SolucionesTol2\N101C',num2str(c) ,'s0.5.csv']));
    estabfrac(101,c,u,2*pi/0.8,1,0.5);
frame = getframe(gcf);
writeVideo(v,frame);
end
close(v);
