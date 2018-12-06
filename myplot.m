function [mplot]=myplot(x,V,F,N,SIG,l,v,Energy)%,%n,im,f,map,fig)
% % fig = figure();
% % создание первого пустого кадра
% set(fig,'Position',[350,200,700,700]);
% frame = getframe(fig);
    Energy_v=Energy*ones(N,1);

     x_coord(1:N) = x(1,1:N);
     y_coord(1:N) = x(2,1:N);
     subplot(2,2,1);
     
     scatter(x_coord,y_coord);
   %  axis([-l 2*l -1 2*l ]);
%          hold on;
     detF=ones(1,N);
     for i = 1:N
           detF(1,i)=det(F(1:2,1:2,i));
     end
       
        tri=delaunay(x_coord,y_coord);
        subplot(2,2,2);
        trisurf(tri,x_coord,y_coord,detF(1,1:N));
%          xlim([0, 2*l]);
%          ylim([0, 2*l]);
%          zlim([0.99, 1.01]);
%          hold on;
         
        subplot(2,2,3);
        trisurf(tri,x_coord,y_coord,v(1,1:N));
%          xlim([0, 2*l]);
%          ylim([0, 2*l]);
%          zlim([-0.8, 0.8]);
%       hold on;
       errSIG=zeros(1,N);
       for i = 1:N
               errSIG(i)=SIG(1,2,i);
       end
       
        subplot(2,2,4);
        trisurf(tri,x_coord,y_coord, Energy_v);  
%          xlim([0, 2*l]);
%          ylim([0, 2*l]);
%          zlim([-1, 1]);
         pause(0.000001);
%         hold on;
      
%       
%      
%         frame = getframe(fig);
%         [im,map] = rgb2ind(frame.cdata,4);
%         imwrite(im,map,'animation3.gif','DelayTime',0.1,'WriteMode','Append');

mplot=0;