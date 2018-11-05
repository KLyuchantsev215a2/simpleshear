function [mplot]=myplot(x,V,F,N)

     x_coord(1:N) = x(1,1:N);
     y_coord(1:N) = x(2,1:N);
     subplot(2,2,1);
     scatter(x_coord,y_coord);
       
        
     detF=ones(1,N);
     for i = 1:N
           detF(1,i)=det(F(1:2,1:2,i));
     end
       
        tri=delaunay(x_coord,y_coord);
        subplot(2,2,2);
        trisurf(tri,x_coord,y_coord,detF(1,1:N));
        
        subplot(2,2,3);
        trisurf(tri,x_coord,y_coord,V(1,1:N));
        pause(0.0000001);
mplot=0;