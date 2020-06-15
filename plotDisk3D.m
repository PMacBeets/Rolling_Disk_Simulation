function plotDisk3D(center,normal,diam,g)
% rotate by angle gamma


theta1=g:0.1:g+pi;
theta2=g+pi:0.1:g+2*pi;


v=null(normal);

% These points are for drawing a disk with z pointing perpendicular to
% plane - therefor normal x and z spots are swapped
point1=repmat(center',1,size(theta1,2))+diam*(v(:,1)*cos(theta1)+v(:,2)*sin(theta1));
point2=repmat(center',1,size(theta2,2))+diam*(v(:,1)*cos(theta2)+v(:,2)*sin(theta2));
fill3(point1(1,:),point1(2,:),point1(3,:),'r-'); hold on
fill3(point2(1,:),point2(2,:),point2(3,:),'y-')


xlabel('x')
ylabel('y')
zlabel('z')
end