function d = dif(u,h)

d = [-3*u(:,1)+4*u(:,2)-u(:,3),...
    u(:,3:end)-u(:,1:end-2),...
    u(:,end-2)-4*u(:,end-1)+3*u(:,end)]/(2*h);