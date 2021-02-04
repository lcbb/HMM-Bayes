function draw_std(means,sigmas,colors)
    [~,num] = size(means);
    
    hold all
    
    for i=1:num
        plot(means(1,i),means(2,i),'kx');
        circle(means(1,i),means(2,i),sigmas(i),colors(i,:));
        circle(means(1,i),means(2,i),2*sigmas(i),colors(i,:));
    end
end

function circle(x,y,r,color)
    %x and y are the coordinates of the center of the circle
    %r is the radius of the circle

    ang=0:0.001:2*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    plot(x+xp,y+yp,'Color',color);
end