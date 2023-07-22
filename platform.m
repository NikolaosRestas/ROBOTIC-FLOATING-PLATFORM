function [] = platform(pinakas)
  
  figure;
  grid on;
  
  h(1) = line([-16.82 , 27.64],[17.5 , 0]); %making the platform in 0 degrees
  h(2) = line([27.64 , -16.82],[0,-17.5]);
  h(3) = line([-16.82, -16.82],[-17.5, 17.5]);
  
  %h(1) = line([-17.36 , 27.64],[-17.5 , 0]); %making the platform in 0 degrees [X]->[Y]
  %h(2) = line([27.64 , -17.36],[0,17.5]);
  %h(3) = line([-17.36, -17.36],[17.7, -17.5]);
  
  
  %h(1) = line([0,-17.5],[27.64 , -17.36]); %making the platform in 90 degrees
  %h(2) = line([0,17.5],[27.64,-17.36]);
  %h(3) = line([-17.5, 17.5],[-17.36,-17.36]);
  
  %h(1) = line([-27.64,17.36],[0,17.5]); %making the platform in 180 degrees
  %h(2) = line([-27.64 , 17.36],[0,-17.5]);
  %h(3) = line([17.36, 17.36],[-17.5, 17.5]);
 
  
  model = hgtransform('Parent',gca); % making the object
  set(h,'Parent',model);
  axis(gca,"equal");
  axis(gca,[-100,600],[-100,600]);
  xlim([-100, 600]);
  ylim([-100, 600]);
  
  for t = 1:1:size(pinakas,2)
    dx = [pinakas(1,t); pinakas(2,t); 0]; % z-axis is always 0 (platform)
    translate = [eye(3) dx; zeros(1, 3) 1];
    theta = pinakas(3, t);
    R = [cosd(theta) -sind(theta) 0; 
         sind(theta) cosd(theta) 0;
         0 0 1];
    rotate = [R zeros(3, 1); zeros(1, 3) 1];
    set(model,'Matrix',translate*rotate);%* rotate);
    F(t) = getframe(gcf);
    drawnow;
  endfor
  
endfunction

function R = rotation(theta)
  
  R = [cosd(theta) -sind(theta) 0;
      sind(theta) cosd(theta) 0;
      0 0 1];
endfunction
