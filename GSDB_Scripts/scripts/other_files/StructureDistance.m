function [ dist] = StructureDistance(xyz,n)
%% This finds the culidean distance between ppoints of the cordintes
% xyz_cord = The cordinates of final structure
% n = The number of regions
% Extract distanc efrom the Model
dist = zeros(n,n);
 for i=1:n
     
         for j=1:n           
               %Calculayte the Euclidean distance between distance from
               %structure and the one from genome location
              dist(i,j)=norm(xyz(i,:)-xyz(j,:)); % find the distance from i to rest(i+1 to n)
            
         end   
     
 end
 
end


