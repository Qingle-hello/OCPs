function [mass_local,stiff_local,xg,wg]=local_matrix(L, M, r)
%LOCAL_MATRICX Summary of this function goes here
%   Detailed explanation goes here 

h=L/M; 

switch r
    case 1 
        xg=[-1,1];
        wg=[1,1];
        mass_local = h/2*diag(wg);
        stiff_local= zeros(r+1);
        a=xg(1);b=xg(2);
        phi{1}=@(x) x-x+1/(a-b);
        phi{2}=@(x) x-x+1/(b-a);
        for i=1:2
           for j=1:2
              stiff_local(i,j)=phi{i}(xg).*phi{j}(xg)*wg'*2/h;
           end
        end 
        
    case 2 
        xg=[-1,0,1];
        wg=[1/3,4/3,1/3];
        mass_local = h/2*diag(wg);
        stiff_local= zeros(r+1);
        a=xg(1);b=xg(2);c=xg(3);
        phi{1}=@(x) ((x-c)+(x-b))/((a-b)*(a-c));  
        phi{2}=@(x) ((x-a)+(x-c))/((b-a)*(b-c));  
        phi{3}=@(x) ((x-a)+(x-b))/((c-a)*(c-b)); 
        for i=1:3
           for j=1:3
              stiff_local(i,j)=phi{i}(xg).*phi{j}(xg)*wg'*2/h;
           end
        end         
        
    
    
    
    case 3
xg=[-1,-0.447213595499957939282,0.447213595499957939282,1];
wg=[1/6,5/6,5/6,1/6];
mass_local = h/2*diag(wg);
stiff_local= zeros(r+1);
a=xg(1);b=xg(2);c=xg(3);d=xg(4);
phi{1}=@(x) ((x-b).*(x-c)+(x-b).*(x-d)+(x-c).*(x-d))/((a-b)*(a-c)*(a-d));  
phi{2}=@(x) ((x-a).*(x-c)+(x-a).*(x-d)+(x-c).*(x-d))/((b-a)*(b-c)*(b-d));  
phi{3}=@(x) ((x-a).*(x-b)+(x-a).*(x-d)+(x-b).*(x-d))/((c-a)*(c-b)*(c-d));  
phi{4}=@(x) ((x-b).*(x-c)+(x-a).*(x-c)+(x-b).*(x-a))/((d-a)*(d-b)*(d-c));  
for i=1:4
    for j=1:4
    stiff_local(i,j)=phi{i}(xg).*phi{j}(xg)*wg'*2/h;
    end
end 

    case 4
xg=[-1,-0.6546536707079771437983,0,0.654653670707977143798,1];
wg=[0.1,0.544444444444444444444,0.7111111111111111111111,0.544444444444444444444,0.1];
mass_local = h/2*diag(wg);
stiff_local= zeros(r+1);
a=xg(1);b=xg(2);c=xg(3);d=xg(4);e=xg(5);
phi{1}=@(x) ((x-c).*(x-d).*(x-e) +(x-b).*(x-d).*(x-e) +(x-b).*(x-c).*(x-e) +(x-b).*(x-c).*(x-d)) /((a-b)*(a-c)*(a-d)*(a-e));  
phi{2}=@(x) ((x-c).*(x-d).*(x-e) +(x-a).*(x-d).*(x-e) +(x-a).*(x-c).*(x-e) +(x-a).*(x-c).*(x-d)) /((b-a)*(b-c)*(b-d)*(b-e));  
phi{3}=@(x) ((x-b).*(x-d).*(x-e) +(x-a).*(x-d).*(x-e) +(x-a).*(x-b).*(x-e) +(x-a).*(x-b).*(x-d)) /((c-a)*(c-b)*(c-d)*(c-e));   
phi{4}=@(x) ((x-b).*(x-c).*(x-e) +(x-a).*(x-c).*(x-e) +(x-a).*(x-b).*(x-e) +(x-a).*(x-b).*(x-c)) /((d-a)*(d-b)*(d-c)*(d-e)); 
phi{5}=@(x) ((x-b).*(x-c).*(x-d) +(x-a).*(x-c).*(x-d) +(x-a).*(x-b).*(x-d) +(x-a).*(x-b).*(x-c)) /((e-a)*(e-b)*(e-c)*(e-d)); 
for i=1:5
    for j=1:5
    stiff_local(i,j)=phi{i}(xg).*phi{j}(xg)*wg'*2/h;
    end
end 

end

