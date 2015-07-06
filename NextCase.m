function [] = NextCase(d2,d3,a,b,m,C,epsilon,limit,start,acc)
%Coded by Michael Gallaugher for Dr. M. Wang NSERC Summer 2013
%Runge Kutta method for 3 equations g1, g2 and g3 and potential u

%Set initial values
    
     tt(1)=0;
    z1(1)=0;
    z3(1)=a;
    z5(1)=b;
    z7(1)=m;
    
    t=start;
    h=acc;
    
    %Set up Series Solution
    
    intz1=t+((1/24)*(-a^2*b^2*d2*epsilon-a^2*b^2*d3*epsilon+2*a^2*b^2*epsilon*m+2*C*a^2*b^2+2*a^2*b^2*epsilon-2*a^2*d3^2-2*b^2*d2^2+2*a^2*d3+2*b^2*d2)/(b^2*a^2))*t^3+((1/1920)*(a^4*b^4*d2^2*epsilon^2+5*a^4*b^4*d2*d3*epsilon^2-10*a^4*b^4*d2*epsilon^2*m+4*a^4*b^4*d3^2*epsilon^2-16*a^4*b^4*d3*epsilon^2*m+16*a^4*b^4*epsilon^2*m^2-10*C*a^4*b^4*d2*epsilon-16*C*a^4*b^4*d3*epsilon+32*C*a^4*b^4*epsilon*m-a^4*b^4*d2*epsilon^2-7*a^4*b^4*d3*epsilon^2+20*a^4*b^4*epsilon^2*m+16*C^2*a^4*b^4+20*C*a^4*b^4*epsilon+4*a^4*b^4*epsilon^2+10*a^4*b^2*d2*d3^2*epsilon+16*a^4*b^2*d3^3*epsilon-32*a^4*b^2*d3^2*epsilon*m+4*a^2*b^4*d2^3*epsilon+10*a^2*b^4*d2^2*d3*epsilon-20*a^2*b^4*d2^2*epsilon*m-32*C*a^4*b^2*d3^2-20*C*a^2*b^4*d2^2-10*a^4*b^2*d2*d3*epsilon-24*a^4*b^2*d3^2*epsilon+32*a^4*b^2*d3*epsilon*m-10*a^2*b^4*d2*d3*epsilon+20*a^2*b^4*d2*epsilon*m+32*C*a^4*b^2*d3+20*C*a^2*b^4*d2+8*a^4*b^2*d3*epsilon+16*a^4*d3^4-4*a^2*b^4*d2*epsilon+20*a^2*b^2*d2^2*d3^2+4*b^4*d2^4-20*a^4*d3^3-20*a^2*b^2*d2^2*d3-20*a^2*b^2*d2*d3^2+4*b^4*d2^3-8*a^4*d3^2+20*a^2*b^2*d2*d3-20*b^4*d2^2+12*a^4*d3+12*b^4*d2)/(b^4*a^4))*t^5;
    intz2=1+3*((1/24)*(-a^2*b^2*d2*epsilon-a^2*b^2*d3*epsilon+2*a^2*b^2*epsilon*m+2*C*a^2*b^2+2*a^2*b^2*epsilon-2*a^2*d3^2-2*b^2*d2^2+2*a^2*d3+2*b^2*d2)/(b^2*a^2))*t^2+5*((1/1920)*(a^4*b^4*d2^2*epsilon^2+5*a^4*b^4*d2*d3*epsilon^2-10*a^4*b^4*d2*epsilon^2*m+4*a^4*b^4*d3^2*epsilon^2-16*a^4*b^4*d3*epsilon^2*m+16*a^4*b^4*epsilon^2*m^2-10*C*a^4*b^4*d2*epsilon-16*C*a^4*b^4*d3*epsilon+32*C*a^4*b^4*epsilon*m-a^4*b^4*d2*epsilon^2-7*a^4*b^4*d3*epsilon^2+20*a^4*b^4*epsilon^2*m+16*C^2*a^4*b^4+20*C*a^4*b^4*epsilon+4*a^4*b^4*epsilon^2+10*a^4*b^2*d2*d3^2*epsilon+16*a^4*b^2*d3^3*epsilon-32*a^4*b^2*d3^2*epsilon*m+4*a^2*b^4*d2^3*epsilon+10*a^2*b^4*d2^2*d3*epsilon-20*a^2*b^4*d2^2*epsilon*m-32*C*a^4*b^2*d3^2-20*C*a^2*b^4*d2^2-10*a^4*b^2*d2*d3*epsilon-24*a^4*b^2*d3^2*epsilon+32*a^4*b^2*d3*epsilon*m-10*a^2*b^4*d2*d3*epsilon+20*a^2*b^4*d2*epsilon*m+32*C*a^4*b^2*d3+20*C*a^2*b^4*d2+8*a^4*b^2*d3*epsilon+16*a^4*d3^4-4*a^2*b^4*d2*epsilon+20*a^2*b^2*d2^2*d3^2+4*b^4*d2^4-20*a^4*d3^3-20*a^2*b^2*d2^2*d3-20*a^2*b^2*d2*d3^2+4*b^4*d2^3-8*a^4*d3^2+20*a^2*b^2*d2*d3-20*b^4*d2^2+12*a^4*d3+12*b^4*d2)/(b^4*a^4))*t^4;
    intz3=a+((1/8)*(a^2*epsilon+2*d2-2)/a)*t^2+((1/384)*(-a^4*b^2*d2*epsilon^2-a^4*b^2*d3*epsilon^2+2*a^4*b^2*epsilon^2*m+2*C*a^4*b^2*epsilon+2*a^4*b^2*epsilon^2-2*a^4*d3^2*epsilon-4*a^2*b^2*d2^2*epsilon-2*a^2*b^2*d2*d3*epsilon+4*a^2*b^2*d2*epsilon*m+4*C*a^2*b^2*d2+2*a^4*d3*epsilon+8*a^2*b^2*d2*epsilon+2*a^2*b^2*d3*epsilon-4*a^2*b^2*epsilon*m-4*C*a^2*b^2-4*a^2*b^2*epsilon-4*a^2*d2*d3^2-4*b^2*d2^3+4*a^2*d2*d3+4*a^2*d3^2+8*b^2*d2^2-4*a^2*d3-4*b^2*d2)/(b^2*a^3))*t^4;
    intz4=2*((1/8)*(a^2*epsilon+2*d2-2)/a)*t+4*((1/384)*(-a^4*b^2*d2*epsilon^2-a^4*b^2*d3*epsilon^2+2*a^4*b^2*epsilon^2*m+2*C*a^4*b^2*epsilon+2*a^4*b^2*epsilon^2-2*a^4*d3^2*epsilon-4*a^2*b^2*d2^2*epsilon-2*a^2*b^2*d2*d3*epsilon+4*a^2*b^2*d2*epsilon*m+4*C*a^2*b^2*d2+2*a^4*d3*epsilon+8*a^2*b^2*d2*epsilon+2*a^2*b^2*d3*epsilon-4*a^2*b^2*epsilon*m-4*C*a^2*b^2-4*a^2*b^2*epsilon-4*a^2*d2*d3^2-4*b^2*d2^3+4*a^2*d2*d3+4*a^2*d3^2+8*b^2*d2^2-4*a^2*d3-4*b^2*d2)/(b^2*a^3))*t^3;
    intz5=b+((1/8)*(b^2*epsilon+2*d3-2)/b)*t^2+((1/384)*(-a^2*b^4*d2*epsilon^2-a^2*b^4*d3*epsilon^2+2*a^2*b^4*epsilon^2*m+2*C*a^2*b^4*epsilon+2*a^2*b^4*epsilon^2-2*a^2*b^2*d2*d3*epsilon-4*a^2*b^2*d3^2*epsilon+4*a^2*b^2*d3*epsilon*m-2*b^4*d2^2*epsilon+4*C*a^2*b^2*d3+2*a^2*b^2*d2*epsilon+8*a^2*b^2*d3*epsilon-4*a^2*b^2*epsilon*m+2*b^4*d2*epsilon-4*C*a^2*b^2-4*a^2*b^2*epsilon-4*a^2*d3^3-4*b^2*d2^2*d3+8*a^2*d3^2+4*b^2*d2^2+4*b^2*d2*d3-4*a^2*d3-4*b^2*d2)/(b^3*a^2))*t^4;
    intz6=2*((1/8)*(b^2*epsilon+2*d3-2)/b)*t+4*((1/384)*(-a^2*b^4*d2*epsilon^2-a^2*b^4*d3*epsilon^2+2*a^2*b^4*epsilon^2*m+2*C*a^2*b^4*epsilon+2*a^2*b^4*epsilon^2-2*a^2*b^2*d2*d3*epsilon-4*a^2*b^2*d3^2*epsilon+4*a^2*b^2*d3*epsilon*m-2*b^4*d2^2*epsilon+4*C*a^2*b^2*d3+2*a^2*b^2*d2*epsilon+8*a^2*b^2*d3*epsilon-4*a^2*b^2*epsilon*m+2*b^4*d2*epsilon-4*C*a^2*b^2-4*a^2*b^2*epsilon-4*a^2*d3^3-4*b^2*d2^2*d3+8*a^2*d3^2+4*b^2*d2^2+4*b^2*d2*d3-4*a^2*d3-4*b^2*d2)/(b^3*a^2))*t^3;
    intz7=m+((1/4)*epsilon*m+(1/4)*C)*t^2+((1/192)*(-a^2*b^2*d2*epsilon^2*m-a^2*b^2*d3*epsilon^2*m+2*a^2*b^2*epsilon^2*m^2-C*a^2*b^2*d2*epsilon-C*a^2*b^2*d3*epsilon+4*C*a^2*b^2*epsilon*m+2*a^2*b^2*epsilon^2*m+2*C^2*a^2*b^2+2*C*a^2*b^2*epsilon-2*a^2*d3^2*epsilon*m-2*b^2*d2^2*epsilon*m-2*C*a^2*d3^2-2*C*b^2*d2^2+2*a^2*d3*epsilon*m+2*b^2*d2*epsilon*m+2*C*a^2*d3+2*C*b^2*d2)/(b^2*a^2))*t^4;
    intz8=2*((1/4)*epsilon*m+(1/4)*C)*t+4*((1/192)*(-a^2*b^2*d2*epsilon^2*m-a^2*b^2*d3*epsilon^2*m+2*a^2*b^2*epsilon^2*m^2-C*a^2*b^2*d2*epsilon-C*a^2*b^2*d3*epsilon+4*C*a^2*b^2*epsilon*m+2*a^2*b^2*epsilon^2*m+2*C^2*a^2*b^2+2*C*a^2*b^2*epsilon-2*a^2*d3^2*epsilon*m-2*b^2*d2^2*epsilon*m-2*C*a^2*d3^2-2*C*b^2*d2^2+2*a^2*d3*epsilon*m+2*b^2*d2*epsilon*m+2*C*a^2*d3+2*C*b^2*d2)/(b^2*a^2))*t^3;
    
    
    tt(2)=start;
    z1(2)=intz1;
    z3(2)=intz3;
    z5(2)=intz5;
    z7(2)=intz7;
    y=[intz1,intz2,intz3,intz4,intz5,intz6,intz7,intz8];
    ynew=y;
    i=3;
    
    %Apply Runge-Kutta
    
    while ((ynew(1)>0)&(ynew(3)>0))
        
        ytemp=y;
        K1=func2(ytemp,d2,d3,epsilon,C);
        ytemp=y+0.5*h*K1;
        K2=func2(ytemp,d2,d3,epsilon,C);
        ytemp=y+0.5*K2*h;
        K3=func2(ytemp,d2,d3,epsilon,C);
        ytemp=y+K3*h;
        K4=func2(ytemp,d2,d3,epsilon,C);
        ynew=y+(1/6)*h*(K1+2*K2+2*K3+K4);
        y=ynew;
        
       
         
          
          
         t=t+h
          
          if (ynew(5)<=0)
             
              break
          end
          
          if (ynew(1)<=0)
              break
          end
          
          if (ynew(3)<=0)
              break
          end
          
          if (t>limit);
              
            
              
              break
          end
          rest1=ynew(1)
          
          
       
      
       
      z1(i)=ynew(1);
          z3(i)=ynew(3);
          z5(i)=ynew(5);
          z7(i)=ynew(7);
          tt(i)=t;
          
          i=i+1;
           
      
      
          
         
    end
  
     
    
    %Plot The Graph
    plot(tt,z1,'blue',tt,z3,'red',tt,z5,'green',tt,z7,'m')
    text(tt(i-1000),ynew(1),'g1')
    text(tt(i-1000),ynew(3),'g2')
    text(tt(i-1000),ynew(5),'g3')
    text(tt(i-1000),ynew(7),'U')
    title(str)
    
end

function [ yp ] = func2( y,d2,d3,epsilon,C)
%Calculates the values of the equations for the 3 equations and potential u
yp(1)=y(2);
yp(2)=-d2*y(2)*y(4)/y(3)-d3*y(2)*y(6)/y(5)+y(8)*y(2)+(epsilon/2)*y(2);
yp(3)=y(4);
yp(4)=(1-d2)*y(4)^2/y(3)-y(2)*y(4)/y(1)-d3*y(6)*y(4)/y(5)+y(8)*y(4)+(d2-1)/y(3)+(epsilon/2)*y(3);
yp(5)=y(6);
yp(6)=(1-d3)*y(6)^2/y(5)-y(2)*y(6)/y(1)-d2*y(4)*y(6)/y(3)+y(8)*y(6)+(d3-1)/y(5)+(epsilon/2)*y(5);
yp(7)=y(8);
yp(8)=-y(8)*(y(2)/y(1)+d2*y(4)/y(3)+d3*y(6)/y(5))+y(8)^2+epsilon*y(7)+C;

end
     