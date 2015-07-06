
     function [] = experiment1(d1,d2,A2,A3,intz3,intz5,C,epsilon,limit,start,acc,type)

     
     %Coded by Michael Gallaugher for Dr. M. Wang NSERC Summer 2013
     %Runge Kutta algorithm for 2 functions (f and h) and potential u
     
     %d1, d2, A2, A3, C and epsilon are the constants in the equations
     %intz3 = the initial condition of f2 or h
     %intz5 = the initial condition of u
     %limit = the cut off point for t
     %start = the point at which to evaluate the series solution
     %acc = the desired accuracy of the Runge Kutta algorithm
     %note* func 3 is the basically a function that evaluates the
     %differential equations at the Runge Kutta Stages


    
    
        
     
     
     
    
    
    %Set the initial values
    
    tt(1)=0;
    z1(1)=0;
    z3(1)=intz3;
    z5(1)=intz5;
    
    
    %Set up and evaluate the series solution at the specified start time
    %The series solution was found using maple. 
    
    t=start;
    intz1=t+(-(1/12)*(-d1*epsilon*intz3^2+d2*epsilon*intz3^2-2*epsilon*intz3^2*intz5-2*C*intz3^2-epsilon*intz3^2+2*A2)/(intz3^2*(d1+1)*d1))*t^3+((1/480)*(d1^3*d2*epsilon^2*intz3^4-2*d1^2*d2^2*epsilon^2*intz3^4+28*d1^2*d2*epsilon^2*intz3^4*intz5+13*d1*d2^3*epsilon^2*intz3^4-52*d1*d2^2*epsilon^2*intz3^4*intz5+52*d1*d2*epsilon^2*intz3^4*intz5^2+28*C*d1^2*d2*epsilon*intz3^4-52*C*d1*d2^2*epsilon*intz3^4+104*C*d1*d2*epsilon*intz3^4*intz5+5*d1^2*d2*epsilon^2*intz3^4-20*d1*d2^2*epsilon^2*intz3^4+40*d1*d2*epsilon^2*intz3^4*intz5+3*d2^3*epsilon^2*intz3^4-12*d2^2*epsilon^2*intz3^4*intz5+12*d2*epsilon^2*intz3^4*intz5^2+52*C^2*d1*d2*intz3^4+40*C*d1*d2*epsilon*intz3^4-12*C*d2^2*epsilon*intz3^4+24*C*d2*epsilon*intz3^4*intz5+7*d1*d2*epsilon^2*intz3^4-6*d2^2*epsilon^2*intz3^4+12*d2*epsilon^2*intz3^4*intz5+20*A2*d1^2*d2*epsilon*intz3^2+52*A2*d1*d2^2*epsilon*intz3^2-104*A2*d1*d2*epsilon*intz3^2*intz5+12*C^2*d2*intz3^4+12*C*d2*epsilon*intz3^4+3*d2*epsilon^2*intz3^4-104*A2*C*d1*d2*intz3^2-40*A2*d1*d2*epsilon*intz3^2+12*A2*d2^2*epsilon*intz3^2-24*A2*d2*epsilon*intz3^2*intz5-24*A2*C*d2*intz3^2-12*A2*d2*epsilon*intz3^2+144*A3*d1^3*d2+48*A2^2*d1^2+52*A2^2*d1*d2+288*A3*d1^2*d2+12*A2^2*d2+144*A3*d1*d2)/(intz3^4*(d1+1)^2*(d1+3)*d1^2*d2))*t^5;
    intz2=1+3*(-(1/12)*(-d1*epsilon*intz3^2+d2*epsilon*intz3^2-2*epsilon*intz3^2*intz5-2*C*intz3^2-epsilon*intz3^2+2*A2)/(intz3^2*(d1+1)*d1))*t^2+5*((1/480)*(d1^3*d2*epsilon^2*intz3^4-2*d1^2*d2^2*epsilon^2*intz3^4+28*d1^2*d2*epsilon^2*intz3^4*intz5+13*d1*d2^3*epsilon^2*intz3^4-52*d1*d2^2*epsilon^2*intz3^4*intz5+52*d1*d2*epsilon^2*intz3^4*intz5^2+28*C*d1^2*d2*epsilon*intz3^4-52*C*d1*d2^2*epsilon*intz3^4+104*C*d1*d2*epsilon*intz3^4*intz5+5*d1^2*d2*epsilon^2*intz3^4-20*d1*d2^2*epsilon^2*intz3^4+40*d1*d2*epsilon^2*intz3^4*intz5+3*d2^3*epsilon^2*intz3^4-12*d2^2*epsilon^2*intz3^4*intz5+12*d2*epsilon^2*intz3^4*intz5^2+52*C^2*d1*d2*intz3^4+40*C*d1*d2*epsilon*intz3^4-12*C*d2^2*epsilon*intz3^4+24*C*d2*epsilon*intz3^4*intz5+7*d1*d2*epsilon^2*intz3^4-6*d2^2*epsilon^2*intz3^4+12*d2*epsilon^2*intz3^4*intz5+20*A2*d1^2*d2*epsilon*intz3^2+52*A2*d1*d2^2*epsilon*intz3^2-104*A2*d1*d2*epsilon*intz3^2*intz5+12*C^2*d2*intz3^4+12*C*d2*epsilon*intz3^4+3*d2*epsilon^2*intz3^4-104*A2*C*d1*d2*intz3^2-40*A2*d1*d2*epsilon*intz3^2+12*A2*d2^2*epsilon*intz3^2-24*A2*d2*epsilon*intz3^2*intz5-24*A2*C*d2*intz3^2-12*A2*d2*epsilon*intz3^2+144*A3*d1^3*d2+48*A2^2*d1^2+52*A2^2*d1*d2+288*A3*d1^2*d2+12*A2^2*d2+144*A3*d1*d2)/(intz3^4*(d1+1)^2*(d1+3)*d1^2*d2))*t^4;
    intz3=intz3+((1/4)*(d2*epsilon*intz3^2+2*A2)/(d2*intz3*(d1+1)))*t^2+(-(1/96)*(-d1*d2^2*epsilon^2*intz3^4+4*d2^3*epsilon^2*intz3^4-8*d2^2*epsilon^2*intz3^4*intz5-8*C*d2^2*epsilon*intz3^4-7*d2^2*epsilon^2*intz3^4+4*A2*d1*d2*epsilon*intz3^2+16*A2*d2^2*epsilon*intz3^2-16*A2*d2*epsilon*intz3^2*intz5-16*A2*C*d2*intz3^2-20*A2*d2*epsilon*intz3^2+48*A3*d1^2*d2+12*A2^2*d1+16*A2^2*d2+96*A3*d1*d2-12*A2^2+48*A3*d2)/(intz3^3*(d1+1)^2*d2^2*(d1+3)))*t^4+((1/5760)*(d1^3*d2^3*epsilon^3*intz3^6+28*d1^2*d2^4*epsilon^3*intz3^6+88*d1^2*d2^3*epsilon^3*intz3^6*intz5+88*d1*d2^5*epsilon^3*intz3^6-352*d1*d2^4*epsilon^3*intz3^6*intz5+352*d1*d2^3*epsilon^3*intz3^6*intz5^2+88*C*d1^2*d2^3*epsilon^2*intz3^6-352*C*d1*d2^4*epsilon^2*intz3^6+704*C*d1*d2^3*epsilon^2*intz3^6*intz5-20*d1^2*d2^3*epsilon^3*intz3^6-180*d1*d2^4*epsilon^3*intz3^6+600*d1*d2^3*epsilon^3*intz3^6*intz5+8*d2^5*epsilon^3*intz3^6-32*d2^4*epsilon^3*intz3^6*intz5+32*d2^3*epsilon^3*intz3^6*intz5^2+352*C^2*d1*d2^3*epsilon*intz3^6+600*C*d1*d2^3*epsilon^2*intz3^6-32*C*d2^4*epsilon^2*intz3^6+64*C*d2^3*epsilon^2*intz3^6*intz5+107*d1*d2^3*epsilon^3*intz3^6-16*d2^4*epsilon^3*intz3^6+32*d2^3*epsilon^3*intz3^6*intz5+62*A2*d1^3*d2^2*epsilon^2*intz3^4+376*A2*d1^2*d2^3*epsilon^2*intz3^4-64*A2*d1^2*d2^2*epsilon^2*intz3^4*intz5+528*A2*d1*d2^4*epsilon^2*intz3^4-1408*A2*d1*d2^3*epsilon^2*intz3^4*intz5+704*A2*d1*d2^2*epsilon^2*intz3^4*intz5^2+32*C^2*d2^3*epsilon*intz3^6+32*C*d2^3*epsilon^2*intz3^6+8*d2^3*epsilon^3*intz3^6-64*A2*C*d1^2*d2^2*epsilon*intz3^4-1408*A2*C*d1*d2^3*epsilon*intz3^4+1408*A2*C*d1*d2^2*epsilon*intz3^4*intz5-280*A2*d1^2*d2^2*epsilon^2*intz3^4-840*A2*d1*d2^3*epsilon^2*intz3^4+1920*A2*d1*d2^2*epsilon^2*intz3^4*intz5+48*A2*d2^4*epsilon^2*intz3^4-128*A2*d2^3*epsilon^2*intz3^4*intz5+64*A2*d2^2*epsilon^2*intz3^4*intz5^2+704*A2*C^2*d1*d2^2*intz3^4+1920*A2*C*d1*d2^2*epsilon*intz3^4-128*A2*C*d2^3*epsilon*intz3^4+128*A2*C*d2^2*epsilon*intz3^4*intz5+394*A2*d1*d2^2*epsilon^2*intz3^4-64*A2*d2^3*epsilon^2*intz3^4+64*A2*d2^2*epsilon^2*intz3^4*intz5+1200*A3*d1^4*d2^2*epsilon*intz3^2+1344*A3*d1^3*d2^3*epsilon*intz3^2-1920*A3*d1^3*d2^2*epsilon*intz3^2*intz5+300*A2^2*d1^3*d2*epsilon*intz3^2+1168*A2^2*d1^2*d2^2*epsilon*intz3^2-480*A2^2*d1^2*d2*epsilon*intz3^2*intz5+1056*A2^2*d1*d2^3*epsilon*intz3^2-1408*A2^2*d1*d2^2*epsilon*intz3^2*intz5+64*A2*C^2*d2^2*intz3^4+64*A2*C*d2^2*epsilon*intz3^4+16*A2*d2^2*epsilon^2*intz3^4-1920*A3*C*d1^3*d2^2*intz3^2+3600*A3*d1^3*d2^2*epsilon*intz3^2+3648*A3*d1^2*d2^3*epsilon*intz3^2-5760*A3*d1^2*d2^2*epsilon*intz3^2*intz5-480*A2^2*C*d1^2*d2*intz3^2-1408*A2^2*C*d1*d2^2*intz3^2-720*A2^2*d1^2*d2*epsilon*intz3^2-1200*A2^2*d1*d2^2*epsilon*intz3^2+1440*A2^2*d1*d2*epsilon*intz3^2*intz5+96*A2^2*d2^3*epsilon*intz3^2-128*A2^2*d2^2*epsilon*intz3^2*intz5-5760*A3*C*d1^2*d2^2*intz3^2+2640*A3*d1^2*d2^2*epsilon*intz3^2+3264*A3*d1*d2^3*epsilon*intz3^2-5760*A3*d1*d2^2*epsilon*intz3^2*intz5+1440*A2^2*C*d1*d2*intz3^2-128*A2^2*C*d2^2*intz3^2+420*A2^2*d1*d2*epsilon*intz3^2-64*A2^2*d2^2*epsilon*intz3^2+3360*A2*A3*d1^4*d2+2688*A2*A3*d1^3*d2^2-5760*A3*C*d1*d2^2*intz3^2-720*A3*d1*d2^2*epsilon*intz3^2+960*A3*d2^3*epsilon*intz3^2-1920*A3*d2^2*epsilon*intz3^2*intz5+360*A2^3*d1^3+1056*A2^3*d1^2*d2+704*A2^3*d1*d2^2+12000*A2*A3*d1^3*d2+7296*A2*A3*d1^2*d2^2-1920*A3*C*d2^2*intz3^2-960*A3*d2^2*epsilon*intz3^2-480*A2^3*d1^2-480*A2^3*d1*d2+64*A2^3*d2^2+13920*A2*A3*d1^2*d2+6528*A2*A3*d1*d2^2+120*A2^3*d1+5280*A2*A3*d1*d2+1920*A2*A3*d2^2)/(d1*intz3^5*(d1+1)^3*d2^3*(d1+3)*(d1+5)))*t^6;
    intz4=2*((1/4)*(d2*epsilon*intz3^2+2*A2)/(d2*intz3*(d1+1)))*t+4*(-(1/96)*(-d1*d2^2*epsilon^2*intz3^4+4*d2^3*epsilon^2*intz3^4-8*d2^2*epsilon^2*intz3^4*intz5-8*C*d2^2*epsilon*intz3^4-7*d2^2*epsilon^2*intz3^4+4*A2*d1*d2*epsilon*intz3^2+16*A2*d2^2*epsilon*intz3^2-16*A2*d2*epsilon*intz3^2*intz5-16*A2*C*d2*intz3^2-20*A2*d2*epsilon*intz3^2+48*A3*d1^2*d2+12*A2^2*d1+16*A2^2*d2+96*A3*d1*d2-12*A2^2+48*A3*d2)/(intz3^3*(d1+1)^2*d2^2*(d1+3)))*t^3+6*((1/5760)*(d1^3*d2^3*epsilon^3*intz3^6+28*d1^2*d2^4*epsilon^3*intz3^6+88*d1^2*d2^3*epsilon^3*intz3^6*intz5+88*d1*d2^5*epsilon^3*intz3^6-352*d1*d2^4*epsilon^3*intz3^6*intz5+352*d1*d2^3*epsilon^3*intz3^6*intz5^2+88*C*d1^2*d2^3*epsilon^2*intz3^6-352*C*d1*d2^4*epsilon^2*intz3^6+704*C*d1*d2^3*epsilon^2*intz3^6*intz5-20*d1^2*d2^3*epsilon^3*intz3^6-180*d1*d2^4*epsilon^3*intz3^6+600*d1*d2^3*epsilon^3*intz3^6*intz5+8*d2^5*epsilon^3*intz3^6-32*d2^4*epsilon^3*intz3^6*intz5+32*d2^3*epsilon^3*intz3^6*intz5^2+352*C^2*d1*d2^3*epsilon*intz3^6+600*C*d1*d2^3*epsilon^2*intz3^6-32*C*d2^4*epsilon^2*intz3^6+64*C*d2^3*epsilon^2*intz3^6*intz5+107*d1*d2^3*epsilon^3*intz3^6-16*d2^4*epsilon^3*intz3^6+32*d2^3*epsilon^3*intz3^6*intz5+62*A2*d1^3*d2^2*epsilon^2*intz3^4+376*A2*d1^2*d2^3*epsilon^2*intz3^4-64*A2*d1^2*d2^2*epsilon^2*intz3^4*intz5+528*A2*d1*d2^4*epsilon^2*intz3^4-1408*A2*d1*d2^3*epsilon^2*intz3^4*intz5+704*A2*d1*d2^2*epsilon^2*intz3^4*intz5^2+32*C^2*d2^3*epsilon*intz3^6+32*C*d2^3*epsilon^2*intz3^6+8*d2^3*epsilon^3*intz3^6-64*A2*C*d1^2*d2^2*epsilon*intz3^4-1408*A2*C*d1*d2^3*epsilon*intz3^4+1408*A2*C*d1*d2^2*epsilon*intz3^4*intz5-280*A2*d1^2*d2^2*epsilon^2*intz3^4-840*A2*d1*d2^3*epsilon^2*intz3^4+1920*A2*d1*d2^2*epsilon^2*intz3^4*intz5+48*A2*d2^4*epsilon^2*intz3^4-128*A2*d2^3*epsilon^2*intz3^4*intz5+64*A2*d2^2*epsilon^2*intz3^4*intz5^2+704*A2*C^2*d1*d2^2*intz3^4+1920*A2*C*d1*d2^2*epsilon*intz3^4-128*A2*C*d2^3*epsilon*intz3^4+128*A2*C*d2^2*epsilon*intz3^4*intz5+394*A2*d1*d2^2*epsilon^2*intz3^4-64*A2*d2^3*epsilon^2*intz3^4+64*A2*d2^2*epsilon^2*intz3^4*intz5+1200*A3*d1^4*d2^2*epsilon*intz3^2+1344*A3*d1^3*d2^3*epsilon*intz3^2-1920*A3*d1^3*d2^2*epsilon*intz3^2*intz5+300*A2^2*d1^3*d2*epsilon*intz3^2+1168*A2^2*d1^2*d2^2*epsilon*intz3^2-480*A2^2*d1^2*d2*epsilon*intz3^2*intz5+1056*A2^2*d1*d2^3*epsilon*intz3^2-1408*A2^2*d1*d2^2*epsilon*intz3^2*intz5+64*A2*C^2*d2^2*intz3^4+64*A2*C*d2^2*epsilon*intz3^4+16*A2*d2^2*epsilon^2*intz3^4-1920*A3*C*d1^3*d2^2*intz3^2+3600*A3*d1^3*d2^2*epsilon*intz3^2+3648*A3*d1^2*d2^3*epsilon*intz3^2-5760*A3*d1^2*d2^2*epsilon*intz3^2*intz5-480*A2^2*C*d1^2*d2*intz3^2-1408*A2^2*C*d1*d2^2*intz3^2-720*A2^2*d1^2*d2*epsilon*intz3^2-1200*A2^2*d1*d2^2*epsilon*intz3^2+1440*A2^2*d1*d2*epsilon*intz3^2*intz5+96*A2^2*d2^3*epsilon*intz3^2-128*A2^2*d2^2*epsilon*intz3^2*intz5-5760*A3*C*d1^2*d2^2*intz3^2+2640*A3*d1^2*d2^2*epsilon*intz3^2+3264*A3*d1*d2^3*epsilon*intz3^2-5760*A3*d1*d2^2*epsilon*intz3^2*intz5+1440*A2^2*C*d1*d2*intz3^2-128*A2^2*C*d2^2*intz3^2+420*A2^2*d1*d2*epsilon*intz3^2-64*A2^2*d2^2*epsilon*intz3^2+3360*A2*A3*d1^4*d2+2688*A2*A3*d1^3*d2^2-5760*A3*C*d1*d2^2*intz3^2-720*A3*d1*d2^2*epsilon*intz3^2+960*A3*d2^3*epsilon*intz3^2-1920*A3*d2^2*epsilon*intz3^2*intz5+360*A2^3*d1^3+1056*A2^3*d1^2*d2+704*A2^3*d1*d2^2+12000*A2*A3*d1^3*d2+7296*A2*A3*d1^2*d2^2-1920*A3*C*d2^2*intz3^2-960*A3*d2^2*epsilon*intz3^2-480*A2^3*d1^2-480*A2^3*d1*d2+64*A2^3*d2^2+13920*A2*A3*d1^2*d2+6528*A2*A3*d1*d2^2+120*A2^3*d1+5280*A2*A3*d1*d2+1920*A2*A3*d2^2)/(d1*intz3^5*(d1+1)^3*d2^3*(d1+3)*(d1+5)))*t^5;
    intz5=intz5+((1/2)*(epsilon*intz5+C)/(d1+1))*t^2+(-(1/12)*(-d1*epsilon^2*intz3^2*intz5+d2*epsilon^2*intz3^2*intz5-2*epsilon^2*intz3^2*intz5^2-C*d1*epsilon*intz3^2+C*d2*epsilon*intz3^2-4*C*epsilon*intz3^2*intz5-epsilon^2*intz3^2*intz5-2*C^2*intz3^2-C*epsilon*intz3^2+2*A2*epsilon*intz5+2*A2*C)/(intz3^2*(d1+1)^2*(d1+3)))*t^4+((1/360)*(2*d1^3*d2*epsilon^3*intz3^4*intz5-4*d1^2*d2^2*epsilon^3*intz3^4*intz5+26*d1^2*d2*epsilon^3*intz3^4*intz5^2+11*d1*d2^3*epsilon^3*intz3^4*intz5-44*d1*d2^2*epsilon^3*intz3^4*intz5^2+44*d1*d2*epsilon^3*intz3^4*intz5^3+2*C*d1^3*d2*epsilon^2*intz3^4-4*C*d1^2*d2^2*epsilon^2*intz3^4+52*C*d1^2*d2*epsilon^2*intz3^4*intz5+11*C*d1*d2^3*epsilon^2*intz3^4-88*C*d1*d2^2*epsilon^2*intz3^4*intz5+132*C*d1*d2*epsilon^2*intz3^4*intz5^2+5*d1^2*d2*epsilon^3*intz3^4*intz5+30*d1*d2*epsilon^3*intz3^4*intz5^2+d2^3*epsilon^3*intz3^4*intz5-4*d2^2*epsilon^3*intz3^4*intz5^2+4*d2*epsilon^3*intz3^4*intz5^3+26*C^2*d1^2*d2*epsilon*intz3^4-44*C^2*d1*d2^2*epsilon*intz3^4+132*C^2*d1*d2*epsilon*intz3^4*intz5+5*C*d1^2*d2*epsilon^2*intz3^4+60*C*d1*d2*epsilon^2*intz3^4*intz5+C*d2^3*epsilon^2*intz3^4-8*C*d2^2*epsilon^2*intz3^4*intz5+12*C*d2*epsilon^2*intz3^4*intz5^2+4*d1*d2*epsilon^3*intz3^4*intz5-2*d2^2*epsilon^3*intz3^4*intz5+4*d2*epsilon^3*intz3^4*intz5^2+10*A2*d1^2*d2*epsilon^2*intz3^2*intz5+44*A2*d1*d2^2*epsilon^2*intz3^2*intz5-88*A2*d1*d2*epsilon^2*intz3^2*intz5^2+44*C^3*d1*d2*intz3^4+30*C^2*d1*d2*epsilon*intz3^4-4*C^2*d2^2*epsilon*intz3^4+12*C^2*d2*epsilon*intz3^4*intz5+4*C*d1*d2*epsilon^2*intz3^4-2*C*d2^2*epsilon^2*intz3^4+8*C*d2*epsilon^2*intz3^4*intz5+d2*epsilon^3*intz3^4*intz5+10*A2*C*d1^2*d2*epsilon*intz3^2+44*A2*C*d1*d2^2*epsilon*intz3^2-176*A2*C*d1*d2*epsilon*intz3^2*intz5+30*A2*d1*d2*epsilon^2*intz3^2*intz5+4*A2*d2^2*epsilon^2*intz3^2*intz5-8*A2*d2*epsilon^2*intz3^2*intz5^2+4*C^3*d2*intz3^4+4*C^2*d2*epsilon*intz3^4+C*d2*epsilon^2*intz3^4-88*A2*C^2*d1*d2*intz3^2+30*A2*C*d1*d2*epsilon*intz3^2+4*A2*C*d2^2*epsilon*intz3^2-16*A2*C*d2*epsilon*intz3^2*intz5-4*A2*d2*epsilon^2*intz3^2*intz5+48*A3*d1^3*d2*epsilon*intz5+36*A2^2*d1^2*epsilon*intz5+44*A2^2*d1*d2*epsilon*intz5-8*A2*C^2*d2*intz3^2-4*A2*C*d2*epsilon*intz3^2+48*A3*C*d1^3*d2+96*A3*d1^2*d2*epsilon*intz5+36*A2^2*C*d1^2+44*A2^2*C*d1*d2+60*A2^2*d1*epsilon*intz5+4*A2^2*d2*epsilon*intz5+96*A3*C*d1^2*d2+48*A3*d1*d2*epsilon*intz5+60*A2^2*C*d1+4*A2^2*C*d2+48*A3*C*d1*d2)/(d1*intz3^4*(d1+1)^3*(d1+3)*d2*(d1+5)))*t^6;
    intz6=2*((1/2)*(epsilon*intz5+C)/(d1+1))*t+4*(-(1/12)*(-d1*epsilon^2*intz3^2*intz5+d2*epsilon^2*intz3^2*intz5-2*epsilon^2*intz3^2*intz5^2-C*d1*epsilon*intz3^2+C*d2*epsilon*intz3^2-4*C*epsilon*intz3^2*intz5-epsilon^2*intz3^2*intz5-2*C^2*intz3^2-C*epsilon*intz3^2+2*A2*epsilon*intz5+2*A2*C)/(intz3^2*(d1+1)^2*(d1+3)))*t^3+6*((1/360)*(2*d1^3*d2*epsilon^3*intz3^4*intz5-4*d1^2*d2^2*epsilon^3*intz3^4*intz5+26*d1^2*d2*epsilon^3*intz3^4*intz5^2+11*d1*d2^3*epsilon^3*intz3^4*intz5-44*d1*d2^2*epsilon^3*intz3^4*intz5^2+44*d1*d2*epsilon^3*intz3^4*intz5^3+2*C*d1^3*d2*epsilon^2*intz3^4-4*C*d1^2*d2^2*epsilon^2*intz3^4+52*C*d1^2*d2*epsilon^2*intz3^4*intz5+11*C*d1*d2^3*epsilon^2*intz3^4-88*C*d1*d2^2*epsilon^2*intz3^4*intz5+132*C*d1*d2*epsilon^2*intz3^4*intz5^2+5*d1^2*d2*epsilon^3*intz3^4*intz5+30*d1*d2*epsilon^3*intz3^4*intz5^2+d2^3*epsilon^3*intz3^4*intz5-4*d2^2*epsilon^3*intz3^4*intz5^2+4*d2*epsilon^3*intz3^4*intz5^3+26*C^2*d1^2*d2*epsilon*intz3^4-44*C^2*d1*d2^2*epsilon*intz3^4+132*C^2*d1*d2*epsilon*intz3^4*intz5+5*C*d1^2*d2*epsilon^2*intz3^4+60*C*d1*d2*epsilon^2*intz3^4*intz5+C*d2^3*epsilon^2*intz3^4-8*C*d2^2*epsilon^2*intz3^4*intz5+12*C*d2*epsilon^2*intz3^4*intz5^2+4*d1*d2*epsilon^3*intz3^4*intz5-2*d2^2*epsilon^3*intz3^4*intz5+4*d2*epsilon^3*intz3^4*intz5^2+10*A2*d1^2*d2*epsilon^2*intz3^2*intz5+44*A2*d1*d2^2*epsilon^2*intz3^2*intz5-88*A2*d1*d2*epsilon^2*intz3^2*intz5^2+44*C^3*d1*d2*intz3^4+30*C^2*d1*d2*epsilon*intz3^4-4*C^2*d2^2*epsilon*intz3^4+12*C^2*d2*epsilon*intz3^4*intz5+4*C*d1*d2*epsilon^2*intz3^4-2*C*d2^2*epsilon^2*intz3^4+8*C*d2*epsilon^2*intz3^4*intz5+d2*epsilon^3*intz3^4*intz5+10*A2*C*d1^2*d2*epsilon*intz3^2+44*A2*C*d1*d2^2*epsilon*intz3^2-176*A2*C*d1*d2*epsilon*intz3^2*intz5+30*A2*d1*d2*epsilon^2*intz3^2*intz5+4*A2*d2^2*epsilon^2*intz3^2*intz5-8*A2*d2*epsilon^2*intz3^2*intz5^2+4*C^3*d2*intz3^4+4*C^2*d2*epsilon*intz3^4+C*d2*epsilon^2*intz3^4-88*A2*C^2*d1*d2*intz3^2+30*A2*C*d1*d2*epsilon*intz3^2+4*A2*C*d2^2*epsilon*intz3^2-16*A2*C*d2*epsilon*intz3^2*intz5-4*A2*d2*epsilon^2*intz3^2*intz5+48*A3*d1^3*d2*epsilon*intz5+36*A2^2*d1^2*epsilon*intz5+44*A2^2*d1*d2*epsilon*intz5-8*A2*C^2*d2*intz3^2-4*A2*C*d2*epsilon*intz3^2+48*A3*C*d1^3*d2+96*A3*d1^2*d2*epsilon*intz5+36*A2^2*C*d1^2+44*A2^2*C*d1*d2+60*A2^2*d1*epsilon*intz5+4*A2^2*d2*epsilon*intz5+96*A3*C*d1^2*d2+48*A3*d1*d2*epsilon*intz5+60*A2^2*C*d1+4*A2^2*C*d2+48*A3*C*d1*d2)/(d1*intz3^4*(d1+1)^3*(d1+3)*d2*(d1+5)))*t^5;
    if type==-1
        intz1=-0.7878157749e-3*t^7+0.5105707396e-2*t^5-0.3604430210e-1*t^3+t;
        intz2=-0.5514710424e-2*t^6+0.2552853698e-1*t^4-.1081329063*t^2+1;
        intz3=2+1/8*t^2-0.006159018880*t^4+0.0005976764885*t^6-0.00007145740183*t^8;
        intz4=(1/4)*t-0.2463607552e-1*t^3+0.3586058931e-2*t^5-0.5716592146e-3*t^7;
        intz5=-1.567468375+.3918670937*t^2-0.7062287949e-2*t^4+0.6669195741e-3*t^6-0.7717976962e-4*t^8;
        intz6=.7837341874*t-0.2824915180e-1*t^3+0.4001517445e-2*t^5-0.6174381570e-3*t^7;
    end
    %Set accuracy
    h=acc;
    
    %Set values before running the Runge Kutta algorithm
    tt(2)=t;
    z1(2)=intz1;
    z3(2)=intz3;
    z5(2)=intz5;
    y=[intz1,intz2,intz3,intz4,intz5,intz6];
    ynew=y;
    i=3;
    
    %Run the Runge Kutta algorithm
    while ((ynew(1)>0)&(ynew(3)>0))
        
        ytemp=y;
        K1=func3(ytemp,d1,d2,epsilon,A2,A3,C);
        ytemp=y+0.5*h*K1;
        K2=func3(ytemp,d1,d2,epsilon,A2,A3,C);
        ytemp=y+0.5*K2*h;
        K3=func3(ytemp,d1,d2,epsilon,A2,A3,C);
        ytemp=y+K3*h;
        K4=func3(ytemp,d1,d2,epsilon,A2,A3,C);
        ynew=y+(1/8)*h*(K1+3*K2+3*K3+K4);
        y=ynew;
        
       
         
          
          
         t=t+h
         
         %Test that the values still meet the requirements
          if (ynew(1)<=0)
             
              break
          end
          
          %Stop when t is greater than the limit specified
          if (t>limit)
              
            
              
              break
          end
          
         
       
      %Set the new values
          z1(i)=ynew(1);
          z3(i)=ynew(3);
          z5(i)=ynew(5);
          tt(i)=t;
          
          i=i+1;
           
      
      
          
         
    end
    
 %Plot the graphs
 str=sprintf('d1=%g d2=%g A2=%g A3=%g hbar=%g ubar=%g C=%g Epsilon=%g',d1,d2,A2,A3,z3(1),z5(1),C,epsilon);
    plot(tt,z1,'blue',tt,z3,'red',tt,z5,'green')
    text(tt(i-1000),ynew(1),'g1')
    text(tt(i-1000),ynew(3),'g2')
    text(tt(i-1000),ynew(5),'u')
   title(str)
    xlabel('t')
    
     end

     
     function [ yp ] = func3( y,d1,d2,epsilon,A2,A3,C0)
     
     


yp(1)=y(2);
yp(2)=(-(d1-1)*(y(2)^2))/(y(1))-(d2*y(2)*y(4))/(y(3))+(d1-1)/(y(1))+(((A3)/(d1))*y(1)^3)/(y(3)^4)+(y(2)*y(6))+(epsilon/2)*y(1);
yp(3)=y(4);
yp(4)=(-(d2-1)*(y(4)^(2)))/(y(3))-(d1*y(2)*y(4))/(y(1))+((A2/d2)/y(3))-((A3/d2)*2*y(1)^(2))/(y(3)^(3))+y(4)*y(6)+(epsilon/2)*y(3);
yp(5)=y(6);
yp(6)=-y(6)*((d1*y(2))/(y(1))+(d2*y(4))/(y(3)))+y(6)^(2)+(epsilon*y(5))+C0;

end

    
        