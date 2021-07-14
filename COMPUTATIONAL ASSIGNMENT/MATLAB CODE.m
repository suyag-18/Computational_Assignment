%Computational Assignment

%System Values %Napthalene
Tc = 784 ;
Pc = 40.0;
w = 0.302;

%EOS %Soave Redlich Kwong
Ru = 0.0821 ;
sigma=1;
epsilon=0;
si=0.42748;
ohm=0.08664;
a = (0.42748*(Ru^2)*(Tc^2))/Pc ;
b = (0.08664*Ru*Tc)/Pc ;

%taking approx guesses
%Then analyzing the isotherms to get a close values
%which would give real roots

 guesses=[38.0934,38.4976,38.8946,39.2964,39.6548,39.9999];
 
 %taking a vector to store all the saturation points
 SatP=zeros(1,12);
 SatV=zeros(1,12);
 
        %loop for getting saturation points and isotherm at different temp
        
        for T=779:1:784        
            i = T-778;
    	hold on

        
        %plotting isotherm
        
        alpha= ((1 + (0.48508 + 1.574*w - 0.176*(w^2))*(1 - (T/Tc)^0.5))^2) ;
        aalpha = a*((1 + (0.48508 + 1.574*w - 0.176*(w^2))*(1 - (T/Tc)^0.5))^2) ;
        
        syms v;
        title('PV isotherm for Napthalene',' ');
        fplot((((Ru*T)/(v-b)) - ((aalpha)/(v*(v+b)))) , [0 2]);
        
        xlim([0.24 1.15])
    	ylim([37.5 41])
        
        xlabel("Volume(L)");
        ylabel("Pressure(Atm)");
        
    %code for finding saturation points
    
Psat=guesses(i);
Psat=Psat+0.0001;
 %finding roots
 a0=1;
 a1=-(Ru*T/Psat);
 a2=-(Ru*(T/Psat)*b - (aalpha/Psat) +(b^2));
 a3= -aalpha*(b/Psat);
 
 k=[a0 a1 a2 a3];
 
 vol=roots(k);
 vol=sort(vol);
 
 Vl=vol(1);
 Vv=vol(3);
 
 if(imag(Vv))~=0 
     Vv=Vl;
 end
 

 
      %departure function calculation
      depl= depfl(Vl,Psat,T,aalpha,b,Ru);
      depg= depfg(Vv,Psat,T,aalpha,b,Ru);
 
      %getting better approx of Psat using departure function
       
        temp = abs(depl-depg);
      if(temp<1e-2)
          
          SatP(2*i-1)=Psat;
          SatP(2*i)=Psat;
          
      end
      
      
          
          Ptemp=Psat+0.01;
         
          a0=1;
            a1=-(Ru*T/Ptemp);
            a2=-(Ru*(T/Ptemp)*b - (aalpha/Ptemp) +(b^2));
            a3= -aalpha*(b/Ptemp);
 
            k=[a0 a1 a2 a3];
 
            voltemp=roots(k);
            voltemp=sort(voltemp);
 
            Vltemp=voltemp(1);
            Vgtemp=voltemp(3);
 
                if(imag(Vgtemp))~=0 
                        Vgtemp=Vltemp;
                end
          depltemp= depfl(Vltemp,Psat,T,aalpha,b,Ru);
          depgtemp= depfg(Vgtemp,Psat,T,aalpha,b,Ru);
         
          temp1=abs(depltemp-depgtemp);
          if(temp1<temp)
              while(abs(depg-depl)>1e-2)
                  Psat=Psat+0.01;
          a0=1;
            a1=-(Ru*T/Psat);
            a2=-(Ru*(T/Psat)*b - (aalpha/Psat) +(b^2));
            a3= -aalpha*(b/Psat);
 
            k=[a0 a1 a2 a3];
 
            vol=roots(k);
            vol=sort(vol);
 
            Vl=vol(1);
            Vv=vol(3);
 
                if(imag(Vv))~=0 
                        Vv=Vl;
                end
          depl= depfl(Vl,Psat,T,aalpha,b,Ru);
          depg= depfg(Vv,Psat,T,aalpha,b,Ru);
                  
              end
          end
       if(temp1>temp)
              while(abs(depg-depl)>1e-2)
                  Psat=Psat-0.01;
          a0=1;
            a1=-(Ru*T/Psat);
            a2=-(Ru*(T/Psat)*b - (aalpha/Psat) +(b^2));
            a3= -aalpha*(b/Psat);
 
            k=[a0 a1 a2 a3];
 
            vol=roots(k);
            vol=sort(vol);
 
            Vl=vol(1);
            Vv=vol(3);
 
                if(imag(Vv))~=0 
                        Vv=Vl;
                end
          depl= depfl(Vl,Psat,T,aalpha,b,Ru);
          depg= depfg(Vv,Psat,T,aalpha,b,Ru);
                  
              end
       end
        
    SatP(2*i-1)=Psat;
    SatP(2*i)=Psat;
    SatV(2*i-1)=Vl;
    SatV(2*i)=Vv;

        end

%making last vector null
SatV(12) = [];  
SatP(12) = [];

%spline and saturation points plotting
y=SatP;
x=SatV;

plot(SatV,SatP,'o');

xx = 0.42:0.01:0.78;
hold on
yy = spline(x,y,xx);
plot(x,y,'o',xx,yy,'-.','Linewidth',1.5);

%isotherm for T>Tc
 for T=785:1:786     
        alpha= ((1 + (0.48508 + 1.574*w - 0.176*(w^2))*(1 - (T/Tc)^0.5))^2) ;
        aalpha = a*((1 + (0.48508 + 1.574*w - 0.176*(w^2))*(1 - (T/Tc)^0.5))^2) ;
        syms v;
        fplot((((Ru*T)/(v-b)) - ((aalpha)/(v*(v+b)))) , [0 2]);
 end

 %departure functions      
function dep=depfl(Vl,Psat,T,aalpha,b,Ru)
        Roh=1/Vl;
        Z=Psat*Vl/(Ru*T);
        dep=-log(1-Roh*b)-(aalpha/(b*Ru*T)*log(1+Roh*b))+Z-1-log(Z);
end

function dep=depfg(Vg,Psat,T,aalpha,b,Ru)
        Roh=1/Vg;        
        Z=Psat*Vg/(Ru*T);
        dep=-log(1-Roh*b)-(aalpha/(b*Ru*T)*log(1+Roh*b))+Z-1-log(Z);
end       