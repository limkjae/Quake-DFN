# ********************************************************************   
# *****                                                          *****   
# *****    DISPLACEMENT AND STRAIN AT DEPTH                      *****  
# *****    DUE TO BURIED FINITE FAULT IN A SEMIINFINITE MEDIUM   *****  
# *****                         CODED BY  Y.OKADA ... SEP 1991   *****  
# *****                         REVISED   Y.OKADA ... NOV 1991   *****  
# *****           CONVERTED TO JULIA                             *****
# *****           WITH MINOR FIX          K. IM   ... OCT 2022   *****
# *****                                                          *****  
# ********************************************************************  

function DCCON0_V(ALPHA,DIP)

  F0 = 0 #zeros(N_CELL);
  F1 = 1 #ones(N_CELL);
  F2 = 2 #ones(N_CELL).*2.0;
  PI2 = 6.283185307179586; #ones(N_CELL).*6.283185307179586;
  EPS = 1.0e-6;#ones(N_CELL).*1.0e-6;
  
  ALP1=(F1-ALPHA)/F2;
  ALP2= ALPHA/F2;
  ALP3=(F1-ALPHA)/ALPHA;
  ALP4= F1-ALPHA;
  ALP5= ALPHA;
  
  P18=PI2/360.0; 
  SD=sin(DIP*P18);    
  CD=cos(DIP*P18); 
  
  c1 = Int8(abs(CD) < EPS);
  c2 = Int8(abs(CD) >= EPS);
  s1 = Int8(SD > F0);
  s2 = Int8(SD == F0);
  s3 = Int8(SD < F0);
  
  CD = F0*c1 + CD*c2;
  SD = c1*(F1*s1 + SD*s2 + (-1.0)*F1*s3) + c2*SD;      
  SDSD=SD*SD;     
  CDCD=CD*CD;      
  SDCD=SD*CD;   
  #global S2D=F2.*SDCD;     
  #global C2D=CDCD-SDSD;

  return ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD


  
end


function DCCON2_V(XI,ET,Q,SD,CD,N_CELL)
F0 = 0
F1 = 1
F2 = 2.0;
EPS = 0.000001;

c1 = Int8.(abs.(XI) .< EPS);
c2 = Int8.(abs.(XI) .>= EPS);
XI = F0.*c1 + XI.*c2;

c1 = Int8.(abs.(ET) .< EPS);
c2 = Int8.(abs.(ET) .>= EPS);
ET = F0.*c1 + ET.*c2;

c1 = Int8.(abs.(Q) .< EPS);
c2 = Int8.(abs.(Q) .>= EPS);
Q = F0.*c1 + Q.*c2;

XI2 = XI .* XI;
ET2 = ET .* ET;
Q2 = Q .* Q;
R2 = XI2 .+ ET2 .+ Q2;
R = sqrt.(R2);

c1 = Int8.(R .== F0);
c1_sum = sum(c1);
if c1_sum > 0
    return
end
R3 = R .* R2;
R5 = R3 .* R2;
Y = ET .* CD + Q .* SD;
D = ET .* SD - Q .* CD;

c1 = Int8.(Q .== F0);
c2 = Int8.(Q .!= F0);

TT = c1 .* F0 .+ c2 .* atan.(XI .* ET ./ (Q .* R));

c1 = Int8.(XI .< F0); 
c2 = Int8.(Q .== F0); 
c3 = Int8.(ET .== F0);
c4 = c1 .* c2 .* c3;
# c5 = zeros(N_CELL,1); 
c5 = 1 .- c4;
RXI = R .+ XI;
########### KJ Revision for Stability!!! ###########
##### Added if c5==0 due to occational NaN ######### 
ALX = zeros(N_CELL)
X11 = zeros(N_CELL)
X32 = zeros(N_CELL)
for i = 1:N_CELL
  if c5[i] == 0
    ALX[i] = (-log(R[i]-XI[i])).*c4[i] ;
    X11[i] = F0.*c4[i] ;
    X32[i] = F0.*c4[i] ;
  else
  ALX[i] = (-log(R[i]-XI[i])).*c4[i] + log(RXI[i]).*c5[i];
  X11[i] = F0.*c4[i] + (F1./(R[i].*RXI[i])).*c5[i];
  X32[i] = F0.*c4[i] + ((R[i]+RXI[i]).*X11[i].*X11[i]./R[i]) .*c5[i];
  end
end


#if sum(isnan.(  (F1/(R*RXI))*c5  ))>0;  
#  println("NaN Here!!!!!    ",(R.*RXI) ,"   ",c5); 
#end

c1 = Int8.(ET .< F0)
c2 = Int8.(Q .== F0)
c3 = Int8.(XI .== F0)
c4 = c1.*c2.*c3;
c5 = 1.0 .- c4 ;
RET = R .+ ET;

########### KJ Revision for Stability!!! ###########
##### Added if c5==0 due to occational NaN ######### 
ALE = zeros(N_CELL)
Y11 = zeros(N_CELL)
Y32 = zeros(N_CELL)
for i = 1:N_CELL
  if c5[i] == 0
    ALE[i] = (-log(R[i]-ET[i])).*c4[i] ;
    Y11[i] = F0.*c4[i] ;
    Y32[i] = F0.*c4[i] ;
  else
    ALE[i] = (-log(R[i]-ET[i])).*c4[i] + log(RET[i]).*c5[i];
    Y11[i] = F0.*c4[i] + (F1./(R[i].*RET[i])).*c5[i];
    Y32[i] = F0.*c4[i] + ((R[i]+RET[i]).*Y11[i].*Y11[i]./R[i]).*c5[i];
  end
end

EY=SD./R-Y.*Q./R3;  
EZ=CD./R+D.*Q./R3;   
FY=D./R3+XI2.*Y32.*SD;  
FZ=Y./R3+XI2.*Y32.*CD;  
GY=F2.*X11.*SD-Y.*Q.*X32;
GZ=F2.*X11.*CD+D.*Q.*X32;    
HY=D.*Q.*X32+XI.*Q.*Y32.*SD;      
HZ=Y.*Q.*X32+XI.*Q.*Y32.*CD;     
return  XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ
end



function UA_V(XI,ET,Q,DISL1,DISL2,DISL3, ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD,
XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ, N_CELL)

F0 = 0.0
F2 = 2.0 
PI2 =6.283185307179586 


du1 = zeros(N_CELL,12)
du2 = zeros(N_CELL,12)
du3 = zeros(N_CELL,12) 

U=zeros(N_CELL,12);

XY=XI.*Y11;
QX=Q .*X11;
QY=Q .*Y11;

c1 = Int8(DISL1 != F0);
du1[:,1]=    TT./F2 +ALP2.*XI.*QY;
du1[:,2]=           ALP2.*Q./R;
du1[:,3]= ALP1.*ALE -ALP2.*Q.*QY;
du1[:,4]=-ALP1.*QY  -ALP2.*XI2.*Q.*Y32;
du1[:,5]=          -ALP2.*XI.*Q./R3;
du1[:,6]= ALP1.*XY  +ALP2.*XI.*Q2.*Y32;
du1[:,7]= ALP1.*XY.*SD        +ALP2.*XI.*FY+D./F2.*X11;
du1[:,8]=                    ALP2.*EY;
du1[:,9]= ALP1.*(CD./R+QY.*SD) -ALP2.*Q.*FY;
du1[:,10]= ALP1.*XY.*CD        +ALP2.*XI.*FZ+Y./F2.*X11;
du1[:,11]=                    ALP2.*EZ;
du1[:,12]=-ALP1.*(SD./R-QY.*CD) -ALP2.*Q.*FZ;
U=U + (ones(N_CELL,12)*DISL1/PI2) .* du1 .* (ones(N_CELL,12) * c1) 

c2 =  Int8(DISL2 != F0);
du2[:,1]=           ALP2.*Q./R;
du2[:,2]=    TT./F2 +ALP2.*ET.*QX;
du2[:,3]= ALP1.*ALX -ALP2.*Q.*QX;
du2[:,4]=        -ALP2.*XI.*Q./R3;
du2[:,5]= -QY./F2 -ALP2.*ET.*Q./R3;
du2[:,6]= ALP1./R +ALP2.*Q2./R3;
du2[:,7]=                      ALP2.*EY;
du2[:,8]= ALP1.*D.*X11+XY./F2.*SD +ALP2.*ET.*GY;
du2[:,9]= ALP1.*Y.*X11          -ALP2.*Q.*GY;
du2[:,10]=                      ALP2.*EZ;
du2[:,11]= ALP1.*Y.*X11+XY./F2.*CD +ALP2.*ET.*GZ;
du2[:,12]=-ALP1.*D.*X11          -ALP2.*Q.*GZ;
U=U + (ones(N_CELL,12)*DISL2/PI2) .* du2 .* (ones(N_CELL,12) * c2) 


c3 = Int8(DISL3 != F0);
du3[:,1]=-ALP1.*ALE -ALP2.*Q.*QY;
du3[:,2]=-ALP1.*ALX -ALP2.*Q.*QX;
du3[:,3]=    TT./F2 -ALP2.*(ET.*QX+XI.*QY);
du3[:,4]=-ALP1.*XY  +ALP2.*XI.*Q2.*Y32;
du3[:,5]=-ALP1./R   +ALP2.*Q2./R3;
du3[:,6]=-ALP1.*QY  -ALP2.*Q.*Q2.*Y32;
du3[:,7]=-ALP1.*(CD./R+QY.*SD)  -ALP2.*Q.*FY;
du3[:,8]=-ALP1.*Y.*X11         -ALP2.*Q.*GY;
du3[:,9]= ALP1.*(D.*X11+XY.*SD) +ALP2.*Q.*HY;
du3[:,10]= ALP1.*(SD./R-QY.*CD)  -ALP2.*Q.*FZ;
du3[:,11]= ALP1.*D.*X11         -ALP2.*Q.*GZ;
du3[:,12]= ALP1.*(Y.*X11+XY.*CD) +ALP2.*Q.*HZ;
U=U + (ones(N_CELL,12)*DISL3/PI2) .* du3 .* (ones(N_CELL,12) * c3) 


return U
end



function UB_V(XI,ET,Q,DISL1,DISL2,DISL3, ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD,
XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ,N_CELL)

F0 = 0.0
F1 = 1.0 
F2 = 2.0 
PI2 =6.283185307179586 

DU = zeros(N_CELL,12)  

#-----------------

RD = R + D
D11 = F1 ./ (R.*RD)
AJ2 = XI .*Y ./ RD .* D11
AJ5 =-(D .+ Y.*Y ./ RD) .* D11

c1 = Int8.(CD .!= F0)
c2 = Int8.(CD .== F0)
s1 = Int8.(XI .== F0)
s2 = Int8.(XI .!= F0)

tempCD = CD;
tempCDCD = CDCD;
CD = c1.*CD + c2.*1.0e-12;
CDCD = c1.*CDCD + c2.*1.0e-12;

X = sqrt.(XI2 .+ Q2);

RD2 = RD.*RD;
AI4 = c1.*(s1.*F0 + s2.*(F1./CDCD.*( XI./RD.*SDCD+
      F2.*atan.((ET.*(X+Q.*CD)+X.*(R+X).*SD)./(XI.*(R+X).*CD)))))+
      c2.*(XI.*Y./RD2./F2);
AI3 = c1.*((Y.*CD./RD-ALE+SD.*log.(RD))./CDCD)+
      c2.*((ET./RD+Y.*Q./RD2-ALE)./F2);
AK1 = c1.*(XI.*(D11-Y11.*SD)./CD)+c2.*(XI.*Q./RD.*D11);
AK3 = c1.*((Q.*Y11-Y.*D11)./CD)+c2.*(SD./RD.*(XI2.*D11.-F1));
AJ3 = c1.*((AK1-AJ2.*SD)./CD)+c2.*(-XI./RD2.*(Q2.*D11.-F1./F2));
AJ6 = c1.*((AK3-AJ5.*SD)./CD)+c2.*(-Y./RD2.*(XI2.*D11.-F1./F2));

CD = tempCD;
CDCD = tempCDCD;

XY=XI.*Y11;
AI1=-XI./RD.*CD-AI4.*SD;
AI2= log.(RD)+AI3.*SD;
AK2= F1./R+AK3.*SD;
AK4= XY.*CD-AK1.*SD;
AJ1= AJ5.*CD-AJ6.*SD;
AJ4=-XY-AJ2.*CD+AJ3.*SD;

U=zeros(N_CELL,12);

QX=Q.*X11;
QY=Q.*Y11;

c1 = Int8(DISL1 != F0);
DU[:,1]=-XI.*QY-TT -ALP3.*AI1.*SD;
DU[:,2]=-Q./R      +ALP3.*Y./RD.*SD;
DU[:,3]= Q.*QY     -ALP3.*AI2.*SD;
DU[:,4]= XI2.*Q.*Y32 -ALP3.*AJ1.*SD;
DU[:,5]= XI.*Q./R3   -ALP3.*AJ2.*SD;
DU[:,6]=-XI.*Q2.*Y32 -ALP3.*AJ3.*SD;
DU[:,7]=-XI.*FY-D.*X11 +ALP3.*(XY+AJ4).*SD;
DU[:,8]=-EY          +ALP3.*(F1./R+AJ5).*SD;
DU[:,9]= Q.*FY        -ALP3.*(QY-AJ6).*SD;
DU[:,10]=-XI.*FZ-Y.*X11 +ALP3.*AK1.*SD;
DU[:,11]=-EZ          +ALP3.*Y.*D11.*SD;
DU[:,12]= Q.*FZ        +ALP3.*AK2.*SD;
U=U + (ones(N_CELL,12)*DISL1/PI2) .* DU .* (ones(N_CELL,12) * c1) 


c2 = Float64( DISL2 != F0);
DU[:,1]=-Q./R      +ALP3.*AI3.*SDCD;
DU[:,2]=-ET.*QX-TT -ALP3.*XI./RD.*SDCD;
DU[:,3]= Q.*QX     +ALP3.*AI4.*SDCD;
DU[:,4]= XI.*Q./R3     +ALP3.*AJ4.*SDCD;
DU[:,5]= ET.*Q./R3+QY  +ALP3.*AJ5.*SDCD;
DU[:,6]=-Q2./R3       +ALP3.*AJ6.*SDCD;
DU[:,7]=-EY          +ALP3.*AJ1.*SDCD;
DU[:,8]=-ET.*GY-XY.*SD +ALP3.*AJ2.*SDCD;
DU[:,9]= Q.*GY        +ALP3.*AJ3.*SDCD;
DU[:,10]=-EZ          -ALP3.*AK3.*SDCD;
DU[:,11]=-ET.*GZ-XY.*CD -ALP3.*XI.*D11.*SDCD;
DU[:,12]= Q.*GZ        -ALP3.*AK4.*SDCD;
U=U + (ones(N_CELL,12)*DISL2/PI2) .* DU .* (ones(N_CELL,12) * c2)


c3 =  Float64(DISL3 != F0);
DU[:,1]= Q.*QY           -ALP3.*AI3.*SDSD;
DU[:,2]= Q.*QX           +ALP3.*XI./RD.*SDSD;
DU[:,3]= ET.*QX+XI.*QY-TT -ALP3.*AI4.*SDSD;
DU[:,4]=-XI.*Q2.*Y32 -ALP3.*AJ4.*SDSD;
DU[:,5]=-Q2./R3     -ALP3.*AJ5.*SDSD;
DU[:,6]= Q.*Q2.*Y32  -ALP3.*AJ6.*SDSD;
DU[:,7]= Q.*FY -ALP3.*AJ1.*SDSD;
DU[:,8]= Q.*GY -ALP3.*AJ2.*SDSD;
DU[:,9]=-Q.*HY -ALP3.*AJ3.*SDSD;
DU[:,10]= Q.*FZ +ALP3.*AK3.*SDSD;
DU[:,11]= Q.*GZ +ALP3.*XI.*D11.*SDSD;
DU[:,12]=-Q.*HZ +ALP3.*AK4.*SDSD;
U=U + (ones(N_CELL,12)*DISL3/PI2) .* DU .* (ones(N_CELL,12) * c3)


return U
end


function UC_V(XI,ET,Q,DISL1,DISL2,DISL3, Z, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD,
XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32,  N_CELL)

F0 = 0.0
F1 = 1.0 
F2 = 2.0 
F3 = 3.0
PI2 =6.283185307179586 

DU = zeros(N_CELL,12)  

C=D+Z;  
X53=(8.0.*R2+9.0.*R.*XI+F3.*XI2).*X11.*X11.*X11./R2;
Y53=(8.0.*R2+9.0.*R.*ET+F3.*ET2).*Y11.*Y11.*Y11./R2;
H=Q.*CD-Z;
Z32=SD./R3-H.*Y32;
Z53=F3.*SD./R5-H.*Y53;
Y0=Y11-XI2.*Y32;
Z0=Z32-XI2.*Z53;
PPY=CD./R3+Q.*Y32.*SD;
PPZ=SD./R3-Q.*Y32.*CD;
QQ=Z.*Y32+Z32+Z0;
QQY=F3.*C.*D./R5-QQ.*SD;
QQZ=F3.*C.*Y./R5-QQ.*CD+Q.*Y32;
XY=XI.*Y11;
QY=Q.*Y11;
QR=F3.*Q./R5;
CDR=(C+D)./R3;
YY0=Y./R3-Y0.*CD;

U=zeros(N_CELL,12);



c1 = Float64(DISL1 != F0);
DU[:,1]= ALP4.*XY.*CD           -ALP5.*XI.*Q.*Z32;
DU[:,2]= ALP4.*(CD./R+F2.*QY.*SD) -ALP5.*C.*Q./R3;
DU[:,3]= ALP4.*QY.*CD           -ALP5.*(C.*ET./R3-Z.*Y11+XI2.*Z32);
DU[:,4]= ALP4.*Y0.*CD                  -ALP5.*Q.*Z0;
DU[:,5]=-ALP4.*XI.*(CD./R3+F2.*Q.*Y32.*SD) +ALP5.*C.*XI.*QR;
DU[:,6]=-ALP4.*XI.*Q.*Y32.*CD            +ALP5.*XI.*(F3.*C.*ET./R5-QQ);
DU[:,7]=-ALP4.*XI.*PPY.*CD    -ALP5.*XI.*QQY;
DU[:,8]= ALP4.*F2.*(D./R3-Y0.*SD).*SD-Y./R3.*CD -ALP5.*(CDR.*SD-ET./R3-C.*Y.*QR);
DU[:,9]=-ALP4.*Q./R3+YY0.*SD  +ALP5.*(CDR.*CD+C.*D.*QR-(Y0.*CD+Q.*Z0).*SD);
DU[:,10]= ALP4.*XI.*PPZ.*CD    -ALP5.*XI.*QQZ;
DU[:,11]= ALP4.*F2.*(Y./R3-Y0.*CD).*SD+D./R3.*CD -ALP5.*(CDR.*CD+C.*D.*QR);
DU[:,12]= YY0.*CD -ALP5.*(CDR.*SD-C.*Y.*QR-Y0.*SDSD+Q.*Z0.*CD);
U=U + (ones(N_CELL,12)*DISL1/PI2) .* DU .* (ones(N_CELL,12) * c1) 

c2 = Float64(DISL2 != F0);
DU[:,1]= ALP4.*CD./R -QY.*SD -ALP5.*C.*Q./R3;
DU[:,2]= ALP4.*Y.*X11       -ALP5.*C.*ET.*Q.*X32;
DU[:,3]=     -D.*X11-XY.*SD -ALP5.*C.*(X11-Q2.*X32);
DU[:,4]=-ALP4.*XI./R3.*CD +ALP5.*C.*XI.*QR +XI.*Q.*Y32.*SD;
DU[:,5]=-ALP4.*Y./R3     +ALP5.*C.*ET.*QR;
DU[:,6]=    D./R3-Y0.*SD +ALP5.*C./R3.*(F1 .- F3.*Q2./R2);
DU[:,7]=-ALP4.*ET./R3+Y0.*SDSD -ALP5.*(CDR.*SD-C.*Y.*QR);
DU[:,8]= ALP4.*(X11-Y.*Y.*X32) -ALP5.*C.*((D+F2.*Q.*CD).*X32-Y.*ET.*Q.*X53);
DU[:,9]=  XI.*PPY.*SD+Y.*D.*X32 +ALP5.*C.*((Y+F2.*Q.*SD).*X32-Y.*Q2.*X53);
DU[:,10]=      -Q./R3+Y0.*SDCD -ALP5.*(CDR.*CD+C.*D.*QR);
DU[:,11]= ALP4.*Y.*D.*X32       -ALP5.*C.*((Y-F2.*Q.*SD).*X32+D.*ET.*Q.*X53);
DU[:,12]=-XI.*PPZ.*SD+X11-D.*D.*X32-ALP5.*C.*((D-F2.*Q.*CD).*X32-D.*Q2.*X53);
U=U + (ones(N_CELL,12)*DISL2/PI2) .* DU .* (ones(N_CELL,12) * c2) 

c3 = Float64(DISL3 != F0);
DU[:,1]=-ALP4.*(SD./R+QY.*CD)   -ALP5.*(Z.*Y11-Q2.*Z32);
DU[:,2]= ALP4.*F2.*XY.*SD+D.*X11 -ALP5.*C.*(X11-Q2.*X32);
DU[:,3]= ALP4.*(Y.*X11+XY.*CD)  +ALP5.*Q.*(C.*ET.*X32+XI.*Z32);
DU[:,4]= ALP4.*XI./R3.*SD+XI.*Q.*Y32.*CD+ALP5.*XI.*(F3.*C.*ET./R5 .- F2.*Z32-Z0);
DU[:,5]= ALP4.*F2.*Y0.*SD-D./R3 +ALP5.*C./R3.*(F1 .- F3.*Q2./R2);
DU[:,6]=-ALP4.*YY0           -ALP5.*(C.*ET.*QR-Q.*Z0);
DU[:,7]= ALP4.*(Q./R3+Y0.*SDCD)   +ALP5.*(Z./R3.*CD+C.*D.*QR-Q.*Z0.*SD);
DU[:,8]=-ALP4.*F2.*XI.*PPY.*SD-Y.*D.*X32 +ALP5.*C.*((Y+F2.*Q.*SD).*X32-Y.*Q2.*X53);
DU[:,9]=-ALP4.*(XI.*PPY.*CD-X11+Y.*Y.*X32) + ALP5.*(C.*((D+F2.*Q.*CD).*X32-Y.*ET.*Q.*X53)+XI.*QQY);
DU[:,10]=  -ET./R3+Y0.*CDCD -ALP5.*(Z./R3.*SD-C.*Y.*QR-Y0.*SDSD+Q.*Z0.*CD);
DU[:,11]= ALP4.*F2.*XI.*PPZ.*SD-X11+D.*D.*X32 -ALP5.*C.*((D-F2.*Q.*CD).*X32-D.*Q2.*X53);
DU[:,12]= ALP4.*(XI.*PPZ.*CD+Y.*D.*X32) +ALP5.*(C.*((Y-F2.*Q.*SD).*X32+D.*ET.*Q.*X53)+XI.*QQZ);
U=U + (ones(N_CELL,12)*DISL3/PI2) .* DU .* (ones(N_CELL,12) * c3) 

   
return U

end



function Okada_DC3D_Vector(ALPHA, X,Y,Z,DEPTH,DIP_Angle, AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3)



  XOriginal=copy(X)
  YOriginal=copy(Y)
  ZOriginal=copy(Z)
  
  # from here function begins
  
  N_CELL= length(X);
  if maximum(Z) .> 0.0
      println(" ** POSITIVE Z WAS GIVEN IN SUB-DC3D");
  end
  U    = zeros(N_CELL,12);
  DUA  = zeros(N_CELL,12);
  DUB  = zeros(N_CELL,12);
  DUC  = zeros(N_CELL,12);
  IRET = zeros(N_CELL);
  
  DIP = DIP_Angle;
  
  # DCCON0 function here
  
  ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD = DCCON0_V(ALPHA,DIP)
  


  # C======================================                                 05080000
  # C=====  REAL-SOURCE CONTRIBUTION  =====                                 05090000
  # C======================================                                 05100000
  D = DEPTH .+ Z;
  P = Y .* CD + D .* SD;
  Q = Y .* SD - D .* CD;
  JXI = zeros(N_CELL);
  JET = zeros(N_CELL);
  aa = (X .+ AL1) .* (X .- AL2);
  cneg = Int8.(aa .<= 0.0);
  JXI = JXI .+ cneg;
  bb = (P .+ AW1) .* (P .- AW2);
  cneg = Int8.(bb .<= 0.0); 
  JET = JET .+ cneg;
  DD1 = DISL1;
  DD2 = DISL2;
  DD3 = DISL3;
  
  K=0
  J=0
  ET=0
  DU=zeros(N_CELL,12)

  for K= 1:2 #2
    if K==1
        ET = P .+ AW1;
    end
    if K==2
       ET = P .- AW2;
    end
    for J=1:2 #2
        if J==1
          XI = X .+ AL1;
        end
        if J==2
          XI = X .- AL2;
        end
        XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ = DCCON2_V(XI,ET,Q,SD,CD, N_CELL)
        cjxi1 = Int8.(JXI .== 1);
        cjet1 = Int8.(JET .== 1);
        cq1   = Int8.(abs.(Q)   .<= 1.0e-12);
        cet1  = Int8.(abs.(ET)  .<= 1.0e-12);
        cxi1  = Int8.(abs.(XI)  .<= 1.0e-12);
        cc1 = cjxi1 .* cq1 .* cet1; 
        cc3 = cjet1 .* cq1 .* cxi1
        cc0 = Int8.((cc1 .+ cc3) .>= 1);
        IRET = IRET .+ cc0;
  
        DUA = UA_V(XI,ET,Q,DD1,DD2,DD3, ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD,
        XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ, N_CELL);
  
        for I=1:3:10
          DU[:,I]  =-DUA[:,I];
          DU[:,I+1]=-DUA[:,I+1].*CD+DUA[:,I+2].*SD;
          DU[:,I+2]=-DUA[:,I+1].*SD-DUA[:,I+2].*CD;
          if I<10.0
            continue;
          end
          DU[:,I]  =-DU[:,I];
          DU[:,I+1]=-DU[:,I+1];
          DU[:,I+2]=-DU[:,I+2];
        end
  
        
        if (J+K)!=3
          U[:,1:12]=U[:,1:12]+DU[:,1:12];
        end
  
        if (J+K)==3
          U[:,1:12]=U[:,1:12]-DU[:,1:12];
        end
      end
  end
      

  X=XOriginal
  Y=YOriginal
  Z=ZOriginal
  
  D = DEPTH .- Z;
  P = Y .* CD .+ D.*SD;
  Q = Y .* SD .- D.*CD;
  
  JET = 1;
  
  cc = (P .+ AW1) .* (P .- AW2);
  
  c1 = Int8.(cc .<= 0.0);
  JET = c1 .* JET;
          
  for K= 1:2 #2
    if K==1
      ET = P .+ AW1;
    end
    if K==2
      ET = P .- AW2;
    end
    for J=1:2 #2
      if J==1
        XI = X .+ AL1;
      end
      if J==2
        XI = X .- AL2;
      end
      XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ = DCCON2_V(XI,ET,Q,SD,CD,N_CELL);
      
      DUA = UA_V(XI,ET,Q,DD1,DD2,DD3, ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD,
      XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ,N_CELL);

      DUB = UB_V(XI,ET,Q,DD1,DD2,DD3, ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD,
      XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ,N_CELL);

      DUC = UC_V(XI,ET,Q,DD1,DD2,DD3, Z, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD,
      XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, N_CELL);
  
      for I=1:3:10
        DU[:,I]=DUA[:,I]+DUB[:,I]+Z.*DUC[:,I];
        DU[:,I+1]=(DUA[:,I+1]+DUB[:,I+1]+Z.*DUC[:,I+1]).*CD -(DUA[:,I+2]+DUB[:,I+2]+Z.*DUC[:,I+2]).*SD;
        DU[:,I+2]=(DUA[:,I+1]+DUB[:,I+1]-Z.*DUC[:,I+1]).*SD +(DUA[:,I+2]+DUB[:,I+2]-Z.*DUC[:,I+2]).*CD;
        if I<10.0
          continue;
        end
        DU[:,10]=DU[:,10]+DUC[:,1];
        DU[:,11]=DU[:,11]+DUC[:,2].*CD-DUC[:,3].*SD;
        DU[:,12]=DU[:,12]-DUC[:,2].*SD-DUC[:,3].*CD;
      end
  
  
      if (J+K)!=3
        U[:,1:12]=U[:,1:12]+DU[:,1:12];
      end
      if (J+K)==3
        U[:,1:12]=U[:,1:12]-DU[:,1:12];
      end
  
      
    end
  end
  
  UX=U[:,1];
  UY=U[:,2];
  UZ=U[:,3];
  UXX=U[:,4];
  UYX=U[:,5];
  UZX=U[:,6];
  UXY=U[:,7];
  UYY=U[:,8];
  UZY=U[:,9];
  UXZ=U[:,10];
  UYZ=U[:,11];
  UZZ=U[:,12];
  cc5 = Int8.(IRET .>= 1);
  IRET = cc5;


  return UX, UY, UZ, UXX, UYX, UZX, UXY, UYY, UZY, UXZ, UYZ, UZZ, IRET
  #return U

end








#################################################################################
#################################################################################
##################### Same code with Single Receiver ############################
#################################################################################
#################################################################################

#=


function DCCON0(ALPHA,DIP)

  F0 = 0 #zeros(N_CELL,1);
  F1 = 1 #ones(N_CELL,1);
  F2 = 2 #ones(N_CELL,1).*2.0;
  PI2 = 6.283185307179586; #ones(N_CELL,1).*6.283185307179586;
  EPS = 1.0e-6;#ones(N_CELL,1).*1.0e-6;
  
  ALP1=(F1.-ALPHA)./F2;
  ALP2= ALPHA./F2;
  ALP3=(F1.-ALPHA)./ALPHA;
  ALP4= F1.-ALPHA;
  ALP5= ALPHA;
  
  P18=PI2./360.0; 
  SD=sin(DIP.*P18);    
  CD=cos(DIP.*P18); 
  
  c1 = Float64(abs(CD) < EPS);
  c2 = Float64(abs(CD) >= EPS);
  s1 = Float64(SD > F0);
  s2 = Float64(SD == F0);
  s3 = Float64(SD < F0);
  
  CD = F0.*c1 + CD.*c2;
  SD = c1.*(F1.*s1 + SD.*s2 + (-1.0).*F1.*s3) + c2.*SD;      
  SDSD=SD.*SD;     
  CDCD=CD.*CD;      
  SDCD=SD.*CD;   
  #global S2D=F2.*SDCD;     
  #global C2D=CDCD-SDSD;

  return ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD


  
end


function DCCON2(XI,ET,Q,SD,CD)
F0 = 0; #zeros(N_CELL,1,'double');
F1 = 1; #ones(N_CELL,1,'double');
F2 = 2; #ones(N_CELL,1,'double').*2.0;
EPS = 0.000001; #ones(N_CELL,1,'double').*0.000001;

c1 = Float64(abs(XI) < EPS);
c2 = Float64(abs(XI) >= EPS);
XI = F0.*c1 + XI.*c2;

c1 = Float64(abs(ET) < EPS);
c2 = Float64(abs(ET) >= EPS);
ET = F0.*c1 + ET.*c2;

c1 = Float64(abs(Q) < EPS);
c2 = Float64(abs(Q) >= EPS);
Q = F0.*c1 + Q.*c2;

XI2=XI.*XI;
ET2=ET.*ET;
Q2=Q.*Q;
R2=XI2+ET2+Q2;
R =sqrt(R2);

c1 = R==F0;
c1_sum = c1;
if c1_sum > 0
    return
end
R3=R .*R2;
R5=R3.*R2;
Y =ET.*CD+Q.*SD;
D =ET.*SD-Q.*CD;

c1 = Float64(Q == F0);
c2 = Float64(Q != F0);
s1 = Float64(Q.*R == F0);
s2 = Float64(Q.*R != F0);


TT = c1.*F0 + c2.*atan(XI.*ET./(Q.*R));

c1 = Float64(XI < F0); 
c2 = Float64(Q == F0); 
c3 = Float64(ET == F0);
c4 = c1.*c2.*c3;
c5 = 0; 
c5 = (c5 - c4)+1.0;
RXI=R+XI;

########### KJ Revision for Stability!!! ###########
##### Added if c5==0 due to occational NaN ######### 
if c5==0    
  ALX = (-log(R-XI)).*c4 ;
  X11 = F0.*c4 ;
  X32 = F0.*c4 ;
else
ALX = (-log(R-XI)).*c4 + log(RXI).*c5;
X11 = F0.*c4 + (F1./(R.*RXI)).*c5;
X32 = F0.*c4 + ((R+RXI).*X11.*X11./R) .*c5;
end

#if sum(isnan.(  (F1/(R*RXI))*c5  ))>0;  
#  println("NaN Here!!!!!    ",(R.*RXI) ,"   ",c5); 
#end

c1 = ET < F0
c2 = Q == F0
c3 = XI == F0
c4 = c1.*c2.*c3;
c5 = 0 
c5 = (c5 - c4)+1.0;
RET=R+ET;

########### KJ Revision for Stability!!! ###########
##### Added if c5==0 due to occational NaN ######### 
if c5==0   
ALE = (-log(R-ET)).*c4 ;
Y11 = F0.*c4 ;
Y32 = F0.*c4 ;
else
ALE = (-log(R-ET)).*c4 + log(RET).*c5;
Y11 = F0.*c4 + (F1./(R.*RET)).*c5;
Y32 = F0.*c4 + ((R+RET).*Y11.*Y11./R).*c5;
end


EY=SD./R-Y.*Q./R3;  
EZ=CD./R+D.*Q./R3;   
FY=D./R3+XI2.*Y32.*SD;  
FZ=Y./R3+XI2.*Y32.*CD;  
GY=F2.*X11.*SD-Y.*Q.*X32;
GZ=F2.*X11.*CD+D.*Q.*X32;    
HY=D.*Q.*X32+XI.*Q.*Y32.*SD;      
HZ=Y.*Q.*X32+XI.*Q.*Y32.*CD;     
return  XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ
end



function UA(XI,ET,Q,DISL1,DISL2,DISL3, ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD,
XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ)

F0 = 0.0
F2 = 2.0 
PI2 =6.283185307179586 

DU = zeros(1,12)
du1 = zeros(1,12)
du2 = zeros(1,12)
du3 = zeros(1,12) 

U=zeros(1,12);

XY=XI.*Y11;
QX=Q .*X11;
QY=Q .*Y11;

c1 = Float64(DISL1 != F0);
du1[1,1]=    TT./F2 +ALP2.*XI.*QY;
du1[1,2]=           ALP2.*Q./R;
du1[1,3]= ALP1.*ALE -ALP2.*Q.*QY;
du1[1,4]=-ALP1.*QY  -ALP2.*XI2.*Q.*Y32;
du1[1,5]=          -ALP2.*XI.*Q./R3;
du1[1,6]= ALP1.*XY  +ALP2.*XI.*Q2.*Y32;
du1[1,7]= ALP1.*XY.*SD        +ALP2.*XI.*FY+D./F2.*X11;
du1[1,8]=                    ALP2.*EY;
du1[1,9]= ALP1.*(CD./R+QY.*SD) -ALP2.*Q.*FY;
du1[1,10]= ALP1.*XY.*CD        +ALP2.*XI.*FZ+Y./F2.*X11;
du1[1,11]=                    ALP2.*EZ;
du1[1,12]=-ALP1.*(SD./R-QY.*CD) -ALP2.*Q.*FZ;
U=U + (ones(1,12)*DISL1/PI2) .* du1 .* (ones(1,12) * c1) 

c2 =  Float64(DISL2 != F0);
du2[1,1]=           ALP2.*Q./R;
du2[1,2]=    TT./F2 +ALP2.*ET.*QX;
du2[1,3]= ALP1.*ALX -ALP2.*Q.*QX;
du2[1,4]=        -ALP2.*XI.*Q./R3;
du2[1,5]= -QY./F2 -ALP2.*ET.*Q./R3;
du2[1,6]= ALP1./R +ALP2.*Q2./R3;
du2[1,7]=                      ALP2.*EY;
du2[1,8]= ALP1.*D.*X11+XY./F2.*SD +ALP2.*ET.*GY;
du2[1,9]= ALP1.*Y.*X11          -ALP2.*Q.*GY;
du2[1,10]=                      ALP2.*EZ;
du2[1,11]= ALP1.*Y.*X11+XY./F2.*CD +ALP2.*ET.*GZ;
du2[1,12]=-ALP1.*D.*X11          -ALP2.*Q.*GZ;
U=U + (ones(1,12)*DISL2/PI2) .* du2 .* (ones(1,12) * c2) 


c3 = Float64(DISL3 != F0);
du3[1,1]=-ALP1.*ALE -ALP2.*Q.*QY;
du3[1,2]=-ALP1.*ALX -ALP2.*Q.*QX;
du3[1,3]=    TT./F2 -ALP2.*(ET.*QX+XI.*QY);
du3[1,4]=-ALP1.*XY  +ALP2.*XI.*Q2.*Y32;
du3[1,5]=-ALP1./R   +ALP2.*Q2./R3;
du3[1,6]=-ALP1.*QY  -ALP2.*Q.*Q2.*Y32;
du3[1,7]=-ALP1.*(CD./R+QY.*SD)  -ALP2.*Q.*FY;
du3[1,8]=-ALP1.*Y.*X11         -ALP2.*Q.*GY;
du3[1,9]= ALP1.*(D.*X11+XY.*SD) +ALP2.*Q.*HY;
du3[1,10]= ALP1.*(SD./R-QY.*CD)  -ALP2.*Q.*FZ;
du3[1,11]= ALP1.*D.*X11         -ALP2.*Q.*GZ;
du3[1,12]= ALP1.*(Y.*X11+XY.*CD) +ALP2.*Q.*HZ;
U=U + (ones(1,12)*DISL3/PI2) .* du3 .* (ones(1,12) * c3) 


return U
end



function UB(XI,ET,Q,DISL1,DISL2,DISL3, ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD,
XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ)

F0 = 0.0
F1 = 1.0 
F2 = 2.0 
PI2 =6.283185307179586 

DU = zeros(1,12)  

#-----------------

RD=R+D;
D11=F1./(R.*RD);
AJ2=XI.*Y./RD.*D11;
AJ5=-(D+Y.*Y./RD).*D11;

c1 = Float64(CD != F0);
c2 = Float64( CD == F0);
s1 = Float64( XI == F0);
s2 = Float64( XI != F0);


tempCD = CD;
tempCDCD = CDCD;
CD = c1.*CD + c2.*1.0e-12;
CDCD = c1.*CDCD + c2.*1.0e-12;

X=sqrt(XI2+Q2);

RD2=RD.*RD;
AI4 = c1.*(s1.*F0 + s2.*(F1./CDCD.*( XI./RD.*SDCD+
      F2.*atan((ET.*(X+Q.*CD)+X.*(R+X).*SD)./(XI.*(R+X).*CD)))))+
      c2.*(XI.*Y./RD2./F2);
AI3 = c1.*((Y.*CD./RD-ALE+SD.*log(RD))./CDCD)+
      c2.*((ET./RD+Y.*Q./RD2-ALE)./F2);
AK1 = c1.*(XI.*(D11-Y11.*SD)./CD)+c2.*(XI.*Q./RD.*D11);
AK3 = c1.*((Q.*Y11-Y.*D11)./CD)+c2.*(SD./RD.*(XI2.*D11-F1));
AJ3 = c1.*((AK1-AJ2.*SD)./CD)+c2.*(-XI./RD2.*(Q2.*D11-F1./F2));
AJ6 = c1.*((AK3-AJ5.*SD)./CD)+c2.*(-Y./RD2.*(XI2.*D11-F1./F2));

CD = tempCD;
CDCD = tempCDCD;

XY=XI.*Y11;
AI1=-XI./RD.*CD-AI4.*SD;
AI2= log(RD)+AI3.*SD;
AK2= F1./R+AK3.*SD;
AK4= XY.*CD-AK1.*SD;
AJ1= AJ5.*CD-AJ6.*SD;
AJ4=-XY-AJ2.*CD+AJ3.*SD;

U=zeros(1,12);

QX=Q.*X11;
QY=Q.*Y11;

c1 = Float64(DISL1 != F0);
DU[1,1]=-XI.*QY-TT -ALP3.*AI1.*SD;
DU[1,2]=-Q./R      +ALP3.*Y./RD.*SD;
DU[1,3]= Q.*QY     -ALP3.*AI2.*SD;
DU[1,4]= XI2.*Q.*Y32 -ALP3.*AJ1.*SD;
DU[1,5]= XI.*Q./R3   -ALP3.*AJ2.*SD;
DU[1,6]=-XI.*Q2.*Y32 -ALP3.*AJ3.*SD;
DU[1,7]=-XI.*FY-D.*X11 +ALP3.*(XY+AJ4).*SD;
DU[1,8]=-EY          +ALP3.*(F1./R+AJ5).*SD;
DU[1,9]= Q.*FY        -ALP3.*(QY-AJ6).*SD;
DU[1,10]=-XI.*FZ-Y.*X11 +ALP3.*AK1.*SD;
DU[1,11]=-EZ          +ALP3.*Y.*D11.*SD;
DU[1,12]= Q.*FZ        +ALP3.*AK2.*SD;
U=U + (ones(1,12)*DISL1/PI2) .* DU .* (ones(1,12) * c1) 


c2 = Float64( DISL2 != F0);
DU[1,1]=-Q./R      +ALP3.*AI3.*SDCD;
DU[1,2]=-ET.*QX-TT -ALP3.*XI./RD.*SDCD;
DU[1,3]= Q.*QX     +ALP3.*AI4.*SDCD;
DU[1,4]= XI.*Q./R3     +ALP3.*AJ4.*SDCD;
DU[1,5]= ET.*Q./R3+QY  +ALP3.*AJ5.*SDCD;
DU[1,6]=-Q2./R3       +ALP3.*AJ6.*SDCD;
DU[1,7]=-EY          +ALP3.*AJ1.*SDCD;
DU[1,8]=-ET.*GY-XY.*SD +ALP3.*AJ2.*SDCD;
DU[1,9]= Q.*GY        +ALP3.*AJ3.*SDCD;
DU[1,10]=-EZ          -ALP3.*AK3.*SDCD;
DU[1,11]=-ET.*GZ-XY.*CD -ALP3.*XI.*D11.*SDCD;
DU[1,12]= Q.*GZ        -ALP3.*AK4.*SDCD;
U=U + (ones(1,12)*DISL2/PI2) .* DU .* (ones(1,12) * c2)


c3 =  Float64(DISL3 != F0);
DU[1,1]= Q.*QY           -ALP3.*AI3.*SDSD;
DU[1,2]= Q.*QX           +ALP3.*XI./RD.*SDSD;
DU[1,3]= ET.*QX+XI.*QY-TT -ALP3.*AI4.*SDSD;
DU[1,4]=-XI.*Q2.*Y32 -ALP3.*AJ4.*SDSD;
DU[1,5]=-Q2./R3     -ALP3.*AJ5.*SDSD;
DU[1,6]= Q.*Q2.*Y32  -ALP3.*AJ6.*SDSD;
DU[1,7]= Q.*FY -ALP3.*AJ1.*SDSD;
DU[1,8]= Q.*GY -ALP3.*AJ2.*SDSD;
DU[1,9]=-Q.*HY -ALP3.*AJ3.*SDSD;
DU[1,10]= Q.*FZ +ALP3.*AK3.*SDSD;
DU[1,11]= Q.*GZ +ALP3.*XI.*D11.*SDSD;
DU[1,12]=-Q.*HZ +ALP3.*AK4.*SDSD;
U=U + (ones(1,12)*DISL3/PI2) .* DU .* (ones(1,12) * c3)


return U
end


function UC(XI,ET,Q,DISL1,DISL2,DISL3, Z, ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD,
XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ)

F0 = 0.0
F1 = 1.0 
F2 = 2.0 
F3 = 3.0
PI2 =6.283185307179586 

DU = zeros(1,12)  

C=D+Z;  
X53=(8.0.*R2+9.0.*R.*XI+F3.*XI2).*X11.*X11.*X11./R2;
Y53=(8.0.*R2+9.0.*R.*ET+F3.*ET2).*Y11.*Y11.*Y11./R2;
H=Q.*CD-Z;
Z32=SD./R3-H.*Y32;
Z53=F3.*SD./R5-H.*Y53;
Y0=Y11-XI2.*Y32;
Z0=Z32-XI2.*Z53;
PPY=CD./R3+Q.*Y32.*SD;
PPZ=SD./R3-Q.*Y32.*CD;
QQ=Z.*Y32+Z32+Z0;
QQY=F3.*C.*D./R5-QQ.*SD;
QQZ=F3.*C.*Y./R5-QQ.*CD+Q.*Y32;
XY=XI.*Y11;
QX=Q.*X11;
QY=Q.*Y11;
QR=F3.*Q./R5;
CQX=C.*Q.*X53;
CDR=(C+D)./R3;
YY0=Y./R3-Y0.*CD;

U=zeros(1,12);



c1 = Float64(DISL1 != F0);
DU[1,1]= ALP4.*XY.*CD           -ALP5.*XI.*Q.*Z32;
DU[1,2]= ALP4.*(CD./R+F2.*QY.*SD) -ALP5.*C.*Q./R3;
DU[1,3]= ALP4.*QY.*CD           -ALP5.*(C.*ET./R3-Z.*Y11+XI2.*Z32);
DU[1,4]= ALP4.*Y0.*CD                  -ALP5.*Q.*Z0;
DU[1,5]=-ALP4.*XI.*(CD./R3+F2.*Q.*Y32.*SD) +ALP5.*C.*XI.*QR;
DU[1,6]=-ALP4.*XI.*Q.*Y32.*CD            +ALP5.*XI.*(F3.*C.*ET./R5-QQ);
DU[1,7]=-ALP4.*XI.*PPY.*CD    -ALP5.*XI.*QQY;
DU[1,8]= ALP4.*F2.*(D./R3-Y0.*SD).*SD-Y./R3.*CD -ALP5.*(CDR.*SD-ET./R3-C.*Y.*QR);
DU[1,9]=-ALP4.*Q./R3+YY0.*SD  +ALP5.*(CDR.*CD+C.*D.*QR-(Y0.*CD+Q.*Z0).*SD);
DU[1,10]= ALP4.*XI.*PPZ.*CD    -ALP5.*XI.*QQZ;
DU[1,11]= ALP4.*F2.*(Y./R3-Y0.*CD).*SD+D./R3.*CD -ALP5.*(CDR.*CD+C.*D.*QR);
DU[1,12]= YY0.*CD -ALP5.*(CDR.*SD-C.*Y.*QR-Y0.*SDSD+Q.*Z0.*CD);
U=U + (ones(1,12)*DISL1/PI2) .* DU .* (ones(1,12) * c1) 

c2 = Float64(DISL2 != F0);
DU[1,1]= ALP4.*CD./R -QY.*SD -ALP5.*C.*Q./R3;
DU[1,2]= ALP4.*Y.*X11       -ALP5.*C.*ET.*Q.*X32;
DU[1,3]=     -D.*X11-XY.*SD -ALP5.*C.*(X11-Q2.*X32);
DU[1,4]=-ALP4.*XI./R3.*CD +ALP5.*C.*XI.*QR +XI.*Q.*Y32.*SD;
DU[1,5]=-ALP4.*Y./R3     +ALP5.*C.*ET.*QR;
DU[1,6]=    D./R3-Y0.*SD +ALP5.*C./R3.*(F1-F3.*Q2./R2);
DU[1,7]=-ALP4.*ET./R3+Y0.*SDSD -ALP5.*(CDR.*SD-C.*Y.*QR);
DU[1,8]= ALP4.*(X11-Y.*Y.*X32) -ALP5.*C.*((D+F2.*Q.*CD).*X32-Y.*ET.*Q.*X53);
DU[1,9]=  XI.*PPY.*SD+Y.*D.*X32 +ALP5.*C.*((Y+F2.*Q.*SD).*X32-Y.*Q2.*X53);
DU[1,10]=      -Q./R3+Y0.*SDCD -ALP5.*(CDR.*CD+C.*D.*QR);
DU[1,11]= ALP4.*Y.*D.*X32       -ALP5.*C.*((Y-F2.*Q.*SD).*X32+D.*ET.*Q.*X53);
DU[1,12]=-XI.*PPZ.*SD+X11-D.*D.*X32-ALP5.*C.*((D-F2.*Q.*CD).*X32-D.*Q2.*X53);
U=U + (ones(1,12)*DISL2/PI2) .* DU .* (ones(1,12) * c2) 

c3 = Float64(DISL3 != F0);
DU[1,1]=-ALP4.*(SD./R+QY.*CD)   -ALP5.*(Z.*Y11-Q2.*Z32);
DU[1,2]= ALP4.*F2.*XY.*SD+D.*X11 -ALP5.*C.*(X11-Q2.*X32);
DU[1,3]= ALP4.*(Y.*X11+XY.*CD)  +ALP5.*Q.*(C.*ET.*X32+XI.*Z32);
DU[1,4]= ALP4.*XI./R3.*SD+XI.*Q.*Y32.*CD+ALP5.*XI.*(F3.*C.*ET./R5-F2.*Z32-Z0);
DU[1,5]= ALP4.*F2.*Y0.*SD-D./R3 +ALP5.*C./R3.*(F1-F3.*Q2./R2);
DU[1,6]=-ALP4.*YY0           -ALP5.*(C.*ET.*QR-Q.*Z0);
DU[1,7]= ALP4.*(Q./R3+Y0.*SDCD)   +ALP5.*(Z./R3.*CD+C.*D.*QR-Q.*Z0.*SD);
DU[1,8]=-ALP4.*F2.*XI.*PPY.*SD-Y.*D.*X32 +ALP5.*C.*((Y+F2.*Q.*SD).*X32-Y.*Q2.*X53);
DU[1,9]=-ALP4.*(XI.*PPY.*CD-X11+Y.*Y.*X32) + ALP5.*(C.*((D+F2.*Q.*CD).*X32-Y.*ET.*Q.*X53)+XI.*QQY);
DU[1,10]=  -ET./R3+Y0.*CDCD -ALP5.*(Z./R3.*SD-C.*Y.*QR-Y0.*SDSD+Q.*Z0.*CD);
DU[1,11]= ALP4.*F2.*XI.*PPZ.*SD-X11+D.*D.*X32 -ALP5.*C.*((D-F2.*Q.*CD).*X32-D.*Q2.*X53);
DU[1,12]= ALP4.*(XI.*PPZ.*CD+Y.*D.*X32) +ALP5.*(C.*((Y-F2.*Q.*SD).*X32+D.*ET.*Q.*X53)+XI.*QQZ);
U=U + (ones(1,12)*DISL3/PI2) .* DU .* (ones(1,12) * c3) 

   
return U

end



function Okada_DC3D(ALPHA, X,Y,Z,DEPTH,DIP_Angle, AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3)




  XOriginal=copy(X)
  YOriginal=copy(Y)
  ZOriginal=copy(Z)
  
  # from here function begins
  
  N_CELL= length(X);
  F0 = zeros(N_CELL,1);
  if Z>0.0
      println(" ** POSITIVE Z WAS GIVEN IN SUB-DC3D");
  end
  U    = zeros(N_CELL,12);
  DUA  = zeros(N_CELL,12);
  DUB  = zeros(N_CELL,12);
  DUC  = zeros(N_CELL,12);
  IRET = 0; #zeros(N_CELL,1);
  
  AALPHA=ALPHA;
  DIP=DIP_Angle;
  
  # DCCON0 function here
  
  ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD = DCCON0(ALPHA,DIP)
  
  # C======================================                                 05080000
  # C=====  REAL-SOURCE CONTRIBUTION  =====                                 05090000
  # C======================================                                 05100000
  D=DEPTH+Z;
  P=Y.*CD+D.*SD;
  Q=Y.*SD-D.*CD;
  JXI =0 # zeros(Int8, N_CELL,1);
  JET =0 # zeros(Int8, N_CELL,1);
  aa = (X+AL1).*(X-AL2);
  cneg = Int8(aa <= 0.0);
  JXI = JXI + cneg;
  jxi_sum = JXI;
  bb = (P+AW1).*(P-AW2);
  cneg = Int8(bb <= 0.0); 
  JET = JET + cneg;
  jet_sum = JET;
  DD1=DISL1;
  DD2=DISL2;
  DD3=DISL3;
  
  
  K=0
  J=0
  ET=0
  DU=zeros(1,12)
  for K= 1:2 #2
    if K==1
        ET=P+AW1;
    end
    if K==2
       ET=P-AW2;
    end
    for J=1:2 #2
        if J==1
          XI=X+AL1;
        end
        if J==2
          XI=X-AL2;
        end
        XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ = DCCON2(XI,ET,Q,SD,CD);
        cjxi1 = Float64(JXI == 1);
        cjxi2 = Float64(JXI != 1);
        cjet1 = Float64(JET == 1);
        cjet2 = Float64(JET != 1);
        cq1   = Float64(abs(Q)   <= 1.0e-12);
        cq2   = Float64(abs(Q)   > 1.0e-12);
        cet1  = Float64(abs(ET)  <= 1.0e-12);
        cet2  = Float64(abs(ET)  > 1.0e-12);
        cxi1  = Float64(abs(XI)  <= 1.0e-12);
        cxi2  = Float64(abs(XI)  > 1.0e-12);
        cc1 = cjxi1.*cq1.*cet1; cc2 = (cc1 - 1.0).*(-1.0);
        cc3 = cjet1.*cq1.*cxi1; cc4 = (cc3 - 1.0).*(-1.0);
        cc0 = Float64((cc1 + cc3) >= 1);
        cc5 = Float64((cc1 + cc3) < 1);
        IRET = IRET + cc0;
  
        DUA = UA(XI,ET,Q,DD1,DD2,DD3, ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD,
        XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ);
  
        for I=1:3:10
          DU[1,I]  =-DUA[1,I];
          DU[1,I+1]=-DUA[1,I+1].*CD+DUA[1,I+2].*SD;
          DU[1,I+2]=-DUA[1,I+1].*SD-DUA[1,I+2].*CD;
          if I<10.0
            continue;
          end
          DU[1,I]  =-DU[1,I];
          DU[1,I+1]=-DU[1,I+1];
          DU[1,I+2]=-DU[1,I+2];
        end
  
        
        if (J+K)!=3
          U[1,1:12]=U[1,1:12]+DU[1,1:12];
        end
  
        if (J+K)==3
          U[1,1:12]=U[1,1:12]-DU[1,1:12];
        end
      end
  end
  
  
  X=XOriginal
  Y=YOriginal
  Z=ZOriginal
  
  ZZ=Z;
  D=DEPTH-Z;
  P=Y.*CD+D.*SD;
  Q=Y.*SD-D.*CD;
  
  JET = 1;
  
  cc=(P+AW1).*(P-AW2);
  
  c1 = Int8(cc <= 0.0);
  c2 = Int8(cc >  0.0);
  JET = c1*JET;
  
  
  
  
  for K= 1:2 #2
    if K==1
      ET=P+AW1;
    end
    if K==2
      ET=P-AW2;
    end
    for J=1:2 #2
      if J==1
        XI=X+AL1;
      end
      if J==2
        XI=X-AL2;
      end
      XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ = DCCON2(XI,ET,Q,SD,CD);
      
      DUA = UA(XI,ET,Q,DD1,DD2,DD3, ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD,
      XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ);
      DUB = UB(XI,ET,Q,DD1,DD2,DD3, ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD,
      XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ);
      DUC = UC(XI,ET,Q,DD1,DD2,DD3, Z, ALP1, ALP2, ALP3, ALP4, ALP5, SD, CD, SDSD, CDCD, SDCD,
      XI2, ET2, Q2, R2, R, R3, R5, Y, D, TT, ALX, X11, X32, ALE, Y11, Y32, EY,EZ,FY,FZ,GY,GZ,HY,HZ);
  
      for I=1:3:10
        DU[1,I]=DUA[1,I]+DUB[1,I]+Z.*DUC[1,I];
        DU[1,I+1]=(DUA[1,I+1]+DUB[1,I+1]+Z.*DUC[1,I+1]).*CD -(DUA[1,I+2]+DUB[1,I+2]+Z.*DUC[1,I+2]).*SD;
        DU[1,I+2]=(DUA[1,I+1]+DUB[1,I+1]-Z.*DUC[1,I+1]).*SD +(DUA[1,I+2]+DUB[1,I+2]-Z.*DUC[1,I+2]).*CD;
        if I<10.0
          continue;
        end
        DU[1,10]=DU[1,10]+DUC[1,1];
        DU[1,11]=DU[1,11]+DUC[1,2].*CD-DUC[1,3].*SD;
        DU[1,12]=DU[1,12]-DUC[1,2].*SD-DUC[1,3].*CD;
      end
  
  
      if (J+K)!=3
        U[1,1:12]=U[1,1:12]+DU[1,1:12];
      end
      if (J+K)==3
        U[1,1:12]=U[1,1:12]-DU[1,1:12];
      end
  
      
    end
  end
  
  
  UX=U[1,1];
  UY=U[1,2];
  UZ=U[1,3];
  UXX=U[1,4];
  UYX=U[1,5];
  UZX=U[1,6];
  UXY=U[1,7];
  UYY=U[1,8];
  UZY=U[1,9];
  UXZ=U[1,10];
  UYZ=U[1,11];
  UZZ=U[1,12];
  cc5 = Int8(IRET >= 1);
  IRET = cc5;



  return UX, UY, UZ, UXX, UYX, UZX, UXY, UYY, UZY, UXZ, UYZ, UZZ, IRET
  #return U

end




=#