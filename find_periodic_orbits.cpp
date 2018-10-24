// Find periodic orbits in a non axisymmetric potential
// Sofia G. Gallego, 2017
// Bar potential used in Gallego & Cuadra 2017

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <sstream>
#include <string>
#include <math.h>

double accx(double x, double y, double z, double u, double v, double w);
double accy(double x, double y, double z, double u, double v, double w);
double accz(double x, double y, double z, double u, double v, double w);
double accr(double r, double G, double Mb0, double ab);
double dpot(double x, double y, double z, int coord);
void runge_kutta(double& x,double& y, double& z,double& u, double& v, double& w, double h);
double P(double c2,double c3,double b20,double b22,double x,double y,double z,double r);
double dP_dx(double c3,double b20,double b22,double x,double y,double z,double r);
double dP_dy(double c3,double b20,double b22,double x,double y,double z,double r);
double dP_dz(double c3,double b20,double b22,double x,double y,double z,double r);
double Y(double b20,double b22,double x,double y,double z,double r);
double dY_dx(double b20,double b22,double x,double y,double z,double r);
double dY_dy(double b20,double b22,double x,double y,double z,double r);
double dY_dz(double b20,double b22,double x,double y,double z,double r);
double f1(double x,double y,double z);
double df1_dx(double x,double y,double z,double r);
double df1_dy(double x,double y,double z,double r);
double df1_dz(double x,double y,double z,double r);
double f2(double x,double y);
double df2_dx(double x,double y);
double df2_dy(double x,double y);
void acc(double x,double y, double z,double u,double v,double w,double& ax,double& ay,double& az);
void dpot(double x, double y,double z,double& dx,double& dy,double& dz);
double dPr(double r[], double ct,double ct2,double ct4,double ct6,double c2p,double c4p,double c6p,double invr16);
double dPt(double r[],double st,double ct2,double ct4,double ct6,double c2p,double c4p,double c6p,double invr13);
double dPp(double r[], double ct,double ct2,double ct4,double ct6,double s2p,double s4p,double s6p,double invr13);

using namespace std;

int main()
{
  double h=0.00001;// 0.00001 timestep of 10kyr aprox (9785yr) t=1 aprox 1Gyr
  double dt=0.,x=0.,y=0.,z=0.,u=0.,v=0.,w=0.,x0=0.,y0=0.,z0=0.,u0=0.,v0=0.,w0=0.,ya=0.,xa=0.,za=0.,ua=0.,va=0.,wa=0.,yaa=0.,xaa=0.,zaa=0.,uaa=0.,vaa=0.,waa=0.;
  int i=0,iold=0,j=0;
  double b, ub;
  int ran,tries=0;
  ifstream fin("periodic-sols.dat");
  ofstream per2("total_orbits.dat");
  double precision=4.e-3,dif=0.,difa=9.e9;
  bool end=false;
  double period=0;
  int N=4000,k=0,NN,n=0;
  double xx[N],vv[N],pp[N],bb[N],uu[N],hh=0,zz=0;
  double pmin=99999999, temp, div;
  srand (time(NULL));
  
  x0=.18;
  y0=0.;
  z0=0.;
  dt=0.;
  u0=0.;
  w0=0;
  
  if(fin.good()){
    
    cout << "Periodic solutions file already exists\n";
    string aaa;
    getline(fin,aaa);
    cout << aaa << endl;
    while(!fin.eof()){
      fin >> xx[j] >> bb[j] >> vv[j] >> uu[j] >> pp[j];
      j++;
    }
    NN=j;
  }
  
  else{
    ofstream per("periodic-sols.dat");
    per << "#semi-major_axis semi-minor_axis vely_a velx_b period\n";
    
    v0=80;//guess for initial velocity in the y axis
    
    for(x0=.05;x0<=.3;x0+=.0001){
      x=x0; y=y0; z=z0; u=u0; v=v0; w=w0;
      
      while(!end){
        i++;
        xaa=xa; yaa=ya; uaa=ua; xa=x; ya=y; za=z; ua=u; va=v; wa=w;
        
        runge_kutta(x,y,z,u,v,w,h);

        if(x*xaa<0 && y>0){// end of quarter of an orbit
          b=y;
          ub=(u+ua)/2.;
        }

        if(((y)*(yaa)<0 && x<0)){ // end of half an orbit
          //-fabs((xa-x)*ya/(ya-y))-fabs(xa) is the interpolation of x in the y axis between 
          //the timesteps close to y axis (one of them could have y very close to 0 or exactly 0)
          dif=fabs(x0-fabs((xa-x)*ya/(ya-y))-fabs(xa))/x0;
          //dif=fabs(x0+x)/x0;//simple interpolation  
          period=2.*(i-1.+ya/(ya-y))*h;//also interpolation          
          
          if(dif<=precision){
            
            printf("Periodic orbit found: a %f v0 %f period %f Myr i %d u %f ua %f\n",x0,v0,period*978.5,i,u,ua);
            
            xx[j]=x0;vv[j]=v0;pp[j]=period;bb[j]=b;uu[j]=abs(ub);
            if(pmin>period) pmin=period;
            printf("Minimum period: %f\n",pmin*978.5);
            j++;
            end=true;
          }
          
          else{			
            v0=fabs(-v0-(x+x0));
            tries++;
            
            if(rand()%1000==0) printf("v0 %f dif %f i %i x %f x0 %f period %f\n",v0,dif,i,abs(x),x0,period*978.5);
            if(tries>1e4){
              v0=v0-1+2*rand()/RAND_MAX;
              printf("Too many tries v0 %f dif %f i %i x %f x0 %f\n",v0,dif,i,abs(x),x0);
              tries=0;
            } 
            x=x0;
            y=y0;
            z=z0;
            u=u0;
            v=v0;
            w=w0;
            
          }
          
          difa=dif;	
          
          i=0;
        }// end of half an orbit
      }// end while(end==0) 
      end=false;
    }// end for(x0=.18-.03,i=0;x0<=.18+.03;x0+=.005,i++)
    
    NN=j;
    
    for(k=0;k<NN;k++){
      
      if(pp[k]>1.5*pmin){
        cout << "Wrong period " << pp[k]*978.5;
        div=double(int(pp[k]/pmin));
        pp[k]=pp[k]/div;
        //if((k!=0)&&(pp[k]<pp[k-1])) pp[k]=pp[k]/double(int(pp[k-1]/pp[k]));
        cout << " corrected " << pp[k]*978.5 << " division " << div << endl;
        pmin=pp[k];
      }
      per << xx[k] << " " << bb[k] << " " << vv[k] << " " << uu[k] << " " << pp[k] << endl;
      
    }
    cout << "Writing periodic-sols file\n";
    per.close();
  }//end else
  cout << "end correcting periodic orbs... " << NN << endl;
  
  //Now create total orbits with x0, period and v0 info
  
  for(k=0;k<NN-1;k++){
    
    N=int(xx[k]*1e5);
    hh=pp[k]*1.2/float(N);
    cout << "Number of steps " << N << " orbit " << k << " of " << NN << " timestep " << hh*978.5 << " period " << pp[k]*978.5 << " Myr" << endl;
    x=xx[k];y=0.;z=0.;u=0.;v=vv[k];w=0.;
    
    per2 << x << " " << y << " " << u << " " << v << " " << xx[k] << endl;
    
    for(j=1;j<N/4;j++){
      runge_kutta(x,y,z,u,v,w,hh);
      per2 << x << " " << y << " " << u << " " << v << " " << xx[k] << endl;
      if(x<=0||y<0) j=N; 
    }
  }
  cout << "end detecting periodic orbs... " << endl;
  
  per2.close();
  
  return 0;  
  
}


void runge_kutta(double& x,double& y, double& z,double& u, double& v, double& w, double h){
  
  double k1=0.,k2=0.,k3=0.,k4=0.;
  double l1=0.,l2=0.,l3=0.,l4=0.;
  double m1=0.,m2=0.,m3=0.,m4=0.;
  double n1=0.,n2=0.,n3=0.,n4=0.;
  double o1=0.,o2=0.,o3=0.,o4=0.;
  double p1=0.,p2=0.,p3=0.,p4=0.;
  double ax,ay,az;
  
  acc(x,y,z,u,v,w,ax,ay,az);
  
  k1=ax;l1=ay;o1=az;
  m1=u;n1=v;p1=w;
  acc(x+h*m1/2.,y+h*n1/2.,z+p1/2.,u+k1*h/2.,v+l1*h/2.,w+o1*h/2.,ax,ay,az);
  k2=ax;l2=ay;o2=az;
  m2=u+k1*h/2.;n2=v+l1*h/2.;p2=w+o1*h/2.;
  acc(x+h*m2/2.,y+h*n2/2.,z+p2/2.,u+k2*h/2.,v+l2*h/2.,w+o2*h/2.,ax,ay,az);
  k3=ax;l3=ay;o3=az;
  m3=u+k2*h/2.;n3=v+l2*h/2.;p3=w+o2*h/2.;
  acc(x+h*m3,y+h*n3,z+h*p3,u+k3*h,v+l3*h,w+o3*h,ax,ay,az);
  k4=ax;l4=ay;o4=az;
  m4=u+k3*h;n4=v+l3*h;p4=w+o3*h;
  
  u+=h*(k1+2.*k2+2.*k3+k4)/6.;
  v+=h*(l1+2.*l2+2.*l3+l4)/6.;
  w+=h*(o1+2.*o2+2.*o3+o4)/6.;
  x+=h*(m1+2.*m2+2.*m3+m4)/6.;
  y+=h*(n1+2.*n2+2.*n3+n4)/6.;
  z+=h*(p1+2.*p2+2.*p3+p4)/6.;
  
  return;
}


void acc(double x,double y, double z,double u,double v,double w,double& ax,double& ay,double& az){
  
  double omega=63.; // bar pattern speed from Kim et al. 2011, 63 km s^-1 kpc^-1
  double dx,dy,dz;
  
  dpot(x,y,z,dx,dy,dz);

  ax=-dx+2.*omega*v+omega*omega*x;
  ay=-dy-2.*omega*u+omega*omega*y;
  az=-dz;
  if(rand()%1000000000==0){
    printf("x %f y %f z %f\n",x,y,z);
    printf("dx %f dy %f dz %f\n",dx,dy,dz);
    printf("ax %f ay %f az %f\n",ax,ay,az);
  }
  return;
}

double accr(double r, double G, double M, double a){
  double temp;
  //temp=(-4.*r*r*(a-.5*r)*log(-r)+4.*r*r*(a-.5*r)*log(a)-5.*r*r*(r+a))/(r*pow(r+a,4.))+a/pow(r+a,3.);
  //return 2.*G*M*a*temp;
  return G*M/r/pow(r+a,2.);
}

void dpot(double x, double y,double z,double& dx,double& dy,double& dz){
  double r0=0.1,rho0=4.,a=0.25,b22=0.1,b20=0.3,G=4.307e4; // for the bar potential
  double c1=4*3.141592*G*rho0*pow(r0,2-a);
  double c2=1/a/(1.+a);
  double c3=1./(2.-a)/(3.+a);
  double dpot=0,r=0,theta=0;
  int model=0;
  double Mb0=1.88, ab=4.56, acc; // values taken from Portail 16?

  if(model==0){
    double r;
    r=sqrt(x*x+y*y+z*z);
    dx=c1*pow(r,a)*(a*pow(r,-2.)*x*P(c2,c3,b20,b22,x,y,z,r)+dP_dx(c3,b20,b22,x,y,z,r));
    dy=c1*pow(r,a)*(a*pow(r,-2.)*y*P(c2,c3,b20,b22,x,y,z,r)+dP_dy(c3,b20,b22,x,y,z,r));
    dz=c1*pow(r,a)*(a*pow(r,-2.)*z*P(c2,c3,b20,b22,x,y,z,r)+dP_dz(c3,b20,b22,x,y,z,r));
  }

  if(model==3){
    //Zhao v2.0 to match Portail+16 model
    a=0.8;b22=0.1;b20=.6;
    double r;
    r=sqrt(x*x+y*y+z*z);
    dx=c1*pow(r,a)*(a*pow(r,-2.)*x*P(c2,c3,b20,b22,x,y,z,r)+dP_dx(c3,b20,b22,x,y,z,r));
    dy=c1*pow(r,a)*(a*pow(r,-2.)*y*P(c2,c3,b20,b22,x,y,z,r)+dP_dy(c3,b20,b22,x,y,z,r));
    dz=c1*pow(r,a)*(a*pow(r,-2.)*z*P(c2,c3,b20,b22,x,y,z,r)+dP_dz(c3,b20,b22,x,y,z,r));
  }
  
  if(model==1){
    acc=accr(r,G,Mb0,ab);
    dx=acc*x;
    dy=acc*y;
    dz=0;
  }
  if(model==2){
    double dr,dt,dp;
    double r[14],invr,invr16,invr13,p;
    r[0]=sqrt(x*x+y*y+z*z);
    invr=1./r[0];
    invr13=pow(r[0]+1.,-13.);
    invr16=pow(r[0]+1.,-16.);
    
    int i;
    for(i=1;i<14;i++) r[i]=r[i-1]*r[0];
    double ct,st,ct2,ct4,ct6,cp,sp,c2p,c4p,c6p,s2p,s4p,s6p;
    ct=y/sqrt(x*x+y*y);
    st=x/sqrt(x*x+y*y);
    ct2=ct*ct;
    ct4=ct2*ct2;
    ct6=ct4*ct2;
    cp=z*invr;
    c2p=2*cp*cp-1;
    p=acos(cp);
    c4p=cos(4*p);
    c6p=cos(6*p);
    sp=sin(p);
    s2p=sin(2*p);
    s4p=sin(4*p);
    s6p=sin(6*p);
    
    dr=dPr(r,ct,ct2,ct4,ct6,c2p,c4p,c6p,invr16);
    dt=dPt(r,st,ct2,ct4,ct6,c2p,c4p,c6p,invr13);
    dp=dPp(r,ct,ct2,ct4,ct6,c2p,c4p,c6p,invr13);
    //printf("dr %f dt %f dp %f (1+r)^-13 %e (1+r)^-16 %e\nr %f p %f\n",dr,dt,dp,invr13,invr16,r[0],p);
    dx=G*Mb0*(dr*ct*sp-dt*st*invr/sp+dp*ct*cp*invr);
    dy=G*Mb0*(dr*st*sp+dt*ct*invr/sp+dp*st*cp*invr);
    if(z==0) dz=0; else dz=G*Mb0*(dr*cp-dp*sp*invr);
  }
}


double dPr(double r[], double ct,double ct2,double ct4,double ct6,double c2p,double c4p,double c6p,double invr16){
  double a;
  a=((2983.365*r[7]+3409.56000024852*r[6]-2130.97499950296*r[5]-2557.16999975148*r[4])*ct6+(-6238.15500*r[7]-2856.*r[8]-871.5*r[9]+9750.82500*r[5]-1069.32*r[6]-420.*r[2]
    -598.5*r[3]+6387.99000*r[4])*ct4+(22.32*r[11]+153.672*r[10]+1412.832*r[9]+443.064*r[2]+8.16*r[1]+2.148*r[0]+374.832000*r[3]-5129.598*r[4]-3774.36*r[6]-9887.085*r[5]
    +3510.87900*r[7]+3747.01200*r[8])*ct2+2267.23500*r[5]+223.668000*r[3]-8.1*r[1]-22.32*r[11]-153.672000*r[10]+1434.12000*r[6]-256.089000*r[7]-891.012000*r[8]-541.332000*r[9]
    +1298.77800*r[4]-23.0640000*r[2]-2.148*r[0])*c2p+((779.624999850887*r[5]+935.549999925444*r[4]-1247.47456*r[6]-1091.475*r[7])*ct6+(2541.73500*r[7]-215.04*r[8]-96.6*r[9]
    -575.925000*r[5]+3750.60000*r[6]-114.240000*r[2]-303.240000*r[3]-1855.35000*r[4])*ct4+(193.2*r[9]+228.48*r[2]+606.48*r[3]+904.050000*r[4]-3759.*r[6]-1187.02500*r[5]
   -1809.04500*r[7]+430.080000*r[8])*ct2+983.325000*r[5]-303.240000*r[3]+1255.80000*r[6]+358.785000*r[7]-215.040000*r[8]-96.6*r[9]+15.75*r[4]-114.240000*r[2])*c4p
    +((-51.9749999857988*r[5]-62.37*r[4]+83.1671006*r[6]+72.765*r[7])*ct6+(155.924999950296*r[5]+187.109999975148*r[4]-249.4824852*r[6]-218.295*r[7])*ct4
    +(-155.924999950296*r[5]-187.109999975148*r[4]+249.480000024852*r[6]+218.295*r[7])*ct2+62.37*r[4]+51.9749999857988*r[5]-83.1671006*r[6]-72.765*r[7])*c6p+(422.946563*r[5]
    +507.535875*r[4]-676.714500*r[6]-592.125187*r[7])*ct6+(2129.13531200001*r[7]+2235.73875*r[8]+678.365625*r[9]-5560.98593800001*r[5]-2262.03250*r[6]+318.4825*r[2]+436.323125*r[3]
   -2609.35062500001*r[4])*ct4+(-22.9635*r[11]-183.438000*r[10]-1201.32075*r[9]-84.8430000*r[2]+26.7375000*r[1]+.723*r[0]+232.391250*r[3]+2921.918625*r[4]+2348.83650*r[6]
    +5343.82968799999*r[5]-2401.869562*r[7]-3026.39250*r[8])*ct2+3965.97068837264*r[5]+1405.55462512918*r[3]+172.826500010249*r[1]+.853*r[13]+18.5543451*r[12]+155.663500006809*r[11]
    +711.734000048858*r[10]+5404.38250056856*r[6]+5671.43843862905*r[7]+4187.53575043265*r[8]+2105.23212518211*r[9]+2535.70212524520*r[4]+588.708500043227*r[2]+34.0170000010054*r[0]
    +3.095;
  //printf("a dPr %f\n",a);
  return a*invr16;
}

double dPt(double r[],double st,double ct2,double ct4,double ct6,double c2p,double c4p,double c6p,double invr13){
  double a;
  a=14.88*r[1]*c6p*((171.8528226*ct4*r[3]+(28.22580645*r[1]-46.85483871*r[5]-112.3387097*r[4]-232.9475807*r[3]+37.82258065*r[2])*ct2-.1443548387-16.04435484*r[2]+r[7]
    +42.61209678*r[5]+6.163709678*r[6]+85.78306452*r[4]+81.02620968*r[3]-16.31209678*r[1]-.7024193549*r[0])*c2p+(-62.87298388*ct4*r[3]+(7.677419355*r[1]+95.09274194*r[3]
    +17.83870968*r[2]-5.193548387*r[5]-7.903225807*r[4])*ct2+7.903225807*r[4]+5.193548387*r[5]-7.677419355*r[1]-32.21975807*r[3]-17.83870968*r[2])*c4p+(4.191532258*ct4*r[3]
    -8.383064517*ct2*r[3]+4.191532258*r[3])*c6p-34.10859375*ct4*r[3]+(-27.73891129*r[2]-21.40339382*r[1]+76.21144154*r[3]+88.01041667*r[4]+36.47127016*r[5])*ct2-37.23014113*r[5]
    -61.09740424*r[3]+1.296673387*r[1]-7.192741936*r[6]-1.028830645*r[7]-73.88326613*r[4]-10.55302419*r[2]-1.311290323*r[0]-0.4858870968e-1)*st;
  return a*invr13;
}

double dPp(double r[], double ct,double ct2,double ct4,double ct6,double s2p,double s4p,double s6p,double invr13){
  double a;
  a=(14.88*((57.28427419*ct6*r[3]+(18.91129032*r[2]+14.11290322*r[1]-23.42741935*r[5]-56.16935483*r[4]-116.4737903*r[3])*ct4+(r[7]+6.16370967761637*r[6]+42.6120967742612*r[5]
    +85.7830645185783*r[4]+81.0262096781026*r[3]-16.0443548416044*r[2]-16.3120967716312*r[1]-.702419354770242*r[0]-.144354838714435)*ct2+.144354838714435-1.*r[7]-6.16370967761637*r[6]
    -29.6137096729614*r[4]-21.8366935521837*r[3]-2.86693548428669*r[2]+2.19919354821992*r[1]+.702419354770242*r[0]-19.1846774219185*r[5])*s2p+(-41.91532258*ct6*r[3]+(95.09274192*r[3]
    +17.83870968*r[2]+7.677419354*r[1]-5.193548386*r[5]-7.903225805*r[4])*ct4+(-15.35483871*r[1]+10.38709677*r[5]+15.80645161*r[4]-64.43951612*r[3]-35.67741935*r[2])*ct2+17.83870968*r[2]
    -5.193548386*r[5]-7.903225805*r[4]+11.26209677*r[3]+7.677419354*r[1])*s4p+(12.57459677*ct2*r[3]-12.57459677*ct4*r[3]+4.191532258*ct6*r[3]-4.191532258*r[3])*s6p))*r[1];
  return a*invr13;
}


double P(double c2,double c3,double b20,double b22,double x,double y,double z,double r){
  return c2+c3*Y(b20,b22,x,y,z,r);
}

double dP_dx(double c3,double b20,double b22,double x,double y,double z,double r){
  return c3*dY_dx(b20,b22,x,y,z,r);    
}

double dP_dy(double c3,double b20,double b22,double x,double y,double z,double r){
  return c3*dY_dy(b20,b22,x,y,z,r);
}

double dP_dz(double c3,double b20,double b22,double x,double y,double z,double r){
  
  return c3*dY_dz(b20,b22,x,y,z,r);
}


double Y(double b20,double b22,double x,double y,double z,double r){
  return -b20/2.*(3.*f1(x,y,z)-1.)-3.*b22*(f1(x,y,z)-1.)*f2(x,y);
}

double dY_dx(double b20,double b22,double x,double y,double z,double r){
  return -3.*(df1_dx(x,y,z,r)*(b20/2.+b22*f2(x,y))+b22*(f1(x,y,z)-1.)*df2_dx(x,y));
}

double dY_dy(double b20,double b22,double x,double y,double z,double r){
  return -3.*(df1_dy(x,y,z,r)*(b20/2.+b22*f2(x,y))+b22*(f1(x,y,z)-1.)*df2_dy(x,y));
}

double dY_dz(double b20,double b22,double x,double y,double z,double r){
  return -3.*df1_dz(x,y,z,r)*(b20/2.+b22*f2(x,y));
}


double f1(double x,double y,double z){
  
  if(z==0.) return 0.;
  else return pow(cos(atan2(sqrt(x*x+y*y),z)),2.);
  
}

double df1_dx(double x,double y,double z,double r){
  if(z==0.) return 0.;
  else return -2.*x*z*z/pow(r,4.);
}

double df1_dy(double x,double y,double z,double r){
  if(z==0.) return 0.;
  else return -2.*y*z*z/pow(r,4.);
}

double df1_dz(double x,double y,double z,double r){
  if(z==0.) return 0.;
  else return 2.*z*(x*x+y*y)/pow(r,4.);
}


double f2(double x,double y){
  //if(x==0.) return -1.;
  return cos(2.*atan2(y,x));
}

double df2_dx(double x,double y){
  //if(x==0.) return 0.;
  return 2.*y*sin(2.*atan2(y,x))/(x*x+y*y);
}

double df2_dy(double x,double y){
  //if(x==0.) return 0.;
  return -2.*x*sin(2.*atan2(y,x))/(x*x+y*y);
}




