// SimulatePct.cpp in /pctsimulatorTest4

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <string.h>
#include <time.h>
#include <fstream>
#include <sstream>
using namespace std;
float PI=3.1415926; 

class ConfigObj {
public:
  ConfigObj(bool d3,int xp,int yp,int zp,float bw,float bl,
   float bh,int h,int tpa,int sa,int as,string sfn,int po1,
   int po2) : 
   dimensional3(d3),xpartitions(xp),ypartitions(yp),
   zpartitions(zp),boxwidth(bw),boxlength(bl),
   boxheight(bh),histories(h),totalprojectionangles(tpa),
   startangle(sa),anglespace(as),simulationfoldername(sfn),
   pathoption(po1),phantomoption(po2) {};
  bool dimensional3; 
  int xpartitions;
  int ypartitions;
  int zpartitions;
  float boxwidth;
  float boxlength;
  float boxheight;
  int histories;
  int totalprojectionangles;
  int startangle;
  int anglespace;
  string simulationfoldername;
  int pathoption;
  int phantomoption;
};
ConfigObj ReadConfigObjFromFile() {
  fstream infile;
  infile.open("configfile",ios::in);
  bool d3;
  int xp,yp,zp,h,tpa,sa,as,po1,po2;
  float bw,bl,bh;
  string sfn,tempS; 
  infile>>tempS>>d3;
  infile>>tempS>>xp;
  infile>>tempS>>yp;
  infile>>tempS>>zp;
  infile>>tempS>>bw;
  infile>>tempS>>bl;
  infile>>tempS>>bh;  
  infile>>tempS>>h;
  infile>>tempS>>tpa;
  infile>>tempS>>sa;
  infile>>tempS>>as;
  infile>>tempS>>sfn;
  infile>>tempS>>po1;
  infile>>tempS>>po2;
  infile.close();
  return ConfigObj(d3,xp,yp,zp,bw,bl,bh,h,tpa,sa,as,sfn,po1,
   po2);
}

string ReadEntryStringFromConfigFile(int entryId) {
  fstream infile;
  infile.open("configfile",ios::in);
  string temp1,temp2;
  for (int i=0; i<entryId; i++) 
    infile>>temp1>>temp2; //move to correct entry
  infile>>temp1>>temp2;
  return temp2; 
}

/***********************************PHANTOM******************************/

class Ellipse {
 public:
  Ellipse(float xxc,float yyc,float aa,float bb) :
    xc(xxc),yc(yyc),a(aa),b(bb) {}; 
  float xc;
  float yc;
  float a;
  float b;
};
// NonuniformEllipseObject phantom:
// a vector of ellipses for a convex object;
// the last element of v must be the boundary ellipse;
// if an ellipse is inside another ellipse than the 
// smaller one (on the inside) must come before the 
// bigger one in the E vector
class neo {
public:
  neo(vector<Ellipse> e,vector<float> r) : E(e),rsp(r) {};
  vector<Ellipse> E;
  vector<float> rsp;
};
neo CreateDefaultNeo(float s) {
  vector<Ellipse> E;
  vector<float> rsp;
  Ellipse e1(s*0,s*0,s*70,s*90);
  Ellipse e2(s*0,s*85,s*10,s*2.5);
  Ellipse e3(s*0,s*0,s*60,s*80);  
  Ellipse e4(s*20,s*0,s*10,s*20);  
  Ellipse e5(s*-20,s*0,s*10,s*20);
  E.push_back(e5);  
  E.push_back(e4);
  E.push_back(e3);
  E.push_back(e2);
  E.push_back(e1);
  rsp.push_back(0.9); // ventricles
  rsp.push_back(0.9); // ventricles
  rsp.push_back(1.04); // brain tissue
  rsp.push_back(0.0); // air pocket 
  rsp.push_back(1.6); // skull  
  return neo(E,rsp);  
}
class LPoint {
 public:
  LPoint() {
    x=0;
    y=0;
    l=0;   
  }
  LPoint(float xx,float yy,int ll) :
    x(xx),y(yy),l(ll) {};
  float x;
  float y;
  int l;
};
float calculateTriangleArea(LPoint A,LPoint B,
 LPoint C) {
  float output=fabs((A.x*(B.y-C.y)+B.x*(C.y-A.y)
   +C.x*(A.y-B.y))/2.0);
  return output;
}
bool equalf(float a,float b) {
  float tol=0.00001; 
  if (a<b+tol && a>b-tol)
    return true;
  else
    return false;
}
bool IsPointInEllipse(float xc,float yc,float a,float b,
 float k,float h) {
  float yup=h+b*sqrt(1.0-pow(xc-k,2)/pow(a,2));
  float ylo=h-b*sqrt(1.0-pow(xc-k,2)/pow(a,2));
  if (yc<=yup && yc>=ylo) 
    return true;
  else 
    return false;
}
float FindRSPatPoint(float x,float y,neo phant) {
  float RSP=0;
  float xcPhan,ycPhan,aPhan,bPhan;
  for (int i=0; i<phant.E.size(); i++) {
    xcPhan=phant.E[i].xc;
    ycPhan=phant.E[i].yc;
    aPhan=phant.E[i].a;
    bPhan=phant.E[i].b;
    if (IsPointInEllipse(x,y,aPhan,bPhan,xcPhan,ycPhan)) {
      RSP=phant.rsp[i];
      return RSP;
    }
  }
  return RSP;
}
float FindRSPatPoint_DefaultNeo(float xc,float yc) {
  float RSP;
  if (IsPointInEllipse(xc,yc,10,20,-20,0) || 
   IsPointInEllipse(xc,yc,10,20,20,0))
    RSP=0.9;
  else if (IsPointInEllipse(xc,yc,60,80,0,0)) {
    RSP=1.04;
  }
  else if (IsPointInEllipse(xc,yc,10,2.5,0,85))
    RSP=0;
  else if (IsPointInEllipse(xc,yc,70,90,0,0))
    RSP=1.6;
  else 
    RSP=0;
  return RSP;
}
// method==0 --> centerPointMethod      vs
// method==1 --> Corner Point Averaging
float* CreateNeoSlice(neo phantom,int ypartitions,int xpartitions,
  float boxLength,float boxWidth,int method) {
  int sizeXtrue=xpartitions*ypartitions*sizeof(float);
  float *Xtrue=(float*)malloc(sizeXtrue);
  float yshift=boxLength/2.0;
  float xshift=boxWidth/2.0;
  float ystep=boxLength/ypartitions;
  float xstep=boxWidth/xpartitions;
  float x0,x1,x2,x3,x4,y0,y1,y2,y3,y4;
  float xc,yc,RSP,temp0,temp1,temp2,temp3,temp4;
  float xdot1l0,xdot2l0,xdot1l2,xdot2l2,
        ydot1l1,ydot2l1,ydot1l3,ydot2l3;
  float base,rbase,height,rheight,area1,area2;
  float a,b;
  float voxelArea=xstep*ystep;
  int counter1=0;
  if (method==0) { //centerPointMethod
    for (int i=0; i<ypartitions; i++) {
      for (int j=0; j<xpartitions; j++) {
        x1=j*xstep-xshift;
        x2=x1+xstep;
        y1=yshift-i*ystep;
        y2=y1-ystep;
        xc=(x1+x2)/2.0;
        yc=(y1+y2)/2.0;
        RSP=FindRSPatPoint(xc,yc,phantom);
        Xtrue[i*xpartitions+j]=RSP;
      }
    }
  }
  else if (method==1) { //corner point averaging method
    for (int i=0; i<ypartitions; i++) {
      for (int j=0; j<xpartitions; j++) {
        x1=j*xstep-xshift; 
        y1=yshift-i*ystep;
        x2=x1+xstep; 
        y2=y1;
        x3=x2; 
        y3=y1-ystep;
        x4=x1; 
        y4=y3;
        temp1=FindRSPatPoint(x1,y1,phantom);
        temp2=FindRSPatPoint(x2,y2,phantom);
        temp3=FindRSPatPoint(x3,y3,phantom);
        temp4=FindRSPatPoint(x4,y4,phantom);
        RSP=(temp1+temp2+temp3+temp4)/4.0;
        Xtrue[i*xpartitions+j]=RSP; 
      }      
    }
  }
  else { //weighted area method (needs more testing)
    for (int i=0; i<ypartitions; i++) {
      for (int j=0; j<xpartitions; j++) {
        cout<< i << " " << j << endl;
        x0=j*xstep-xshift; 
        y0=yshift-i*ystep;
        x1=x0+xstep; 
        y1=y0;
        x2=x1; 
        y2=y0-ystep;
        x3=x0; 
        y3=y2;
        temp0=FindRSPatPoint(x0,y0,phantom);
        temp1=FindRSPatPoint(x1,y1,phantom);
        temp2=FindRSPatPoint(x2,y2,phantom);
        temp3=FindRSPatPoint(x3,y3,phantom);
        if (equalf(temp0,temp1) && equalf(temp0,temp2) 
         && equalf(temp0,temp3)) {
          RSP=temp0;
          Xtrue[i*xpartitions+j]=RSP;
        }
        else {
          vector<LPoint> iP; //P is a temporary vector  
          // to hold the intersection points of the 
          // ellipses and the voxels.  It should be
          // of size 2, unless more than 1 ellipse 
          // intersect a single voxel, however we'll 
          // keep the voxel size small enough to avoid
          // this.  So if the size is not 2, use center 
          // point method                              
          for (int k=0; k<phantom.E.size(); k++) {
            xc=phantom.E[k].xc;
            yc=phantom.E[k].yc;
            a=phantom.E[k].a;
            b=phantom.E[k].b;
            xdot1l0=xc+a*sqrt(1-pow(y0-yc,2)/pow(b,2));    
            xdot2l0=xc-a*sqrt(1-pow(y0-yc,2)/pow(b,2));    
            xdot1l2=xc+a*sqrt(1-pow(y2-yc,2)/pow(b,2));    
            xdot2l2=xc-a*sqrt(1-pow(y2-yc,2)/pow(b,2));    
            ydot1l1=yc+b*sqrt(1-pow(x1-xc,2)/pow(a,2));    
            ydot2l1=xc-b*sqrt(1-pow(x1-xc,2)/pow(a,2));    
            ydot1l3=xc+b*sqrt(1-pow(x3-xc,2)/pow(a,2));    
            ydot2l3=xc-b*sqrt(1-pow(x3-xc,2)/pow(a,2));    
            if (xdot1l0<x1 && xdot1l0>x0 && 
             1-pow(y0-yc,2)/pow(b,2)>0)
              iP.push_back(LPoint(xdot1l0,y0,0));
            if (xdot2l0<x1 && xdot2l0>x0 && 
             1-pow(y0-yc,2)/pow(b,2)>0)
              iP.push_back(LPoint(xdot2l0,y0,0));
            if (xdot1l2<x1 && xdot1l2>x0 && 
             1-pow(y2-yc,2)/pow(b,2)>0)
              iP.push_back(LPoint(xdot1l2,y2,2));
            if (xdot2l2<x1 && xdot2l2>x0 && 
             1-pow(y2-yc,2)/pow(b,2)>0)
              iP.push_back(LPoint(xdot2l2,y2,2));
            if (ydot1l1<y0 && ydot1l1>y3 && 
             1-pow(x1-xc,2)/pow(a,2)>0)
              iP.push_back(LPoint(x1,ydot1l1,1));
            if (ydot2l1<y0 && ydot2l1>y3 && 
             1-pow(x1-xc,2)/pow(a,2)>0)
              iP.push_back(LPoint(x1,ydot2l1,1));
            if (ydot1l3<y0 && ydot1l3>y3 && 
             1-pow(x3-xc,2)/pow(a,2)>0)
              iP.push_back(LPoint(x0,ydot1l3,3));
            if (ydot2l3<y0 && ydot2l3>y3 && 
             1-pow(x3-xc,2)/pow(a,2)>0)
              iP.push_back(LPoint(x0,ydot2l3,3));
          }
          if (iP.size()!=2) { //center-point-method
            xc=(x0+x1)/(float)2;
            yc=(y0+y3)/(float)2;
            RSP=FindRSPatPoint(xc,yc,phantom);
            Xtrue[i*xpartitions+j]=RSP;
            //cout<<i<<" "<<j<<endl;            
          }          
          else {
            vector<LPoint> cP; //corner points
            cP.push_back(LPoint(x0,y0,-1));
            cP.push_back(LPoint(x1,y1,-1));
            cP.push_back(LPoint(x2,y2,-1));
            cP.push_back(LPoint(x3,y3,-1));            
            vector<LPoint> sP; //shape points
            if (iP[0].l>iP[1].l) {
              LPoint tempLPoint=iP[0];
              iP[0]=iP[1];
              iP[1]=tempLPoint; 
            }
            int la=iP[0].l;
            int lb=iP[1].l;
            sP.push_back(iP[0]);
            for (int ind=la+1; ind<=lb; ind++)
              sP.push_back(cP[ind]);
            sP.push_back(iP[1]);
            if (sP.size()<3) { //center-point-method
              xc=(x0+x1)/(float)2;
              yc=(y0+y3)/(float)2;
              RSP=FindRSPatPoint(xc,yc,phantom);
              Xtrue[i*xpartitions+j]=RSP;                          
              counter1++; 
            }                   
            else {     
              float area1=0;
              int tempIndex=1;               
              for (int s=0; s<sP.size()-2; s++) {
                LPoint tempA=sP[0];
                LPoint tempB=sP[tempIndex];
                LPoint tempC=sP[tempIndex+1];              
                area1+=calculateTriangleArea(tempA,tempB,tempC);
                tempIndex++;
              }              
              area2=voxelArea-area1;
              RSP=(1.0/voxelArea)*(area1*temp2+area2*temp4);
              Xtrue[i*xpartitions+j]=RSP;
            }
          }
        }
      }
    }
  }
  return Xtrue;
}
float* CreateXtrueNeo3d(neo MyNeo,int zparts,int yparts,int xparts,
 float boxLength,float boxWidth,int method) {
  float* Neo2d=CreateNeoSlice(MyNeo,yparts,xparts,boxLength,boxWidth,method);
  float* XtrueNeo3d=(float*)malloc(zparts*yparts*xparts*sizeof(float));
  for (int k=0; k<zparts; k++) {
    for (int ii=0; ii<yparts*xparts; ii++) {
      XtrueNeo3d[k*yparts*xparts+ii]=Neo2d[ii];
    }
  }
  free(Neo2d);
  return XtrueNeo3d;
}
void WriteNeoToFile(float* Xtrue,int zs,int ys,int xs,string filename) {
  fstream outfile;
  outfile.open(filename.c_str(),ios::out);
  outfile<<zs<<' '<<ys<<' '<<xs<<endl;
  for (int ii=0; ii<zs*ys*xs; ii++)
    outfile<<Xtrue[ii]<<' ';
  outfile.close();
}

/************************************************************************/

/***********************************PROTONPATH*******************************/

class PointVec {
 public:
  PointVec(float xx,float yy) : x(xx), y(yy) {};
  PointVec() {};
  void set(float xx,float yy) {
    x=xx;
    y=yy;
  }
  float x;
  float y;
};

class PointVec3 {
public:
  PointVec3(float xx,float yy,float zz) : x(xx),y(yy),z(zz) {};
  float x;
  float y;
  float z;
};

float normEu(PointVec P1,PointVec P2) {
  return sqrt(pow(P2.x-P1.x,2)+pow(P2.y-P1.y,2));
}

// each PointVec3  (varLat, cov, VarAng)
vector<PointVec3> scatteringParameters() {
  vector<PointVec3> output;
  output.push_back(PointVec3(0,0,0));
  output.push_back(PointVec3(0.00112,  0.0001686, 3.397*pow(10,-5)));
  output.push_back(PointVec3(0.009335, 0.0007052, 7.154*pow(10,-5)));  
  output.push_back(PointVec3(0.0324,   0.001638,  0.0001117));
  output.push_back(PointVec3(0.07861,  0.002992,  0.0001542));    
  output.push_back(PointVec3(0.1567,   0.004793,  0.0001994));
  output.push_back(PointVec3(0.2761,   0.007067,  0.0002472));        
  output.push_back(PointVec3(0.4466,   0.009843,  0.0002979));        
  output.push_back(PointVec3(0.6786,   0.01315,   0.0003519));   
  output.push_back(PointVec3(0.9833,   0.01703,   0.0004094));             
  output.push_back(PointVec3(1.372,    0.0215,    0.0004709));             
  output.push_back(PointVec3(1.859,    0.02663,   0.0005368));             
  output.push_back(PointVec3(2.456,    0.03245,   0.0006078));             
  output.push_back(PointVec3(3.178,    0.03902,   0.0006847));            
  output.push_back(PointVec3(4.041,    0.0464,    0.0007683));             
  output.push_back(PointVec3(5.063,    0.05467,   0.0008599));             
  output.push_back(PointVec3(6.261,    0.06392,   0.0009611));            
  output.push_back(PointVec3(7.658,    0.07425,   0.001074));            
  output.push_back(PointVec3(9.275,    0.08579,   0.001201));            
  output.push_back(PointVec3(11.14,    0.09871,   0.001347));             
  output.push_back(PointVec3(13.28,    0.1132,    0.001518));
  return output;             
}

// return uniform random number from [0,1]
float ranf() {
  return rand()%1000001/(float)1000000;
}

// return uniform rv from [-125,125]
float MyUniformRand() {
  float rv=rand()%1000001/(float)1000000;
  rv=rv*250.0; //scale rv to be in [0,250]
  rv=rv-125.0; //shift rv to be in [-125,125]
  return rv;
}

//for vertical displacement
float MyUniformRand2() {
  float rv=rand()%1000001/(float)1000000;
  rv=rv*100.0; //scale rv to be in [0,250]
  rv=rv-50.0; //shift rv to be in [-125,125]
  return rv;
}

// Marsaglia Polar Method ~ modified from Box-Muller
// to take 2 uniform rvs and transform them into 
// 2 independent standard normal rvs
PointVec generate2RandStdNorm() {
  float u,v,x,y,s,t;
  s=1.0;
  while (s>=1.0 or equalf(s,0)) {
    u=2.0*ranf()-1.0;
    v=2.0*ranf()-1.0;
    s=u*u+v*v;
  }
  t=sqrt((-2.0*log(s))/s);
  x=u*t;
  y=v*t;
  return PointVec(x,y);
}

// hard coded calculation of determinant
float Determinant3x3(float* A) {
  float d=A[0]*(A[4]*A[8]+A[5]*A[7]) 
         -A[1]*(A[3]*A[8]+A[5]*A[6])
         +A[2]*(A[3]*A[7]+A[4]*A[6]);
  return d;
}

// solve linear 3x3 system of equations Ax=b
float* CramerSolve3x3(float* A,float* b) {
  float* output=(float*)malloc(3*sizeof(float));
  float detA=Determinant3x3(A);
  if (detA==0) {
    cout<<"error: zero determinant for J"<<endl;
    output[0]=0;
    output[1]=0;
    output[2]=0;   
    return output;
  }
  float* tempM=(float*)malloc(9*sizeof(float));
  for (int i=0; i<9; i++)
    tempM[i]=A[i];
  //calculate the x of x=(x,y,z)
  tempM[0]=b[0];
  tempM[3]=b[1];
  tempM[6]=b[2];
  float x=Determinant3x3(tempM)/detA;
  tempM[0]=A[0];
  tempM[3]=A[3];
  tempM[6]=A[6];
  //calculate y
  tempM[1]=b[0];
  tempM[4]=b[1];
  tempM[7]=b[2];
  float y=Determinant3x3(tempM)/detA;
  tempM[1]=A[1];
  tempM[4]=A[4];
  tempM[7]=A[7];
  //calculate z
  tempM[2]=b[0];
  tempM[5]=b[1];
  tempM[8]=b[2];
  float z=Determinant3x3(tempM)/detA;
  output[0]=x;
  output[1]=y;
  output[2]=z;
  return output;
}

// F is vector function for nonlinear system 
// of equations, this function returns -F 
// which is needed in newton method
// v = (var1,cov,var2) 
float* MyF(float*x, float* v) {
  float* F=(float*)malloc(3*sizeof(float));
  F[0]=-(x[0]*x[0]+x[1]*x[1]-v[0]);
  F[1]=-(x[1]*x[1]+x[2]*x[2]-v[2]);
  F[2]=-(x[0]*x[1]+x[1]*x[2]-v[1]);
  return F;
}

// hard coded calculation of the jacobian 
// for the nonlinear system of equation to
// solve for the entries in matrix to 
// multiply to the vector of std norms to 
// get the joint norms
float* MyJacobian(float* x) {
  float* J=(float*)malloc(9*sizeof(float));
  J[0]=2*x[0]; 
  J[1]=2*x[1]; 
  J[2]=0;
  J[3]=0; 
  J[4]=2*x[1];
  J[5]=2*x[2];
  J[6]=x[1];
  J[7]=x[0]+x[2];
  J[8]=x[1];
  return J; 
}

// J~jacobian of F, v is needed to calculate F
// x~solution vector
float* NewtonNonlinearSysMethod(float* v) {
  float *J,*F,*s; 
  float* x=(float*)malloc(3*sizeof(float));
  float* xprevious=(float*)malloc(3*sizeof(float));
  x[0]=10; x[1]=0.1; x[2]=0.01; //intitial x guess 
  for (int k=0; k<20; k++) { //20 iterations
    J=MyJacobian(x); // 3x3 jacobian
    F=MyF(x,v); //returns 3x1 vector
    //solve Js=-f for s using Cramer
    s=CramerSolve3x3(J,F); 
    x[0]=x[0]+s[0];
    x[1]=x[1]+s[1];
    x[2]=x[2]+s[2];  
    free(J); free(F); free(s);     
  }  
  return x;
}

PointVec generateJointRandNorm(float var1,float cov,float var2) {
  PointVec U=generate2RandStdNorm();;
  float* v=(float*)malloc(3*sizeof(float)); 
  v[0]=var1; v[1]=cov; v[2]=var2;
  float* temp=NewtonNonlinearSysMethod(v);
  PointVec X;
  X.x=temp[0]*U.x+temp[1]*U.y;
  X.y=temp[1]*U.x+temp[2]*U.y;
  free(v); free(temp);
  return X;
}

// m~slopeLine,b~yinterceptLine,xa~lowerXboundBox,
// ya~lowerYboundBox,xb~UpperXboundBox,yb~upperYboundBox
// function returns the 2 points of interesection 
// of the line and the box
vector<PointVec> LineBoxHits(float m,float b,float xa,
 float xb,float ya,float yb) {
  // calculate intersection of proton entry line with 
  // reconstruction box
  float xdot,ydot;
  vector<PointVec> v; //entry and exit of the box
  xdot=(1.0/m)*(yb-b);
  if (xa<=xdot and xdot<=xb) //hit top line
    v.push_back(PointVec(xdot,yb));      
  xdot=(1.0/m)*(ya-b);
  if (xa<=xdot and xdot<=xb) //hit bottom line
    v.push_back(PointVec(xdot,ya));      
  ydot=m*xb+b;
  if (ya<=ydot and ydot<=yb) //hit right line
    v.push_back(PointVec(xb,ydot));
  ydot=m*xa+b;
  if (ya<=ydot and ydot<=yb) //hit left line
    v.push_back(PointVec(xa,ydot));
  if (v.size()!=2) 
    cout<<"number of interesections="<<v.size()<<endl;
  return v;
}
// Point (x0,y0) is for the upper left corner of each pixel
float CalculateX0(float x,float xstep) {
  float tempDivide=x/xstep;
  float remainPercent=tempDivide-floor(tempDivide);
  float xmove=remainPercent*xstep;
  float x0=x-xmove;   
  return x0;
}
float CalculateY0(float y,float ystep) {
  float tempDivide=y/ystep;
  float remainPercent=tempDivide-floor(tempDivide);
  float ymove=(1-remainPercent)*ystep;
  float y0=y+ymove;     
  return y0;
}

class History2d {
public:
  History2d(vector<int> c,float w) 
   : colInd(c),wepl(w) {};
  vector<int> colInd;
  float wepl;
};

// KnownHull indicates that we intersected the proton path with the ellipse boundary of the neo phantom
History2d* generateProtonPath2d_KnownHull(float theta,float* Xtrue,int yparts,int xparts,float boxLength,float boxWidth,neo MyNeo,int pathOption,int sliceIndex) {
  int startCol=sliceIndex*yparts*xparts;
  vector<int> colInd,nullvec;
  float wepl=0;
  int col,localCol;
  PointVec Pin,Pt,Pout;
  float thetaOut;
  float a=MyNeo.E[MyNeo.E.size()-1].a; // semi major xaxis length of outer ellipse
  float b=MyNeo.E[MyNeo.E.size()-1].b; // semi minor yaxis ...
  float ystep=boxLength/(float)yparts;
  float xstep=boxWidth/(float)xparts;
  // the shifts to center the recon space 
  float yshift=boxLength/2; 
  float xshift=boxWidth/2;
  float gantryRad=3000;
  vector<PointVec3> scatParam=scatteringParameters();
  PointVec P1(gantryRad*cos(theta),gantryRad*sin(theta));
  PointVec v1(cos(theta),sin(theta));
  PointVec v2(-sin(theta),cos(theta));
  float varLat,covariance,varAng;
  float d1=MyUniformRand();  
  //float d1=0; /**************random********/
  //calculate intersection of proton entry line with ellipse
  PointVec P3(P1.x+d1*v2.x,P1.y+d1*v2.y);
  float m=v1.y/v1.x; //slope of proton entry line
  float B=P3.y-m*P3.x; //y-intercept of proton entry line
  //to calculate the entry point Pin:
  //the intersection test of line and ellipse --> solve for x in
  //the quadratic equation:  alpha*x^2+beta*x+gamma=0
  float alpha=1/pow(a,2)+pow(m/b,2); 
  float beta=2*m*B/pow(b,2);
  float gamma=pow(B/b,2)-1;
  float disc=pow(beta,2)-4*alpha*gamma; //discriminant
  if ( disc<0 or equalf(disc,0) ) { 
    // proton path missed the object -> return blank history
    History2d* T=new History2d(nullvec,0);
    return T;
  }
  else {
    float x1=(-beta+sqrt(disc))/(2*alpha);
    float x2=(-beta-sqrt(disc))/(2*alpha);
    float y1=m*x1+B;
    float y2=m*x2+B;
    float y11=b*sqrt(1-pow(x1/a,2));
    float y22=-b*sqrt(1-pow(x2/a,2));
    float norm1=normEu(PointVec(x1,y1),P3);
    float norm2=normEu(PointVec(x2,y2),P3);
    if (norm1<norm2) { 
      Pin.set(x1,y1);
      Pt.set(x2,y2); 
    }
    else {
      Pin.set(x2,y2);
      Pt.set(x1,y1);
    }
    //depth is the distance proton traveled through object 
    //assuming it moves in a straight line
    float depth=normEu(Pin,Pt);
    //the lookup table is indexed by cm so divide by 10
    //to convert mm to cm
    int index=ceil(depth/10);
    if (index>20) index=20; //avoid seg faults
    varLat=scatParam[index].x;
    covariance=scatParam[index].y;
    varAng=scatParam[index].z;
    PointVec R2=generateJointRandNorm(varLat,covariance,varAng);
    float d2=R2.x; //exiting lateral displacement 
    float psi=R2.y; //exiting angular displacement
    //generate temporary line to intersect with ellipse to find Pout
    PointVec P4(Pt.x+d2*v2.x,Pt.y+d2*v2.y);
    float m2=m; //slope
    float B2=P4.y-m2*P4.x; //intercept
    //line-ellipse interesection quadratic coefficients
    float alpha2=1/pow(a,2)+pow(m2/b,2); 
    float beta2=2*m2*B2/pow(b,2);
    float gamma2=pow(B2/b,2)-1;
    float disc2=pow(beta2,2)-4*alpha2*gamma2; //discriminant
    if (disc2<0) { // very rare circumstance 
      History2d* T=new History2d(nullvec,0);
      return T;
    }
    else if (equalf(disc2,0)) {
      float x21=-beta2/(2*alpha2);
      float y21=m2*x21+B2;
      Pout.set(x21,y21);
    }
    else {
      float x21=(-beta2+sqrt(disc2))/(2*alpha2);
      float x22=(-beta2-sqrt(disc2))/(2*alpha2);
      float y21=m2*x21+B2;
      float y22=m2*x22+B2;
      float norm21=normEu(PointVec(x21,y21),P4);
      float norm22=normEu(PointVec(x22,y22),P4);
      if (norm21<norm22) 
        Pout.set(x21,y21);        
      else 
        Pout.set(x22,y22);      
    }
    //exiting angle relative to xy axis
    thetaOut=theta+psi; 
    //calculate upper left corner coordinates of the entry voxel
    //and exit voxel to find xmin and ymin
    float xin0=CalculateX0(Pin.x,xstep);
    float xout0=CalculateX0(Pout.x,xstep);
    float xmin=min(xin0,xout0);
    float xmax=max(xin0,xout0);
    if (pathOption==1) { //straight-line between Pin and Pout
      float ms=(Pout.y-Pin.y)/(Pout.x-Pin.x); //slope
      float Bs=Pout.y-ms*Pout.x; //y-intercept
      for (float xx=xmin; xx<xmax; xx+=xstep) {
        float ys=ms*xx+Bs;
        float y0=CalculateY0(ys,ystep);
        int i=(int)((yshift-y0)/ystep);
        int j=(int)((xshift+xx)/xstep);
        //if the slope of the line is greater than 1 in abs val
        //then the line will intersect more voxels above i (if
        //positive slope) before the next xx so need to mark all
        //those voxels by calculating y at the next xx
        float ysnext=ms*(xx+xstep)+Bs;
        float y0next=CalculateY0(ysnext,ystep);
        int inext=(int)((yshift-y0next)/ystep);
        if (i>=0 and inext>=0) { 
          for (int ii=min(i,inext); ii<=max(i,inext); ii++) {
            localCol=ii*xparts+j;
            if (0<localCol and localCol<xparts*yparts) {
              col=startCol+localCol;
              colInd.push_back(col);
              wepl+=Xtrue[localCol];   
            }            
          }
        }
      } 
    }
    else if (pathOption=2) { //cubic-spline
      float thetaIn=theta;
      float c=thetaIn*(Pout.x-Pin.x)-(Pout.y-Pin.x);
      float d=-thetaOut*(Pout.x-Pin.x)-(Pout.y-Pin.y);
      for (float xx=xmin; xx<xmax; xx+=xstep) {
        float t=(xx-Pin.x)/(Pout.x-Pin.x);
        float q=(1-t)*Pin.y+t*Pout.y+t*(1-t)*(c*(1-t)+d*t);
        float y0=CalculateY0(q,ystep);
        int i=(int)((yshift-y0)/ystep);
        int j=(int)((xshift+xx)/xstep);
        float tnext=(xx+xstep-Pin.x)/(Pout.x-Pin.x);
        float qnext=(1-tnext)*Pin.y+tnext*Pout.y+tnext*
         (1-tnext)*(c*(1-tnext)+d*tnext);
        float y0next=CalculateY0(qnext,ystep);
        int inext=(int)((yshift-y0next)/ystep);
        if (i>=0 and inext>=0) { 
          for (int ii=min(i,inext); ii<=max(i,inext); ii++) {
            localCol=ii*xparts+j;
            if (0<=localCol and localCol<xparts*yparts) {
              col=startCol+localCol;
              colInd.push_back(col);
              wepl+=Xtrue[localCol];   
            }
          }
        }
      } 
    } 
  } 
  History2d* T=new History2d(colInd,wepl);
  return T;
}

// 2dHistories -> no proton scattering in the vertical direction
void SimulatePct_KnownHull_2dHistories(neo MyNeo,int histories,
 int zparts,int yparts,int xparts,float boxLength,
 float boxWidth,int totalProjAngles,float startAngle,
 float angleSpace,int pathOption,int phantomOption,
 string MatrixFolderName,string NeoXtrueFilename) {
  fstream outfileX;
  outfileX.open(NeoXtrueFilename.c_str(),ios::out); //clobber old
  float* Xtrue=CreateNeoSlice(MyNeo,yparts,xparts,boxLength,
   boxWidth,phantomOption);
  for (int i=0; i<yparts*xparts; i++) 
    outfileX<<Xtrue[i]<<' ';
  outfileX<<endl;
  outfileX.close();  
  int rowsPerAngle=ceil(histories/(float)totalProjAngles);
  int rowsPerAnglePerSlice=ceil(rowsPerAngle/(float)zparts);
  float endAngle=startAngle+totalProjAngles*angleSpace;
  for (float theta=startAngle; theta<endAngle; theta+=angleSpace) {
    fstream outfile;
    int tempTheta=(int)theta;
    ostringstream o;
    o<<tempTheta;
    string num=o.str();
    string filenameMat=MatrixFolderName+"AcsrB_"+num+".txt";
    outfile.open(filenameMat.c_str(),ios::out|ios::app);
    float thetaRad=theta*PI/180.0;
    for (int sliceId=0; sliceId<zparts; sliceId++ ) {
      for (int localRow=0; localRow<rowsPerAnglePerSlice; localRow++) {
        History2d* T=generateProtonPath2d_KnownHull(thetaRad,Xtrue,
         yparts,xparts,boxLength,boxWidth,MyNeo,pathOption,sliceId);
        if (T->colInd.size()!=0) {
          outfile<<tempTheta<<' ';
          outfile<<localRow<<' ';
          outfile<<T->wepl<<' ';
          outfile<<T->colInd.size()<<' ';
          for (int j=0; j<T->colInd.size(); j++)
            outfile<<T->colInd[j]<<' ';
          outfile<<endl;
          //if (localRow%100000==0) 
          //cout<<"angle="<<theta<<" || localRow="<<localRow<<endl;  
        }
        else {
          localRow--;
          //cout<<"localRow--"<<localRow<<endl;
        }
        delete T;
      }
    }
    cout<<"file finished for angle "<<theta<<endl;
    outfile.close();
  }
}

// File Format:
// rowIndex angle B hits colInd[0] colInd[1] ... colInd[hits-1] endl
// startAngle,angleSpace in degrees
void SimulatePct2d(neo MyNeo,int histories,int yparts,int xparts,
 float boxLength,float boxWidth,int totalProjAngles,float startAngle,
 float angleSpace,int pathOption,int phantomOption,
 string MatrixFilename,string NeoXtrueFilename) {
  int sliceIndex=0; // 2d sim has only 1 slice
  fstream outfileX;
  outfileX.open(NeoXtrueFilename.c_str(),ios::out); //clobber old
  float* Xtrue=CreateNeoSlice(MyNeo,yparts,xparts,boxLength,
   boxWidth,phantomOption);
  for (int i=0; i<yparts*xparts; i++) 
    outfileX<<Xtrue[i]<<' ';
  outfileX<<endl;
  outfileX.close();  
  fstream outfile;
  outfile.open(MatrixFilename.c_str(),ios::out|ios::app);
  int rowIndex=0;
  int rowsPerAngle=ceil(histories/(float)totalProjAngles);
  float endAngle=startAngle+totalProjAngles*angleSpace;
  for (float theta=startAngle; theta<endAngle; theta+=angleSpace) {
    float thetaRad=theta*PI/180.0;
    if (equalf(thetaRad,0.0) or equalf(thetaRad,PI) or 
     equalf(thetaRad,PI/2) or equalf(thetaRad,3/2*PI))
      thetaRad+=1*PI/180.0;  // to prevent division by zero (this will be fixed) 
    for (int localRow=0; localRow<rowsPerAngle; localRow++) {
      if (rowIndex<histories) {
        History2d* T=generateProtonPath2d_KnownHull(thetaRad,
         Xtrue,yparts,xparts,boxLength,boxWidth,MyNeo,
         pathOption,sliceIndex);
        if (T->colInd.size()!=0) {
          outfile<<rowIndex<<' '<<thetaRad*180/PI<<' ';
          outfile<<T->wepl<<' '<<T->colInd.size()<<' ';
          for (int j=0; j<T->colInd.size(); j++)
            outfile<<T->colInd[j]<<' ';
          outfile<<endl;
          rowIndex++; 
        }
        else {
          localRow--;  
        }
      }
    }
    cout<<"finished projecting virtual protons for angle ";
    cout<<thetaRad*180/PI<<endl;
  }
  outfile.close();
}

int main() {

  ConfigObj C1=ReadConfigObjFromFile();
  bool threeDim=C1.dimensional3; // 1
  int xparts=C1.xpartitions; // 160
  int yparts=C1.ypartitions;  // 200
  int zparts=C1.zpartitions; // 200
  float boxWidth=C1.boxwidth; // 160
  float boxLength=C1.boxlength; // 200
  float boxHeight=C1.boxheight; // 200
  int histories=C1.histories; // 64000000
  int voxels=C1.xpartitions*C1.ypartitions; //6400000
  int totalProjAngles=C1.totalprojectionangles; // 180
  int startAngle=C1.startangle; // 1
  int angleSpace=C1.anglespace; // 2
  string SimName=C1.simulationfoldername;
  int pathOption=C1.pathoption; // 0 (straight-line)
  int phantomOption=C1.phantomoption; // 1 (corner-points)
  string SimulationFolderName="SimulationData/"+SimName;
  string MatrixFolderName=SimulationFolderName+"/MatrixData/";
  string NeoXtrueFilename=SimulationFolderName+"/NEOXtrue.txt";
  neo MyNeo=CreateDefaultNeo(1.0); // 1.0 is scale factor 

  if (!threeDim) { // for 2d simulation
    string MatrixFilename=MatrixFolderName+"AcsrB.txt";
    SimulatePct2d(MyNeo,histories,yparts,xparts,boxLength,
     boxWidth,totalProjAngles,startAngle,angleSpace,
     pathOption,phantomOption,MatrixFilename,
     NeoXtrueFilename);
  }
  else {
    SimulatePct_KnownHull_2dHistories(MyNeo,histories,zparts,yparts,xparts,boxLength,boxWidth,totalProjAngles,startAngle,angleSpace,pathOption,phantomOption,MatrixFolderName,NeoXtrueFilename);

  }
  
}


