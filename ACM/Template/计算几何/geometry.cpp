#include<cmath>
#include<cstdio>
#include<vector>
#include<cstdlib>
#include<cstring>
#include<iostream>
#include<algorithm>
using namespace std;
const double Pi=acos(-1.0);
const double eps=1e-8;
const double R=6875.0/2;
struct P
{
    double x,y;
    P(double xx,double yy)
    {
        x=xx,y=yy;
    }
    P()
    {
        x=0,y=0;
    }
};
P O=P(0,0);
void Read(P &a)
{
    scanf("%lf%lf",&a.x,&a.y);
}
typedef P V;
V operator+(P A,P B)
{
    return V(A.x+B.x,A.y+B.y);
}
V operator-(P A,P B)
{
    return V(A.x-B.x,A.y-B.y);
}
V operator*(V A,double p)
{
    return V(A.x*p,A.y*p);
}
V operator*(double p,V A)
{
    return V(A.x*p,A.y*p);
}
V operator/(V A,double p)
{
    return V(A.x/p,A.y/p);
}
V operator/(double p,V A)
{
    return V(A.x/p,A.y/p);
}
double Abs(double x)
{
    return x>=0?x:-x;
}
int dcmp(double x)
{
    if(Abs(x)<eps)return 0;
    return x<0?-1:1;
}
double Min(double x,double y)
{
    if(dcmp(x-y)>=0)return y;
    return x;
}
double Max(double x,double y)
{
    if(dcmp(x-y)>=0)return x;
    return y;
}
double Sqr(double x)
{
    return x*x;
}
bool operator==(const P &a,const P &b)
{
    return dcmp(a.x-b.x)==0&&dcmp(a.y-b.y)==0;
}
double Dot(V a,V b)//点积
{
    return a.x*b.x+a.y*b.y;
}
double Cross(V a,V b)//叉积
{
    return a.x*b.y-a.y*b.x;
}
double Length(V a)//长度
{
    return sqrt(Dot(a,a));
}
V Normal(V A)//单位法向量
{
    double L=Length(A);
    return V(-A.y/L,A.x/L);
}
double Angle(V A,V B)//夹角
{
    double w=Dot(A,B)/(Length(A)*Length(B));
    //if(dcmp(w-1)==0)w=1;
    //if(dcmp(w+1)==0)w=-1;
    return acos(w);
}
double Area2(P A,P B,P C)//三角形面积两倍
{
    return Cross(B-A,C-A);
}
V Rotate(V a,double rad)//向量旋转逆时针rad,范围：0-2*Pi
{
    return V(a.x*cos(rad)-a.y*sin(rad),a.x*sin(rad)+a.y*cos(rad));
}
double Atan2(V A)//与x轴所成的夹角，范围[-Pi,Pi]
{
    return atan2(A.y,A.x);
}
int Quadrant(P A)//判断点在哪个象限
{
    if(dcmp(A.x)>=0)
    {
        if(dcmp(A.y)>=0)return 1;
        return 4;
    }
    else
    {
        if(dcmp(A.y)>=0)return 2;
        return 3;
    }
}
bool cmp1(P A,P B)//极角排序
{
    int u=Quadrant(A),v=Quadrant(B);
    if(u!=v)
        return u<v;
    if(dcmp(Cross(A-O,B-O))<0)
        return 1;
    return dcmp(Length(A-B))<0;
}
P LineLine(P p,V v,P q,V w)//计算两直线的交点（计算前应保证有交点,即Cross（v，w）！=0）
{
    V u=p-q;
    double t=Cross(w,u)/Cross(v,w);
    return p+v*t;
}
bool OnSegment(P a,P b,P p)//判断一点是否在线段上
{
    if(p==a||p==b)return 1;
    return dcmp(Cross(a-p,b-p))==0&&dcmp(Dot(a-p,b-p))<0;
}
bool SegmentIntersection(P a1,P a2,P b1,P b2)//线段与线段是否相交   规范相交（不会共线，没有端点交）
{
    //if(OnSegment(a1,b1,b2)||OnSegment(a2,b1,b2)||OnSegment(b1,a1,a2)||OnSegment(b2,a1,a2))return 1;
    double c1=Cross(a2-a1,b1-a1);
    double c2=Cross(a2-a1,b2-a1);
    double c3=Cross(b2-b1,a1-b1);
    double c4=Cross(b2-b1,a2-b1);
    return dcmp(c1)*dcmp(c2)<0&&dcmp(c3)*dcmp(c4)<0;
}
int LineSegment(P a1,P a2,P b1,P b2)//线段与直线相交判断,a1、a2为直线
{
    double x=Cross(a2-a1,b1-a1);
    double y=Cross(a2-a1,b2-a1);
    if(dcmp(x*y)<0)return 1;//不包含端点
    if(dcmp(x*y)==0)return 2;//包含端点
    return 0;
}
double PToSegment(P p,P a,P b)//点p到线段ab的距离
{
    if(a==b)
        return Length(p-a);
    V v1=b-a,v2=p-a,v3=p-b;
    if(dcmp(Dot(v1,v2))<0)
        return Length(v2);
    if(dcmp(Dot(v1,v3))>0)
        return Length(v3);
    return Abs(Area2(p,a,b)/Length(v1));
}
double PToLine(P p,P a,P b)//点到直线的距离
{
    if(a==b)
        return Length(p-a);
    return Abs(Area2(p,a,b))/Length(a-b);
}
double SegmentToSegment(P a1,P a2,P b1,P b2)//线段到线段的距离
{
    if(SegmentIntersection(a1,a2,b1,b2))return 0;
    return Min(Min(PToSegment(a1,b1,b2),PToSegment(a2,b1,b2)),Min(PToSegment(b1,a1,a2),PToSegment(b2,a1,a2)));
}
P PLineProjection(P p,P a,P b)//点p在直线上的投影点
{
    V v=b-a;
    return a+v*(Dot(v,p-a)/Dot(v,v));
}
//-----------------------------------------------------------------------------
struct Circle
{
    P c;
    double r;
    Circle(P cc,double rr)
    {
        c=cc,r=rr;
    }
    Circle()
    {
        c.x=c.y=r=0;
    }
    P p(double a)//a的范围0-2Pi
    {
        return P(c.x+r*cos(a),c.y+r*sin(a));
    }
    void GetIntersection(Circle &c2, vector<double> &rad)
    {
        V vc = c2.c - c;
        double d=Length(vc);
        if(dcmp(d)==0) return;
        if(dcmp(r+c2.r-d)<0) return;
        if(dcmp(Abs(r-c2.r)-d)>0) return;
        double a=Atan2(vc);
        double da=acos((r*r+d*d-c2.r*c2.r)/(2*r*d));
        rad.push_back((a-da));
        rad.push_back((a+da));
    }
};
double Torad(double deg)//角度转化成弧度
{
    return deg/180*Pi;
}
void GetCoor(double r,double lat,double lng,double& x,double& y,double& z)//纬度经度转化成空间坐标
{
    lat=Torad(lat);
    lng=Torad(lng);
    x=r*cos(lat)*cos(lng);
    y=r*cos(lat)*sin(lng);
    z=r*sin(lat);
}
double Getdistance(double x1,double y1,double z1,double x2,double y2,double z2,double r)//球面最短距离（走大圆）
{
    double d=sqrt(Sqr(x1-x2)+Sqr(y1-y2)+Sqr(z1-z2));
    return asin(d/2/r)*2*r;
}
int GetTagents(P p,Circle C,V *v)//过点p到圆C的切线，v[i]是第i条切线的向量，返回切线条数
{
    V u=C.c-p;
    double dist=Length(u);
    if(dcmp(dist-C.r)<0)
        return 0;
    if(dcmp(dist-C.r)==0)
    {
        v[0]=Rotate(u,Pi/2);
        return 1;
    }
    double ang=asin(C.r/dist);
    v[0]=Rotate(u,-ang);
    v[1]=Rotate(u,ang);
    return 2;
}
int LineCircle(P p1,P p2,Circle C,P &ans1,P &ans2)//直线和圆交点
{
    double a=p2.x-p1.x,b=p1.x-C.c.x,c=p2.y-p1.y,d=p1.y-C.c.y;
    double e=a*a+c*c,f=2*(a*b+c*d),g=b*b+d*d-C.r*C.r;
    double del=f*f-4*e*g;
    double t1,t2;
    if(dcmp(del)<0)return 0;
    if(dcmp(del)==0)
    {
        t1=-f/(2*e);
        ans1=p1+t1*(p2-p1);
        return 1;
    }
    t1=(-f-sqrt(del))/(2*e);
    t2=(-f+sqrt(del))/(2*e);
    ans1=p1+t1*(p2-p1);
    ans2=p1+t2*(p2-p1);
    return 2;
}
int CircleCircle(Circle C1,Circle C2,P &ans1,P &ans2)//圆和圆的交点
{
    double d=Length(C1.c-C2.c);
    if(dcmp(d)==0)
    {
        if(dcmp(C1.r-C2.r)==0)return -2;//两圆重合
        return 0;//包含关系
    }
    if(dcmp(C1.r+C2.r-d)<0)return 0;//相离
    if(dcmp(Abs(C1.r-C2.r)-d)>0)return 0;//内含
    double a=Atan2(C2.c-C1.c);
    double da=acos((Sqr(C1.r)+Sqr(d)-Sqr(C2.r))/(2*C1.r*d));
    ans1=C1.p(a+da),ans2=C1.p(a-da);
    if(ans1==ans2)return 1;
    return 2;
}
double AreaCircleCircle(Circle C1,Circle C2)//两圆的相交面积
{
    P ans1,ans2;
    double d=Length(C1.c-C2.c);
    if(dcmp(d-C1.r-C2.r)>=0)
        return 0;
    int k=CircleCircle(C1,C2,ans1,ans2);
    if(k<=1)
    {
        if(dcmp(C1.r>C2.r))
            return Pi*Sqr(C2.r);
        return Pi*Sqr(C1.r);
    }
    double ans=0;
    double angle1=Angle(ans1-C1.c,C2.c-C1.c);
    ans+=Sqr(C1.r)*angle1;
    double angle2=Angle(ans1-C2.c,C1.c-C2.c);
    ans+=Sqr(C2.r)*angle2;
    ans=ans-Abs(Area2(C1.c,C2.c,ans1)/2)-Abs(Area2(C1.c,C2.c,ans2)/2);
    return ans;
}
Circle OutCircle(P p1,P p2,P p3)//三角形外接圆,保证不共线
{
    double bx=p2.x-p1.x;
    double by=p2.y-p1.y;
    double cx=p3.x-p1.x;
    double cy=p3.y-p1.y;
    double d=2*(bx*cy-cx*by);
    double ccx=(cy*(bx*bx+by*by)-by*(cx*cx+cy*cy))/d+p1.x;
    double ccy=(bx*(cx*cx+cy*cy)-cx*(bx*bx+by*by))/d+p1.y;
    P p=P(ccx,ccy);
    return Circle(p,Length(p1-p));
}
Circle InCircle(P p1,P p2,P p3)//三角形内切圆，保证不共线
{
    double a=Length(p2-p3);
    double b=Length(p3-p1);
    double c=Length(p1-p2);
    P p=(p1*a+p2*b+p3*c)/(a+b+c);
    return Circle(p,PToLine(p,p1,p2));
}
P Chuixin(P a,P b,P c)//三角形的垂心
{
    double dx=2*(a.x-b.x),dy=2*(a.y-b.y);
    double ex=2*(b.x-c.x),ey=2*(b.y-c.y);
    double dz=a.x*a.x-b.x*b.x+a.y*a.y-b.y*b.y;
    double ez=b.x*b.x-c.x*c.x+b.y*b.y-c.y*c.y;
    double k=(dz*ex-ez*dx)/(ex*dy-dx*ey);
    double h=(dz*ey-ez*dy)/(ey*dx-dy*ex);
    double anx=a.x+b.x+c.x-h-h;
    double any=a.y+b.y+c.y-k-k;
    return P(anx,any);
}
bool cmp2(P a,P b)
{
    if(dcmp(a.x-b.x)==0)
        return dcmp(a.y-b.y)<=0;
    return dcmp(a.x-b.x)<0;
}
void Convex(P *p,int n,P *q,int &m)//计算凸包。输入点集p，输出点集q和个数m
{
    sort(p,p+n,cmp2);
    //n=unique(p,p+n)-p;//去掉重复的点
    int cnt=0;
    for(int i=0;i<n;i++)
    {
        while(cnt>1&&dcmp(Cross(q[cnt-1]-q[cnt-2],p[i]-q[cnt-2]))<=0)//如果点能在凸包边上，去掉等号
            cnt--;
        q[cnt++]=p[i];
    }
    int num=cnt;
    for(int i=n-2;i>=0;i--)
    {
        while(cnt>num&&dcmp(Cross(q[cnt-1]-q[cnt-2],p[i]-q[cnt-2]))<=0)
            cnt--;
        q[cnt++]=p[i];
    }
    if(n>1)cnt--;
    m=cnt;
}
double PolygonArea(P *p,int n)//多边形有向面积
{
    double area=0;
    for(int i=1;i<n-1;i++)
    {
        area+=Cross(p[i]-p[0],p[i+1]-p[0]);
    }
    return area/2;
}
double PolygonC(P *p,int n)//多边形周长
{
    double ans=0;
    p[n]=p[0];
    for(int i=0;i<n;i++)
        ans+=Length(p[i]-p[i+1]);
    return ans;
}
int InPolygon(P *p,int n,P a)//判断一点是否在多边形内部
{
    p[n]=p[0];
    for(int i=0;i<n;i++)
    {
        if(OnSegment(p[i],p[i+1],a))
            return 1;
        double s=Cross(p[i+1]-p[i],a-p[i]);
        if(dcmp(s)<0)
            return 0;
    }
    return 2;
}
/*多边形的重心公式：
每个三角形重心：cx = x1 + x2 + x3；cy同理。
每个三角形面积：s =  ( (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1) ) / 2；
多边形重心：cx = (∑ cx[i]*s[i]) / （3*∑s[i]）;  cy = (∑ cy[i]*s[i] ) / （3*∑s[i]）;
*/
double MinRectangle(P *p,int n)//外接矩形最小面积
{
    double ans=-1;
    p[n]=p[0];
    int u=0,v=0,w=0;
    for(int i=0;i<n;i++)
    {
        if(dcmp(p[i].x-p[u].x)>0)u=i;
        if(dcmp(p[i].y-p[w].y)<0)w=i;
        if(dcmp(p[i].y-p[v].y)>0)v=i;
    }
    p[n]=p[0];p[n+1]=p[1];
    for(int i=0;i<n;i++)
    {
        while(dcmp(PToLine(p[u],p[i],p[i+1])-PToLine(p[u+1],p[i],p[i+1]))<=0)
            u=(u+1)%n;
        while(dcmp(Angle(p[v+1]-p[v],p[i+1]-p[i])-Pi/2)>=0)
            v=(v+1)%n;
        while(dcmp(Angle(p[w+1]-p[w],p[i+1]-p[i])-Pi/2)<=0)
            w=(w+1)%n;
        double S=PToLine(p[u],p[i],p[i+1])*Length(PLineProjection(p[v],p[i],p[i+1])-PLineProjection(p[w],p[i],p[i+1]));
        if(ans==-1||dcmp(S-ans)<0)
            ans=S;
    }
    return ans;
}
double TriangleCircle1(P a,P b,Circle C)
{
    P q[10];
    double d1=Length(C.c-a),d2=Length(C.c-b);
    int w=dcmp(Cross(a-C.c,b-C.c));
    double angle=Abs(Angle(a-C.c,b-C.c));
    if(w==0)
        return 0;
    if(dcmp(d1-C.r)<=0&&dcmp(d2-C.r)<=0)
        return Cross(a-C.c,b-C.c)/2;
    if(dcmp(d1-C.r)>=0&&dcmp(d2-C.r)>=0)
    {
        int k=LineCircle(a,b,C,q[0],q[1]);
        if(k<2||!OnSegment(a,b,q[0])||!OnSegment(a,b,q[1]))
        {
//            LineCircle(a,C.c,C,q[0],q[1]);
//                if(!OnSegment(a,C.c,q[0]))swap(q[0],q[1]);
//            LineCircle(b,C.c,C,q[0],q[1]);
//                if(!OnSegment(b,C.c,q[1]))swap(q[0],q[1]);
            return Sqr(C.r)*angle/2*w;
        }
        else
        {
            double ans=Sqr(C.r)*angle/2;
            angle=Abs(Angle(q[0]-C.c,q[1]-C.c));
            ans=ans-(Sqr(C.r)*angle/2-Abs(Area2(C.c,q[0],q[1]))/2);
            return ans*w;
        }
    }
    if(dcmp(d2-C.r)>=0)swap(a,b);
    LineCircle(a,b,C,q[0],q[1]);
    if(!OnSegment(a,b,q[0]))swap(q[0],q[1]);
    angle=Abs(Angle(q[0]-C.c,a-C.c));
    return w*(Sqr(C.r)*angle/2+Abs(Area2(q[0],b,C.c)/2));
}
double PolygonCircle2(P *p,int n,Circle C)//多边形与圆的面积交
{
    p[n]=p[0];
    double ans=0;
    for(int i=0;i<n;i++)
    {
        ans+=TriangleCircle1(p[i],p[i+1],C);
    }
    return Abs(ans);
}
struct Line
{
    P p;
    V v;
    double ang;
    Line(){}
    Line(P p,V v):p(p),v(v){ang=atan2(v.y,v.x);}
    bool operator<(const Line &L)const
    {
        return dcmp(ang-L.ang)<0;
    }
};
bool OnLeft(Line L,P p)//判断一点是否在一条直线的左边
{
    return dcmp(Cross(L.v,p-L.p))>=0;
}
P GetIntersection(Line a,Line b)//计算两直线交点
{
    V u=a.p-b.p;
    double t=Cross(b.v,u)/Cross(a.v,b.v);
    return a.p+t*a.v;
}
int HPI(Line *L,int n,P *poly)
{
    sort(L,L+n);
    int fi=0,la=0;
    P *p=new P[n];
    Line *q=new Line[n];
    q[0]=L[0];
    for(int i=1;i<n;i++)
    {
        while(fi<la&&!OnLeft(L[i],p[la-1]))la--;
        while(fi<la&&!OnLeft(L[i],p[fi]))fi++;
        q[++la]=L[i];
        if(dcmp(Abs(Cross(q[la].v,q[la-1].v))-eps)<0)//平行
        {
            la--;
            if(OnLeft(q[la],L[i].p))q[la]=L[i];
        }
        if(fi<la)p[la-1]=GetIntersection(q[la-1],q[la]);
    }
    while(fi<la&&!OnLeft(q[fi],p[la-1]))la--;
    if(la-fi<=1)return 0;//空集
    p[la]=GetIntersection(q[la],q[fi]);//首尾两个半平面交点
    int m=0;
    for(int i=fi;i<=la;i++)
        poly[m++]=p[i];
    return m;
}
//---------------------------------------------------------------------
struct P3
{
    double x,y,z;
    P3(double x=0,double y=0,double z=0):x(x),y(y),z(z){}
};
typedef P3 V3;
V3 operator+(V3 A,V3 B)
{
    return V3(A.x+B.x,A.y+B.y,A.z+B.z);
}
V3 operator-(P3 A,P3 B)
{
    return V3(A.x-B.x,A.y-B.y,A.z-B.z);
}
V3 operator*(V3 A,double p)
{
    return V3(A.x*p,A.y*p,A.z*p);
}
V3 operator/(V3 A,double p)
{
    return V3(A.x/p,A.y/p,A.z/p);
}
double Dot3(V3 A,V3 B)
{
    return A.x*B.x+A.y*B.y+A.z*B.z;
}
double Length3(V3 A)
{
    return sqrt(Dot3(A,A));
}
double Angle3(V3 A,V3 B)//两向量的夹角
{
    return acos(Dot3(A,B)/(Length3(A)*Length3(B)));
}
double PToPlane3(P3 p,P3 p0,V3 n)//p到p0-n的距离
{
    return Abs(Dot3(p-p0,n))/Length3(n);//如果去掉绝对值，得到的是有向距离；
}
V3 UnitLength(V3 A)//单位长度
{
    return A/Length3(A);
}
P3 GetProjection(P3 p,P3 p0,V3 n)//p在p0-n上的投影点
{
    double d=Dot3(p-p0,n)/Length3(n);//直线与平面的交点
    return p-UnitLength(n)*d;
}
P3 LinePlaneIntersection(P3 p1,P3 p2,P3 p0,V3 n)//直线与平面交点
{
    V3 v=p2-p1;
    double t=Dot3(n,p0-p1)/Dot3(n,v);//判断分母是否为0
    return p1+v*t;
}
V3 Cross3(V3 A,V3 B)
{
    return V3(A.y*B.z-A.z*B.y,A.z*B.x-A.x*B.z,A.x*B.y-A.y*B.x);
}
double P3ToLine(P3 p,P3 a,P3 b)//p到直线ab的距离
{
    V3 v1=b-a,v2=p-a;
    return Length3(Cross3(v1,v2))/Length3(v1);
}
double LineLine3(P3 p,V3 v,P3 q,V3 w)
{
    V3 n=Cross3(v,w);
    return Dot3(n,p-q)/Length3(n);
}
double Area3(P3 A,P3 B,P3 C)//空间中的三角形面积的2倍
{
    return Length3(Cross3(B-A,C-A));
}
#define N 70
struct Face
{
    int a,b,c;                                                        //表示凸包一个面上三个点的编号
    bool ok;                                                          //表示该面是否属于最终凸包中的面
};
struct CH3D
{
    int n;                                                                   //初始顶点数
    P3 p[N];                                                           //初始顶点
    int num;                                                                 //凸包表面的三角形数
    Face F[8*N];
    int g[N][N];                                                       //凸包表面的三角形
    double vlen(P3 a)                            //向量长度
    {
        return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
    }
    P3 cross(const P3 &a, const P3 &b, const P3 &c)             //叉乘
    {
         return P3((b.y-a.y)*(c.z-a.z)-(b.z-a.z)*(c.y-a.y),-((b.x-a.x)*(c.z-a.z)
             -(b.z-a.z)*(c.x-a.x)),(b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x));
    }
    double area(P3 a,P3 b,P3 c)                                   //三角形面积*2
    {
          return vlen(Cross3(b-a,c-a));
    }
    double volume(P3 a,P3 b,P3 c,P3 d)                        //四面体有向体积*6
    {
          return Dot3(Cross3(b-a,c-a),d-a);
    }
    double dblcmp(P3 &q,Face &f)                                       //正:点在面同向
    {
          P3 m=p[f.b]-p[f.a];
          P3 n=p[f.c]-p[f.a];
          P3 t=q-p[f.a];
          return Dot3(Cross3(m,n),t);
    }
    void deal(int x,int a,int b)
    {
        int f=g[a][b];
        Face add;
        if(F[f].ok)
        {
             if(dcmp(dblcmp(p[x],F[f]))>0)
                 dfs(x,f);
             else
             {
                 add.a=b;
                 add.b=a;
                 add.c=x;
                 add.ok=1;
                 g[x][b]=g[a][x]=g[b][a]=num;
                 F[num++]=add;
             }
        }
    }
    void dfs(int p,int now)
    {
        F[now].ok=0;
        deal(p,F[now].b,F[now].a);
        deal(p,F[now].c,F[now].b);
        deal(p,F[now].a,F[now].c);
    }
    bool same(int s,int t)
    {
        P3 &a=p[F[s].a];
        P3 &b=p[F[s].b];
        P3 &c=p[F[s].c];
        return dcmp(volume(a,b,c,p[F[t].a]))==0&&dcmp(volume(a,b,c,p[F[t].b]))==0&&dcmp(volume(a,b,c,p[F[t].c]))==0;
    }
    void solve()                                                         //构建三维凸包
    {
        int i,j,tmp;
        Face add;
        bool flag=true;
        num=0;
        if(n<4)
           return;
        for(i=1;i<n;i++)                                              //此段是为了保证前四个点不共面,若以保证,则可去掉
        {
            if(dcmp(vlen(p[0]-p[i]))>0)
            {
                swap(p[1],p[i]);
                flag=false;
                break;
            }
        }
        if(flag)
            return;
        flag=true;
        for(i=2;i<n;i++)                                             //使前三点不共线
        {
            if(dcmp(vlen(Cross3(p[0]-p[1],p[1]-p[i])))>0)
            {
                swap(p[2],p[i]);
                flag=false;
                break;
            }
        }
        if(flag)
            return;
        flag=true;
        for(i=3;i<n;i++)                                            //使前四点不共面
        {
            if(dcmp(Dot3(Cross3(p[0]-p[1],p[1]-p[2]),p[0]-p[i]))!=0)
            {
                swap(p[3],p[i]);
                flag=false;
                break;
            }
        }
        if(flag)
            return;
        for(i=0;i<4;i++)
        {
            add.a=(i+1)%4;
            add.b=(i+2)%4;
            add.c=(i+3)%4;
            add.ok=true;
            if(dcmp(dblcmp(p[i],add))>0)
               swap(add.b,add.c);
            g[add.a][add.b]=g[add.b][add.c]=g[add.c][add.a]=num;
            F[num++]=add;
        }
        for(i=4;i<n;i++)
        {
            for(j=0;j<num;j++)
            {
                if(F[j].ok&&dcmp(dblcmp(p[i],F[j]))>0)
                {
                    dfs(i,j);
                    break;
                }
            }
        }
        tmp=num;
        for(i=num=0;i<tmp;i++)
          if(F[i].ok)
          {
                 F[num++]=F[i];
          }
    }
    double area() //表面积
    {
        double res=0.0;
        if(n==3)
        {
           P3 q=cross(p[0],p[1],p[2]);
           res=vlen(q)/2.0;
           return res;
        }
        for(int i=0;i<num;i++)
          res+=area(p[F[i].a],p[F[i].b],p[F[i].c]);
        return res/2.0;
    }
    double volume()                                                  //体积
    {
        double res=0.0;
        P3 tmp(0,0,0);
        for(int i=0;i<num;i++)
         res+=volume(tmp,p[F[i].a],p[F[i].b],p[F[i].c]);
        return Abs(res/6.0);
    }
    int triangle()                                                  //表面三角形个数
    {
        return num;
    }
    int polygon()                                                   //表面多边形个数
    {
       int i,j,res,flag;
       for(i=res=0;i<num;i++)
       {
            flag=1;
            for(j=0;j<i;j++)
             if(same(i,j))
             {
                  flag=0;
                  break;
             }
            res+=flag;
       }
       return res;
    }
    P3 getcent()//求凸包质心
    {
       P3 ans(0,0,0),temp=p[F[0].a];
       double v = 0.0,t2;
       for(int i=0;i<num;i++)
       {
           if(F[i].ok)
           {
               P3 p1=p[F[i].a],p2=p[F[i].b],p3=p[F[i].c];
               t2 = volume(temp,p1,p2,p3)/6.0;//体积大于0，也就是说，点 temp 不在这个面上
               if(t2>0)
               {
                   ans.x+=(p1.x+p2.x+p3.x+temp.x)*t2;
                   ans.y+=(p1.y+p2.y+p3.y+temp.y)*t2;
                   ans.z+=(p1.z+p2.z+p3.z+temp.z)*t2;
                   v+=t2;
               }
           }
       }
       ans.x/=(4*v);
       ans.y/=(4*v);
       ans.z/=(4*v);
       return ans;
    }
    double fun(P3 fuck)
    {//点到凸包上的最近距离（枚举每个面到这个点的距离）
       double mi=99999999;
       for(int i=0;i<num;i++)
       {
           if(F[i].ok==true)
           {
               P3 p1=p[F[i].a],p2=p[F[i].b],p3=p[F[i].c];
               double a=(p2.y-p1.y)*(p3.z-p1.z)-(p2.z-p1.z)*(p3.y-p1.y);
               double b=(p2.z-p1.z)*(p3.x-p1.x)-(p2.x-p1.x)*(p3.z-p1.z);
               double c=(p2.x-p1.x)*(p3.y-p1.y)-(p2.y-p1.y)*(p3.x-p1.x);
               double d=-(a*p1.x+b*p1.y+c*p1.z);
               double temp=Abs(a*fuck.x+b*fuck.y+c*fuck.z+d)/sqrt(a*a+b*b+c*c);
               if(dcmp(temp-mi)<0)mi=temp;
           }
       }
       return mi;
    }
};
P3 GetCent(Face *face,int n)
{

}
//struct Face
//{
//    int v[3];
//    bool ok;
//    V3 Normal(P3 *p)const
//    {
//        return Cross3(p[v[1]]-p[v[0]],p[v[2]]-p[v[0]]);
//    }
//    int Cansee(P3 *p,int i)const
//    {
//        return Dot3(p[i]-p[v[0]],Normal(p))>0?1:0;
//    }
//};
//bool vis[N][N];
//vector<Face> CH3D(P3 *p,int n)//三维凸包
//{
//    vector<Face> cur;
//    cur.push_back((Face){0,1,2});
//    cur.push_back((Face){2,1,0});
//    for(int i=3;i<n;i++)
//    {
//        if(p[i].x==0&&p[i].y==2&&p[i].z==0)puts("fdsa");
//        vector<Face> next;
//        int num=cur.size();
//        for(int j=0;j<num;j++)
//        {
//            Face &f=cur[j];
//            int res=f.Cansee(p,i);
//            if(!res)next.push_back(f);
//            for(int k=0;k<3;k++)
//                vis[f.v[k]][f.v[(k+1)%3]]=res;
//        }
//        num=cur.size();
//        for(int j=0;j<num;j++)
//            for(int k=0;k<3;k++)
//        {
//            int a=cur[j].v[k],b=cur[j].v[(k+1)%3];
//            if(vis[a][b]!=vis[b][a]&&vis[a][b])
//            next.push_back((Face){a,b,i});
//        }
//        cur=next;
//    }
//    return cur;
//}
CH3D A,B;
int main()
{
//    freopen("test1.in","r",stdin);
//    freopen("test1.out","w",stdout);
    //printf("%lf\n",LineLine3(P3(0,0,0),V3(1,0,0),P3(0,0,2),V3(0,1,0)));
    int T,n;
    P3 p[35][3],c[35];
    double r[35];
    V3 v[35];
    scanf("%d",&T);
    while(T--)
    {
        scanf("%d",&n);
        for(int i=0;i<n;i++)
            for(int j=0;j<3;j++)
                scanf("%lf%lf%lf",&p[i][j].x,&p[i][j].y,&p[i][j].z);
        for(int i=0;i<n;i++)
            c[i]=p[i][0],r[i]=Length3(p[i][0]-p[i][1]),v[i]=Cross3(p[i][1]-p[i][0],p[i][2]-p[i][0]);
        double mi=-1;
        bool is=0;
        for(int i=0;i<n&&!is;i++)
            for(int j=i+1;j<n;j++)
        {
            V3 w=Cross3(v[i],v[j]);
            double d=Dot3(w,c[i]-c[j])/Length3(w);
            d=Abs(d);
            d=d-r[i]-r[j];
            if(d<=0)
            {
                is=1;
                break;
            }
            else if(dcmp(d-mi)<0||mi==-1)mi=d;
        }
        if(is)puts("Lucky");
        else printf("%.2f\n",mi);
    }
    return 0;
}
/*









*/

