#include <iostream>
#include <cmath> // for calculating square root. without it this task would be just ridiculous.


struct Point
{
    float x, y, z;
};
struct Sphere
{
    Point p;
    float r;
};

float determ(float a[4][4],int n) {
    int  p, i, j, k, h ;
    float det=0, temp[4][4];
    if(n==1) {
        return a[0][0];
    } else if(n==2) {
        det=(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
        return det;
    } else {
        for(p=0;p<n;p++) {
            h = 0;
            k = 0;
            for(i=1;i<n;i++) {
                for( j=0;j<n;j++) {
                    if(j==p) {
                        continue;
                    }
                    temp[h][k] = a[i][j];
                    k++;
                    if(k==n-1) {
                        h++;
                        k = 0;
                    }
                }
            }
            det=det+a[0][p]*powf(-1,p)*determ(temp,n-1);
        }
        return det;
    }
}


Sphere calcSphere(Point* SPoints)
{
    float t[4][4] = {{SPoints[0].x,SPoints[0].y,SPoints[0].z,1},
                     {SPoints[1].x,SPoints[1].y,SPoints[1].z,1},
                     {SPoints[2].x,SPoints[2].y,SPoints[2].z,1},
                     {SPoints[3].x,SPoints[3].y,SPoints[3].z,1}};
    float T = determ(t,4);

    float d[4][4] = {{-SPoints[0].x*SPoints[0].x+SPoints[0].y*SPoints[0].y+SPoints[0].z*SPoints[0].z,SPoints[0].y,SPoints[0].z,1},
                     {-SPoints[1].x*SPoints[1].x+SPoints[0].y*SPoints[1].y+SPoints[1].z*SPoints[1].z,SPoints[1].y,SPoints[1].z,1},
                     {-SPoints[2].x*SPoints[2].x+SPoints[0].y*SPoints[2].y+SPoints[2].z*SPoints[2].z,SPoints[2].y,SPoints[2].z,1},
                     {-SPoints[3].x*SPoints[3].x+SPoints[0].y*SPoints[3].y+SPoints[3].z*SPoints[3].z,SPoints[3].y,SPoints[3].z,1}};
    float D = determ(d,4)/T;

    float e[4][4] = {{SPoints[0].x,-SPoints[0].x*SPoints[0].x+SPoints[0].y*SPoints[0].y+SPoints[0].z*SPoints[0].z,SPoints[0].z,1},
                     {SPoints[1].x,-SPoints[1].x*SPoints[1].x+SPoints[0].y*SPoints[1].y+SPoints[1].z*SPoints[1].z,SPoints[1].z,1},
                     {SPoints[2].x,-SPoints[2].x*SPoints[2].x+SPoints[0].y*SPoints[2].y+SPoints[2].z*SPoints[2].z,SPoints[2].z,1},
                     {SPoints[3].x,-SPoints[3].x*SPoints[3].x+SPoints[0].y*SPoints[3].y+SPoints[3].z*SPoints[3].z,SPoints[3].z,1}};
    float E = determ(e,4)/T;

    float f[4][4] = {{SPoints[0].x,SPoints[0].y,-SPoints[0].x*SPoints[0].x+SPoints[0].y*SPoints[0].y+SPoints[0].z*SPoints[0].z,1},
                     {SPoints[1].x,SPoints[1].y,-SPoints[1].x*SPoints[1].x+SPoints[0].y*SPoints[1].y+SPoints[1].z*SPoints[1].z,1},
                     {SPoints[2].x,SPoints[2].y,-SPoints[2].x*SPoints[2].x+SPoints[0].y*SPoints[2].y+SPoints[2].z*SPoints[2].z,1},
                     {SPoints[3].x,SPoints[3].y,-SPoints[3].x*SPoints[3].x+SPoints[0].y*SPoints[3].y+SPoints[3].z*SPoints[3].z,1}};
    float F = determ(f,4)/T;

    float g[4][4] = {{SPoints[0].x,SPoints[0].y,SPoints[0].z,-SPoints[0].x*SPoints[0].x+SPoints[0].y*SPoints[0].y+SPoints[0].z*SPoints[0].z},
                     {SPoints[1].x,SPoints[1].y,SPoints[1].z,-SPoints[1].x*SPoints[1].x+SPoints[0].y*SPoints[1].y+SPoints[1].z*SPoints[1].z},
                     {SPoints[2].x,SPoints[2].y,SPoints[2].z,-SPoints[2].x*SPoints[2].x+SPoints[0].y*SPoints[2].y+SPoints[2].z*SPoints[2].z},
                     {SPoints[3].x,SPoints[3].y,SPoints[3].z,-SPoints[3].x*SPoints[3].x+SPoints[0].y*SPoints[3].y+SPoints[3].z*SPoints[3].z}};
    float G = determ(g,4)/T;

    Sphere S;
    S.p.x = -D/2;
    S.p.y = -E/2;
    S.p.z = -F/2;
    S.r = 1/2*sqrtf(D*D+E*E+F*F-4G);
    return S;
}

Point midpoint(Point p0, Point p1)
{
    Point p;
    p.x = p0.x+p1.x/2;
    p.y = p0.y+p1.y/2;
    p.z = p0.z+p1.z/2;
    return p;
}

float distance(Point p0, Point p1)
{
    return sqrtf((p0.x-p1.x)*(p0.x-p1.x)+(p0.y-p1.y)*(p0.y-p1.y)+(p0.z-p1.z)*(p0.z-p1.z));
}

Sphere* sed(const Point* points, Point* sPoints, uint32_t numPoints, uint32_t numSPoints)
{
    Sphere* S;

    if(numPoints == 0 || numSPoints >= 4)
    {
        if (numSPoints == 1)
        {
            Point p0 = sPoints[0];
            S->p = p0;
            S->r = 0;
            return S;
        }
        else if (numSPoints == 2)
        {
            Point p0 = sPoints[0];
            Point p1 = sPoints[1];
            Point center = midpoint(p0,p1);
            float diameter = distance(p0, p1);
            S->p = center;
            S->r = diameter / 2.0f;
            return S;
        }
        else if (numSPoints == 3)
        {
            //circle
        }
    }
}

Sphere calculateBoundingSphere(const Point* points, uint32_t numPoints)
{
    Sphere result;
    {

    }
    return result;
}


int main() {
    std::cout << "Hello, World!" << std::endl;
    return 0;
}