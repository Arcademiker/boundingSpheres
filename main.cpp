#include <iostream>
#include <cmath> // for calculating square root and pow nothing else. I hope that was allowed.
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
using namespace std;

struct Point
{
    float x, y, z;

    Point operator-(Point subP){
        Point retP;
        retP.x = subP.x-this->x;
        retP.y = subP.y-this->y;
        retP.z = subP.z-this->z;
        return retP;
    }
    Point operator+(Point subP){
        Point retP;
        retP.x = subP.x+this->x;
        retP.y = subP.y+this->y;
        retP.z = subP.z+this->z;
        return retP;
    }
    float operator*(Point subP){
        float x = subP.x*this->x;
        float y = subP.y*this->y;
        float z = subP.z*this->z;
        return x+y+z;
    }
    Point operator*(float scalar){
        Point retP;
        retP.x = scalar*this->x;
        retP.y = scalar*this->y;
        retP.z = scalar*this->z;
        return retP;
    }
};

struct Sphere
{
    Point p;
    float r;
    bool exist;
};

float determ(float a[4][4],int n)
{
    int  p, i, j, k, h ;
    float det=0.0f, temp[4][4];
    if(n==1)
    {
        return a[0][0];
    }
    else if(n==2)
    {
        det=(a[0][0]*a[1][1]-a[0][1]*a[1][0]);
        return det;
    }
    else
    {
        for(p=0;p<n;p++)
        {
            h = 0;
            k = 0;
            for(i=1;i<n;i++)
            {
                for( j=0;j<n;j++)
                {
                    if(j==p)
                    {
                        continue;
                    }
                    temp[h][k] = a[i][j];
                    k++;
                    if(k==n-1)
                    {
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

float distance(Point p0, Point p1)
{
    return sqrtf((p0.x-p1.x)*(p0.x-p1.x)+(p0.y-p1.y)*(p0.y-p1.y)+(p0.z-p1.z)*(p0.z-p1.z));
}



Sphere calcCircle(Point* SPoints)
{
    Sphere S;

    float Cx = SPoints[1].x-SPoints[0].x;
    float Cy = SPoints[1].y-SPoints[0].y;
    float Cz = SPoints[1].z-SPoints[0].z;
    float Bx = SPoints[2].x-SPoints[0].x;
    float By = SPoints[2].y-SPoints[0].y;
    float Bz = SPoints[2].z-SPoints[0].z;
    float B2 = SPoints[0].x*SPoints[0].x-SPoints[2].x*SPoints[2].x+SPoints[0].y*SPoints[0].y-SPoints[2].y*SPoints[2].y+SPoints[0].z*SPoints[0].z-SPoints[2].z*SPoints[2].z;
    float C2 = SPoints[0].x*SPoints[0].x-SPoints[1].x*SPoints[1].x+SPoints[0].y*SPoints[0].y-SPoints[1].y*SPoints[1].y+SPoints[0].z*SPoints[0].z-SPoints[1].z*SPoints[1].z;

    float CByz = Cy*Bz-Cz*By;
    float CBxz = Cx*Bz-Cz*Bx;
    float CBxy = Cx*By-Cy*Bx;
    float ZZ1 = -(Bz-Cz*Bx/Cx)/(By-Cy*Bx/Cx);
    float Z01 = -(B2-Bx/Cx*C2)/(2.0f*(By-Cy*Bx/Cx));
    float ZZ2 = -(ZZ1*Cy+Cz)/Cx;
    float Z02 = -(2.0f*Z01*Cy+C2)/(2.0f*Cx);

    S.p.z = -((Z02-SPoints[0].x)*CByz-(Z01-SPoints[0].y)*CBxz-SPoints[0].z*CBxy)/(ZZ2*CByz-ZZ1*CBxz+CBxy);
    S.p.x = ZZ2*S.p.z + Z02;
    S.p.y = ZZ1*S.p.z + Z01;

    S.r = distance(S.p,SPoints[0]);

    std::cout << "circle center: " << S.p.x << " " << S.p.y << " " << S.p.z << std::endl;
    std::cout << "with point0: " << SPoints[0].x << " " << SPoints[0].y << " " << SPoints[0].z << " distance: " << distance(S.p, SPoints[0] ) << std::endl;
    std::cout << "with point1: " << SPoints[1].x << " " << SPoints[1].y << " " << SPoints[1].z << " distance: " << distance(S.p, SPoints[1] ) << std::endl;
    std::cout << "with point2: " << SPoints[2].x << " " << SPoints[2].y << " " << SPoints[2].z << " distance: " << distance(S.p, SPoints[2] ) << std::endl;
    std::cout << "radius: " << S.r << std::endl;
    // the plane of the three point:
    float a = (SPoints[1].y-SPoints[0].y)*(SPoints[2].z-SPoints[0].z)-(SPoints[2].y-SPoints[0].y)*(SPoints[1].z-SPoints[0].z);
    float b = (SPoints[1].z-SPoints[0].z)*(SPoints[2].x-SPoints[0].x)-(SPoints[2].z-SPoints[0].z)*(SPoints[1].x-SPoints[0].x);
    float c = (SPoints[1].x-SPoints[0].x)*(SPoints[2].y-SPoints[0].y)-(SPoints[2].x-SPoints[0].x)*(SPoints[1].y-SPoints[0].y);
    float n = -(a*SPoints[1].x+b*SPoints[1].y+c*SPoints[1].z);
    std::cout << ">>>>plane test: 0? " << a*S.p.x+b*S.p.y+c*S.p.z+n << std::endl << std::endl;

    return S;
}


// solve:
//
//    | ( x^2 +  y^2 +  z^2)  x   y   z   1  |
//    |                                      |
//    | (x1^2 + y1^2 + z1^2)  x1  y1  z1  1  |
//    |                                      |
//    | (x2^2 + y2^2 + z2^2)  x2  y2  z2  1  | = 0
//    |                                      |
//    | (x3^2 + y3^2 + z3^2)  x3  y3  z3  1  |
//    |                                      |
//    | (x4^2 + y4^2 + z4^2)  x4  y4  z4  1  |
Sphere calcSphere(Point* SPoints)
{
    Sphere S;
    float t[4][4] = {{SPoints[0].x,SPoints[0].y,SPoints[0].z,1},
                     {SPoints[1].x,SPoints[1].y,SPoints[1].z,1},
                     {SPoints[2].x,SPoints[2].y,SPoints[2].z,1},
                     {SPoints[3].x,SPoints[3].y,SPoints[3].z,1}};
    float T = determ(t,4);
    if (T!=0)
    {
        float d[4][4] = {{-(SPoints[0].x*SPoints[0].x+SPoints[0].y*SPoints[0].y+SPoints[0].z*SPoints[0].z),SPoints[0].y,SPoints[0].z,1},
                         {-(SPoints[1].x*SPoints[1].x+SPoints[1].y*SPoints[1].y+SPoints[1].z*SPoints[1].z),SPoints[1].y,SPoints[1].z,1},
                         {-(SPoints[2].x*SPoints[2].x+SPoints[2].y*SPoints[2].y+SPoints[2].z*SPoints[2].z),SPoints[2].y,SPoints[2].z,1},
                         {-(SPoints[3].x*SPoints[3].x+SPoints[3].y*SPoints[3].y+SPoints[3].z*SPoints[3].z),SPoints[3].y,SPoints[3].z,1}};
        float D = determ(d,4)/T;

        float e[4][4] = {{SPoints[0].x,-(SPoints[0].x*SPoints[0].x+SPoints[0].y*SPoints[0].y+SPoints[0].z*SPoints[0].z),SPoints[0].z,1},
                         {SPoints[1].x,-(SPoints[1].x*SPoints[1].x+SPoints[1].y*SPoints[1].y+SPoints[1].z*SPoints[1].z),SPoints[1].z,1},
                         {SPoints[2].x,-(SPoints[2].x*SPoints[2].x+SPoints[2].y*SPoints[2].y+SPoints[2].z*SPoints[2].z),SPoints[2].z,1},
                         {SPoints[3].x,-(SPoints[3].x*SPoints[3].x+SPoints[3].y*SPoints[3].y+SPoints[3].z*SPoints[3].z),SPoints[3].z,1}};
        float E = determ(e,4)/T;

        float f[4][4] = {{SPoints[0].x,SPoints[0].y,-(SPoints[0].x*SPoints[0].x+SPoints[0].y*SPoints[0].y+SPoints[0].z*SPoints[0].z),1},
                         {SPoints[1].x,SPoints[1].y,-(SPoints[1].x*SPoints[1].x+SPoints[1].y*SPoints[1].y+SPoints[1].z*SPoints[1].z),1},
                         {SPoints[2].x,SPoints[2].y,-(SPoints[2].x*SPoints[2].x+SPoints[2].y*SPoints[2].y+SPoints[2].z*SPoints[2].z),1},
                         {SPoints[3].x,SPoints[3].y,-(SPoints[3].x*SPoints[3].x+SPoints[3].y*SPoints[3].y+SPoints[3].z*SPoints[3].z),1}};
        float F = determ(f,4)/T;

        float g[4][4] = {{SPoints[0].x,SPoints[0].y,SPoints[0].z,-(SPoints[0].x*SPoints[0].x+SPoints[0].y*SPoints[0].y+SPoints[0].z*SPoints[0].z)},
                         {SPoints[1].x,SPoints[1].y,SPoints[1].z,-(SPoints[1].x*SPoints[1].x+SPoints[1].y*SPoints[1].y+SPoints[1].z*SPoints[1].z)},
                         {SPoints[2].x,SPoints[2].y,SPoints[2].z,-(SPoints[2].x*SPoints[2].x+SPoints[2].y*SPoints[2].y+SPoints[2].z*SPoints[2].z)},
                         {SPoints[3].x,SPoints[3].y,SPoints[3].z,-(SPoints[3].x*SPoints[3].x+SPoints[3].y*SPoints[3].y+SPoints[3].z*SPoints[3].z)}};
        float G = determ(g,4)/T;

        S.p.x = -D / 2.0f;
        S.p.y = -E / 2.0f;
        S.p.z = -F / 2.0f;
        S.r = 1.0f/2.0f * sqrtf(D * D + E * E + F * F - 4.0f * G);

        S.exist= true;

        std::cout << "sphere center: " << S.p.x << " " << S.p.y << " " << S.p.z << std::endl;
        std::cout << "with point0: " << SPoints[0].x << " " << SPoints[0].y << " " << SPoints[0].z << " distance: " << distance(S.p, SPoints[0] ) <<std::endl;
        std::cout << "with point1: " << SPoints[1].x << " " << SPoints[1].y << " " << SPoints[1].z << " distance: " << distance(S.p, SPoints[1] ) <<std::endl;
        std::cout << "with point2: " << SPoints[2].x << " " << SPoints[2].y << " " << SPoints[2].z << " distance: " << distance(S.p, SPoints[2] ) <<std::endl;
        std::cout << "with point3: " << SPoints[3].x << " " << SPoints[3].y << " " << SPoints[3].z << " distance: " << distance(S.p, SPoints[3] ) <<std::endl;
        std::cout << "radius: " << S.r << std::endl  << std::endl;
    }
    else
    {
        S.p.x = 0.0f;
        S.p.y = 0.0f;
        S.p.z = 0.0f;
        S.r = 0.0f;
        S.exist= false;
    }
    return S;
}

Point midpoint(Point p0, Point p1)
{
    Point p;
    p.x = (p0.x+p1.x)/2.0f;
    p.y = (p0.y+p1.y)/2.0f;
    p.z = (p0.z+p1.z)/2.0f;
    return p;
}

Sphere bound(Point* sPoints, int numSPoints)
{
    Sphere S;

    if (numSPoints == 1)
    {
        Point p0 = sPoints[0];
        S.p = p0;
        S.r = 0.0f;
        S.exist = true;
        std::cout << "p0: (" << p0.x << "," << p0.y << "," << p0.z << ")" << std::endl;
        std::cout << "sphere center (1): " << S.p.x << " " << S.p.y << " " << S.p.z << std::endl;
        std::cout << "radius: " << S.r << std::endl << std::endl;
        return S;
    }
    else if (numSPoints == 2)
    {
        Point p0 = sPoints[0];
        Point p1 = sPoints[1];
        Point center = midpoint(p0,p1);
        float diameter = distance(p0, p1);
        S.p = center;
        S.r = diameter / 2.0f;
        S.exist = true;
        std::cout << "p0: (" << p0.x << "," << p0.y << "," << p0.z << ")" << std::endl;
        std::cout << "p1: (" << p1.x << "," << p1.y << "," << p1.z << ")" << std::endl;
        std::cout << "SC: (" << S.p.x << ", " << S.p.y << ", " << S.p.z << ")" << std::endl;
        std::cout << "radius: " << S.r << std::endl << std::endl;
        return S;
    }
    else if (numSPoints == 3)
    {
        S = calcCircle(sPoints);
        return S;
    }
    else if(numSPoints == 4)
    {
        S = calcSphere(sPoints);
        return S;
    }
    else
    {
        S.exist= false;
        return S; //undefined
    }

}

Sphere sed(Point* points, Point* sPoints, uint32_t numPoints, uint32_t numSPoints, float sphereSize, int r)
{
    Sphere S;
    if(numPoints > 0)
    {
        Sphere D;
        Sphere E;
        Point newPoints[numPoints - 1];
        Point newSPoints[numSPoints + 1];
        int p = rand() % numPoints;

        for (int i = 0; i < numPoints - 1; ++i)
        {
            if (i < p)
            {
                newPoints[i] = points[i];
            }
            else
            {
                newPoints[i] = points[i + 1];
            }
        }

        if(numSPoints > 0)
        {
            D = bound(sPoints, numSPoints);
            // is point p in sphere D?
            if (distance(points[p], D.p) <= D.r) {
                std::cout << "depth: " << r << ", " << numSPoints << " rim, new point inside: " << "to test: "
                          << numPoints << " radius: " << D.r << std::endl << std::endl;
                E = sed(newPoints, sPoints, numPoints - 1, numSPoints, sphereSize, ++r);

                if(E.exist) {
                    return E;
                }
                else
                {
                    return D;
                }
            }
        }



        float tmpD;
        float maxD = 0;
        int pMax = 0;

        newSPoints[0] = points[p];
        for (int i = 0; i < numSPoints; ++i)
        {
            newSPoints[i + 1] = sPoints[i];
            tmpD = distance(sPoints[i], points[p]);
            if (tmpD > maxD)
            {
                maxD = tmpD;
                pMax = i;
            }
        }

        if(numSPoints > 0)
        {
            newSPoints[pMax + 1] = newSPoints[1];
            newSPoints[1] = sPoints[pMax];
        }

        D = bound(newSPoints,2);
        //else if(numSPoints+1 > 1)
        //{
        //    D = sed(newPoints, newSPoints, numPoints - 1, 2, sphereSize, r++);
        //}

        //inside?
        bool oldSPin = true;
        for (int i = 2; i < numSPoints+1; ++i)
        {
            if (distance(newSPoints[i], D.p)/2.0f > D.r)
            {
                oldSPin = false;
            }
        }
        if (oldSPin && numSPoints+1 > 1)
        {
            if(numSPoints+1 > 2)
            {
                std::cout << "points: " << numPoints + numSPoints+1 - 3 << "numSP: " << numSPoints+1 << std::endl;
                Point resetPoints[numPoints + numSPoints+1 - 3];
                for (int i = 0; i < numPoints + numSPoints+1 - 3; ++i)
                {
                    if (i < numPoints - 1)
                    {
                        resetPoints[i] = points[i];
                        std::cout << "point " << i << ">> (" << resetPoints[i].x << "," << resetPoints[i].y << ","<< resetPoints[i].z << ")" << std::endl;
                    }
                    else
                    {
                        resetPoints[i] = newSPoints[i - numPoints + 3];
                        std::cout << "point " << i << " in sP " << i - numPoints + 3 << " >> (" << resetPoints[i].x << "," << resetPoints[i].y << ","<< resetPoints[i].z << ")" << std::endl;
                    }
                }

                std::cout << "depth: " << r << ", 2 points rim: " << "to test: " << numPoints << " radius: " << D.r
                          << std::endl << std::endl;
                E = sed(resetPoints, newSPoints, numPoints + numSPoints+1 - 3, 2, sphereSize, ++r);


                if(E.exist) {
                    return E;
                }
                else
                {
                    return D;
                }
            }
            else
            {
                std::cout << "depth: " << r << ", II points rim: " << "to test: " << numPoints << " radius: " << D.r
                          << std::endl << std::endl;
                E = sed(newPoints, newSPoints, numPoints - 1, 2, sphereSize, ++r);


                if(E.exist) {
                    return E;
                }
                else
                {
                    return D;
                }
            }
        }
        else
        {
            D = bound(newSPoints,numSPoints + 1);
            std::cout << "depth: "<< r << ", " << numSPoints+1 <<" points rim: " << "to test: " << numPoints << " radius: " << D.r <<  std::endl  << std::endl;
            E = sed(newPoints, newSPoints, numPoints - 1, numSPoints + 1, sphereSize, ++r);
            if(E.exist) {
                return E;
            }
            else
            {
                return D;
            }
        }



    }
    S.exist= false;
    return S; //undefined
}

Sphere calculateBoundingSphere(const Point* points, uint32_t numPoints)
{
    Sphere result;
    {
        srand (time(NULL));
        Point pointSet[numPoints];
        Point sPoints[0];
        // copy Points because Points must not be changed
        for (int i = 0; i < numPoints; ++i)
        {
            pointSet[i] = points[i];
        }
        result = sed(pointSet,sPoints,numPoints,0,0,0);
    }
    return result;
}


int main()
{
    Point points[10];
    points[0].x=1;
    points[0].y=3;
    points[0].z=2;

    points[1].x=10;
    points[1].y=1;
    points[1].z=1;

    points[2].x=4;
    points[2].y=2;
    points[2].z=1;

    points[3].x=6;
    points[3].y=3;
    points[3].z=4;

    points[4].x=3;
    points[4].y=11;
    points[4].z=9;

    points[5].x=11;
    points[5].y=7;
    points[5].z=2;

    points[6].x=14;
    points[6].y=9;
    points[6].z=9;

    points[7].x=14;
    points[7].y=21;
    points[7].z=1;

    points[8].x=1;
    points[8].y=0;
    points[8].z=-9;

    points[9].x=-14;
    points[9].y=9;
    points[9].z=5;
    Sphere result = calculateBoundingSphere(points,10);

    std::cout << "bounding sphere center: " <<  result.p.x << " " <<  result.p.y << " " <<  result.p.z << std::endl;
    std::cout << "radius: " <<  result.r << std::endl;
    for(int i = 0; i < 10; ++i)
    {
        std::cout << "point " << i << " is";
        if(distance(points[i],result.p) <= result.r)
        {
            std::cout << " with distance " << distance(points[i],result.p) << "       inside with coordinates: " << points[i].x << " " << points[i].y << " " << points[i].z << std::endl;
        }
        else
        {
            std::cout << " with distance " << distance(points[i],result.p) << " !!NOT inside!! with coordinates: " << points[i].x << " " << points[i].y << " " << points[i].z << std::endl;
        }
    }
    return 0;
}

