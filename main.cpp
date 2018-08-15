#include <iostream>
#include <cmath> // for calculating square root and pow nothing else. I hope that was allowed.
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
using namespace std;

struct Point
{
    float x, y, z;
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
//    | 0                     a   b   c   n  |  // plane where the three points lie in
Sphere calcCircle(Point* SPoints)
{
    Sphere S;
    // the plane of the three point:
    float a = (SPoints[1].y-SPoints[0].y)*(SPoints[2].z-SPoints[0].z)-(SPoints[2].y-SPoints[0].y)*(SPoints[1].z-SPoints[0].z);
    float b = (SPoints[1].z-SPoints[0].z)*(SPoints[2].x-SPoints[0].x)-(SPoints[2].z-SPoints[0].z)*(SPoints[1].x-SPoints[0].x);
    float c = (SPoints[1].x-SPoints[0].x)*(SPoints[2].y-SPoints[0].y)-(SPoints[2].x-SPoints[0].x)*(SPoints[1].y-SPoints[0].y);
    float n = -(a*SPoints[0].x+b*SPoints[0].y+c*SPoints[0].z);
    //std::cout << "plane test: 0?" << a*SPoints[1].x+b*SPoints[1].y+c*SPoints[1].z+n << std::endl;


    float t[4][4] = {{SPoints[0].x,SPoints[0].y,SPoints[0].z,1},
                     {SPoints[1].x,SPoints[1].y,SPoints[1].z,1},
                     {SPoints[2].x,SPoints[2].y,SPoints[2].z,1},
                     {a, b, c, n}};
    float T = determ(t,4);
    if (T!=0)
    {
        float d[4][4] = {{-(SPoints[0].x*SPoints[0].x+SPoints[0].y*SPoints[0].y+SPoints[0].z*SPoints[0].z),SPoints[0].y,SPoints[0].z,1},
                         {-(SPoints[1].x*SPoints[1].x+SPoints[0].y*SPoints[1].y+SPoints[1].z*SPoints[1].z),SPoints[1].y,SPoints[1].z,1},
                         {-(SPoints[2].x*SPoints[2].x+SPoints[0].y*SPoints[2].y+SPoints[2].z*SPoints[2].z),SPoints[2].y,SPoints[2].z,1},
                         {0, b, c, n}};
        float D = determ(d,4)/T;

        float e[4][4] = {{SPoints[0].x,-(SPoints[0].x*SPoints[0].x+SPoints[0].y*SPoints[0].y+SPoints[0].z*SPoints[0].z),SPoints[0].z,1},
                         {SPoints[1].x,-(SPoints[1].x*SPoints[1].x+SPoints[0].y*SPoints[1].y+SPoints[1].z*SPoints[1].z),SPoints[1].z,1},
                         {SPoints[2].x,-(SPoints[2].x*SPoints[2].x+SPoints[0].y*SPoints[2].y+SPoints[2].z*SPoints[2].z),SPoints[2].z,1},
                         {0, a, c, n}};
        float E = determ(e,4)/T;

        float f[4][4] = {{SPoints[0].x,SPoints[0].y,-(SPoints[0].x*SPoints[0].x+SPoints[0].y*SPoints[0].y+SPoints[0].z*SPoints[0].z),1},
                         {SPoints[1].x,SPoints[1].y,-(SPoints[1].x*SPoints[1].x+SPoints[0].y*SPoints[1].y+SPoints[1].z*SPoints[1].z),1},
                         {SPoints[2].x,SPoints[2].y,-(SPoints[2].x*SPoints[2].x+SPoints[0].y*SPoints[2].y+SPoints[2].z*SPoints[2].z),1},
                         {0, a, b, n}};
        float F = determ(f,4)/T;

        float g[4][4] = {{SPoints[0].x,SPoints[0].y,SPoints[0].z,-(SPoints[0].x*SPoints[0].x+SPoints[0].y*SPoints[0].y+SPoints[0].z*SPoints[0].z)},
                         {SPoints[1].x,SPoints[1].y,SPoints[1].z,-(SPoints[1].x*SPoints[1].x+SPoints[0].y*SPoints[1].y+SPoints[1].z*SPoints[1].z)},
                         {SPoints[2].x,SPoints[2].y,SPoints[2].z,-(SPoints[2].x*SPoints[2].x+SPoints[0].y*SPoints[2].y+SPoints[2].z*SPoints[2].z)},
                         {0, a, b, c}};
        float G = determ(g,4)/T;

        S.p.x = -D / 2.0f;
        S.p.y = -E / 2.0f;
        S.p.z = -F / 2.0f;
        S.r = 1.0f / 2.0f * sqrtf(D * D + E * E + F * F - 4.0f * G);

        S.exist= true;

        std::cout << "circle center: " << S.p.x << " " << S.p.y << " " << S.p.z << std::endl;
        std::cout << "with point0: " << SPoints[0].x << " " << SPoints[0].y << " " << SPoints[0].z << " distance: " << distance(S.p, SPoints[0] ) << std::endl;
        std::cout << "with point1: " << SPoints[1].x << " " << SPoints[1].y << " " << SPoints[1].z << " distance: " << distance(S.p, SPoints[1] ) << std::endl;
        std::cout << "with point2: " << SPoints[2].x << " " << SPoints[2].y << " " << SPoints[2].z << " distance: " << distance(S.p, SPoints[2] ) << std::endl;
        std::cout << "radius: " << S.r << std::endl;
        std::cout << ">>>>plane test: 0? " << a*S.p.x+b*S.p.y+c*S.p.z+n << std::endl;
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
                         {-(SPoints[1].x*SPoints[1].x+SPoints[0].y*SPoints[1].y+SPoints[1].z*SPoints[1].z),SPoints[1].y,SPoints[1].z,1},
                         {-(SPoints[2].x*SPoints[2].x+SPoints[0].y*SPoints[2].y+SPoints[2].z*SPoints[2].z),SPoints[2].y,SPoints[2].z,1},
                         {-(SPoints[3].x*SPoints[3].x+SPoints[0].y*SPoints[3].y+SPoints[3].z*SPoints[3].z),SPoints[3].y,SPoints[3].z,1}};
        float D = determ(d,4)/T;

        float e[4][4] = {{SPoints[0].x,-(SPoints[0].x*SPoints[0].x+SPoints[0].y*SPoints[0].y+SPoints[0].z*SPoints[0].z),SPoints[0].z,1},
                         {SPoints[1].x,-(SPoints[1].x*SPoints[1].x+SPoints[0].y*SPoints[1].y+SPoints[1].z*SPoints[1].z),SPoints[1].z,1},
                         {SPoints[2].x,-(SPoints[2].x*SPoints[2].x+SPoints[0].y*SPoints[2].y+SPoints[2].z*SPoints[2].z),SPoints[2].z,1},
                         {SPoints[3].x,-(SPoints[3].x*SPoints[3].x+SPoints[0].y*SPoints[3].y+SPoints[3].z*SPoints[3].z),SPoints[3].z,1}};
        float E = determ(e,4)/T;

        float f[4][4] = {{SPoints[0].x,SPoints[0].y,-(SPoints[0].x*SPoints[0].x+SPoints[0].y*SPoints[0].y+SPoints[0].z*SPoints[0].z),1},
                         {SPoints[1].x,SPoints[1].y,-(SPoints[1].x*SPoints[1].x+SPoints[0].y*SPoints[1].y+SPoints[1].z*SPoints[1].z),1},
                         {SPoints[2].x,SPoints[2].y,-(SPoints[2].x*SPoints[2].x+SPoints[0].y*SPoints[2].y+SPoints[2].z*SPoints[2].z),1},
                         {SPoints[3].x,SPoints[3].y,-(SPoints[3].x*SPoints[3].x+SPoints[0].y*SPoints[3].y+SPoints[3].z*SPoints[3].z),1}};
        float F = determ(f,4)/T;

        float g[4][4] = {{SPoints[0].x,SPoints[0].y,SPoints[0].z,-(SPoints[0].x*SPoints[0].x+SPoints[0].y*SPoints[0].y+SPoints[0].z*SPoints[0].z)},
                         {SPoints[1].x,SPoints[1].y,SPoints[1].z,-(SPoints[1].x*SPoints[1].x+SPoints[0].y*SPoints[1].y+SPoints[1].z*SPoints[1].z)},
                         {SPoints[2].x,SPoints[2].y,SPoints[2].z,-(SPoints[2].x*SPoints[2].x+SPoints[0].y*SPoints[2].y+SPoints[2].z*SPoints[2].z)},
                         {SPoints[3].x,SPoints[3].y,SPoints[3].z,-(SPoints[3].x*SPoints[3].x+SPoints[0].y*SPoints[3].y+SPoints[3].z*SPoints[3].z)}};
        float G = determ(g,4)/T;

        S.p.x = -D / 2.0f;
        S.p.y = -E / 2.0f;
        S.p.z = -F / 2.0f;
        S.r = 1.0f / 2.0f * sqrtf(D * D + E * E + F * F - 4.0f * G);

        S.exist= true;

        std::cout << "sphere center: " << S.p.x << " " << S.p.y << " " << S.p.z << std::endl;
        std::cout << "with point0: " << SPoints[0].x << " " << SPoints[0].y << " " << SPoints[0].z << " distance: " << distance(S.p, SPoints[0] ) <<std::endl;
        std::cout << "with point1: " << SPoints[1].x << " " << SPoints[1].y << " " << SPoints[1].z << " distance: " << distance(S.p, SPoints[1] ) <<std::endl;
        std::cout << "with point2: " << SPoints[2].x << " " << SPoints[2].y << " " << SPoints[2].z << " distance: " << distance(S.p, SPoints[2] ) <<std::endl;
        std::cout << "with point3: " << SPoints[3].x << " " << SPoints[3].y << " " << SPoints[3].z << " distance: " << distance(S.p, SPoints[3] ) <<std::endl;
        std::cout << "radius: " << S.r << std::endl;
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

Sphere sed(Point* points, Point* sPoints, uint32_t numPoints, uint32_t numSPoints, float sphereSize, int opt)
{
    Sphere S;

    if(numPoints == 0 || numSPoints >= 4)
    {
        if (numSPoints == 1)
        {
            Point p0 = sPoints[0];
            S.p = p0;
            S.r = 0.0f;
            S.exist = true;
            std::cout << "sphere center (1): " << S.p.x << " " << S.p.y << " " << S.p.z << std::endl;
            std::cout << "radius: " << S.r << std::endl;
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
            std::cout << "sphere center (2): " << S.p.x << " " << S.p.y << " " << S.p.z << std::endl;
            std::cout << "radius: " << S.r << std::endl;
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

    if(numPoints > 0) {
        int p;
        Sphere D;
        Point newPoints[numPoints - 1];
        Point newSPoints[numSPoints + 1];
        if (opt > 0 && numSPoints > opt)
        {
            p = rand() % (numSPoints-opt);
            D = sed(points, sPoints, numPoints , numSPoints, sphereSize, opt--);
            // is point p in sphere D?
            if (D.exist && distance(sPoints[p+opt], D.p) <= D.r)
            {
                return D;
            }

            for (int i = 0; i < numPoints - 1; ++i)
            {
                newPoints[i] = points[i];
            }

            newSPoints[0] = sPoints[p];
            for (int i = 0; i < numSPoints; ++i)
            {
                if (i < p)
                {
                    newSPoints[i + 1] = sPoints[i];
                }
                else
                {
                    newSPoints[i] = sPoints[i];
                }
            }

        }
        else
        {

            p = rand() % numPoints;
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

            D = sed(newPoints, sPoints, numPoints - 1, numSPoints, sphereSize, opt--);
            // is point p in sphere D?
            if (D.exist && distance(points[p], D.p) <= D.r)
            {
                return D;
            }

            newSPoints[0] = points[p];
            for (int i = 0; i < numSPoints; ++i)
            {
                newSPoints[i + 1] = sPoints[i];
            }

        }

        //looking for a smaller alternative
        if(numSPoints > 2)
        {
            float newRadius = distance(sPoints[0],sPoints[1]);
            if (newRadius < D.r && newRadius > sphereSize)
            {
                Sphere E = sed(points, newSPoints, numPoints, 2, newRadius, numSPoints);
                if(E.exist)
                {
                    return E;
                }
            }

        }
        else
        {
            if(numSPoints > 4)
            {
                Sphere E = sed(points, newSPoints, numPoints, 2, sphereSize, numSPoints);
                if(E.exist)
                {
                    return E;
                }
            }
            else
            {
                D = sed(newPoints, newSPoints, numPoints - 1, numSPoints + 1, D.r, opt--);
                return D;
            }
        }

    }
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