#include <iostream>
#include <cmath> // for calculating square root and pow
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
using namespace std;

// a struct to store 3D Points with the ability to perform basic mathematical vector operations
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

// solves 4x4 determinant:
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

// calc Euclidean distance:
float distance(Point p0, Point p1)
{
    return sqrtf((p0.x-p1.x)*(p0.x-p1.x)+(p0.y-p1.y)*(p0.y-p1.y)+(p0.z-p1.z)*(p0.z-p1.z));
}


// calc circle from three points. circle becomes aequator of new sphere.
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

    /*
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
    // is the center of the circle on the Plane with the 3 points:
    //std::cout << ">>>>plane test: 0? " << a*S.p.x+b*S.p.y+c*S.p.z+n << std::endl << std::endl;
    */

    return S;
}


// calculate sphere from 4 points by solving:
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

        /*
        std::cout << "sphere center: " << S.p.x << " " << S.p.y << " " << S.p.z << std::endl;
        std::cout << "with point0: " << SPoints[0].x << " " << SPoints[0].y << " " << SPoints[0].z << " distance: " << distance(S.p, SPoints[0] ) <<std::endl;
        std::cout << "with point1: " << SPoints[1].x << " " << SPoints[1].y << " " << SPoints[1].z << " distance: " << distance(S.p, SPoints[1] ) <<std::endl;
        std::cout << "with point2: " << SPoints[2].x << " " << SPoints[2].y << " " << SPoints[2].z << " distance: " << distance(S.p, SPoints[2] ) <<std::endl;
        std::cout << "with point3: " << SPoints[3].x << " " << SPoints[3].y << " " << SPoints[3].z << " distance: " << distance(S.p, SPoints[3] ) <<std::endl;
        std::cout << "radius: " << S.r << std::endl  << std::endl;
         */
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

//point between two points
Point midpoint(Point p0, Point p1)
{
    Point p;
    p.x = (p0.x+p1.x)/2.0f;
    p.y = (p0.y+p1.y)/2.0f;
    p.z = (p0.z+p1.z)/2.0f;
    return p;
}

//bound a sphere with the set of sPoints (sphere Points) at its surface
Sphere bound(Point* sPoints, int numSPoints)
{
    Sphere S;

    if (numSPoints == 1)
    {
        Point p0 = sPoints[0];
        S.p = p0;
        S.r = 0.0f;
        S.exist = true;
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

//-D: current Sphere  (E for experimental Sphere)
//-P (ponits): Points to test
//-R (sPonints): possibly Points candidate for Sphere creation
//
////bounding sphere algorithm
//sed (Point* P, Point* R, int t, Sphere& D) {
//if (R > 3)
//  E = sed( P, R, 2, D) // recursive creation of experimental Sphere with two Points
//  if (E smaller then D) // (D was created with 4 Point)
//      D = E
//else
//  Point p = get random Point out of P
//  if(p inside Sphere D)
//      remove p from set P
//      E = sed( P, R, t, D) // recursive creation of experimental Sphere with Points in R
//      if (E.exist)
//         return E
//      else
//         return D
//  else
//      add p to set of Points R
//      D = bound( R ) // create Sphere with Points in set R
//      E = sed( P, R, t+1, D) // recursive creation of experimental Sphere with Points in R + new Point p
//      if (E.exist)
//         return E
//      else
//         return D
// }
Sphere sed(Point* points, Point* sPoints, uint32_t numPoints, uint32_t numSPoints, int inside, int r, Sphere* D)
{

    std::cout << "Points to test: " << numPoints-inside << " recursive depth: " << r << std::endl;
    if(numPoints-inside > 0)
    {
        Sphere E;
        Point newPoints[numPoints - 1];
        Point newSPoints[numSPoints + 1];

        //remove random Point p from set points (P)
        int p = (rand() % (numPoints-inside))+inside;

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

        // if sPoint (R) contains at least one Point try to test wheather p lies within the Sphere D
        if(numSPoints > 0)
        {
            // if sPoints (R) has more then 3 Points try to find a smaller Sphere created out of 2 of the sPoints.
            if(numSPoints > 3)
            {
                int resetsize = numPoints-1 + numSPoints - 2;
                Point resetPoints[resetsize];
                for (int i = 0; i < resetsize; ++i)
                {
                    if (i < numPoints - 1)
                    {
                        resetPoints[i] = newPoints[i];
                    }
                    else
                    {
                        resetPoints[i] = newSPoints[i - numPoints + 2];
                    }
                }

                E = sed(resetPoints, sPoints, numPoints, 2, 0, ++r, D);
                if(E.r < D->r)
                {
                    D = &E;
                }
            }

            // is point p in sphere D?
            if (D->exist && distance(points[p], D->p) <= D->r +0.0001)
            {
                Point tmp = points[inside];
                points[inside] = points[p];
                points[p] = tmp;
                E = sed(points, sPoints, numPoints, numSPoints, ++inside, ++r, D);

                // if there exists a Sphere E which includes all Points then return it or else return the old Sphere D
                if(E.exist)
                {
                    return E;
                }
                else
                {
                    return *D;
                }
            }
        }
        // point p is not in sphere D!

        float tmpD;
        float maxD = 0;
        int pMax = 0;

        // add point p to sPoints
        // find out which old sPoint has the maximum distance to p
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

        
        // place the point with the maximal distance to p at position 2 in the sPoint set
        if(numSPoints > 0)
        {
            newSPoints[pMax + 1] = newSPoints[1];
            newSPoints[1] = sPoints[pMax];
        }

        // try to find Sphere which contains all sPoints (R) and is defined by p and the sPoint with maximum distance to p
        if(numSPoints+1>1)
        {
            *D = bound(newSPoints, 2);
        }
        bool oldSPin = true;
        bool oldSPin3 = false;
        for (int i = 2; i < numSPoints + 1; ++i)
        {
            if (distance(newSPoints[i], D->p) > D->r)
            {
                oldSPin = false;
            }
        }

        // if D doesn't contain all sPoints,
        // try to find Sphere which contains all sPoints (R) and is defined by p and two other points in sPoints
        if (!oldSPin && numSPoints+1 > 2)
        {
            oldSPin3 = true;
            // try all combinations of p and 2 other points in sPoints
            for(int i = 2; i < numSPoints +1 ; ++i)
            {
                if(i>2)
                {
                    Point tmp = newSPoints[2];
                    newSPoints[2] = newSPoints[i];
                    newSPoints[i] = tmp;
                }
                *D = bound(newSPoints, 3);
                for (int i = 3; i < numSPoints + 1; ++i)
                {
                    if (distance(newSPoints[i], D->p) > D->r)
                    {
                        oldSPin3 = false;
                        oldSPin = true;
                    }
                }
                if(oldSPin3)
                {
                    break;
                }

            }
        }

        // if the solution p with one other Points in sPoints (R) works and we have more then one sPoints
        // try to add other points in the set points (P) to the bounding Sphere and check if they are inside
        // by calling this function sed(...) recursivly
        if (oldSPin && numSPoints+1 > 1)
        {
            if(numSPoints+1 > 2 && (!oldSPin3 || numSPoints+1 > 3))
            {
                int newNumS = 2;
                if (oldSPin3) {newNumS = 3;}
                int resetsize = numPoints-1 + numSPoints+1 - newNumS;
                Point resetPoints[resetsize];
                for (int i = 0; i < resetsize; ++i)
                {
                    if (i < numPoints - 1)
                    {
                        resetPoints[i] = newPoints[i];
                    }
                    else
                    {
                        resetPoints[i] = newSPoints[i - numPoints + newNumS+1];
                    }
                }

                E = sed(resetPoints, newSPoints, resetsize, newNumS, 0, ++r, D);


                if(E.exist)
                {
                    return E;
                }
                else
                {
                    return *D;
                }
            }
            else
            {
                // if a 3 Points solution for creating the bounding Sphere exists, then call sed(...) recursivly                

                if(oldSPin3)
                {
                    E = sed(newPoints, newSPoints, numPoints - 1, 3, 0, ++r, D);
                }
                else
                {
                    E = sed(newPoints, newSPoints, numPoints - 1, 2, 0, ++r, D);
                }



                if(E.exist) {
                    return E;
                }
                else
                {
                    return *D;
                }
            }
        }
        else
        {
            // if all previosly branches fail, then create a Sphere
            // with all current sPoints including Point P as Sphere defining Points
            *D = bound(newSPoints,numSPoints + 1);
            E = sed(newPoints, newSPoints, numPoints - 1, numSPoints + 1, 0, ++r, D);
            if(E.exist) {
                return E;
            }
            else
            {
                return *D;
            }
        }



    }
    D->exist= false;
    return *D; //undefined
}

// Interface for calling the recursive function sed(...)
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
        result = sed(pointSet,sPoints,numPoints,0,0,0, &result);
    }
    return result;
}


int main()
{
    // add 20 3D Points to test the algorithm
    Point points[20];
    points[0].x=1;
    points[0].y=5;
    points[0].z=4;

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
    points[9].y=-11;
    points[9].z=-7;

    points[10].x=1;
    points[10].y=3;
    points[10].z=2;

    points[11].x=10;
    points[11].y=22;
    points[11].z=11;

    points[12].x=14;
    points[12].y=2;
    points[12].z=21;

    points[13].x=16;
    points[13].y=3;
    points[13].z=14;

    points[14].x=23;
    points[14].y=41;
    points[14].z=19;

    points[15].x=11;
    points[15].y=17;
    points[15].z=12;

    points[16].x=75;
    points[16].y=-9;
    points[16].z=29;

    points[17].x=24.5;
    points[17].y=121;
    points[17].z=21;

    points[18].x=1;
    points[18].y=10;
    points[18].z=19;

    points[19].x=3;
    points[19].y=24;
    points[19].z=8;
    Sphere result = calculateBoundingSphere(points,20);

    std::cout << "bounding sphere center: " <<  "at Point ("<<result.p.x << "," <<  result.p.y << "," <<  result.p.z << ")" << std::endl;
    std::cout << "with radius: " <<  result.r << std::endl;
    for(int i = 0; i < 20; ++i)
    {
        std::cout << "point " << i << " is";
        if(distance(points[i],result.p) < result.r+0.0001)
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

