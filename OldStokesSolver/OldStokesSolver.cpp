// SphericalHarmonics.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "MobilityProblem.h"
#include <array>
#include <iomanip>
#include <time.h>
#include <fstream>
#include<type_traits>


void testMobilityWilson3SpherePerp(double distance, double trueMUF, double trueMOmegaF, int numterms, std::fstream& f);


int main()
{

    const double PI = 4.0 * atan(1.0);

    /*
        //testInsideTranslation({ 1.0 , 1.0 , 1.0 }, 3);
        VReal vr;
        VImag vi;
        WReal wr;
        WImag wi;
        XReal xr;
        XImag xi;


        std::cout << "VReal\n";
        std::cout << "n = " << 0 << "\nm = " << 0 << "\n";
        vr.reset(0, 0);
        testOperators(&vr);


        for(int n = 1; n < 6; n++)
            for (int m = 0; m <= n; m++)
            {

                std::cout << "VReal\n";
                std::cout << "n = " << n << "\nm = " << m << "\n";
                vr.reset(m, n);
                testOperators(&vr);
                if (m > 0)
                {
                    std::cout << "VImag\n";
                    std::cout << "n = " << n << "\nm = " << m << "\n";
                    vi.reset(m, n);
                    testOperators(&vi);
                }
                std::cout << "WReal\n";
                std::cout << "n = " << n << "\nm = " << m << "\n";
                wr.reset(m, n);
                testOperators(&wr);

                if (m > 0)
                {
                    std::cout << "WImag\n";
                    std::cout << "n = " << n << "\nm = " << m << "\n";
                    wi.reset(m, n);
                    testOperators(&wi);
                }
                std::cout << "XReal\n";
                std::cout << "n = " << n << "\nm = " << m << "\n";
                xr.reset(m, n);
                testOperators(&xr);
                if (m > 0)
                {
                    std::cout << "XImag\n";
                    std::cout << "n = " << n << "\nm = " << m << "\n";
                    xi.reset(m, n);
                    testOperators(&xi);
                }
            }

     */

     //testTranslation({ 0.0 , 0.0 , 1.0 });
 /*
     std::array<int , 4> legendreindices = { 3 , 5 , 7 ,9 };

     RectCoord direction;
     RectCoord axis;
     std::cout << "Enter the three coordinates(rectangular) of the direction of flow. The direction will be normalized: ";
     std::cin >> direction.x >> direction.y >> direction.z;

     std::cout << "Enter the three coordinates(rectangular) of the axis of rotation. The direction will be normalized: ";
     std::cin >> axis.x >> axis.y >> axis.z;

     std::array<double, 8> rs = {0.01, 0.05, 0.1, 1 , 4 , 8 ,1000.0,1000000.0};
     for (double r : rs)
         for(int i : legendreindices)
         {
             std::cout << "r = " << r << "\n";
             std::cout << "Legendre node " << i << " of 16, Trapezoidal node 5 of 16.\n\n";
              SurfaceCoord s(PI / 2.0 * (GLnodes[i] + 1), 2.0 * PI * (double)5 / (double)NUMTRAPNODES);
              testOperatorsPointwise(SphereCoord(1.0, s),direction,axis, r);
         }

    */
    //testMobilityWilson3Sphere(4.0);



    //twoSphereMobility(100.0, Fs, Ts);

    //SolveBIEWithKnown(5);

 /*
    SingleParticle solution1, solution2, solutionSum;

    solution1 = testMobilitySingle(Fs[0], Ts[0]);
    solution2 = testMobilitySingle(Fs[1], Ts[1]);
    solutionSum = testMobilitySingle(coef1 *Fs[0] + coef2*Fs[1], coef1 * Ts[0] + coef2 *Ts[1]);

    std::cout << (solutionSum - solution1 * coef1 - solution2 * coef2).norm() << std::endl;;
   */
   //testMobilitySingle(F1, 0.0);




   //checkDestructable();


    std::fstream outputdata("output_data8.txt");
    std::array<double, 14> rs = { 3.0 , 2.9 , 2.8 , 2.7 , 2.6 , 2.5 , 2.4 , 2.3, 2.25 , 2.2,2.15,2.1,2.05 , 2.01 };
    std::array<double, 14> MUFs = { 1.53416156 , 1.55461831 ,  1.57664602, 1.60039150, 1.62599574 , 1.65356957 ,
                                   1.68314045,1.71452988,1.73072508, 1.74703222, 1.76311004, 1.77826951, 1.79070892,1.79223228 };
    std::array<double, 14> MOmegaFs = { 0.140187608,0.149081378,0.158630515,0.168801210,0.179481457,0.190411103,
                                        0.201041187 ,0.210224103 ,0.213577673 ,0.215412007 , 0.214756269 , 0.209731929 , 0.195478184 , 0.159607490 };

    std::time_t t0 = time(0);
    for (int i = 0; i < 14; i++)
    {
        testMobilityWilson3SpherePerp(rs[i], MUFs[i], MOmegaFs[i], 8, outputdata);
        std::cout << "s = " << rs[i] << "\n";
        std::cout << " total time elapsed: " << time(0) - t0 << "\n";
    }

    outputdata.close();

    //testMobilityWilson3SpherePerp(2.05, 1.79070892, 0.195478184,8);
     //SolveMobilitySingleSphere(RectCoord(0.0, 0.0, 6.0 * PI), RectCoord(), 8);
    return 0;


}

void testMobilityWilson3SpherePerp(double distance, double trueMUF, double trueMOmegaF, int numterms, std::fstream& f)
{
    const double PI = 4.0 * std::atan(1.0);

    RectCoord north(0.0, 1.0, 0.0);
    RectCoord southwest(std::cos(7.0 / 6.0 * PI), std::sin(7.0 / 6.0 * PI), 0.0); //points on an equilateral trangle.
    RectCoord southeast(std::cos(-PI / 6.0), std::sin(-PI / 6.0), 0.0);

    double normalizer = norm(north - southeast); // get the side length

    std::vector<RectCoord> centers;
    centers.push_back(distance / normalizer * north);
    centers.push_back(distance / normalizer * southwest); // create centers with given side length.
    centers.push_back(distance / normalizer * southeast);

    //centers.push_back(RectCoord(0.0, 0.0, 10.0));

    std::vector<RectCoord> forces;
    forces.push_back(RectCoord(0.0, 0.0, 6.0 * PI)); //force on first particle, pointing to the origin.
    forces.push_back(RectCoord(0.0, 0.0, 6.0 * PI));                    //force on second particle = 0.
    forces.push_back(RectCoord(0.0, 0.0, 6.0 * PI));                    //force on third particle = 0.
    //forces.push_back(RectCoord());

    std::vector<RectCoord> torques;
    torques.push_back(RectCoord());
    torques.push_back(RectCoord());                  //zero torque on all particles.
    torques.push_back(RectCoord());
    //torques.push_back(RectCoord());

    time_t time = std::time(0);
    time_t time0 = time;

    ParticleSystem ps(centers, 1);
    //system = std::move(SolveMobilityManySphere(centers, forces, torques, 1,system, 30,1e-5));
   //std::cout << "Time for first solve (N = 1): " << std::time(0) - time << "\n";

    double U3 = 1.0, U3prev = 0.0;

    MathConstants consts(2 * numterms + 1, 2 * numterms + 1);

    ps = ParticleSystem(centers, numterms);
    time = std::time(0);
    ps = SolveMobilityManySphere(centers, forces, torques, numterms, ps, 30, 1e-5);
    f << "Time for solve (N = " << numterms << "): " << std::time(0) - time << "\n";

    RectCoord particle2velocity = Integrate(&ps, 1.0, ps.particles[1].center, consts) / (4.0 * PI);
    RectCoord particle1rotation = rotationCoefficient(discretize(&ps, ps.particles[0].center, consts), ps.particles[0].center, 1.0, consts) * 0.375 / PI;
    RectCoord particle2rotation = rotationCoefficient(discretize(&ps, ps.particles[1].center, consts), ps.particles[1].center, 1.0, consts) * 0.375 / PI;
    RectCoord particle3rotation = rotationCoefficient(discretize(&ps, ps.particles[2].center, consts), ps.particles[2].center, 1.0, consts) * 0.375 / PI;

    f << "MUF = " << std::setprecision(16) << Integrate(&ps, 1.0, ps.particles[0].center, consts) / (4.0 * PI) << "\n";
    f << "MUF error = " << (Integrate(&ps, 1.0, ps.particles[0].center, consts) / (4.0 * PI) - RectCoord(0.0, 0.0, trueMUF)) / trueMUF << "\n";
    f << "MOmegaF = " << norm(particle2rotation) << "\n";
    f << "MOmegaF error = " << (norm(particle2rotation) - trueMOmegaF) / trueMOmegaF << "\n";




    //std::cout << "Force on particle 1 after solving: " << Integrate(&system.traction, 1.0, system.particles[0].center) <<"\n";
    //std::cout << "Force on particle 2 after solving: " << Integrate(&system.traction,  1.0, system.particles[1].center) <<"\n";
    //std::cout << "Force on particle 3 after solving: " << Integrate(&system.traction, 1.0, system.particles[2].center) <<"\n";


}