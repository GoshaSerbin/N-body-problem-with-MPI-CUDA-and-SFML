#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include "mpi.h"
#include <iomanip>

using namespace std;

typedef double T;

#define  MAX(a,b) ((a)>(b)) ? (a) : (b)

const T G = 6.67e-20;
const T eps = 1e-5;

const bool createNewRandomFile = true;
const bool writeTrajFiles = true;

std::random_device rd;
std::mt19937 gen(rd());

#define MAX_VAL 1e13
std::uniform_real_distribution<> distribution(0, MAX_VAL);

T fRand(T fMin, T fMax)
{
    T f = distribution(gen) / MAX_VAL;
    return fMin + f * (fMax - fMin);
}

struct point
{
    
    T coord[3];
    T vel[3];
    T mass;


    point& operator= (const point& p)
    {
        mass = p.mass;
        for (int i = 0; i < 3; ++i) {
            coord[i] = p.coord[i];
            vel[i] = p.vel[i];
        }
        return *this;
    }
};

std::ostream& operator<< (std::ostream& out, const point& point)
{
    out << point.coord[0] << " " << point.coord[1] << " " << point.coord[2];
    return out;
}

int getNum(const char* const fileName)
{
    int n;
    ifstream ifile;
    ifile.open(fileName);
    if (!ifile.is_open())
    {
        cerr << " Error : file with settings is not open !\n";
        return -1;
    }
    ifile >> n;
    ifile.close();
    return n;
}

int readFile(const int n, const string fileName, point* points)
{
    ifstream ifile;
    ifile.open(fileName);
    string str;
    if (!ifile.is_open())
    {
        cerr << " Error : file with settings is not open !\n";
        return -1;
    }
    getline(ifile, str);
    for (int i = 0; i < n; ++i)
    {
        ifile >> points[i].mass;
        for (int j = 0; j < 3; ++j)
            ifile >> points[i].coord[j];

        for (int j = 0; j < 3; ++j)
            ifile >> points[i].vel[j];
    }
    ifile.close();
    return 0;
}

int createRandomFile(const int n, const string fileName)
{
    ofstream ofile;
    ofile.open(fileName);
    string str;
    if (!ofile.is_open())
    {
        cerr << " Error : file with settings is not open !\n";
        return -1;
    }
    ofile << n << endl;
    for (int i = 0; i < n; ++i)
    {
        ofile << fRand(9e9, 10e9) << " ";//mass
        for (int j = 0; j < 3; ++j)
            ofile << fRand(-3, 3) << " "; //coord
        for (int j = 0; j < 3; ++j)
            ofile << fRand(-0.3, 0.3) << " ";//vel
        ofile << endl;
    }
    ofile.close();
    return 0;
}


T* f(point* points, const int N0, const int n, T* accel, T* result) {
                                    //y = {rx,ry,rz,vx,vy,vz}
    for (int i = 0; i < n * 6; i += 6) {
        result[i] = points[N0 + i / 6].vel[0];
        result[i + 1] = points[N0 + i / 6].vel[1];
        result[i + 2] = points[N0 + i / 6].vel[2];
        result[i + 3] = accel[i / 2];
        result[i + 4] = accel[i / 2 + 1];
        result[i + 5] = accel[i / 2 + 2];
    }

    return result;
}


T dnorm(const T* vec1, const T* vec2)
{
    T sum = 0;
    for (int i = 0; i < 3; ++i)
    {
        sum += pow(vec1[i] - vec2[i], 2);
    }
    return sqrt(sum);
}

void calculateAccelerations(const int N, const point* points, const int N0, const int n, T* accel)
{
    for (int pi = 0; pi < n; ++pi) {
        point p = points[N0 + pi];
        for (int dim = 0; dim < 3; ++dim)
        {
            T a = 0;
            for (int i = 0; i < N; ++i)
            {
                T denominator = pow(dnorm(p.coord, points[i].coord), 3);
                a +=  points[i].mass * (p.coord[dim] - points[i].coord[dim]) / (MAX(denominator, eps));
            }

            accel[3 * pi + dim] = -G * a;
        }
    }

}

void show(const point& point)
{
    cout << "\n\np mass = " << point.mass;
    cout << "\np coord = { " << point.coord[0] << " " << point.coord[1] << " " << point.coord[2] <<" }";
    cout << "\np vel = { " << point.vel[0] << " " << point.vel[1] << " " << point.vel[2] << " }\n"<<endl;
}


T simulate(const int N, const int n, point* points, const T tmax, const T tau)
{

    int np, myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    ofstream* ofiles = nullptr; 

    int N0 = myid * (N / np);

    if (writeTrajFiles) {
        ofiles = new ofstream[n];
        for (int i = 0; i < n; ++i)
        {
            ofiles[i].open("my_traj" + to_string(N0 + i + 1) + ".txt");
        }
        ofstream dataFile;
        dataFile.open("data.txt");
        dataFile << setprecision(13) << N << " " << tmax << " " << tau <<std::endl;
        dataFile.close();
    }
    point* ps1 = new point[n];

    T* k1 = new T[n * 6];
    T* k2 = new T[n * 6];
   
    T* accel = new T[n * 3];


    int* displs = nullptr;
    int* recvcounts = nullptr;

    displs = new int[np];
    recvcounts = new int[np];
    for (int i = 0; i < np - 1; ++i) {
        recvcounts[i] = N / np;
        displs[i] = N / np * i;
    }
    displs[np - 1] = N / np * (np - 1); recvcounts[np - 1] = N / np + (N % np);

    MPI_Datatype mpi_point;
    int length[3] = { 3, 3, 1 };
    MPI_Aint disp[3] = { offsetof(point,coord) ,offsetof(point,vel), offsetof(point,mass) };
    MPI_Datatype type[3] = { MPI_DOUBLE, MPI_DOUBLE ,MPI_DOUBLE };
    MPI_Type_create_struct(3, length, disp, type, &mpi_point);
    MPI_Type_commit(&mpi_point);

    MPI_Datatype mpi_point_reduced;
    int len2[2] = { 3,1 };
    MPI_Aint pos2[2] = { offsetof(point,coord), sizeof(point) };
    MPI_Datatype typ2[2] = { MPI_DOUBLE, MPI_UB };
    MPI_Type_create_struct(2, len2, pos2, typ2, &mpi_point_reduced);
    MPI_Type_commit(&mpi_point_reduced);

    double mpi_t = -MPI_Wtime();

    const T tau05 = tau / 2;
    const T tlim = tmax + tau05;

    for (T t = 0; t < tlim; t += tau)
    {

        if (writeTrajFiles ) {
            for (int pi = 0; pi < n; ++pi) {

                ofiles[pi] << setprecision(13) << t << " " << points[N0 + pi] << std::endl;

            }

        }

        calculateAccelerations(N, points, N0, n, accel);

        
        f(points, N0, n, accel, k1);
        for (int pi = 0; pi < n; ++pi)
        {
            ps1[pi] = points[pi+N0];


            for (int j = 0; j < 3; ++j) {                   
                points[pi+N0].coord[j] = ps1[pi].coord[j] + tau05 * k1[6*(pi) + j];
                points[pi+N0].vel[j] = ps1[pi].vel[j] + tau05 * k1[6*(pi)+ 3 + j];
            }

        }

        
  
        MPI_Allgatherv(points + N0, n, mpi_point_reduced, points, recvcounts, displs, mpi_point_reduced, MPI_COMM_WORLD);

        calculateAccelerations(N, points, N0, n, accel);

        f(points, N0, n, accel, k2);


        for (int pi = 0; pi <  n; ++pi) {

            for (int j = 0; j < 3; ++j) {
                points[pi+N0].coord[j] = ps1[pi].coord[j] + tau * k2[6 * (pi) + j];
                points[pi+N0].vel[j] = ps1[pi].vel[j] + tau * k2[6 * (pi) + 3 + j];
            }

        }
        MPI_Allgatherv(points + N0, n, mpi_point_reduced, points, recvcounts, displs, mpi_point_reduced, MPI_COMM_WORLD);

    }
    mpi_t += MPI_Wtime();

    if (writeTrajFiles)
    {
        for (int i = 0; i < n; ++i)
        {
            ofiles[i].close();
        }
    }
        
    delete[] recvcounts;
    delete[] displs;
    delete[] ofiles;
    delete[] k1;
    delete[] k2;
    return mpi_t;
}

T getError(const int N)
{
    T t = 0;
    point p, my_p;
    T error = 0, max_error = 0;
    for (int i = 0; i < N; ++i)
    {
        ifstream ifile;
        ifstream my_ifile;
        ifile.open("traj" + to_string(i + 1) + ".txt");
        my_ifile.open("my_traj" + to_string(i + 1) + ".txt");
        int count = 0;
        while (count < 201) {
            ifile >> t >> p.coord[0] >> p.coord[1] >> p.coord[2];
            my_ifile >> t >> my_p.coord[0] >> my_p.coord[1] >> my_p.coord[2];
            error = MAX(MAX(fabs(p.coord[0] - my_p.coord[0]), fabs(p.coord[1] - my_p.coord[1])), fabs(p.coord[2] - my_p.coord[2]));
            if (error > max_error) {
                max_error = error;
            }
            ++count;
        }
        ifile.close();
        my_ifile.close();
    }
    return max_error;
}

int main(int argc, char** argv)
{

    MPI_Init(&argc, &argv);
    int myid, np;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if (createNewRandomFile && myid == 0) {
        createRandomFile(10000, "Nbody.txt");
    }
    const char* filename = "4body.txt";//4body.txt or Nbody.txt
    int N;
    
    MPI_Status st;

    if (myid == 0)
    {
        N = getNum(filename);
    }
    
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    point* points = new point[N];
    const int n = N / np + (N % np) * (myid == np - 1);


    const T tmax = 10*31536000;
    const T tau = 60*60;
    T t = 0;

    if (myid == 0)
    {           
        readFile(N, filename, points);  
        cout << "N = " << N << endl;
    }


    MPI_Datatype arrtype;
    MPI_Type_contiguous(7, MPI_DOUBLE, &arrtype);
    MPI_Type_commit(&arrtype);

    MPI_Bcast(points, N, arrtype, 0, MPI_COMM_WORLD);




    t = simulate(N, n, points, tmax, tau);

    if (myid == 0) {
        cout << "total time: " << t << endl;
        //if (filename == "4body.txt") 
        //{
        //    cout<< "error = " << getError(4);
        //}
    }
    

    MPI_Finalize();
    delete[] points;
    return 0;
}