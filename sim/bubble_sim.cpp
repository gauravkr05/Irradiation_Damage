#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <cmath>
#include <random>

using namespace std;

const int NX = 100;
const int NY = 100;
const int NZ = 100;

const int N = NX * NY * NZ;

int id(int x, int y, int z)
{
    return x + NX * (y + NY * z);
}

struct Bubble
{
    double x;
    double y;
    double z;
};

int main()
{
    // ------------------------------------------------
    // Number of initial bubbles
    // ------------------------------------------------

    const int nbubbles = 80;

    vector<Bubble> bubbles;

    mt19937 rng(12345);

    uniform_real_distribution<double> distX(5, NX - 5);
    uniform_real_distribution<double> distY(5, NY - 5);
    uniform_real_distribution<double> distZ(5, NZ - 5);

    // Random initial bubble positions

    for(int i = 0; i < nbubbles; i++)
    {
        Bubble b;

        b.x = distX(rng);
        b.y = distY(rng);
        b.z = distZ(rng);

        bubbles.push_back(b);
    }

    // ------------------------------------------------
    // Output
    // ------------------------------------------------

    ofstream out("bubble_ratio.dat");

    out << "# coverage ratio total_volume largest_volume\n";

    // ------------------------------------------------
    // Increase bubble radius
    // ------------------------------------------------

    for(double radius = 2.0; radius <= 15.0; radius += 0.25)
    {
        vector<unsigned char> field(N, 0);

        // ------------------------------------------------
        // Create spherical bubbles
        // ------------------------------------------------

        for(auto b : bubbles)
        {
            int xmin = max(0, (int)(b.x - radius));
            int xmax = min(NX-1, (int)(b.x + radius));

            int ymin = max(0, (int)(b.y - radius));
            int ymax = min(NY-1, (int)(b.y + radius));

            int zmin = max(0, (int)(b.z - radius));
            int zmax = min(NZ-1, (int)(b.z + radius));

            for(int z = zmin; z <= zmax; z++)
            {
                for(int y = ymin; y <= ymax; y++)
                {
                    for(int x = xmin; x <= xmax; x++)
                    {
                        double dx = x - b.x;
                        double dy = y - b.y;
                        double dz = z - b.z;

                        double r2 =
                            dx*dx +
                            dy*dy +
                            dz*dz;

                        if(r2 <= radius*radius)
                        {
                            field[id(x,y,z)] = 1;
                        }
                    }
                }
            }
        }

        // ========================================================
        // Calculate total bubble volume
        // ========================================================

        int totalVolume = 0;

        for(int i = 0; i < N; i++)
        {
            if(field[i])
                totalVolume++;
        }

        // ========================================================
        // Find connected bubbles
        // ========================================================

        vector<unsigned char> visited(N, 0);

        vector<int> volumes;

        int dx[6] = {1,-1,0,0,0,0};
        int dy[6] = {0,0,1,-1,0,0};
        int dz[6] = {0,0,0,0,1,-1};

        for(int z = 0; z < NZ; z++)
        {
            for(int y = 0; y < NY; y++)
            {
                for(int x = 0; x < NX; x++)
                {
                    int start = id(x,y,z);

                    if(!field[start] || visited[start])
                        continue;

                    queue<int> q;

                    q.push(start);
                    visited[start] = 1;

                    int volume = 0;

                    while(!q.empty())
                    {
                        int current = q.front();
                        q.pop();

                        volume++;

                        int cz = current / (NX*NY);
                        int rem = current % (NX*NY);

                        int cy = rem / NX;
                        int cx = rem % NX;

                        for(int k = 0; k < 6; k++)
                        {
                            int xx = cx + dx[k];
                            int yy = cy + dy[k];
                            int zz = cz + dz[k];

                            if(xx < 0 || xx >= NX ||
                               yy < 0 || yy >= NY ||
                               zz < 0 || zz >= NZ)
                                continue;

                            int next = id(xx,yy,zz);

                            if(field[next] &&
                               !visited[next])
                            {
                                visited[next] = 1;
                                q.push(next);
                            }
                        }
                    }

                    volumes.push_back(volume);
                }
            }
        }

        // ========================================================
        // Largest bubble
        // ========================================================

        int largestVolume = 0;

        for(int v : volumes)
        {
            if(v > largestVolume)
                largestVolume = v;
        }

        // ========================================================
        // Coverage
        //
        // Projection on XY plane
        // ========================================================

        int coveredArea = 0;

        for(int y = 0; y < NY; y++)
        {
            for(int x = 0; x < NX; x++)
            {
                bool covered = false;

                for(int z = 0; z < NZ; z++)
                {
                    if(field[id(x,y,z)])
                    {
                        covered = true;
                        break;
                    }
                }

                if(covered)
                    coveredArea++;
            }
        }

        double coverage =
            (double)coveredArea /
            (NX * NY);

        // ========================================================
        // Largest / total
        // ========================================================

        double ratio = 0.0;

        if(totalVolume > 0)
        {
            ratio =
                (double)largestVolume /
                totalVolume;
        }

        // ------------------------------------------------
        // Save
        // ------------------------------------------------

        out
            << coverage << " "
            << ratio << " "
            << totalVolume << " "
            << largestVolume
            << "\n";

        cout
            << "radius = " << radius
            << " coverage = " << coverage
            << " ratio = " << ratio
            << endl;
    }

    out.close();

    cout << "\nFinished.\n";

    return 0;
}