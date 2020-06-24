//为了计算效率，该文件将所有为1的量全部省去，所有的算式均为化简后的结果
#include <iostream>
#include <fstream>
#include <iomanip>
#include <gsl/gsl_rng.h>
#include <cmath>

//时间步，这里如果设为0.1，变化太大，设为0.001，变化太小
#define t 0.01

using namespace std;

//x为坐标，v为速度，f为受力，f1为第n+1个时间步的受力
double x[1000][3], v[1000][3], f[1000][3], f1[1000][3];

int main()
{
    int m = 0, n, a, b;
    double i, j, k, sum = 0.0, s, sx, sy, sz, c;
    ofstream file;//声明输出对象file
    
    gsl_rng* r;
    gsl_rng_default_seed = 3;
    r = gsl_rng_alloc(gsl_rng_mt19937);

    //初始化坐标
    for (i = 1.5; i < 30.0; i = i + 3.0)
    {
        for (j = 1.5; j < 30.0; j = j + 3.0)
        {
            for (k = 1.5; k < 30.0; k = k + 3.0)
            {
                x[m][0] = i;
                x[m][1] = j;
                x[m][2] = k;
                m++;
            }
        }
    }

    //初始化速度并标定
    for (m = 0; m < 1000; m++)
    {
        for (n = 0; n < 3; n++)
        {
            v[m][n] = gsl_rng_uniform(r);
            sum += pow(v[m][n],2);
        }
    }

    for (m = 0; m < 1000; m++)
    {
        for (n = 0; n < 3; n++)
        {
            v[m][n] = sqrt(pow(v[m][n], 2) * 2997.0 / sum);
        }
    }

    //初始化受力
    for (m = 0; m < 1000; m++)
    {
        for (n = 0; n < 3; n++)
        {
            f[m][n] = 0.0;
        }
    }
    
    //输出初始坐标文件
    file.open("data1.gro", ios::out);
    file << "O1" << endl;
    file << setiosflags(ios::right) << setw(5) << 1000 << endl;
    for (m = 0; m < 1000; m++)
    {
        file << setiosflags(ios::right) << setw(5) << m + 1;
        file << "SOL  ";
        file << setiosflags(ios::right) << setw(5) << "OW";
        file << setiosflags(ios::right) << setw(5) << m + 1;
        for (n = 0; n < 3; n++)
        {
            file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(8) << setprecision(3) << x[m][n];
        }
        file << endl;
    }
    file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(10) << setprecision(5) << 30.0;
    file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(10) << setprecision(5) << 30.0;
    file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(10) << setprecision(5) << 30.0 << endl;
    file.close();

    //开始模拟
    for (m = 0; m < 100; m++)
    {
        //速度verlet公式，n+1时刻的位置并更新
        for (a = 0; a < 1000; a++)
        {
            for (b = 0; b < 3; b++)
            {
                x[a][b] = x[a][b] + v[a][b] * t + 0.5 * f[a][b] * pow(t, 2);
            }
        }

        //周期边界
        for (a = 0; a < 1000; a++)
        {
            for (b = 0; b < 3; b++)
            {
                if (x[a][b] < 0.0)  x[a][b] = x[a][b] + 30.0;
                if (x[a][b] > 30.0) x[a][b] = x[a][b] - 30.0;
            }
        }
        
        //n+1时刻的受力
        //老方法，计算量大
        /*
        for (a = 0; a < 1000; a++)
        {
            for (b = 0; b < 1000; b++)
            {
                //最小镜像约定
                sx = x[b][0] - x[a][0];
                sy = x[b][1] - x[a][1];
                sz = x[b][2] - x[a][2];
                if (sx < -15.0) sx = sx + 30.0;
                if (sx > 15.0)  sx = sx - 30.0;
                if (sy < -15.0) sy = sy + 30.0;
                if (sy > 15.0)  sy = sy - 30.0;
                if (sz < -15.0) sz = sz + 30.0;
                if (sz > 15.0)  sz = sz - 30.0;

                //考虑截断半径
                s = sqrt(pow(sx, 2) + pow(sy, 2) + pow(sz, 2));
                if (s < 4.0 && s > 0.0)
                {
                    f1[a][0] += (48 * pow(1.0 / s, 13) - 24 * pow(1.0 / s, 7)) * sx/s;
                    f1[a][1] += (48 * pow(1.0 / s, 13) - 24 * pow(1.0 / s, 7)) * sy/s;
                    f1[a][2] += (48 * pow(1.0 / s, 13) - 24 * pow(1.0 / s, 7)) * sz/s;
                }
            }
        }
        */
        //新，应用牛顿第三定律，重复受力只算一次，来源于the art of molecular dynamics simulation
        for (a = 0; a < 999; a++)
        {
            for (b = a + 1; b < 1000; b++)
            {
                //最小镜像约定
                sx = x[b][0] - x[a][0];
                sy = x[b][1] - x[a][1];
                sz = x[b][2] - x[a][2];
                if (sx < -15.0) sx = sx + 30.0;
                if (sx > 15.0)  sx = sx - 30.0;
                if (sy < -15.0) sy = sy + 30.0;
                if (sy > 15.0)  sy = sy - 30.0;
                if (sz < -15.0) sz = sz + 30.0;
                if (sz > 15.0)  sz = sz - 30.0;

                //考虑截断半径
                s = sqrt(pow(sx, 2) + pow(sy, 2) + pow(sz, 2));
                if (s < 4.0 && s > 0.0)
                {
                    f1[a][0] += (48 * pow(1.0 / s, 13) - 24 * pow(1.0 / s, 7)) * sx / s;
                    f1[a][1] += (48 * pow(1.0 / s, 13) - 24 * pow(1.0 / s, 7)) * sy / s;
                    f1[a][2] += (48 * pow(1.0 / s, 13) - 24 * pow(1.0 / s, 7)) * sz / s;
                    f1[b][0] -= (48 * pow(1.0 / s, 13) - 24 * pow(1.0 / s, 7)) * sx / s;
                    f1[b][1] -= (48 * pow(1.0 / s, 13) - 24 * pow(1.0 / s, 7)) * sy / s;
                    f1[b][2] -= (48 * pow(1.0 / s, 13) - 24 * pow(1.0 / s, 7)) * sz / s;
                }
            }
        }

        //n+1时刻的速度并更新
        for (a = 0; a < 1000; a++)
        {
            for (b = 0; b < 3; b++)
            {
                v[a][b] = v[a][b] + 0.5 * (f1[a][b] + f[a][b]) * t;
            }
        }

        //速度耦合（berendsen热浴）
        sum = 0.0;
        for (a = 0; a < 1000; a++)
        {
            for (b = 0; b < 3; b++)
            {
                sum += pow(v[a][b], 2);
            }
        }
        c = sqrt(1 + t * (2997.0 / sum - 1));
        for (a = 0; a < 1000; a++)
        {
            for (b = 0; b < 3; b++)
            {
                v[m][n] = c * v[m][n];
            }
        }

        //更新受力
        for (a = 0; a < 1000; a++)
        {
            for (b = 0; b < 3; b++)
            {
                f[a][b] = f1[a][b];
            }
        }

        //将f1数组清零
        for (a = 0; a < 1000; a++)
        {
            for (b = 0; b < 3; b++)
            {
                f1[a][b] = 0;
            }
        }
    }

    //输出最终坐标文件
    file.open("data2.gro", ios::out);
    file << "O2" << endl;
    file << setiosflags(ios::right) << setw(5) << 1000 << endl;
    for (m = 0; m < 1000; m++)
    {
        file << setiosflags(ios::right) << setw(5) << m + 1;
        file << "SOL  ";
        file << setiosflags(ios::right) << setw(5) << "OW";
        file << setiosflags(ios::right) << setw(5) << m + 1;
        for (n = 0; n < 3; n++)
        {
            file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(8) << setprecision(3) << x[m][n];
        }
        file << endl;
    }
    file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(10) << setprecision(5) << 30.0;
    file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(10) << setprecision(5) << 30.0;
    file << setiosflags(ios::right) << setiosflags(ios::fixed) << setw(10) << setprecision(5) << 30.0 << endl;
    file.close();
    
    return 0;
}

