#include <omp.h>
#include <time.h>
#include <unistd.h>
#include <stdio.h>

int main(int argc, char **argv)
{
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int counter = 0 ;
        #pragma omp for schedule(dynamic,1)
        for(int i = 0 ; i<30;i++)
        {
            if(tid == 0)
            {
                usleep(10000);
            }
            else
            {
                usleep(30000);
            }
            counter++;
        }
        printf("tid %d, counter %d\n",tid,counter);
    }
    return 0;
}