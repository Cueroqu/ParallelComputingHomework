#include<mpi.h>
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<cstring>
#include<algorithm>

//#define NOW_DEBUG

typedef unsigned DType;
#define MPI_DT_ MPI_INT
const int DEFAULT_N = 40;
const int DEFAULT_MASK = 10;
const int DEFAULT_TIMES = 5;

DType* data;
DType* tmpData;
int* cnt;
int* num;
int* nwi;
int* sendcounts;
int* sdispls;
int* recvcounts;
int* rdispls;
int* spf;
int mask,times,n;
int rank,size,len;

#ifdef NOW_DEBUG
void outSelfLine(DType* data){
    fprintf(stdout,"[PROC %2d]",rank);
    for (int i=0;i<len;++i){
        fprintf(stdout,"%u ",data[i]);
    }
    fprintf(stdout,"\n");
    fflush(stdout);
}
#endif

void init(){
    data = new DType[len];
    tmpData = new DType[len];
    cnt = new int[size*mask];
    num = new int[mask];
    nwi = new int[mask];
    sendcounts = new int[size];
    sdispls = new int[size];
    recvcounts = new int[size];
    rdispls = new int[size];
    spf = new int[size+1];
}
void finz(){
    delete[] data;
    delete[] tmpData;
    delete[] cnt;
    delete[] num;
    delete[] nwi;
    delete[] sendcounts;
    delete[] sdispls;
    delete[] recvcounts;
    delete[] rdispls;
    delete[] spf;
}

void genRandom(){
    srand(time(NULL)+rank);
#ifndef NOW_DEBUG
    int tot = sizeof(DType)*len;
    unsigned char* c = (unsigned char*)data;
    for (int i=0;i<tot;++i){
        c[i] = rand()&255;
    }
#endif
#ifdef NOW_DEBUG
    for (int i=0;i<len;++i){data[i] = rand()%32767;}
#endif
}

void bucketSort(const DType& bv){
    memset(num,0,sizeof(int)*mask);
    for (int i=0;i<len;++i){
        ++num[(data[i]/bv)%mask];
    }
    nwi[0] = 0;
    for (int i=1;i<mask;++i){
        nwi[i] = nwi[i-1]+num[i-1];
    }
    for (int i=0;i<len;++i){
        int tmp = (data[i]/bv)%mask;
        tmpData[nwi[tmp]++] = data[i];
    }
}

void remap(){
    MPI_Allgather(num,mask,MPI_INT,cnt,mask,MPI_INT,MPI_COMM_WORLD);
#ifdef NOW_DEBUG
    if (!rank){
        for (int i=0;i<size;++i){
            for (int j=0;j<mask;++j){
                printf("%d ",cnt[i*mask+j]);
            }
            printf("\n");
        }
    }
    fflush(stdout);
    memset(data,0,sizeof(DType)*len);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    int arrStart = 0;
    int arrEnd = 0;
    int allStart = 0;
    for (int i=0;i<mask;++i){
        spf[0] = allStart;
        for (int j=0;j<size;++j){
            spf[j+1] = spf[j]+cnt[j*mask+i];
        }
        for (int j=0;j<size;++j){
            if (spf[rank]>j*len+len-1 || spf[rank+1]-1<j*len){
                sdispls[j] = arrStart;
                sendcounts[j] = 0;
            } else {
                sdispls[j] = arrStart;
                sendcounts[j] = std::min(spf[rank+1],j*len+len)-std::max(j*len,spf[rank]);
                arrStart += sendcounts[j];
            }
            if (spf[j]>rank*len+len-1 || spf[j+1]-1<rank*len){
                rdispls[j] = arrEnd;
                recvcounts[j] = 0;
            } else {
                rdispls[j] = arrEnd;
                recvcounts[j] = std::min(spf[j+1],rank*len+len)-std::max(rank*len,spf[j]);
                arrEnd += recvcounts[j];
            }
        }

#ifdef NOW_DEBUG
        //if (i==2){
            printf("[PROC %2d i %d sendcounts]",rank,i);
            for (int j=0;j<size;++j){
                printf("%d ",sendcounts[j]);
            }
            printf("\n");
            fflush(stdout);
            MPI_Barrier(MPI_COMM_WORLD);

            printf("[PROC %2d i %d recvcounts]",rank,i);
            for (int j=0;j<size;++j){
                printf("%d ",recvcounts[j]);
            }
            printf("\n");
            fflush(stdout);
            MPI_Barrier(MPI_COMM_WORLD);

            printf("[PROC %2d i %d spf]",rank,i);
            for (int j=0;j<=size;++j){
                printf("%d ",spf[j]);
            }
            printf("\n");
            fflush(stdout);
            MPI_Barrier(MPI_COMM_WORLD);
        //}
#endif

        MPI_Alltoallv(tmpData,sendcounts,sdispls,MPI_DT_,data,recvcounts,rdispls,MPI_DT_,MPI_COMM_WORLD);
/*
#ifdef NOW_DEBUG
        outSelfLine(data);
        MPI_Barrier(MPI_COMM_WORLD);
        if (!rank){
            fprintf(stdout,"[mask = %2d]========\n",i);
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
#endif
*/
        allStart = spf[size];
    }
}

void doRadixSort(){
    DType bv = 1;
    for (int i=0;i<times;++i){
        bucketSort(bv);
        bv *= mask;
        remap();
#ifdef NOW_DEBUG
        outSelfLine(data);
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
        if (!rank){
            fprintf(stdout,"[iter = %2d]========\n",i);
        }
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
#endif

    }
}

void outLine(int id,DType* td){
#ifdef NOW_DEBUG
    fprintf(stdout,"[PROC %2d]",id);
#endif
    for (int j=0;j<len;++j){
        fprintf(stdout,"%u ",td[j]);
    }
    fprintf(stdout,"\n");
}
void output(){
    if (!rank){
        outLine(0,data);
        MPI_Status status;
        for (int i=1;i<size;++i){
            MPI_Recv(tmpData,len,MPI_DT_,i,i,MPI_COMM_WORLD,&status);
            outLine(i,tmpData);
        }
    } else {
        MPI_Send(data,len,MPI_DT_,0,rank,MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef NOW_DEBUG
    printf("[PROC %2d]HERE!!!!\n",rank);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

int main(int argc,char* argv[]){
    n = DEFAULT_N;
    mask = DEFAULT_MASK;
    times = DEFAULT_TIMES;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if (n%size){
        n += size-n%size;
    }
    len = n/size;
    init();
    genRandom();
    doRadixSort();
    output();
    finz();
    MPI_Finalize();
    return 0;
}
