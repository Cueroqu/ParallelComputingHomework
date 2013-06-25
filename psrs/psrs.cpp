#include<mpi.h>
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<cstring>
#include<algorithm>

//#define NOW_DEBUG

typedef unsigned DType;
#define MPI_DT_ MPI_INT
const int DEFAULT_N = 40000000;

DType* data;
DType* tmpData;
DType* divide;
int* cnt;
int* num;
int* nwi;
int* sendcounts;
int* sdispls;
int* recvcounts;
int* rdispls;
int* spf;
int n;
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
    if (size>1){
        divide = new DType[size-1];
    }
    cnt = new int[size*size];
    num = new int[size];
    nwi = new int[size+1];
    sendcounts = new int[size];
    sdispls = new int[size];
    recvcounts = new int[size];
    rdispls = new int[size];
    spf = new int[size+1];
}
void finz(){
    delete[] data;
    if (size>1){
        delete[] divide;
    }
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
    for (int it=0;it<size;++it){
        if (rank==it){
            printf("[PROC %2d] oriEle = ",rank);
            for (int i=0;i<len;++i){
                printf("%u ",data[i]);
            }
            printf("\n");
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
}

void sample(int sampleCount){
    for (int i=0;i<sampleCount;++i){
        int tmp = rand()%(sampleCount-i);
        std::swap(data[i],data[tmp+i]);
    }
    if (!rank){
        DType* allChosen = new DType[sampleCount*size];
        MPI_Gather(data,sampleCount,MPI_DT_,allChosen,sampleCount,MPI_DT_,0,MPI_COMM_WORLD);
        std::sort(allChosen,allChosen+sampleCount*size);
        for (int i=0;i<size-1;++i){
            divide[i] = allChosen[(i+1)*sampleCount];
        }
        MPI_Bcast(divide,size-1,MPI_DT_,0,MPI_COMM_WORLD);
        delete[] allChosen;
    } else {
        MPI_Gather(data,sampleCount,MPI_DT_,NULL,0,MPI_DT_,0,MPI_COMM_WORLD);
        MPI_Bcast(divide,size-1,MPI_DT_,0,MPI_COMM_WORLD);
    }
}
void getDivide(){
    int st = 0;
    while (st<len&&data[st]<divide[0]){
        ++st;
    }
    num[0] = st;
    for (int i=1;i<size-1;++i){
        num[i] = st;
        while (st<len&&data[st]<divide[i]){
            ++st;
        }
        num[i] = st-num[i];
    }
    num[size-1] = len-st;
}
void remap(){
    MPI_Allgather(num,size,MPI_INT,cnt,size,MPI_INT,MPI_COMM_WORLD);
    nwi[0] = 0;
    for (int i=0;i<size;++i){
        nwi[i+1] = 0;
        for (int j=0;j<size;++j){
            nwi[i+1] += cnt[j*size+i];
        }
        nwi[0] = std::max(nwi[0],nwi[i+1]);
    }
    if (!rank){
        tmpData = new DType[nwi[0]];
    } else {
        tmpData = new DType[nwi[rank+1]];
    }
    nwi[0] = 0;
    for (int i=1;i<=size;++i){
        nwi[i] += nwi[i-1];
    }
    for (int i=0;i<size;++i){
        recvcounts[i] = cnt[i*size+rank];
    }
    sdispls[0] = 0;
    rdispls[0] = 0;
    for (int i=1;i<size;++i){
        sdispls[i] = sdispls[i-1]+num[i-1];
        rdispls[i] = rdispls[i-1]+recvcounts[i-1];
    }
/*
#ifdef NOW_DEBUG
    fflush(stdout); MPI_Barrier(MPI_COMM_WORLD);
    if (!rank){
        printf("[PROC %2d] nwi = ",rank);
        for (int j=0;j<=size;++j){
            printf("%d ",nwi[j]);
        }
        printf("\n");
        printf("[PROC %2d] cnt = ",rank);
        for (int j=0;j<size*size;++j){
            printf("%d ",cnt[j]);
        }
        printf("\n");
    }
    fflush(stdout); MPI_Barrier(MPI_COMM_WORLD);
    for (int i=0;i<size;++i){
        if (i==rank){
            printf("[PROC %2d] 4p = ",rank);
            for (int j=0;j<size;++j){
                printf("(%d,%d,%d,%d) ",num[j],sdispls[j],recvcounts[j],rdispls[j]);
            }
            printf("\n");
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
*/
    for (int i=0;i<size;++i){
        MPI_Gatherv(data+sdispls[i],num[i],MPI_DT_,tmpData,recvcounts,rdispls,MPI_DT_,i,MPI_COMM_WORLD);
    }
}
void doPSRS(int sampleCount=-1){
    if (size==1){
        std::sort(data,data+len);
        return;
    }
    if (sampleCount<=0){
        sampleCount = size;
    }
    
    sample(sampleCount);
    std::sort(data,data+len);
    getDivide();
    remap();
    std::sort(tmpData,tmpData+(nwi[rank+1]-nwi[rank]));
}

void outLine(int id,DType* td,int l){
#ifdef NOW_DEBUG
    fprintf(stdout,"[PROC %2d]",id);
#endif
    for (int j=0;j<l;++j){
        fprintf(stdout,"%u ",td[j]);
    }
    fprintf(stdout,"\n");
    fflush(stdout);
}
void output(){
    if (size<=1){
        outLine(0,data,len);
        return;
    }
    if (!rank){
        outLine(0,tmpData,nwi[rank+1]-nwi[rank]);
        MPI_Status status;
        for (int i=1;i<size;++i){
            MPI_Recv(tmpData,nwi[i+1]-nwi[i],MPI_DT_,i,i,MPI_COMM_WORLD,&status);
            outLine(i,tmpData,nwi[i+1]-nwi[i]);
        }
    } else {
        MPI_Send(tmpData,nwi[rank+1]-nwi[rank],MPI_DT_,0,rank,MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef NOW_DEBUG
    printf("[PROC %2d]HERE!!!!\n",rank);
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

int main(int argc,char* argv[]){
    n = DEFAULT_N;
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    if (n%size){
        n += size-n%size;
    }
    len = n/size;
    init();
    genRandom();
    doPSRS();
    output();
    finz();
    MPI_Finalize();
    return 0;
}
