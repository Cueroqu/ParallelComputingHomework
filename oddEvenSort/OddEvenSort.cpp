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
DType* recvData;
int mask,times,n;
int rank,size,len;

void init(){
    data = new DType[len];
    tmpData = new DType[len];
    recvData = new DType[len];
}
void finz(){
    delete[] data;
    delete[] tmpData;
    delete[] recvData;
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

void outLine(int id,DType* td){
#ifdef NOW_DEBUG
    fprintf(stdout,"[PROC %2d]",id);
#endif
    for (int j=0;j<len;++j){
        fprintf(stdout,"%u ",td[j]);
    }
    fprintf(stdout,"\n");
    fflush(stdout);
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
}

bool mergeMax(DType* a,DType* b,DType* c){
    int ia = len-1;
    int ib = len-1;
    for (int i=len-1;i>=0;--i){
        if (a[ia]>=b[ib]){
            c[i] = a[ia--];
        } else {
            c[i] = b[ib--];
        }
    }
    return ib==len-1;
}
bool mergeMin(DType* a,DType* b,DType* c){
    int ia = 0;
    int ib = 0;
    for (int i=0;i<len;++i){
        if (a[ia]<=b[ib]){
            c[i] = a[ia++];
        } else {
            c[i] = b[ib++];
        }
    }
    return ib==0;
}
void doOddEvenSort(){
    std::sort(data,data+len);
    int left = (rank)?rank-1:MPI_PROC_NULL;
    int right = (rank==size-1)?MPI_PROC_NULL:rank+1;
    MPI_Status status;
#ifdef NOW_DEBUG
    output();
#endif
    while (true){
    //for (int itt=0;itt<size;++itt){
        int sorted = 1;
        if (rank%2){//1->0
            if (left!=MPI_PROC_NULL){
                MPI_Sendrecv(data,len,MPI_DT_,left,rank,recvData,len,MPI_DT_,left,left,MPI_COMM_WORLD,&status);
                sorted &= mergeMax(data,recvData,tmpData);
            } else {
                memcpy(tmpData,data,sizeof(DType)*len);
            }
        } else {//0->1
            if (right!=MPI_PROC_NULL){
                MPI_Sendrecv(data,len,MPI_DT_,right,rank,recvData,len,MPI_DT_,right,right,MPI_COMM_WORLD,&status);
                sorted &= mergeMin(data,recvData,tmpData);
            } else {
                memcpy(tmpData,data,sizeof(DType)*len);
            }
        }
#ifdef NOW_DEBUG
        MPI_Barrier(MPI_COMM_WORLD);
        for (int rr=0;rr<size;++rr){
            if (rank==rr){
                printf("[PROC %2d]",rank);
                for (int i=0;i<len;++i){
                    printf("(%d %d %d)",data[i],recvData[i],tmpData[i]);
                }
                printf("\n");
            }
            fflush(stdout);
            MPI_Barrier(MPI_COMM_WORLD);
        }
#endif
        if (rank%2){//1->2
            if (right!=MPI_PROC_NULL){
                MPI_Sendrecv(tmpData,len,MPI_DT_,right,rank,recvData,len,MPI_DT_,right,right,MPI_COMM_WORLD,&status);
                sorted &= mergeMin(tmpData,recvData,data);
            } else {
                memcpy(data,tmpData,sizeof(DType)*len);
            }
        } else {//2->1
            if (left!=MPI_PROC_NULL){
                MPI_Sendrecv(tmpData,len,MPI_DT_,left,rank,recvData,len,MPI_DT_,left,left,MPI_COMM_WORLD,&status);
                sorted &= mergeMax(tmpData,recvData,data);
            } else {
                memcpy(data,tmpData,sizeof(DType)*len);
            }
        }
#ifdef NOW_DEBUG
        MPI_Barrier(MPI_COMM_WORLD);
        for (int rr=0;rr<size;++rr){
            if (rank==rr){
                printf("[PROC %2d]",rank);
                for (int i=0;i<len;++i){
                    printf("(%d %d %d)",data[i],recvData[i],tmpData[i]);
                }
                printf("\n");
            }
            fflush(stdout);
            MPI_Barrier(MPI_COMM_WORLD);
        }
#endif
        int allsorted;
        MPI_Allreduce(&sorted,&allsorted,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
        if (allsorted){
#ifdef NOW_DEBUG
            printf("[PROC %2d]sorted = %d, allsorted = %d\n",rank,sorted,allsorted);
            fflush(stdout);
#endif
            break;
        }
    }
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
    doOddEvenSort();
    output();
    finz();
    MPI_Finalize();
    return 0;
}
