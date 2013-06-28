#include<mpi.h>
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<cstring>
#include<algorithm>
#include<cmath>

//#define NOW_DEBUG

typedef double DType;
#define MPI_DT_ MPI_DOUBLE
const int DEFAULT_N = 4;
const DType ERROR_VALUE = 1e-3;
const char DEFAULT_FILENAME[] = "input.txt";

DType** a;
DType* b;
DType* r;
DType* p;
DType* pn;
DType* ap;
DType* x;
DType rrLast;
DType rrNow;
DType ppNow;

int n;
int rank,size,len;

void init(){
    a = new DType*[len];
    for (int i=0;i<len;++i){
        a[i] = new DType[n];
    }
    b = new DType[n];
    r = new DType[len];
    p = new DType[n];
    pn = new DType[len];
    ap = new DType[len];
    x = new DType[len];
}
void finz(){
    for (int i=0;i<len;++i){
        delete[] a[i];
    }
    delete[] a;
    delete[] b;
    delete[] r;
    delete[] p;
    delete[] pn;
    delete[] ap;
    delete[] x;
}

DType vec_dot(DType* va,DType* vb,int ll){
    DType sub_sum = 0;
    for (int i=0;i<ll;++i){
        sub_sum += va[i]*vb[i];
    }
    DType all_sum;
    MPI_Allreduce(&sub_sum,&all_sum,1,MPI_DT_,MPI_SUM,MPI_COMM_WORLD);
/*
#ifdef NOW_DEBUG
    for (int i=0;i<size;++i){
        if (i==rank){
            printf("[PROC %2d] sub_sum = %lf, all_sum = %lf\n[PROC %2d] v = ",rank,sub_sum,all_sum,rank);
            for (int j=0;j<ll;++j){
                printf("(%lf, %lf) ",va[i],vb[i]);
            }
            printf("\n");
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
#endif
*/
    return all_sum;
}
void mtx_vec_mul(DType* res,DType** mm,DType* vv,int l1,int l2){
    for (int i=0;i<l1;++i){
        res[i] = 0;
        for (int j=0;j<l2;++j){
            res[i] += mm[i][j]*vv[j];
        }
    }
}
void vec_addass(DType* va,DType* vb,DType cc,int ll){
    for (int i=0;i<ll;++i){
        va[i] += cc*vb[i];
    }
}
void vec_mad(DType* vv,DType* va,DType* vb,DType cc,int ll){
    for (int i=0;i<ll;++i){
        vv[i] = va[i]+cc*vb[i];
    }
}

void load_(FILE* fin){
    for (int j=0;j<len;++j){
        for (int k=0;k<n;++k){
            fscanf(fin,"%lf",&a[j][k]);
        }
    }
}
void load(const char* fn=NULL){
    if (rank==size-1){
        FILE* fin = fopen((NULL==fn)?DEFAULT_FILENAME:fn,"r");
        for (int i=0;i<size;++i){
            load_(fin);
            if (i!=rank){
                for (int j=0;j<len;++j){
                    MPI_Send(a[j],n,MPI_DT_,i,i,MPI_COMM_WORLD);
                }
            }
        }
        for (int i=0;i<n;++i){
            fscanf(fin,"%lf",&b[i]);
        }
        fclose(fin);
    } else {
        MPI_Status status;
        for (int i=0;i<len;++i){
            MPI_Recv(a[i],n,MPI_DT_,size-1,rank,MPI_COMM_WORLD,&status);
        }
    }
    MPI_Bcast(b,n,MPI_DT_,size-1,MPI_COMM_WORLD);
    memset(x,0,sizeof(DType)*len);
    memcpy(r,b+len*rank,sizeof(DType)*len);
    memcpy(p,b,sizeof(DType)*n);
    rrNow = vec_dot(r,r,len);
}
#ifdef NOW_DEBUG
void dp_arr(const char* s,DType* arr,int l){
    for (int i=0;i<size;++i){
        if (rank==i){
            printf("[PROC %2d] %s = ",rank,s);
            for (int j=0;j<l;++j){
                printf("%lf ",arr[j]);
            }
            printf("\n");
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}
void dp_1(const char* s,DType res){
    for (int i=0;i<size;++i){
        if (rank==i){
            printf("[PROC %2d] %s = %lf\n",rank,s,res);
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}
#endif
void solve(){
    //int count = 0;
    while (fabs(rrNow)>=ERROR_VALUE){//&&++count<=2){
        mtx_vec_mul(pn,a,p,len,n);
#ifdef NOW_DEBUG
        dp_arr("pn",pn,len);
#endif
        ppNow = vec_dot(pn,p+rank*len,len);
#ifdef NOW_DEBUG
        dp_1("ppNow",ppNow);
#endif
        DType alpha = rrNow/ppNow;
        vec_addass(x,p+rank*len,alpha,len);
#ifdef NOW_DEBUG
        dp_arr("x",x,len);
#endif
        vec_addass(r,pn,-alpha,len);
#ifdef NOW_DEBUG
        dp_arr("r",r,len);
#endif
        rrLast = rrNow;
        rrNow = vec_dot(r,r,len);
#ifdef NOW_DEBUG
        dp_1("rrNow",rrNow);
#endif
        if (fabs(rrNow)<ERROR_VALUE){
            break;
        }
        DType beta = rrNow/rrLast;
#ifdef NOW_DEBUG
        dp_1("beta",beta);
#endif
        vec_mad(pn,r,p+rank*len,beta,len);
#ifdef NOW_DEBUG
        dp_arr("p",pn,len);
#endif
        MPI_Allgather(pn,len,MPI_DT_,p,len,MPI_DT_,MPI_COMM_WORLD);
    }
    MPI_Gather(x,len,MPI_DT_,p,len,MPI_DT_,0,MPI_COMM_WORLD);
    if (!rank){
        printf("x = ");
        for (int i=0;i<n;++i){
            printf("%lf ",p[i]);
        }
        printf("\n");
        fflush(stdout);
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
    load();
    solve();
    finz();
    MPI_Finalize();
    return 0;
}
