#include<mpi.h>
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include<cstring>
#include<algorithm>
#include<cmath>

typedef double DType;
#define MPI_DT_ MPI_DOUBLE
const int DEFAULT_N = 32;
const DType oo = 2e6;
const DType eps = 1e-6;

DType* sameData;
DType* data;
DType* data0;
int* tmpData0;
int n;
int rank,size,len;
DType bestAnswer;
DType tmpBestAnswer;
int curLen;
int otLen;

DType* sampleData;
DType* divide;
int* tmpData1;
int* cnt;
int* num;
int* nwi;
int* sendcounts;
int* sdispls;
int* recvcounts;
int* rdispls;
int* odrX;

void init(){
    sameData = new DType[len+len];
}
void finz(){
    delete[] sameData;
}

DType randamValue(){
    return double(rand())/double(RAND_MAX)*oo;
}
void genRandom(){
    srand(time(NULL)+rank);
    for (int i=0;i<len;++i){
        sameData[i+i] = randamValue();
        sameData[i+i+1] = randamValue();
#ifdef NOW_DEBUG
        sameData[i+i] = rand()%32767;
        sameData[i+i+1] = rand()%32767;
#endif
    }
}

void sample(int sampleCount){
    sampleData = new DType[sampleCount];
    for (int i=0;i<sampleCount;++i){
        int tmp = rand()%(sampleCount-i)+i;
        std::swap(sameData[i+i],sameData[tmp+tmp]);
        std::swap(sameData[i+i+1],sameData[tmp+tmp+1]);
        sampleData[i] = sameData[i+i];
    }
    if (!rank){
        DType* allChosen = new DType[sampleCount*size];
        MPI_Gather(sampleData,sampleCount,MPI_DT_,allChosen,sampleCount,MPI_DT_,0,MPI_COMM_WORLD);
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
    delete[] sampleData;
}
void getDivide(){
    int st = 0;
    while (st<len&&sameData[st+st]<divide[0]){
        ++st;
    }
    num[0] = st;
    for (int i=1;i<size-1;++i){
        num[i] = st;
        while (st<len&&sameData[st+st]<divide[i]){
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
    data = new DType[nwi[rank+1]*2];
    nwi[0] = 0;
    for (int i=1;i<=size;++i){
        nwi[i] += nwi[i-1];
    }
    for (int i=0;i<size;++i){
        recvcounts[i] = cnt[i*size+rank]*2;
    }
    sdispls[0] = 0;
    rdispls[0] = 0;
    for (int i=1;i<size;++i){
        sdispls[i] = sdispls[i-1]+num[i-1];
        rdispls[i] = rdispls[i-1]+recvcounts[i-1];
    }
    for (int i=0;i<size;++i){
        MPI_Gatherv(sameData+sdispls[i]*2,num[i]*2,MPI_DT_,data,recvcounts,rdispls,MPI_DT_,i,MPI_COMM_WORLD);
    }
}

bool cmp(const int& lhs,const int& rhs){
    return data[lhs+lhs]<data[rhs+rhs];
}
inline DType sqr(DType x){
    return x*x;
}
inline DType getDis(DType x1,DType y1,DType x2,DType y2){
    return sqrt(sqr(x1-x2)+sqr(y1-y2));
}
DType go(int left,int right,int* aim,int* tmpa){
    if (right==left){
        aim[left] = odrX[left];
        return oo;
    }
    if (1==right-left){
        if (data[odrX[left]+odrX[left]+1]<data[odrX[right]+odrX[right]+1]){
            aim[left] = odrX[left];
            aim[right]= odrX[right];
        } else {
            aim[left] = odrX[right];
            aim[right]= odrX[left];
        }
        return getDis(data[odrX[left]+odrX[left]],data[odrX[left]+odrX[left]+1],data[odrX[right]+odrX[right]],data[odrX[right]+odrX[right]+1]);
    }
    int mid = (left+right)/2;
    int cur = mid+1;
    DType leftValue = go(left,mid,tmpa,aim);
    DType rightValue = go(cur,right,tmpa,aim);
    if (rightValue<leftValue){
        leftValue = rightValue;
    }
    for (int i=left;i<=mid;++i){
        while ((cur<=right)&&(data[tmpa[cur]+tmpa[cur]+1]<data[tmpa[i]+tmpa[i]+1]-leftValue-eps)){
            ++cur;
        }
        for (int j=cur;j<=right;++j){
            if (data[tmpa[j]+tmpa[j]+1]>data[tmpa[i]+tmpa[i]+1]+leftValue+eps){
                break;
            }
            leftValue = std::min(leftValue,getDis(data[tmpa[i]+tmpa[i]],data[tmpa[i]+tmpa[i]+1],data[tmpa[j]+tmpa[j]],data[tmpa[j]+tmpa[j]+1]));
        }
    }
    cur = mid+1;
    int j = left;
    for (int i=left;i<=right;++i){
        if (j<=mid){
            if (cur<=right){
                if (data[tmpa[j]+tmpa[j]+1]<data[tmpa[cur]+tmpa[cur]+1]){
                    aim[i] = tmpa[j++];
                } else {
                    aim[i] = tmpa[cur++];
                }
            } else {
                aim[i] = tmpa[j++];
            }
        } else {
            aim[i] = tmpa[cur++];
        }
    }
    return leftValue;
}
void dcc(){
    odrX = new int[len];
    tmpData1 = new int[len];
    
    for (int i=0;i<len;++i){
        odrX[i] = i;
    }
    std::sort(odrX,odrX+len,cmp);
    bestAnswer = go(0,len-1,tmpData0,tmpData1);
    delete[] odrX;
    delete[] tmpData1;
    for (int j=0;j<len;++j){
        int i = tmpData0[j];
        data0[j+j] = data[i+i];
        data0[j+j+1] = data[i+i+1];
    }
}

inline int getNodeType(int l,int r,int& ot){
    int ss = 1<<(l-1);
    if (r&(ss-1)){
        return -1;
    }
    ot = r^ss;
    if (ot>=size){
        return -1;
    }
    return r&ss;
}
void selectData(int l){
    int lb = rank;
    int rb = rank+(1<<l)-1;
    bool useL = lb-1>=0;
    bool useR = rb<size-1;
    int j = 0;
    for (int i=0;i<len;++i){
        if ((useL&&data0[i+i]<divide[lb-1]+bestAnswer+eps)||(useR&&data0[i+i]>divide[rb]-bestAnswer-eps)){
            data0[j+j] = data0[i+i];
            data0[j+j+1] = data0[i+i+1];
            ++j;
        }
    }
    curLen = j;
}
void merge(){
    DType* data1 = data+curLen+curLen;
    int j = 0;
    for (int i=0;i<curLen;++i){
        while (j<otLen&&data1[j+j+1]<data0[i+i+1]-bestAnswer-eps){
            ++j;
        }
        for (int k=j;k<otLen;++k){
            if (data1[k+k+1]>data0[i+i+1]+bestAnswer+eps){
                break;
            }
            bestAnswer = std::min(bestAnswer,getDis(data1[k+k],data1[k+k+1],data0[i+i],data0[i+i+1]));
        }
    }
    j = 0;
    int k = 0;
    for (int i=0;i<curLen+otLen;++i){
        if (j<otLen){
            if (k<curLen){
                if (data0[j+j+1]<data1[k+k+1]){
                    data[i+i] = data0[j+j];
                    data[i+i+1] = data0[j+j+1];
                    ++j;
                } else {
                    data[i+i] = data1[k+k];
                    data[i+i+1] = data1[k+k+1];
                    ++k;
                }
            } else {
                data[i+i] = data0[j+j];
                data[i+i+1] = data0[j+j+1];
                ++j;
            }
        } else {
            data[i+i] = data1[k+k];
            data[i+i+1] = data1[k+k+1];
            ++k;
        }
    }
    std::swap(data,data0);
}
void reduceAnswer(int level){
    if (level>0){
        reduceAnswer(level-1);
    }
    tmpBestAnswer = bestAnswer;
    MPI_Allreduce(&tmpBestAnswer,&bestAnswer,1,MPI_DT_,MPI_MIN,MPI_COMM_WORLD);
    int ot;
    int nodeType = getNodeType(level,rank,ot);
    if (1==nodeType){//sender
        selectData(level);
        MPI_Send(&curLen,1,MPI_INT,ot,level,MPI_COMM_WORLD);
        MPI_Send(data0,curLen*2,MPI_DT_,ot,level,MPI_COMM_WORLD);
    } else if (0==nodeType){//receiver
        selectData(level);
        MPI_Status status;
        MPI_Recv(&otLen,1,MPI_INT,ot,level,MPI_COMM_WORLD,&status);
        delete[] data;
        data = new DType[otLen+otLen+curLen+curLen];
        MPI_Recv(data+curLen+curLen,otLen*2,MPI_DT_,ot,level,MPI_COMM_WORLD,&status);
        merge();
    } else {//others
    }
    
}
void solve(){
    if (size==1){
        data = sameData;
        dcc();
        return;
    }
    
    {
        divide = new DType[size-1];
        cnt = new int[size*size];
        num = new int[size];
        nwi = new int[size+1];
        sendcounts = new int[size];
        sdispls = new int[size];
        recvcounts = new int[size];
        rdispls = new int[size];
    }    
    sample(size);
    getDivide();
    remap();
    len = nwi[rank+1]-nwi[rank];
    {
        delete[] cnt;
        delete[] num;
        delete[] nwi;
        delete[] sendcounts;
        delete[] sdispls;
        delete[] recvcounts;
        delete[] rdispls;
    }
    
    tmpData0 = new int[len];
    data0 = new DType[len+len];
    
    dcc();
    int maxLevel = 1;
    while ((1<<maxLevel)<size){
        ++maxLevel;
    }
    reduceAnswer(maxLevel-1);
    tmpBestAnswer = bestAnswer;
    MPI_Allreduce(&tmpBestAnswer,&bestAnswer,1,MPI_DT_,MPI_MIN,MPI_COMM_WORLD);
    
    delete[] tmpData0;
    delete[] divide;
    delete[] data;
    delete[] data0;
    
    if (!rank){
        printf("%.6lf\n",bestAnswer);
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
    genRandom();
    solve();
    finz();
    MPI_Finalize();
    return 0;
}
