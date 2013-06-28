Parallel Sorting by Regular Sampling (PSRS)

Sort n integers.

Algorithm details in:
http://zh.wikipedia.org/wiki/PSRS%E7%AE%97%E6%B3%95

void genRandom(); generate random integers in all processes
void doPSRS(int sampleCount=-1);
    void sample(int sampleCount); for each process of p processes, randomly chooses s=sampleCount elements
    void getDivide(); There're s*p elements chosen, sort and select s,2s,... ,(p-1)s-th elements as partitions and bcast to all processes.
    void remap(); call Alltoall to make each process hold its corresponding elements
