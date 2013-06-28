Closet Pair

Given coordinates of n points in 2-D space, output the distance of the nearest pair of points.

Divide n points in to p parts, s.t. maximum X position of i-th part <= minimum X Position of (i+1)th part
Perform serial O(nlgn) algorithm to find nearest pair.
As serial version, O(lgp) times recursively merge different parts to find nearest pair.

void genRandom(); generate n points with random coordinates
void solve(); main implementation, as described above.
    void sample(int sampleCount); for each process of p processes, randomly chooses s=sampleCount x position.
    void getDivide(); There're s*p elements chosen, sort and select s,2s,... ,(p-1)s-th elements as boundaries and bcast to all processes.
    void remap(); call Alltoall to make ith process hold those between selected (i-1)th and ith boundary element.
    void dcc(); serial version
        void go(); divide and conquer sub routine
    void reduceAnswer(int level); main routine to recursively merge data in all processes
        void selectData(int l); select candidate elements near boundary for each process
        void merge(); merge between two processes
