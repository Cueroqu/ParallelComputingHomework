Odd-Even Sort

Sort n integers.

Algorithm details in:
http://en.wikipedia.org/wiki/Odd-even_sort

Parallize this algorithm by considering processes as array elements in above-mentioned web page.

void genRandom(); generate random integers in all processes
void doOddEvenSort(); do odd exchange and even exchange alternately till sorted
    bool mergeMax(DType* a,DType* b,DType* c); extract from a/b and store last half in c
    bool mergeMin(DType* a,DType* b,DType* c); extract from a/b and store first half in c
