Conjugate Gradient Method

This problem is that given a matrix A(n by n) and vector B(n by 1), output a vector X(n by 1), s.t. AX = B.

Conjugate Gradient Method:
    X = 0
    r[0] = B
    p[1] = r[0]
    for k = 1 to oo, until |rk| <= error :
        alpha[k] = (r[k-1]*r[k-1])/(p[k]*A*p[k])
        x[k] = x[k-1]+alpha[k]*p[k]
        r[k] = r[k-1]-alpha[k]*(A*p[k])
        beta[k+1] = (r[k]*r[k])/(r[k-1]*r[k-1])
        p[k+1] = r[k]+beta[k+1]*p[k]
    return the last x[k]
Suppose p (n%p==0) MPI processes started, each holds n/p rows of A,x[],r[], and the whole B,p[] in order to parallize mtx * vec, vec*vec, vec+vec & vec*scalarValue.

void load(); One process reads A & B from a file and transfer corresponding part to other processes
void solve(); main cg method implementation, as described above.
    DType vec_dot(DType* va,DType* vb,int ll); return va*vb
    void vec_addass(DType* va,DType* vb,DType cc,int ll); va = va+vb*cc;
    void vec_mad(DType* vv,DType* va,DType* vb,DType cc,int ll); vv = va+vb*cc;

