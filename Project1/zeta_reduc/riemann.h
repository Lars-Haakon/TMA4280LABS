//enum Method {ALLREDUCE, RECURSIVEDOUBLE};
typedef enum Method {
    ALLREDUCE = 0,
    RECURSIVEDOUBLE = 1 
} Method;

double riemann(int n, Method method);