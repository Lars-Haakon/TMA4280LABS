typedef enum Method {
    ALLREDUCE = 0,
    RECURSIVEDOUBLE = 1 
} Method;

double arctan(double x, int n, Method method);
double machin(int n, Method method);