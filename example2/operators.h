/*Create Differential Operator A */
inline double A_of(
    global double * u,
    global double * v,
    int n,
    int j,
    int k){
        /* To use as operator instead of multiple, replace 7 with n.*/
        return sigma * (D_of(u,n,j,k))/2.0 - lamb*dt/2.0*(abs2(u,v,7,j,k)-pow(abs2(u,v,7,j,k),1.5))*u[id(n,j,k)];
    }; 
    
    
