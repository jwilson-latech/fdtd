/*Create Differential Operator A */
inline double A_of(
    global double * u,
    global double * v,
    int n,
    int j,
    int k){
        /*double sigma = 1;
        double lamb = 1;
        double dt = 1;*/
        return sigma * (D_of(u,n,j,k))/2.0 - (lamb*dt*abs2(u,v,n,j,k))/2.0*u[id(n,j,k)];
    }; 
    
