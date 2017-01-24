inline int id(
    int n,
    int j,
    int k){
        int id_;
        return id_= j%WIDTH + (k%WIDTH) * WIDTH + n*WIDTH*WIDTH;
    };
    
inline double abs2(
    global double * u,
    global double * v,
    int n,
    int j,
    int k){
       return u[id(n,j,k)]*u[id(n,j,k)] + v[id(n,j,k)]*v[id(n,j,k)];
    };

inline double D_of(
    global double * u,
    int n,
    int j,
    int k){
           return -(u[id(n,j+2,k)]+u[id(n,j,k+2)]+u[id(n,j-2,k)]+u[id(n,j,k-2)]-16.0*u[id(n,j+1,k)]-16.0*u[id(n,j,k+1)]-16.0*u[id(n,j-1,k)]-16.0*u[id(n,j,k-1)]+60.0*u[id(n,j,k)])/12.0;
    };

inline double Dy_of(
    global double * u,
    int n,
    int j,
    int k){
           return  (u[id(n,j,k+1)]+u[id(n,j,k-1)]-2.0*u[id(n,j,k)]);
    };

inline double Dx_of(
    global double * u,
    int n,
    int j,
    int k){
           return (u[id(n,j+1,k)]+u[id(n,j-1,k)]-2.0*u[id(n,j,k)]);
    };

inline double D12y_of(
    global double * u,
    int n,
    int j,
    int k){
           return  -(u[id(n,j,k+2)]+u[id(n,j,k-2)]-16.0*u[id(n,j,k+1)]-16.0*u[id(n,j,k-1)]+30.0*u[id(n,j,k)])/12.0;
    };

inline double D12x_of(
    global double * u,
    int n,
    int j,
    int k){
           return -(u[id(n,j+2,k)]+u[id(n,j-2,k)]-16.0*u[id(n,j+1,k)]-16.0*u[id(n,j-1,k)]+30.0*u[id(n,j,k)])/12.0;
    };
    
