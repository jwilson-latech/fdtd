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
    
inline bool region_total(
    int j,
    int k){
        if ( (j>=0) && (j<WIDTH) && (k>=0) && (k<WIDTH) ){
            return true;
        }
        else { 
            return false;
        }
    };

inline bool region_main(
    int j,
    int k){
        if ( (j>=(bnd1+bnd2)) && (j<(WIDTH-(bnd1+bnd2))) && (k>=(bnd1+bnd2)) && (k<(WIDTH-(bnd1+bnd2))) ){
            return true;
        }
        else { 
            return false;
        }
    };

inline bool region_north1(
    int j,
    int k){
        if ( (k>=(WIDTH-(bnd1+bnd2))) && (k<(WIDTH-bnd1)) && (j>(bnd1-1)) && (j<(WIDTH-bnd1)) ){
            return true;
        }
        else { 
            return false;
        }
    };

inline bool region_west1(
    int j,
    int k){
        if ( (j>(1)) && (j<(bnd1+bnd2)) && (k>(bnd1-1)) && (k<(WIDTH-bnd1)) ){
            return true;
        }
        else { 
            return false;
        }
    };

inline bool region_south1(
    int j,
    int k){
        if ( (k>(1)) && (k<(bnd1+bnd2)) && (j>(bnd1-1)) && (j<(WIDTH-bnd1)) ){
            return true;
        }
        else { 
            return false;
        }
    };

inline bool region_east1(
    int j,
    int k){
        if ( (j>=(WIDTH-(bnd1+bnd2))) && (j<(WIDTH-bnd1)) && (k>(bnd1-1)) && (k<(WIDTH-bnd1)) ){
            return true;
        }
        else { 
            return false;
        }
    };

inline bool region_east(
    int j,
    int k){
        if ( (j>=(WIDTH-bnd1)) && (j<(WIDTH)) && (k>(bnd1-1)) && (k<(WIDTH-bnd1)) ){
            return true;
        }
        else { 
            return false;
        }
    };

inline bool region_north(
    int j,
    int k){
        if ( (k>=(WIDTH-bnd1)) && (k<(WIDTH)) && (j>(bnd1-1)) && (j<(WIDTH-bnd1)) ){
            return true;
        }
        else { 
            return false;
        }
    };

inline bool region_south(
    int j,
    int k){
        if ( (k>=0) && (k<bnd1) && (j>(bnd1-1)) && (j<(WIDTH-bnd1)) ){
            return true;
        }
        else { 
            return false;
        }
    };

inline bool region_west(
    int j,
    int k){
        if ( (j>=0) && (j<bnd1) && (k>(bnd1-1)) && (k<(WIDTH-bnd1)) ){
            return true;
        }
        else { 
            return false;
        }
    };

inline bool region_main1(    
    int j,
    int k){
        if (region_south1(j,k)||region_north1(j,k)||region_west1(j,k)||region_east1(j,k)){
            return true;
        }
        else { 
            return false;
        }
    }

inline bool region_northeast(
    int j,
    int k){
        if ( (j>=(WIDTH-bnd1)) && (j<(WIDTH)) && (k>=(WIDTH-bnd1)) && (k<(WIDTH)) ){
            return true;
        }
        else { 
            return false;
        }
    };
inline bool region_southeast(
    int j,
    int k){
        if ( (j>=(WIDTH-bnd1)) && (j<(WIDTH)) && (k>=0) && (k<bnd1)  ){
            return true;
        }
        else { 
            return false;
        }
    };

inline bool region_northwest(
    int j,
    int k){
        if ( (j>=(0)) && (j<(bnd1))  && (k>=(WIDTH-bnd1)) && (k<(WIDTH)) ){
            return true;
        }
        else { 
            return false;
        }
    };
inline bool region_southwest(
    int j,
    int k){
        if ( (j>=(0)) && (j<(bnd1)) && (k>=0) && (k<bnd1)  ){
            return true;
        }
        else { 
            return false;
        }
    };/*Create Differential Operator A */
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
    
__kernel void operator_A(
        __global double * u, 
        __global double * v,
        __global double * A1_u,
        __global double * A1_v
        ){
    
        int gid_j = get_global_id(1);
        int gid_k = get_global_id(2);
        int j = gid_j % WIDTH;
        int k = gid_k % WIDTH;
        if(region_total(j,k)){
            A1_u[id(1,j,k)] = A_of(u,v,1,j,k);
            A1_v[id(1,j,k)] = A_of(v,u,1,j,k);
        }
        
        };
        
__kernel void gfdtd(
        __global double * u, 
        __global double * v,
        __global double * A1_u,
        __global double * A1_v,
        __global double * A3_u,
        __global double * A3_v
        ){

        int gid_j = get_global_id(1);
        int gid_k = get_global_id(2);
        int j = gid_j % WIDTH;
        int k = gid_k % WIDTH;
        if(region_main(j,k)){
            u[id(2,j,k)] = u[id(0,j,k)] + 2.0*A1_v[id(1,j,k)] - 1.0/12.0*A3_v[id(1,j,k)];
            v[id(2,j,k)] = v[id(0,j,k)] - 2.0*A1_u[id(1,j,k)] + 1.0/12.0*A3_u[id(1,j,k)];
        }
        else if(region_main1(j,k)){
            u[id(2,j,k)] = u[id(0,j,k)] + 2.0*A1_v[id(1,j,k)];
            v[id(2,j,k)] = v[id(0,j,k)] - 2.0*A1_u[id(1,j,k)];
        }
        else{
            u[id(2,j,k)] = u[id(0,j,k)];
            v[id(2,j,k)] = v[id(0,j,k)];
        }

        };
        
__kernel void restore(
        __global double * u, 
        __global double * v
        ){
        int gid_j = get_global_id(1);
        int gid_k = get_global_id(2);
        int j = gid_j % WIDTH;
        int k = gid_k % WIDTH;
        u[id(0,j,k)] = u[id(1,j,k)];
        v[id(0,j,k)] = v[id(1,j,k)];
        barrier(CLK_GLOBAL_MEM_FENCE);
        u[id(1,j,k)] = u[id(2,j,k)];
        v[id(1,j,k)] = v[id(2,j,k)];
        barrier(CLK_GLOBAL_MEM_FENCE);
        u[id(2,j,k)] = 0;
        v[id(2,j,k)] = 0;
        };

__kernel void abc(
    __global double * u, 
    __global double * v){


    int gid_j = get_global_id(1);
    int gid_k = get_global_id(2);
    int j = gid_j % WIDTH;
    int k = gid_k % WIDTH;

    double sigmax0 = (1-cx)/(1+cx);
    double sigmax1 = dt/(1+cx)/dx/dx;
    double sigmax2 = dt/(1+cx)*v2x;
    double sigmax3 = dt/(1+cx)*lamb/4.0;
    double sigmay0 = (1-cy)/(1+cy);
    double sigmay1 = dt/(1+cy)/dy/dy;
    double sigmay2 = dt/(1+cy)*v2y;
    double sigmay3 = dt/(1+cy)*lamb/4.0;
    if (region_main(j,k)||region_main1(j,k)){
        u[id(2,j,k)] = u[id(2,j,k)];
        v[id(2,j,k)] = v[id(2,j,k)];
    }
    else if ( region_east(j,k)){
        
        u[id(2,j,k)] = u[id(0,j-1,k)] 
            + sigmax0*( u[id(0,j,k)] - u[id(2,j-1,k)] )
            + sigmax1*( Dy_of(v,1,j,k) + Dy_of(v,1,j-1,k) )
            + sigmax2*( v[id(1,j,k)] + v[id(1,j-1,k)] )
            - sigmax3*( pow(v[id(1,j,k)] + v[id(1,j-1,k)],2) + pow(u[id(1,j,k)] + u[id(1,j-1,k)],2) )*( v[id(1,j,k)] + v[id(1,j-1,k)] );
        
        v[id(2,j,k)] = v[id(0,j-1,k)] 
            + sigmax0*( v[id(0,j,k)] - v[id(2,j-1,k)] )
            - sigmax1*( Dy_of(u,1,j,k) + Dy_of(u,1,j-1,k) ) 
            - sigmax2*( u[id(1,j,k)] + u[id(1,j-1,k)] )
            + sigmax3*( pow(v[id(1,j,k)] + v[id(1,j-1,k)],2) + pow(u[id(1,j,k)] + u[id(1,j-1,k)],2) )*( u[id(1,j,k)] + u[id(1,j-1,k)] );
                
    }
    else if ( region_north(j,k)||region_northeast(j,k)){
        
        u[id(2,j,k)] = u[id(0,j,k-1)] 
            + sigmay0*( u[id(0,j,k)] - u[id(2,j,k-1)] )
            + sigmay1*( Dx_of(v,1,j,k) + Dx_of(v,1,j,k-1) )
            + sigmay2*( v[id(1,j,k)] + v[id(1,j,k-1)] )
            - sigmay3*( pow(v[id(1,j,k)] + v[id(1,j,k-1)],2) + pow(u[id(1,j,k)] + u[id(1,j,k-1)],2) )*( v[id(1,j,k)] + v[id(1,j,k-1)] );
        
        v[id(2,j,k)] = v[id(0,j,k-1)] 
            + sigmay0*( v[id(0,j,k)] - v[id(2,j,k-1)] )
            - sigmay1*( Dx_of(u,1,j,k) + Dx_of(u,1,j,k-1) ) 
            - sigmay2*( u[id(1,j,k)] + u[id(1,j,k-1)] )
            + sigmay3*( pow(v[id(1,j,k)] + v[id(1,j,k-1)],2) + pow(u[id(1,j,k)] + u[id(1,j,k-1)],2) )*( u[id(1,j,k)] + u[id(1,j,k-1)] );
                
    }

    
    else{
        u[id(2,j,k)] = 0;
        v[id(2,j,k)] = 0;
    }
    };

        
        