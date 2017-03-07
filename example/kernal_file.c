inline int id(
    int n,
    int j,
    int k){
        int id_;
        return id_= j%WIDTH + (k%WIDTH) * WIDTH + n*WIDTH*WIDTH;
    };
    
inline int id2(
    int n,
    int j,
    int k){
        int id_;
        return id_= k%WIDTH + (j%WIDTH) * WIDTH + n*WIDTH*WIDTH;
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
        /* To use as operator instead of multiple, replace 7 with n.*/
        return sigma * (D_of(u,n,j,k))/2.0 - (lamb*dt*abs2(u,v,7,j,k))/2.0*u[id(n,j,k)];
    }; 
    
    
__kernel void operator_A(
        __global double * u, 
        __global double * v,
        int n){
    
        int gid_j = get_global_id(1);
        int gid_k = get_global_id(2);
        int j = gid_j % WIDTH;
        int k = gid_k % WIDTH;
        if(region_total(j,k)){
            u[id(4,j,k)] = A_of(u,v,n-1,j,k);
            v[id(4,j,k)] = A_of(v,u,n-1,j,k);
            barrier(CLK_GLOBAL_MEM_FENCE);
            u[id(5,j,k)] = A_of(u,v,4,j,k);
            v[id(5,j,k)] = A_of(v,u,4,j,k);
            barrier(CLK_GLOBAL_MEM_FENCE);
            u[id(6,j,k)] = A_of(u,v,5,j,k);
            v[id(6,j,k)] = A_of(v,u,5,j,k); 
        }
        };
        
__kernel void operator_mover(
        __global double * u, 
        __global double * v,
        int n){

        int gid_j = get_global_id(1);
        int gid_k = get_global_id(2);
        int j = gid_j % WIDTH;
        int k = gid_k % WIDTH;
        if(region_total(j,k)){
            u[id(7,j,k)] = u[id(n-1,j,k)];
            v[id(7,j,k)] = v[id(n-1,j,k)];
        }
        };
        
__kernel void gfdtd(
        __global double * u, 
        __global double * v,
        int n){

        int gid_j = get_global_id(1);
        int gid_k = get_global_id(2);
        int j = gid_j % WIDTH;
        int k = gid_k % WIDTH;
        if(region_main(j,k)){
            u[id(n,j,k)] = u[id(n-2,j,k)] + 2.0*v[id(4,j,k)] - 1.0/12.0*v[id(6,j,k)];
            v[id(n,j,k)] = v[id(n-2,j,k)] - 2.0*u[id(4,j,k)] + 1.0/12.0*u[id(6,j,k)];
        }
        else if(region_main1(j,k)){
            u[id(n,j,k)] = u[id(n-2,j,k)] + 2.0*v[id(4,j,k)];
            v[id(n,j,k)] = v[id(n-2,j,k)] - 2.0*u[id(4,j,k)];
        }
        else{
            u[id(n,j,k)] = u[id(n-2,j,k)];
            v[id(n,j,k)] = v[id(n-2,j,k)];
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
        
        u[id(0,j,k)] = u[id(2,j,k)];
        v[id(0,j,k)] = v[id(2,j,k)];
        u[id(1,j,k)] = u[id(3,j,k)];
        v[id(1,j,k)] = v[id(3,j,k)];
        
        };
        


__kernel void abc(
    __global double * u, 
    __global double * v,
    int n){


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
    
    int iter;
    for ( iter=0; iter<2; iter=iter+1 ){
        
        if (region_main(j,k)||region_main1(j,k)){
            /* Don't do anything here.*/
        }
        
        else if (region_east(j,k)){
        
            u[id(n,j,k)] = u[id(n-2,j-1,k)] 
                + sigmax0*( u[id(n-2,j,k)] - u[id(n,j-1,k)] )
                + sigmax1*( Dy_of(v,n-1,j,k) + Dy_of(v,n-1,j-1,k) )
                + sigmax2*( v[id(n-1,j,k)] + v[id(n-1,j-1,k)] )
                - sigmax3*( pow(v[id(n-1,j,k)] + v[id(n-1,j-1,k)],2) + pow(u[id(n-1,j,k)] + u[id(n-1,j-1,k)],2) )*( v[id(n-1,j,k)] + v[id(n-1,j-1,k)] );
        
            v[id(n,j,k)] = v[id(n-2,j-1,k)] 
                + sigmax0*( v[id(n-2,j,k)] - v[id(n,j-1,k)] )
                - sigmax1*( Dy_of(u,n-1,j,k) + Dy_of(u,n-1,j-1,k) ) 
                - sigmax2*( u[id(n-1,j,k)] + u[id(n-1,j-1,k)] )
                + sigmax3*( pow(v[id(n-1,j,k)] + v[id(n-1,j-1,k)],2) + pow(u[id(n-1,j,k)] + u[id(n-1,j-1,k)],2) )*( u[id(n-1,j,k)] + u[id(n-1,j-1,k)] );
                
        }
    
        else if (region_north(j,k)){
        
            u[id(n,j,k)] = u[id(n-2,j,k-1)] 
                + sigmay0*( u[id(n-2,j,k)] - u[id(n,j,k-1)] )
                + sigmay1*( Dx_of(v,n-1,j,k) + Dx_of(v,n-1,j,k-1) )
                + sigmay2*( v[id(n-1,j,k)] + v[id(n-1,j,k-1)] )
                - sigmay3*( pow(v[id(n-1,j,k)] + v[id(n-1,j,k-1)],2) + pow(u[id(n-1,j,k)] + u[id(n-1,j,k-1)],2) )*( v[id(n-1,j,k)] + v[id(n-1,j,k-1)] );
        
            v[id(n,j,k)] = v[id(n-2,j,k-1)] 
                + sigmay0*( v[id(n-2,j,k)] - v[id(n,j,k-1)] )
                - sigmay1*( Dx_of(u,n-1,j,k) + Dx_of(u,n-1,j,k-1) ) 
                - sigmay2*( u[id(n-1,j,k)] + u[id(n-1,j,k-1)] )
                + sigmay3*( pow(v[id(n-1,j,k)] + v[id(n-1,j,k-1)],2) + pow(u[id(n-1,j,k)] + u[id(n-1,j,k-1)],2) )*( u[id(n-1,j,k)] + u[id(n-1,j,k-1)] );
                
        }
        else if (region_south(j,k)){

            u[id(n,j,k)] = u[id(n-2,j,k+1)]
                + sigmay0*( u[id(n-2,j,k)] - u[id(n,j,k+1)] )
                + sigmay1*( Dx_of(v,n-1,j,k) + Dx_of(v,n-1,j,k+1) )
                + sigmay2*( v[id(n-1,j,k)] + v[id(n-1,j,k+1)] )
                - sigmay3*( pow(v[id(n-1,j,k)] + v[id(n-1,j,k+1)],2) + pow(u[id(n-1,j,k)] + u[id(n-1,j,k+1)],2) )*( v[id(n-1,j,k)] + v[id(n-1,j,k+1)] );

            v[id(n,j,k)] = v[id(n-2,j,k+1)]
                + sigmay0*( v[id(n-2,j,k)] - v[id(n,j,k+1)] )
                - sigmay1*( Dx_of(u,n-1,j,k) + Dx_of(u,n-1,j,k+1) )
                - sigmay2*( u[id(n-1,j,k)] + u[id(n-1,j,k+1)] )
                + sigmay3*( pow(v[id(n-1,j,k)] + v[id(n-1,j,k+1)],2) + pow(u[id(n-1,j,k)] + u[id(n-1,j,k+1)],2) )*( u[id(n-1,j,k)] + u[id(n-1,j,k+1)] );

        }
        
        else if (region_west(j,k)){

            u[id(n,j,k)] = u[id(n-2,j+1,k)]
                + sigmax0*( u[id(n-2,j,k)] - u[id(n,j+1,k)] )
                + sigmax1*( Dy_of(v,n-1,j,k) + Dy_of(v,n-1,j+1,k) )
                + sigmax2*( v[id(n-1,j,k)] + v[id(n-1,j+1,k)] )
                - sigmax3*( pow(v[id(n-1,j,k)] + v[id(n-1,j+1,k)],2) + pow(u[id(n-1,j,k)] + u[id(n-1,j+1,k)],2) )*( v[id(n-1,j,k)] + v[id(n-1,j+1,k)] );

            v[id(n,j,k)] = v[id(n-2,j+1,k)]
                + sigmax0*( v[id(n-2,j,k)] - v[id(n,j+1,k)] )
                - sigmax1*( Dy_of(u,n-1,j,k) + Dy_of(u,n-1,j+1,k) )
                - sigmax2*( u[id(n-1,j,k)] + u[id(n-1,j+1,k)] )
                + sigmax3*( pow(v[id(n-1,j,k)] + v[id(n-1,j+1,k)],2) + pow(u[id(n-1,j,k)] + u[id(n-1,j+1,k)],2) )*( u[id(n-1,j,k)] + u[id(n-1,j+1,k)] );
        }
        
        else if (region_northeast(j,k)){

            u[id(n,j,k)] = u[id(n-2,j,k-1)] 
                + sigmay0*( u[id(n-2,j,k)] - u[id(n,j,k-1)] )
                + sigmay2*( v[id(n-1,j,k)] + v[id(n-1,j,k-1)] )
                - sigmay3*( pow(v[id(n-1,j,k)] + v[id(n-1,j,k-1)],2) + pow(u[id(n-1,j,k)] + u[id(n-1,j,k-1)],2) )*( v[id(n-1,j,k)] + v[id(n-1,j,k-1)] );
        
            v[id(n,j,k)] = v[id(n-2,j,k-1)] 
                + sigmay0*( v[id(n-2,j,k)] - v[id(n,j,k-1)] )
                - sigmay2*( u[id(n-1,j,k)] + u[id(n-1,j,k-1)] )
                + sigmay3*( pow(v[id(n-1,j,k)] + v[id(n-1,j,k-1)],2) + pow(u[id(n-1,j,k)] + u[id(n-1,j,k-1)],2) )*( u[id(n-1,j,k)] + u[id(n-1,j,k-1)] );
                
        }
        
    
    }
    };
    
    
__kernel void ebc(
    __global double * u, 
    __global double * v,
    int n,
    int thestep){


    int gid_j = get_global_id(1);
    int gid_k = get_global_id(2);
    int j = gid_j % WIDTH;
    int k = gid_k % WIDTH;

    if ( region_north(j,k)||region_east(j,k) ){

        u[id(2,j,k)] =+cos((k2x*j+k2y*k)*dx-k2z0-w2*thestep*dt)/cosh((k1x*j+k1y*k)*dx-k1z0-w1*thestep*dt);
        v[id(2,j,k)] =-sin((k2x*j+k2y*k)*dx-k2z0-w2*thestep*dt)/cosh((k1x*j+k1y*k)*dx-k1z0-w1*thestep*dt);

        u[id(3,j,k)] =+cos((k2x*j+k2y*k)*dx-k2z0-w2*(thestep+1)*dt)/cosh((k1x*j+k1y*k)*dx-k1z0-w1*(thestep+1)*dt);
        v[id(3,j,k)] =-sin((k2x*j+k2y*k)*dx-k2z0-w2*(thestep+1)*dt)/cosh((k1x*j+k1y*k)*dx-k1z0-w1*(thestep+1)*dt);
    }
    else{
    }
    };
    
__kernel void nbc(
    __global double * u, 
    __global double * v,
    int n,
    int thestep){


    int gid_j = get_global_id(1);
    int gid_k = get_global_id(2);
    int j = gid_j % WIDTH;
    int k = gid_k % WIDTH;

    if ( region_north(j,k)||region_east(j,k) ){
        u[id(0,j,k)] =0;
        v[id(0,j,k)] =0;

        u[id(1,j,k)] =0;
        v[id(1,j,k)] =0;
        
        u[id(2,j,k)] =0;
        v[id(2,j,k)] =0;

        u[id(3,j,k)] =0;
        v[id(3,j,k)] =0;
    }
    else{
    }
    };

    
    


        
        