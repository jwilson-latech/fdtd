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

        
        