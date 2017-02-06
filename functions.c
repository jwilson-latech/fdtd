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
                + 0*sigmay1*( Dx_of(v,n-1,j,k) + Dx_of(v,n-1,j,k-1) )
                + sigmay2*( v[id(n-1,j,k)] + v[id(n-1,j,k-1)] )
                - sigmay3*( pow(v[id(n-1,j,k)] + v[id(n-1,j,k-1)],2) + pow(u[id(n-1,j,k)] + u[id(n-1,j,k-1)],2) )*( v[id(n-1,j,k)] + v[id(n-1,j,k-1)] );
        
            v[id(n,j,k)] = v[id(n-2,j,k-1)] 
                + sigmay0*( v[id(n-2,j,k)] - v[id(n,j,k-1)] )
                - 0*sigmay1*( Dx_of(u,n-1,j,k) + Dx_of(u,n-1,j,k-1) ) 
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
    int Z0 = 50;

    if (region_total(j,k)){

        u[id(2,j,k)] =+cos((2*((j+k)*dx-Z0)-12*thestep*dt))/cosh(((j+k)*dx-Z0)-8*thestep*dt);
        v[id(2,j,k)] =-sin((2*((j+k)*dx-Z0)-12*thestep*dt))/cosh(((j+k)*dx-Z0)-8*thestep*dt);

        u[id(3,j,k)] =+cos((2*((j+k)*dx-Z0)-12*thestep*dt-6*thestep*dt))/cosh(((j+k)*dx-Z0)-8*thestep*dt-4*dt);
        v[id(3,j,k)] =-sin((2*((j+k)*dx-Z0)-12*thestep*dt-6*thestep*dt))/cosh(((j+k)*dx-Z0)-8*thestep*dt-4*dt);
    }
    else{
    }
    };
    
    


        
        