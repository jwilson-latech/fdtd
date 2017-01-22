inline int id(
    int j,
    int k){
        int id_;
        return id_= j + (k%WIDTH) * WIDTH;
    };
    
__kernel void solve1(
        __global double * u, 
        __global double * v
        ){
    
        int gid_j = get_global_id(0);
        int gid_k = get_global_id(1);
        int j = gid_j % WIDTH;
        int k = gid_k % WIDTH;
        
        u[id(j,k)] = j+k;
        v[id(j,k)] = j*k;
        
        };
        
__kernel void solve2(
        __global double * u, 
        __global double * v
        ){

        int gid_j = get_global_id(0);
        int gid_k = get_global_id(1);
        int j = gid_j % WIDTH;
        int k = gid_k % WIDTH;

        u[id(j,k)] = j+k;
        v[id(j,k)] = j*k;

        };