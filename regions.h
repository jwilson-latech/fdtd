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
    };