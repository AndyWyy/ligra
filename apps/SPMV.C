#define WEIGHTED 1
#include "ligra.h"
#include "math.h"

struct SPMV_F{
    uintE *vectors;
    uintE *flags;
    intE *result;
    intE round;
    SPMV_F(uintE *_vectors, uintE *_flags, intE* _result, intE _round) :
     vectors(_vectors), flags(_flags), result(_result),round(_round){}

    inline bool update(uintE s, uintE d, intE edgeLen){
         if(flags[d] == UINT_E_MAX){
            flags[d] = s; 
            if(s == round){
                *result += edgeLen * vectors[round];
            }
            return true;
        }
        else return false;
    }

    inline bool updateAtomic(uintE s, uintE d, intE edgeLen){
        if(CAS(&flags[d],UINT_E_MAX,s) == true){
            if(s == round){
                intE res = edgeLen * vectors[round]; 
                writeAdd(result, res);
            }
            return true;
        }
        else return false;
        
    }

    inline bool cond(uintE d){
        return (UINT_E_MAX == flags[d]);
    }

};
template <class vertex>
void Compute(graph<vertex>& GA, commandLine P){
    
    long start = P.getOptionLongValue("-r", 0);
    long n = GA.n;
    uintE *flags = newA(uintE,n);
    parallel_for(long i = 0; i < n; i++) 
        flags[i] = UINT_E_MAX;
    intE * outputs = newA(intE ,n);
    parallel_for(long i = 0; i < n ; i++){
        outputs[i] = 0;
        if(GA.V[i].getOutDegree() != 0){
            for(int j = 0 ; j < GA.V[i].getOutDegree() ; j++){
                if(i == GA.V[i].getOutNeighbor(j)) outputs[i] = GA.V[i].getOutWeight(j);
            }
        }
    }
    uintE *vectors = newA(uintE,n);
    parallel_for(long i = 0; i < n ;i++){
        vectors[i] = 1;
    }
    parallel_for(long i = 0 ; i < n ; i++){
        vertexSubset Frontier(n,i);
        flags[i] = i;
        while(!Frontier.isEmpty()){
            vertexSubset out = edgeMap(GA, Frontier, SPMV_F(vectors,flags,&outputs[i],i));
            Frontier.del();
            Frontier = out;
        }
        parallel_for(long i = 0; i < n; i++)
            flags[i] = UINT_E_MAX;
        Frontier.del();
    }
    for(int i = 0 ; i < n ; i++){
        printf("%d ", outputs[i]);
    }
    free(vectors);
    free(outputs);
}