#include "mex.h"

#define search  if (hashtable[i] == hashval){position[0] = i+1;break;};
#define BLOCKSIZE (4)

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[])
{
    int i, m, blocklimit;
    double *position,*hashtable,*hashsearch,hashval;
        
    hashtable  = mxGetPr(prhs[0]);
    hashsearch = mxGetPr(prhs[1]);
    hashval    = hashsearch[0];
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    position = mxGetPr(plhs[0]);
    position[0] = 0;
    
    m = mxGetNumberOfElements(prhs[0]);
    i=0;
    blocklimit = (m / BLOCKSIZE) * BLOCKSIZE;
    
    /* Simple loop unrolling */
    while (i<blocklimit)
    {
        search
        i++;
        search
        i++;
        search
        i++;
        search
        i++;
    }
    
    if( i < m )
    {
        switch( m - i )
        {
            case 3 : search; i++;
            case 2 : search; i++;
            case 1 : search;
        }
    }
    
      
/*    for (i = 0; i < 25000;i++)
    {
        if (hashtable[i] == hashval)
        {
            position[0] = i+1;
            break;
        }
        i=i+1;
          if (hashtable[i] == hashval)
        {
            position[0] = i+1;
            break;
        }
    }
 */
    
    if (position[0] == 0)
    {
        plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
    }
}
