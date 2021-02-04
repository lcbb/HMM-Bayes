#define _USE_MATH_DEFINES
#include <math.h>
#include "mex.h"

void rescale(double *F, double *F_prev, double *scaling, mwSize num_states);

/* The computational routine */
void hmm_forward(double *obs, double *p_start, double *p_trans, double *mu_emit, double *sigma_emit, 
        mwSize dim, mwSize num_obs, mwSize num_states,
        double *logprob)
{
    
    mwSize obs_idx;
    mwSize state_idx;
    mwSize dim_idx;
    
    double* F_prev = (double*)mxCalloc(num_states,sizeof(double));
    double* F = (double*)mxCalloc(num_states,sizeof(double));
    double* scaling = (double*)mxCalloc(num_obs,sizeof(double));
    
    /*
    double F_prev[num_states];
    double F[num_states];
    double scaling[num_obs];
    */
    
    *logprob = 0;
    
    /* First time point: F(k) = p(y1|k)*pi(k) */
    
    for (state_idx = 0; state_idx < num_states; state_idx++) {
        double exponent = 0;
        
        double obs_sum = 0;
        
        for (dim_idx = 0; dim_idx < dim; dim_idx++) {
            obs_sum += pow(obs[dim_idx] - mu_emit[state_idx*dim+dim_idx],2);
        }
        
        exponent = - obs_sum / (2 * pow(sigma_emit[state_idx],2));
        
        F[state_idx] = p_start[state_idx] / (pow(2*M_PI,(double)dim/2) * pow(sigma_emit[state_idx],(double)dim)) * exp(exponent);
    }
    
    rescale(F,F_prev,scaling,num_states);
    
    /* Other time points: F(k,t) = p(yt|k)*sum_k'(phi(k',k)F(k',t-1)) */
    
    for (obs_idx = 1; obs_idx < num_obs; obs_idx++) {
        for (state_idx = 0; state_idx < num_states; state_idx++) {
            double exponent = 0;

            double obs_sum = 0;

            for (dim_idx = 0; dim_idx < dim; dim_idx++) {
                obs_sum += pow(obs[obs_idx*dim + dim_idx] - mu_emit[state_idx*dim+dim_idx],2);
            }

            exponent = - obs_sum / (2 * pow(sigma_emit[state_idx],2));
            
            int trans_state_idx;
            
            double F_trans = 0;
            
            for (trans_state_idx = 0; trans_state_idx < num_states; trans_state_idx++) {
                F_trans += F_prev[trans_state_idx] * p_trans[state_idx*num_states + trans_state_idx];
            }
            
            F[state_idx] = F_trans / (pow(2*M_PI,(double)dim/2) * pow(sigma_emit[state_idx],(double)dim)) * exp(exponent);
          
        }

        rescale(F,F_prev,&scaling[obs_idx],num_states);
 
    }
     
    
    /*  Final probability: sum_k(F(k,T))   %obs_prob = sum(F(:,T));
     *  Equal to the product of all the scaling factors
     */
    
    for (obs_idx = 0; obs_idx < num_obs; obs_idx++) {
        *logprob += log(scaling[obs_idx]);
    }
    
    mxFree(F_prev);
    mxFree(F);
    mxFree(scaling);
}

void rescale(double *F, double *F_prev, double *scaling, mwSize num_states) {
    mwSize obs_idx;
    mwSize state_idx;
    
    /* Calculate scaling factor */
    
    *scaling = 0;
    
    for (state_idx = 0; state_idx < num_states; state_idx++) {
        *scaling += F[state_idx];
    }
    
    for (state_idx = 0; state_idx < num_states; state_idx++) {
        F_prev[state_idx] = F[state_idx] / *scaling;
    }
}