function [ output_args ] = make_hmm_forward_c( input_args )
    mex hmm_forward_entry.c hmm_forward.c -output hmm_forward
end

