#!/bin/bash

echo "How many time points?"
read t

 /Applications/MATLAB_R2021b.app/bin/matlab -nodisplay -r 'addpath(genpath("/Users/fabian/matlab/")); fn_create_Extra_repmap('$t'); exit; '
