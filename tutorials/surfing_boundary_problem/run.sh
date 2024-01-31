# delta=0.6
l=0.35
# Gc=56e-3
refine=1 # ref2=0.05, 
# sigma_ts=81
# T_K=1573 delta=$delta refine=$refine l=$l Gc=$Gc sigma_ts=$sigma_ts \
filebase=nuc2024_test
# rm -rf .jitcache/
nohup mpiexec -n 1 ~/projects/raccoon/raccoon-opt -i elasticity_ldl.i \
--color off --trap-fpe \
filebase=$filebase refine=$refine l=$l \
Executioner/fixed_point_rel_tol=1e-4 \
Executioner/fixed_point_abs_tol=1e-7 \
Executioner/end_time=0.02 \
> ${filebase}.out 2>&1 &

# nohup mpiexec -n 4 ~/projects/raccoon/raccoon-opt -i elasticity.i \
# --color off --trap-fpe \
# Executioner/fixed_point_rel_tol=1e-4 \
# Executioner/fixed_point_abs_tol=1e-7 \
# Executioner/end_time=0.3 \
# > elasticity.out 2>&1 &

