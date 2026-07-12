# UEL property and state-variable order

## Real properties (`PROPS`, 33 entries)

| Index | Symbol / meaning |
|---:|---|
| 1 | Branch width `b_f` |
| 2 | Straight-segment initial tangential stiffness `k_uu0` |
| 3 | Peak tangential traction `t_U,max` |
| 4 | Tangential displacement jump at peak traction `delta_U,p` |
| 5 | Tangential displacement jump at end of softening `delta_U,f` |
| 6 | Residual tangential traction `t_U,res` |
| 7 | Normal closure stiffness `k_vc` |
| 8 | Normal opening stiffness `k_vt` |
| 9 | Critical opening displacement `delta_V0` |
| 10 | Basic tangential-normal coupling stiffness `k_uv0` |
| 11 | Dimensionless closure enhancement coefficient `alpha_c` |
| 12 | Dimensionless opening degradation coefficient `alpha_t` |
| 13 | Unloading stiffness factor `beta_K` |
| 14 | Peak-strength degradation factor `beta_tau` |
| 15 | Residual-slip coefficient `lambda_r` |
| 16 | Transverse stiffness `k_ww` |
| 17-23 | Node-region multipliers `eta_k_u`, `eta_k_v`, `eta_k_uv`, `eta_t_u`, `eta_delta_v`, `eta_beta_K`, `eta_beta_tau` |
| 24 | Element node-influence weight `omega_n` |
| 25-27 | Local branch axis `e_U` |
| 28-30 | Local interface-normal axis `e_V` |
| 31-33 | Local transverse axis `e_W` |

## State variables (`SVARS`, 9 per Gauss point; 18 per element)

| Local index | Meaning |
|---:|---|
| 1 | Maximum historical absolute tangential displacement jump |
| 2 | Maximum historical absolute normal displacement jump |
| 3 | Unloading reference displacement |
| 4 | Residual tangential displacement jump |
| 5 | Major-excursion index `m` |
| 6 | Loading-state flag (`-1` unloading, `1` reloading below the old maximum, `2` new major loading) |
| 7 | Current `beta_K` |
| 8 | Current `beta_tau` |
| 9 | Element node-influence weight `omega_n` |

## Node ordering

1. matrix side, start;
2. FRP side, start;
3. FRP side, end;
4. matrix side, end.

The displacement jump follows the sign convention used in the manuscript:
`delta = N1(q1-q2) + N2(q4-q3)` in each local direction.
