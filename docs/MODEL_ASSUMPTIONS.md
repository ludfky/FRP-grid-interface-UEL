# Model assumptions that must be checked before publication

This repository is a reference implementation reconstructed from the equations and workflow described in the manuscript. The manuscript does not fully specify every implementation choice needed by executable code. The following choices are made explicitly here and must be compared with the original production model:

1. The mixed-mode traction is evaluated as `t_U = t_U,base + k_uv delta_V` and `t_V = t_V,base + k_uv delta_U`; the derivative of `k_uv(delta_V)` is omitted from the algorithmic matrix.
2. Unloading and reloading use a line of slope `beta_K k_uu0` passing through the local residual-slip intercept and capped by the current degraded envelope.
3. A new major excursion is counted after an unloading reversal when reloading exceeds the previous historical maximum.
4. The code uses two Gauss points and 9 state variables per Gauss point.
5. Element-specific local axes and node-influence weights are supplied through one `*UEL PROPERTY` block per interface element.
6. The preprocessor accepts a user-defined four-node mapping from any placeholder connectivity. For an eight-node cohesive placeholder, the author must supply the node positions that correspond to the four-node line-interface topology.

Do not describe this reconstruction as the exact code used for the published simulations until the original model curves and state histories have been reproduced.
