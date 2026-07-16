# FRP-grid--concrete mixed-mode interface UEL

This repository contains a four-node, three-dimensional zero-thickness interface UEL for Abaqus/Standard and a Python preprocessor for interface-element insertion, local-frame assignment, node-region identification, and element-property allocation.

## Repository contents

- `uel/frp_grid_interface_uel.f90`: Abaqus UEL source.
- `preprocessing/frp_grid_preprocess.py`: input-file preprocessor.
- `examples/`: a minimal placeholder input and JSON configuration.
- `docs/PROPERTY_ORDER.md`: property order, state variables, and node ordering.
- `tests/test_preprocessor.py`: a basic parser/preprocessor test.

## Element topology

The UEL has four nodes and three translational degrees of freedom per node:

1. matrix/concrete side at the element start;
2. FRP-grid side at the element start;
3. FRP-grid side at the element end;
4. matrix/concrete side at the element end.

The paired nodes on the two sides may share the same coordinates. Local axes are therefore calculated by the preprocessor from the branch and interface-normal directions and passed to the UEL as properties.

## Prepare an Abaqus input file

The original `.inp` file must contain a placeholder `*ELEMENT` block with a unique element-set name, for example:

```text
*ELEMENT, TYPE=U1P, ELSET=INTERFACE_PLACEHOLDER
1, 1, 2, 3, 4
```

For a placeholder with more than four nodes, set `node_order` in the JSON configuration to the four local connectivity positions required by the UEL.

Run:

```bash
python preprocessing/frp_grid_preprocess.py \
  examples/single_element_template.inp \
  examples/config_template.json \
  -o examples/single_element_uel.inp
```

The script writes:

- the `*USER ELEMENT` definition;
- one U1 element and `*UEL PROPERTY` block per interface element;
- local `U`, `V`, and `W` axes;
- the node-influence weight `omega_n`;
- an aggregate element set retaining the original placeholder-set name.

## Run Abaqus

```bash
abaqus job=single_element_uel \
  input=examples/single_element_uel.inp \
  user=uel/frp_grid_interface_uel.f90 \
  interactive
```

The exact command and compatible Fortran compiler depend on the Abaqus release and operating system.


## License

MIT License. See [`LICENSE`](LICENSE).
