# Spenso Python API Documentation

A comprehensive Python tensor library for symbolic and numerical tensor computations, with a focus on physics applications.

## Overview

The Spenso Python API provides powerful tools for:
- **Tensor Algebra**: Dense and sparse tensors with flexible data types
- **Symbolic Computation**: Integration with Symbolica for tensors with symbolic expressions
- **Network Operations**: Tensor networks for optimized computation graphs
- **Physics Applications**: Built-in support for HEP tensors (gamma matrices, color structures, etc.)
- **Performance**: Compiled evaluators for high-speed numerical computation

## Installation

```bash
pip install symbolica
# Spenso is available as part of the symbolica.community.spenso module
```

## Quick Start

```python
from symbolica import S
from symbolica.community.spenso import (
    Tensor,
    TensorIndices,
    TensorStructure,
    TensorName,
    Representation,
    Slot,
    TensorNetwork,
    TensorLibrary,
)

# Create a 3D Minkowski representation
lor = Representation.mink(3)

# Create tensor structure with indices
mu = lor("mu")
nu = lor("nu")
indices = TensorIndices(mu, nu)

# Create dense tensor (3x3 identity matrix)
data = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
tensor = Tensor.dense(indices, data)

# Create symbolic tensors
x, y = S("x", "y")
T = TensorName("T")
symbolic_tensor = T(mu, nu)
```

## Core Classes

### Tensor

The main tensor class supporting both dense and sparse storage with numerical or symbolic data.

```python
from symbolica import S
from symbolica.community.spenso import Representation, Tensor, TensorIndices

# Dense tensor creation
euc3 = Representation.euc(3)
structure = TensorIndices(euc3(1), euc3(2))  # 3x3 tensor
data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
tensor = Tensor.dense(structure, data)

# Sparse tensor creation
sparse_tensor = Tensor.sparse(structure, float)
sparse_tensor[0, 0] = 1.0
sparse_tensor[1, 1] = 2.0

# Symbolic tensors
x, y = S("x", "y")
sym_data = [x, y, x * y, x + y]
euc2 = Representation.euc(2)
structure_2x2 = TensorIndices(euc2(1), euc2(2))
sym_tensor = Tensor.dense(structure_2x2, sym_data)

# Element access
tensor[0] = 10.0  # Set element
tensor.to_sparse()  # Convert to sparse storage
```

### Representations

Define the mathematical properties of tensor indices for group theory applications.

```python
from symbolica.community.spenso import Representation

# Standard physics representations
euclidean = Representation.euc(4)  # 4D Euclidean space
minkowski = Representation.mink(4)  # 4D Minkowski space
color_fund = Representation.cof(3)  # SU(3) fundamental
color_adj = Representation.coad(8)  # SU(3) adjoint
bispinor = Representation.bis(4)  # Dirac bispinor

# Custom representations
custom = Representation("MyRep", 5, is_self_dual=False)

# Create slots (representation + index)
mu_slot = euclidean("mu")

# Metric tensors
metric = euclidean.g("mu", "nu")  # Metric tensor g_μν
flat_metric = euclidean.flat("mu", "nu")  # Flat metric η_μν
identity = custom.id("mu", "nu")  # Identity δ_μν
```

### TensorName

Named tensor functions with mathematical properties like symmetry.

```python
from symbolica.community.spenso import TensorName, Representation

# Create tensor names
T = TensorName("T")  # Basic tensor
g = TensorName("g", is_symmetric=True)  # Symmetric tensor
F = TensorName("F", is_antisymmetric=True)  # Antisymmetric tensor

# Predefined physics tensors
gamma = TensorName.gamma  # Dirac gamma matrices
sigma = TensorName.sigma  # Pauli matrices
f_abc = TensorName.f  # SU(N) structure constants

# Use with indices
rep = Representation.mink(3)
mu = rep("mu")
nu = rep("nu")
tensor_expr = T(mu, nu)  # Creates TensorIndices T(μ,ν)
```

### TensorIndices and TensorStructure

Define tensor shapes and index structures.

```python
from symbolica import S
from symbolica.community.spenso import (
    TensorIndices,
    TensorStructure,
    TensorName,
    Representation,
)

# TensorIndices: For tensors with specific abstract indices
rep = Representation.mink(3)
mu = rep("mu")
nu = rep("nu")

# TensorStructure: For tensor templates without specific indices
structure = TensorStructure(rep, rep)

# Named structures
T = TensorName("T")
named_structure = TensorStructure(rep, rep, name=T)

# Convert structure to indices
concrete_indices = structure.index("a", "b")  # Assign indices a, b

# Create symbolic expressions
symbolic_expr = named_structure.symbolic("mu", "nu")  # T(μ,ν)

# With additional arguments
x = S("x")
expr_with_args = named_structure.symbolic(x, ";", "mu", "nu")  # T(x; μ,ν)
```

### TensorNetwork

Computational graphs for optimized tensor operations.

```python
from symbolica import E
from symbolica.community.spenso import (
    TensorNetwork,
    TensorName,
    ExecutionMode,
    Representation,
)

# TensorIndices: For tensors with specific abstract indices
rep = Representation.mink(3)
mu = rep("mu")
nu = rep("nu")
x = E("sin(x)")
T = TensorName("T")
expr = x * T(mu, nu)
network = TensorNetwork(expr)
# Controlled execution
network.execute(mode=ExecutionMode.Scalar)  # Only scalar operations
result = network.result_tensor()
```

### TensorLibrary

Registry for reusable tensor definitions.

```python
from symbolica import S
from symbolica.community.spenso import (
    TensorLibrary,
    LibraryTensor,
    TensorName as N,
    TensorStructure,
    Representation,
)

# Create library
lib = TensorLibrary()

# Register tensors
rep = Representation.euc(2)

sigma_x = S("σ_x")
structure = TensorStructure(rep, rep, name=sigma_x)
pauli_x = LibraryTensor.dense(structure, [0.0, 1.0, 1.0, 0.0])
lib.register(pauli_x)

# Access registered tensors
sigma_structure = lib[sigma_x]

# HEP library with standard physics tensors
hep_lib = TensorLibrary.hep_lib()
gamma_structure = hep_lib[N.gamma]
```

## Evaluation and Performance

### Symbolic Evaluation

```python
# Create symbolic tensor
x, y = S('x','y')
structure = TensorIndices(Representation.euc(2)(1), Representation.euc(2)(1))
tensor = Tensor.dense(structure, [x*y, x+y])

# Create evaluator
evaluator = tensor.evaluator(
    constants={},           # Constant values
    funs={},               # Function definitions
    params=[x, y],         # Variable parameters
    iterations=100,        # Optimization iterations
    n_cores=4             # Parallel cores
)

# Evaluate for multiple parameter sets
inputs = [[1.0, 2.0], [3.0, 4.0]]  # x,y values
results = evaluator.evaluate(inputs)
```

### Compiled Evaluation

```python
# Compile for maximum performance
compiled_evaluator = evaluator.compile(
    function_name="fast_eval",
    filename="tensor_eval.cpp",
    library_name="tensor_lib",
    optimization_level=3
)

# Use compiled evaluator for complex inputs
complex_inputs = [[1.0+2.0j, 3.0+0.0j]]
results = compiled_evaluator.evaluate_complex(complex_inputs)
```

## Physics Applications

### High Energy Physics

```python
from symbolica.community.spenso import TensorLibrary, TensorName, Representation

# Load HEP library
hep_lib = TensorLibrary.hep_lib()

# Standard HEP tensors
gamma = TensorName.gamma  # Dirac gamma matrices
gamma5 = TensorName.gamma5  # γ₅ matrix
sigma = TensorName.sigma  # Pauli matrices
f_abc = TensorName.f  # SU(N) structure constants
t_a = TensorName.t  # SU(N) generators

# Representations
lorentz = Representation.mink(4)  # Lorentz indices
spinor = Representation.bis(4)  # Dirac spinor

# Build physics expressions
mu = lorentz("mu")
nu = lorentz("nu")
alpha = spinor("alpha")
beta = spinor("beta")

# Gamma matrix trace: Tr(γᵘγᵥ)
gamma_trace = gamma(alpha, beta, mu) * gamma(beta, alpha, nu)
```

### Tensor Contractions

```python
from symbolica.community.spenso import (
    TensorLibrary,
    TensorName,
    Representation,
    TensorNetwork,
)

rep = Representation.mink(4)
mu = rep("mu")
nu = rep("nu")
rho = rep("rho")

# Metric contraction: gᵘᵛ gᵤᵥ = 4
g = TensorName.g
contraction = g(mu, nu) * g(mu, nu)  # Repeated indices contract

# Network evaluation
network = TensorNetwork(contraction)
network.execute()
result = network.result_scalar()  # Should give 4 for 4D
```


# Best Practices

1. **Use appropriate storage**: Sparse for mostly-zero tensors, dense for general data
2. **Name your tensors**: Required for symbolic manipulation and library registration
3. **Choose representations carefully**: Self-dual vs dualizable affects contraction rules
4. **Optimize evaluation**: Use compiled evaluators for performance-critical code
5. **Leverage libraries**: Use `TensorLibrary.hep_lib()` for standard physics tensors
6. **Cook your indices**: Indices have to be symbols or numbers. To that end, use the cook_indices function from idenso.
