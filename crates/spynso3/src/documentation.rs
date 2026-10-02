//! Documentation shared by runtime methods and manually specified Python overloads.
//!
//! Most public docs live next to their definitions. These methods require explicit
//! overload metadata, so one literal supplies both PyO3 and the generated stubs.

macro_rules! python_doc {
    ("Tensor.axes") => {
        r###"External axes in the current component and port-position order.

Returns
-------
tuple of Slot or Representation
    Indexed or unresolved axes, respectively. Component coordinates and methods
    taking axis positions use this order. It follows construction and explicit
    permute_axes() calls; structure.axes instead gives the canonical signature.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName
>>> r = Representation.euc(3)
>>> A = TensorName("A")(r("j"), r("i"))
>>> A.axes == (r("j"), r("i"))
True
>>> A.permute_axes([1, 0]).axes == A.axes[::-1]
True"###
    };
    ("Representation.__call__") => {
        r###"Attach an abstract index label to this space.

Parameters
----------
aind : int, str, or Expression
    An index label, not a component coordinate. A string, Python integer, or
    single symbol creates a Slot. Other Expressions retain raw representation
    syntax: supported numeric, tagged named, and scoped labels can be admitted
    directly; other compound payloads require ``intern="indices"`` when
    constructing a tensor.

Returns
-------
Slot or Expression
    A typed index, or symbolic syntax for a compound index payload.

Examples
--------
>>> from symbolica import S
>>> from symbolica.community.tensor import Representation, Slot
>>> space = Representation.euc(3)
>>> isinstance(space(S("i")), Slot)
True"###
    };
    ("TensorName.__call__") => {
        r###"Construct a tensor with scalar arguments and ordered axes.

Parameters
----------
*args : scalar expression, Slot, or Representation
    Scalar arguments first, then labeled Slots or unresolved Representations.
    A compact vector may bind an axis of a generic tensor, producing a
    contraction instead of a stored-data descriptor.

Returns
-------
TensorExpression
    The tensor call, including a scalar tensor when there are no axes.
    Repeated compatible explicit labels are contracted.

Notes
-----
For predefined metrics, Dirac matrices, and color tensors, prefer the
corresponding TensorExpression factory. Fixed-structure names reject direct
calls. Symmetries act on the full argument list; an antisymmetric name
called twice with exactly the same representation argument is zero.

Examples
--------
>>> from symbolica import S
>>> from symbolica.community.tensor import TensorName, Representation
>>> space = Representation.euc(3)
>>> tensor = TensorName("B")(S("x"), 7, space("i"), space)
>>> tensor.rank
2"###
    };
    ("BroadcastFunction.__call__") => {
        r###"Apply this unary function to a scalar or to every tensor component.

Parameters
----------
arg : scalar expression, TensorExpression, Tensor, or TensorNetwork
    Value to transform independently at each component.

Returns
-------
Expression, TensorExpression, or TensorNetwork
    Scalar inputs return Expression; symbolic tensors retain
    TensorExpression; component-bearing inputs return TensorNetwork.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> from symbolica.community.tensor import Tensor
>>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
>>> from symbolica.community.tensor import BroadcastFunction
>>> result = BroadcastFunction.conj()(tensor).to_tensor()
>>> result.shape
(2, 2)"###
    };
    ("Tensor.__getitem__") => {
        r###"Read components in logical row-major order.

Parameters
----------
item : int, slice, or tuple of int or slice
    An integer is a flat position. A tuple gives one selector per axis.
    A flat slice returns a list; coordinate slices return nested lists.
    Negative indices count from the end.

Returns
-------
Expression, float, complex, or list
    A scalar component, or lists for the sliced axes.

Notes
-----
Coordinates follow axes, the current component view; structure.axes is a canonical signature.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> from symbolica.community.tensor import Tensor
>>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
>>> tensor[:, 1]
[2.0, 4.0]
>>> tensor[-1]
4.0"###
    };
    ("Tensor.__setitem__") => {
        r###"Change one stored component in place.

Parameters
----------
item : int or sequence of int
    Flat logical row-major position or one coordinate per logical axis.
    Negative indices count from the end.
value : float, complex, or Expression
    Replacement matching the tensor's component type: float for real
    storage, complex for complex storage, or Expression for symbolic storage.

Returns
-------
None
    The tensor is modified in place.

Notes
-----
Coordinates follow axes, the current component view; structure.axes is a canonical signature.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> from symbolica.community.tensor import Tensor
>>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
>>> tensor[1, 0] = 7.0
>>> tensor[1, 0]
7.0"###
    };
    ("Tensor.__iter__") => {
        r###"Iterate over all components in logical row-major order.

Returns
-------
iterator of Expression, float, or complex
    Includes implicit zeros from sparse storage.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> from symbolica.community.tensor import Tensor
>>> tensor = Tensor.dense(A, [1.0, 2.0, 3.0, 4.0])
>>> list(tensor)
[1.0, 2.0, 3.0, 4.0]"###
    };
    ("TensorExpression.__getitem__") => {
        r###"Convert between flat component positions and coordinates.

Parameters
----------
item : int, sequence of int, or slice
    A nonnegative flat position, coordinates in logical axis order, or
    a slice of flat positions. All axis dimensions must be concrete.

Returns
-------
int, list of int, or list of lists of int
    Coordinates for a flat position, a flat position for coordinates,
    or a coordinate list for a slice. No component values are evaluated.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> A[2]
[1, 0]
>>> A[1, 0]
2"###
    };
    ("dot") => {
        r###"Contract two rank-one tensors into a scalar product.

Parameters
----------
left, right : TensorExpression, Tensor, or TensorNetwork
    Rank-one operands in compatible representation spaces.

Returns
-------
TensorExpression or TensorNetwork
    Symbolic operands return TensorExpression. If either operand carries
    component data, the result is a lazy TensorNetwork.

Notes
-----
The pairing includes the representation metric, for example Minkowski
signs. This is a bilinear product; it does not conjugate either operand.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, dot
>>> space = Representation.euc(3)
>>> p = TensorName.vector("dot_p")(space)
>>> q = TensorName.vector("dot_q")(space)
>>> dot(p, q).is_scalar
True"###
    };
    ("chain") => {
        r###"Compose an ordered product with explicitly labeled endpoints.

Parameters
----------
start_slot, end_slot : Slot
    External input and output indices, compatible with the factors'
    matrix channel. Matching endpoint labels close a trace.
*factors : TensorExpression, scalar expression, Tensor, TensorNetwork, or FactorProjector
    Ordered matrix factors, each with a unique compatible channel.
    Spectator axes remain external. Supply constituent factors separately
    when using a factor projector. Scalar factors multiply the whole chain.
    In the overloads, factor, first, and second name leading entries in
    this same sequence of positional factors.

Returns
-------
TensorExpression or TensorNetwork
    A symbolic chain when all factors are symbolic; a lazy TensorNetwork
    if any factor or projector contains component data.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> from symbolica.community.tensor import chain
>>> product = chain(space("i"), space("j"), A, A)
>>> product.rank
2"###
    };
    ("trace") => {
        r###"Close an ordered product along a representation channel.

Parameters
----------
representation : Representation
    Space over which the matrix indices are traced.
*factors : TensorExpression, scalar expression, Tensor, TensorNetwork, or FactorProjector
    Ordered factors with compatible matrix channels. Additional axes
    are retained as spectator axes. Scalar factors multiply the whole trace.
    In the overloads, factor, first, and second name leading entries in
    this same sequence of positional factors.

Returns
-------
TensorExpression or TensorNetwork
    Symbolic operands produce a symbolic trace. Any component-bearing
    factor produces a lazy TensorNetwork.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> from symbolica.community.tensor import trace
>>> trace(space, A, A).is_scalar
True"###
    };
    ("TensorExpression.__add__") => {
        r###"Add tensors with compatible external interfaces.

Parameters
----------
rhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    Other operand. Component data are retained in a lazy network.

Returns
-------
TensorExpression or TensorNetwork
    Symbolic operands return TensorExpression; component-bearing operands
    return TensorNetwork.

Notes
-----
Scalar zero acts as the additive identity. Other scalars can be added
only to rank-zero tensors.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> result = A + A"###
    };
    ("TensorExpression.__radd__") => {
        r###"Implement reflected addition.

Parameters
----------
lhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    Other operand. Component data are retained in a lazy network.

Returns
-------
TensorExpression or TensorNetwork
    Symbolic operands return TensorExpression; component-bearing operands
    return TensorNetwork.

Notes
-----
Scalar zero acts as the additive identity. Other scalars can be added
only to rank-zero tensors.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> result = A + A"###
    };
    ("TensorExpression.__sub__") => {
        r###"Subtract tensors with compatible external interfaces.

Parameters
----------
rhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    Other operand. Component data are retained in a lazy network.

Returns
-------
TensorExpression or TensorNetwork
    Symbolic operands return TensorExpression; component-bearing operands
    return TensorNetwork.

Notes
-----
Subtraction requires matching tensor axes; a nonzero scalar cannot
be subtracted from a tensor with external axes.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> result = A - A"###
    };
    ("TensorExpression.__rsub__") => {
        r###"Implement reflected subtraction.

Parameters
----------
lhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    Other operand. Component data are retained in a lazy network.

Returns
-------
TensorExpression or TensorNetwork
    Symbolic operands return TensorExpression; component-bearing operands
    return TensorNetwork.

Notes
-----
Subtraction requires matching tensor axes; a nonzero scalar cannot
be subtracted from a tensor with external axes.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> result = A - A"###
    };
    ("TensorExpression.__mul__") => {
        r###"Multiply tensors, contracting unambiguous compatible axes.

Parameters
----------
rhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    Other operand. Component data are retained in a lazy network.

Returns
-------
TensorExpression or TensorNetwork
    Symbolic operands return TensorExpression; component-bearing operands
    return TensorNetwork.

Notes
-----
Matching explicit labels contract. Compatible unresolved axes are
paired to maximize the number of contractions. Among equally complete
pairings, unresolved-unresolved pairs take precedence over unresolved-named
pairs; equally preferred alternatives raise an ambiguity error. Established
matrix channels retain their composition order. Use outer(), contract_ports(),
or compose() to make the intended pairing explicit.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> result = A * 2"###
    };
    ("TensorExpression.__rmul__") => {
        r###"Implement reflected multiplication.

Parameters
----------
lhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    Other operand. Component data are retained in a lazy network.

Returns
-------
TensorExpression or TensorNetwork
    Symbolic operands return TensorExpression; component-bearing operands
    return TensorNetwork.

Notes
-----
Matching explicit labels contract. Compatible unresolved axes are
paired to maximize the number of contractions. Among equally complete
pairings, unresolved-unresolved pairs take precedence over unresolved-named
pairs; equally preferred alternatives raise an ambiguity error. Established
matrix channels retain their composition order. Use outer(), contract_ports(),
or compose() to make the intended pairing explicit.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> result = 2 * A"###
    };
    ("TensorExpression.__truediv__") => {
        r###"Divide tensor components by a scalar expression.

Parameters
----------
rhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    Other operand. Component data are retained in a lazy network.

Returns
-------
TensorExpression or TensorNetwork
    Symbolic operands return TensorExpression; component-bearing operands
    return TensorNetwork.

Notes
-----
The denominator must be scalar. For scalar divided by tensor, the
tensor must also have rank zero.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> result = A / 2"###
    };
    ("TensorExpression.__rtruediv__") => {
        r###"Implement reflected division.

Parameters
----------
lhs : scalar expression, TensorExpression, Tensor, or TensorNetwork
    Other operand. Component data are retained in a lazy network.

Returns
-------
TensorExpression or TensorNetwork
    Symbolic operands return TensorExpression; component-bearing operands
    return TensorNetwork.

Notes
-----
The denominator must be scalar. For scalar divided by tensor, the
tensor must also have rank zero.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> result = 2 / TensorExpression(3)"###
    };
    ("TensorExpression.outer") => {
        r###"Form an outer product without implicit contractions between the operands.

Parameters
----------
rhs : TensorExpression, Tensor, or TensorNetwork
    Tensor whose axes follow this tensor's axes.

Returns
-------
TensorExpression or TensorNetwork
    Symbolic operands return TensorExpression. If rhs carries component
    data, the result is a TensorNetwork.

Notes
-----
Use this when compatible unresolved axes should remain independent.
Choose explicit distinct labels if you need to refer to them separately.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> product = A.outer(A)
>>> product.rank
4"###
    };
    ("TensorExpression.contract_ports") => {
        r###"Contract one chosen pair of axes between two tensors.

Parameters
----------
rhs : TensorExpression, Tensor, or TensorNetwork
    Tensor to contract with this tensor.
left, right : int
    Zero-based axis positions in this tensor and rhs, respectively.
    The representations must be compatible under contraction.

Returns
-------
TensorExpression or TensorNetwork
    A symbolic tensor unless rhs carries component data, in which case
    a lazy TensorNetwork is returned.

Notes
-----
Other axes remain external. Axis numbers refer to ``axes``.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> contracted = A.contract_ports(A, left=1, right=0)
>>> contracted.rank
2"###
    };
    ("TensorExpression.compose") => {
        r###"Multiply two tensors along explicitly selected matrix channels.

Parameters
----------
rhs : TensorExpression, Tensor, or TensorNetwork
    Next factor in the ordered matrix product.
left, right : tuple of int and int
    (input_axis, output_axis) in this tensor and rhs, respectively.
    The left output contracts with the right input. Other axes are
    retained as spectator axes.

Returns
-------
TensorExpression or TensorNetwork
    Symbolic operands return TensorExpression. A concrete rhs produces
    a TensorNetwork that retains its data.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, TensorExpression
>>> space = Representation.euc(2)
>>> A = TensorName("M")(space, space)
>>> product = A.compose(A, left=(0, 1), right=(0, 1))
>>> product.rank
2"###
    };
    ("CanonicalizationError") => {
        r###"Raised when dummy-index or tensor-factor canonicalization cannot be completed.

Notes
-----
Subclass of ValueError. Catch it around the corresponding tensor operation
to handle invalid input or an unsupported symbolic form.

Examples
--------
>>> from symbolica.community.tensor import CanonicalizationError
>>> issubclass(CanonicalizationError, ValueError)
True"###
    };
    ("CookingError") => {
        r###"Raised when a compound index label cannot be encoded by the selected interning mode.

Notes
-----
Subclass of TypeError. Catch it around the corresponding tensor operation
to handle invalid input or an unsupported symbolic form.

Examples
--------
>>> from symbolica.community.tensor import CookingError
>>> issubclass(CookingError, TypeError)
True"###
    };
    ("DiracAdjointError") => {
        r###"Raised when spinor index structure does not define a consistent Dirac adjoint.

Notes
-----
Subclass of ValueError. Catch it around the corresponding tensor operation
to handle invalid input or an unsupported symbolic form.

Examples
--------
>>> from symbolica.community.tensor import DiracAdjointError
>>> issubclass(DiracAdjointError, ValueError)
True"###
    };
    ("NetworkToolingError") => {
        r###"Raised when tensor syntax cannot be interpreted as a supported symbolic network.

Notes
-----
Subclass of ValueError. Catch it around the corresponding tensor operation
to handle invalid input or an unsupported symbolic form.

Examples
--------
>>> from symbolica.community.tensor import NetworkToolingError
>>> issubclass(NetworkToolingError, ValueError)
True"###
    };
    ("module") => {
        r###"Symbolic and component-based tensor manipulation.

Use Representation and TensorName to construct TensorExpression objects;
Tensor stores component data and TensorNetwork executes tensor calculations.
TensorLibrary supplies reusable components, and TensorPattern/TensorRule
provide interface-aware rewriting. Dirac matrices, color tensors, and their
simplifiers are specialized helpers built on these generic tensor operations.

AUTO (also exported as _) leaves a tensor axis unresolved during indexing.
Nc is the registered real Symbolica color-count symbol: its built-in numerical
value is 3. Use your own dimension symbol for formal SU(N) calculations
when that default numerical value is not appropriate.

Examples
--------
>>> from symbolica.community.tensor import Representation, TensorName, AUTO
>>> matrix = TensorName("A")(Representation.euc(3), Representation.euc(3))
>>> matrix(AUTO, "j").rank
2"###
    };
}
