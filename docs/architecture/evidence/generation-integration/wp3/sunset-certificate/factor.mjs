const gcd = (left, right) => {
  let a = left < 0n ? -left : left;
  let b = right < 0n ? -right : right;
  while (b !== 0n) [a, b] = [b, a % b];
  return a;
};

export function rational(numerator, denominator = 1n) {
  if (denominator === 0n) throw new Error('zero rational denominator');
  let n = BigInt(numerator);
  let d = BigInt(denominator);
  if (d < 0n) [n, d] = [-n, -d];
  const factor = gcd(n, d);
  return [n / factor, d / factor];
}

export const add = ([a, b], [c, d]) => rational(a * d + c * b, b * d);
export const multiply = ([a, b], [c, d]) => rational(a * c, b * d);
export const divide = ([a, b], [c, d]) => rational(a * d, b * c);
export const format = ([n, d]) => d === 1n ? String(n) : `${n}/${d}`;
export const same = (left, right) => left[0] === right[0] && left[1] === right[1];

function tokenize(source) {
  const tokens = [];
  let index = 0;
  while (index < source.length) {
    if (/\s/.test(source[index])) {
      index++;
      continue;
    }
    const match = /^(\d+|[A-Za-z_][A-Za-z0-9_]*(?:::(?:\{\}::)?[A-Za-z_][A-Za-z0-9_]*)*|[()+\-*/^,])/.exec(source.slice(index));
    if (!match) throw new Error(`unsupported scalar token at ${index}: ${source.slice(index)}`);
    tokens.push(match[1]);
    index += match[1].length;
  }
  return tokens;
}

export function parse(source) {
  const tokens = tokenize(source);
  let index = 0;
  const peek = () => tokens[index];
  const take = (expected) => {
    const actual = tokens[index++];
    if (actual !== expected) throw new Error(`expected ${expected}, got ${actual}`);
  };
  function primary() {
    const token = peek();
    if (token === '(') {
      index++;
      const nested = expression();
      take(')');
      return nested;
    }
    if (token && /^\d+$/.test(token)) {
      index++;
      return {kind: 'number', value: BigInt(token)};
    }
    if (token && /^[A-Za-z_]/.test(token)) {
      index++;
      take('(');
      const args = [];
      if (peek() !== ')') {
        do {
          args.push(expression());
          if (peek() !== ',') break;
          index++;
        } while (true);
      }
      take(')');
      return {kind: 'call', name: token.split('::').at(-1), namespace: token, args};
    }
    throw new Error(`expected complete scalar atom, got ${token}`);
  }
  function power() {
    const base = primary();
    if (peek() !== '^') return base;
    index++;
    const parenthesized = peek() === '(';
    if (parenthesized) index++;
    let sign = 1n;
    if (peek() === '-' || peek() === '+') {
      if (tokens[index++] === '-') sign = -1n;
    }
    const token = tokens[index++];
    if (!token || !/^\d+$/.test(token)) throw new Error('signed integer exponent required');
    if (parenthesized) take(')');
    if (peek() === '^') throw new Error('chained scalar powers are unsupported');
    return {kind: 'binary', operator: '^', left: base,
      right: {kind: 'number', value: sign * BigInt(token)}};
  }
  function unary() {
    if (peek() === '+') throw new Error('unary plus outside exponent is unsupported');
    if (peek() === '-') {
      index++;
      return {kind: 'unary', operator: '-', value: unary()};
    }
    return power();
  }
  function product() {
    let left = unary();
    while (peek() === '*' || peek() === '/') {
      const operator = tokens[index++];
      left = {kind: 'binary', operator, left, right: unary()};
    }
    return left;
  }
  function expression() {
    let left = product();
    while (peek() === '+' || peek() === '-') {
      const operator = tokens[index++];
      left = {kind: 'binary', operator, left, right: product()};
    }
    return left;
  }
  const ast = expression();
  if (index !== tokens.length) throw new Error(`unconsumed scalar suffix: ${tokens.slice(index).join('')}`);
  return ast;
}

export function evaluate(ast) {
  const heads = [];
  const groups = [];
  function value(node) {
    if (node.kind === 'number') return rational(node.value);
    if (node.kind === 'unary') {
      const result = value(node.value);
      return node.operator === '-' ? rational(-result[0], result[1]) : result;
    }
    if (node.kind === 'binary') {
      const left = value(node.left);
      const right = value(node.right);
      switch (node.operator) {
        case '+': return add(left, right);
        case '-': return add(left, rational(-right[0], right[1]));
        case '*': return multiply(left, right);
        case '/': return divide(left, right);
        case '^': {
          if (right[1] !== 1n || ![-1n, 0n, 1n].includes(right[0])) {
            throw new Error('unsupported scalar power');
          }
          if (left[0] === 0n && right[0] <= 0n) {
            throw new Error('undefined or indeterminate zero power');
          }
          return right[0] === -1n ? divide(rational(1n), left) :
            right[0] === 0n ? rational(1n) : left;
        }
        default: throw new Error(`unsupported scalar operator ${node.operator}`);
      }
    }
    if (node.kind !== 'call') throw new Error('invalid scalar AST');
    const allowed = node.name === 'Group' ? ['Group'] :
      [node.name, `feynkit_generator_factor::${node.name}`,
        `feynkit_generator_factor::{}::${node.name}`];
    if (!allowed.includes(node.namespace)) {
      throw new Error(`unsupported qualified scalar head ${node.namespace}`);
    }
    const values = node.args.map(value);
    if (['AutG', 'InternalFermionLoopSign', 'ExternalFermionOrderingSign',
      'AntiFermionSpinSumSign', 'NumeratorIndependentSymmetryGrouping'].includes(node.name)) {
      if (values.length !== 1 || values[0][1] !== 1n ||
          node.name === 'AutG' && values[0][0] <= 0n ||
          node.name === 'NumeratorIndependentSymmetryGrouping' && values[0][0] <= 0n ||
          ['InternalFermionLoopSign', 'ExternalFermionOrderingSign',
            'AntiFermionSpinSumSign'].includes(node.name) &&
            ![-1n, 1n].includes(values[0][0])) {
        throw new Error(`invalid ${node.name} argument`);
      }
      heads.push({name: node.name, argument: String(values[0][0])});
      return values[0];
    }
    if (['NumeratorDependentGrouping', 'Group'].includes(node.name)) {
      if (values.length !== 3 || values[0][1] !== 1n || values[0][0] < 0n ||
          values[0][0] > BigInt(Number.MAX_SAFE_INTEGER)) {
        throw new Error(`invalid ${node.name} arguments`);
      }
      groups.push({name: node.name, sourceId: Number(values[0][0]),
        ratio: format(values[1]), memberWeight: format(values[2])});
      return multiply(values[1], values[2]);
    }
    throw new Error(`unsupported scalar head ${node.namespace}`);
  }
  return {value: value(ast), heads, groups};
}

export function complete(source) {
  return evaluate(parse(source));
}

export function shape(ast) {
  if (ast.kind === 'number') return ['number', String(ast.value)];
  if (ast.kind === 'unary') return ['unary', ast.operator, shape(ast.value)];
  if (ast.kind === 'binary') {
    if (['+', '*'].includes(ast.operator)) {
      const terms = [];
      const collect = (node) => {
        if (node.kind === 'binary' && node.operator === ast.operator) {
          collect(node.left);
          collect(node.right);
        } else terms.push(shape(node));
      };
      collect(ast);
      return ['commutative', ast.operator, ...terms.sort((a, b) =>
        JSON.stringify(a).localeCompare(JSON.stringify(b)))];
    }
    return ['binary', ast.operator, shape(ast.left), shape(ast.right)];
  }

  if (ast.kind === 'call') return ['call', ast.name, ast.args.map(shape)];
  throw new Error('invalid scalar AST shape');
}
