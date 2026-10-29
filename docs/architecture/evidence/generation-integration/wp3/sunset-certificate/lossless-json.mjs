// Review mutations must retain serialized u128 diagram IDs. Native JSON numbers
// above 2^53 cannot pass through JavaScript Number. Node's source-aware reviver
// and rawJSON preserve their exact decimal spelling, including in subprocesses.
const originalParse = JSON.parse;
JSON.parse = (text, reviver) => originalParse(text, function (name, value, context) {
  let exact = value;
  if (name === 'id' && typeof value === 'number' && !Number.isSafeInteger(value)) {
    if (!/^\d+$/.test(context.source)) throw new Error('noninteger native diagram ID');
    exact = BigInt(context.source);
  }
  return reviver ? reviver.call(this, name, exact, context) : exact;
});
BigInt.prototype.toJSON = function () {return JSON.rawJSON(String(this));};
