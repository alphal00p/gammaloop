import './lossless-json.mjs';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import {spawnSync} from 'node:child_process';
import {fileURLToPath} from 'node:url';

const [historical, current, output] = process.argv.slice(2);
if (!historical || !current || !output) {
  throw new Error('usage: node mutations.mjs <historical-capture> <current-capture> <output.json>');
}
const scripts = path.dirname(fileURLToPath(import.meta.url));
const read = (root, name) => fs.readFileSync(path.join(root, name), 'utf8');
const write = (root, name, value) => fs.writeFileSync(path.join(root, name), value);
const jsonMutation = (name, change) => (root) => {
  const value = JSON.parse(read(root, name));
  change(value);
  write(root, name, JSON.stringify(value, null, 2) + '\n');
};
const textMutation = (name, before, after) => (root) => {
  const value = read(root, name);
  if (!value.includes(before)) throw new Error(`${name}: mutation anchor absent`);
  write(root, name, value.replace(before, after));
};
const mirroredGhostScalar = (root) => {
  const old = 'InternalFermionLoopSign(-1)*ExternalFermionOrderingSign(1)/AutG(2)';
  const changed = `(${old})+7`;
  for (const label of ['all-raw', 'all-grouped']) {
    const witness = `${label}-graphs.tsv`;
    const lines = read(root, witness).trimEnd().split('\n');
    const graph = lines.findIndex((line) => line.startsWith('1\tGL1\t') ||
      line.startsWith('2\tGL1\t'));
    if (graph < 0) throw new Error(`missing mirrored ghost ${label}`);
    const fields = lines[graph].split('\t');
    fields[2] = changed;
    fields[3] = '13/2';
    lines[graph] = fields.join('\t');
    lines[0] = label === 'all-raw' ? '6 | 125/24 = 125/24' :
      lines[0].replace('-19/24+', '149/24+').replace('= -43/24', '= 125/24');
    write(root, witness, lines.join('\n') + '\n');
    const dot = `${label}-dot/processes/amplitudes/empty_empty/default/GL1.dot`;
    const dotSource = read(root, dot);
    write(root, dot, dotSource.replace('InternalFermionLoopSign(-1)";',
      'InternalFermionLoopSign(-1)+7";').replace('overall_factor_evaluated = "-1/2";',
      'overall_factor_evaluated = "13/2";'));
  }
  for (const name of ['all-grouped-native-source.json', 'all-grouped-native-final.json']) {
    const value = JSON.parse(read(root, name));
    value.find((diagram) => diagram.name === 'GL1').overall_factor += '+7';
    write(root, name, JSON.stringify(value, null, 2) + '\n');
  }
  const groups = JSON.parse(read(root, 'all-grouped-native-groups.json'));
  groups.source.find((diagram) => diagram.name === 'GL1').factor = changed;
  groups.groups.find((group) => group.members[0].source_name === 'GL1').members[0].overall_factor = changed;
  write(root, 'all-grouped-native-groups.json', JSON.stringify(groups, null, 2) + '\n');
};
const mirroredLightCoupling = (root) => {
  for (const label of ['all-raw', 'all-grouped', 'nonfactorizable-raw',
    'nonfactorizable-grouped']) {
    for (const graph of ['GL2', 'GL4']) {
      const dot = `${label}-dot/processes/amplitudes/empty_empty/default/${graph}.dot`;
      if (fs.existsSync(path.join(root, dot))) {
        const source = read(root, dot);
        if (!source.includes('UFO::GC_11')) throw new Error(`${dot}: coupling absent`);
        write(root, dot, source.replaceAll('UFO::GC_11', 'UFO::GC_12'));
      }
    }
    if (label.endsWith('grouped')) {
      for (const suffix of ['native-source.json', 'native-final.json']) {
        const name = `${label}-${suffix}`;
        const diagrams = JSON.parse(read(root, name));
        for (const diagram of diagrams.filter((item) => ['GL2', 'GL4'].includes(item.name))) {
          diagram.numerator = diagram.numerator.replaceAll('GC_11', 'GC_12');
          for (const vertex of diagram.vertices) {
            vertex.numerator = vertex.numerator.replaceAll('GC_11', 'GC_12');
          }
        }
        write(root, name, JSON.stringify(diagrams, null, 2) + '\n');
      }
    }
  }
};
const mirroredLightRouting = (root) => {
  for (const label of ['all-raw', 'all-grouped', 'nonfactorizable-raw',
    'nonfactorizable-grouped']) {
    const dot = `${label}-dot/processes/amplitudes/empty_empty/default/GL4.dot`;
    write(root, dot, read(root, dot).replace('lmb_rep="-1*K(1,a___)+K(0,a___)"',
      'lmb_rep="K(1,a___)+K(0,a___)"'));
    if (label.endsWith('grouped')) {
      for (const suffix of ['native-source.json', 'native-final.json']) {
        const name = `${label}-${suffix}`;
        const diagrams = JSON.parse(read(root, name));
        diagrams.find((diagram) => diagram.name === 'GL4')
          .loop_momentum_basis.edge_signatures['2'].loops[1] = 1;
        write(root, name, JSON.stringify(diagrams, null, 2) + '\n');
      }
    }
  }
};
const cases = [
  ['u effective mass', 'model.mjs', jsonMutation('post-import-model.json', (data) => {
    data.particles.find((particle) => particle.name === 'u').mass = 'MT';
  }), 'mass differs'],
  ['closed ghost sign', 'join.mjs', textMutation('all-raw-graphs.tsv',
    'InternalFermionLoopSign(-1)', 'InternalFermionLoopSign(1)'), 'complete independent weight factor differs'],
  ['ghost circuit direction', 'join.mjs', textMutation(
    'all-raw-dot/processes/amplitudes/empty_empty/default/GL1.dot',
    'id=2  lmb_rep', 'id=2 dir=back  lmb_rep'), 'directed Grassmann circuit is not closed'],
  ['unexpected vacuum cut', 'join.mjs', jsonMutation('all-grouped-native-source.json', (data) => {
    data[2].cuts.push({cut: [0]});
  }), 'source 2 does not bind to the raw physical graph'],
  ['local coupling coefficient', 'join.mjs', jsonMutation('all-grouped-native-source.json', (data) => {
    data[2].vertices[0].numerator = data[2].vertices[0].numerator.replace('GC_11', 'GC_12');
  }), 'aggregate differs from complete local product'],
  ['declared spin dimension', 'join.mjs', jsonMutation('all-grouped-native-source.json', (data) => {
    data[2].vertices[0].num_interface[0] = data[2].vertices[0].num_interface[0].replace('bis(4', 'bis(5');
  }), 'CLI payload differs from native owner'],
  ['signed routing coefficient', 'frame.mjs', jsonMutation('all-grouped-native-source.json', (data) => {
    data[2].loop_momentum_basis.edge_signatures['2'].loops[1] = 1;
  }), 'signed affine routing differs'],
  ['selected master', 'join.mjs', jsonMutation('all-grouped-native-groups.json', (data) => {
    data.groups[3].master = 4;
  }), 'final master indices are not a bijection'],
  ['member ratio', 'join.mjs', jsonMutation('all-grouped-native-groups.json', (data) => {
    data.groups[3].members[1].ratio = '-1';
  }), 'scalar head/argument/occurrence contract differs'],
  ['final aggregate payload', 'join.mjs', jsonMutation('all-grouped-native-final.json', (data) => {
    data[3].numerator += '*2';
  }), 'final tensor, routing or frame differs'],
  ['final local payload', 'join.mjs', textMutation(
    'all-grouped-dot/processes/amplitudes/empty_empty/default/GL4.dot',
    'UFO::GC_11', 'UFO::GC_12'), 'CLI payload differs from native owner'],
  ['final rational weight', 'join.mjs', textMutation('all-grouped-graphs.tsv',
    '\t-1\tSome(GroupId(3))', '\t-2\tSome(GroupId(3))'), 'complete final scalar factor differs'],
  ...[
    ['trailing doubled scalar', '*2', 'complete independent weight factor differs'],
    ['trailing additive scalar', '+7', 'complete independent weight factor differs'],
    ['trailing unknown head', '*Bogus(x)', 'expected ('],
    ['squared scalar', '^2', 'unsupported scalar power'],
  ].map(([name, suffix, expected]) => [name, 'join.mjs',
    textMutation('all-raw-graphs.tsv',
      'InternalFermionLoopSign(-1)*ExternalFermionOrderingSign(1)/AutG(2)',
      `(InternalFermionLoopSign(-1)*ExternalFermionOrderingSign(1)/AutG(2))${suffix}`),
    expected]),
  ['inverted automorphism operator', 'join.mjs', textMutation('all-raw-graphs.tsv',
    'ExternalFermionOrderingSign(1)/AutG(2)',
    'ExternalFermionOrderingSign(1)*AutG(2)'), 'complete independent weight factor differs'],
  ['DOT scalar tail', 'join.mjs', textMutation(
    'all-raw-dot/processes/amplitudes/empty_empty/default/GL1.dot',
    'InternalFermionLoopSign(-1)";', 'InternalFermionLoopSign(-1)*2";'),
  'scalar head/argument/occurrence contract differs'],
  ['summary scalar tail', 'join.mjs', textMutation('all-grouped-graphs.tsv',
    ' = -43/24', '+7 = -43/24'), 'independently summed total differs'],
  ['mirrored ghost scalar and totals', 'join.mjs', mirroredGhostScalar,
    'scalar head/argument/occurrence contract differs'],
  ['mirrored light coupling and final tensor', 'join.mjs', mirroredLightCoupling,
    'complete model-derived vertex coefficient differs'],
  ['mirrored signed routing', 'join.mjs', mirroredLightRouting,
    'complete signed vacuum momentum conservation differs'],
];
const results = [];
for (const [name, script, mutate, expected] of cases) {
  const scratch = fs.mkdtempSync(path.join(os.tmpdir(), 'sunset-mutation-'));
  try {
    const cloned = path.join(scratch, 'current');
    fs.cpSync(current, cloned, {recursive: true});
    mutate(cloned);
    const run = spawnSync(process.execPath,
      [path.join(scripts, script), ...(script === 'model.mjs' ? [cloned] : [historical, cloned]),
        path.join(scratch, 'result.json')], {encoding: 'utf8'});
    const detail = `${run.stdout}\n${run.stderr}`;
    if (run.status === 0 || !detail.includes(expected)) {
      throw new Error(`${name}: rejection was absent or wrong: ${detail.slice(0, 800)}`);
    }
    results.push({name, validator: script, expected, exit: run.status});
  } finally {
    fs.rmSync(scratch, {recursive: true, force: true});
  }
}
fs.writeFileSync(output, JSON.stringify({shaGuards: false, results}, null, 2) + '\n');
console.log(`${results.length} semantic mutations rejected without SHA guards`);
