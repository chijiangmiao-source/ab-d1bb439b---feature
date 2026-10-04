// test/erasures.test.js — 擦除（已知不可靠比特）与未知错误联合纠正规则测试。
// 交叉校验要点：
//   - 未填写擦除时结果与现有流程严格一致；
//   - Γ(x)=∏(1+α^(−i)x)，Forney 修正综合症消去擦除贡献，σ 只含未知位置；
//   - 2ν+e ≤ 2t 界限、根证据闭合、翻转后综合症全零三条缺一不可；
//   - 重复/越界/非法/超能力/标记与综合症矛盾时必须拒绝并给出首个原因。
import { test } from 'node:test';
import assert from 'node:assert/strict';
import { analyzeBch, buildGenerator } from '../src/bch.js';
import { makeField, validatePrimitivePolynomial } from '../src/gf.js';

function bMul(a, b) {
  let r = 0n;
  while (b) { if (b & 1n) r ^= a; a <<= 1n; b >>= 1n; }
  return r;
}
function bMod(a, m) {
  const dm = m.toString(2).length - 1;
  let r = a;
  let dr = r === 0n ? -1 : r.toString(2).length - 1;
  while (dr >= dm) { r ^= m << BigInt(dr - dm); dr = r === 0n ? -1 : r.toString(2).length - 1; }
  return r;
}

function setup(m, bits, t) {
  const field = makeField(m, validatePrimitivePolynomial(m, bits));
  const gen = buildGenerator(field, m, t);
  return { field, gen };
}

/** 系统编码，返回 n 位串（最左 x^(n−1)）。 */
function encode(field, g, k, msg) {
  const shifted = BigInt(msg) << BigInt(field.n - k);
  return (shifted ^ bMod(shifted, g)).toString(2).padStart(field.n, '0');
}
function flip(str, positions) {
  const a = str.split('');
  for (const p of positions) a[p] = a[p] === '0' ? '1' : '0';
  return a.join('');
}

test('未填写擦除（省略 / 空数组 / 空串 / 纯空白）结果与现有流程完全一致', () => {
  const { field, gen } = setup(4, '10011', 2);
  const c = encode(field, gen.g, gen.k, 0b1101001);
  const received = flip(c, [3, 10]);
  const base = analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: received });
  for (const erasures of [undefined, [], '', '   ', '，，']) {
    const r = analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: received, erasures });
    assert.deepEqual(r, base, `擦除输入 ${JSON.stringify(erasures)} 改变了无擦除结果`);
    assert.equal(r.erasureLocator.provided, false);
    assert.equal(r.erasureLocator.text, '1');
    assert.equal(r.jointLocator.erasureCount, 0);
    assert.equal(r.jointLocator.capacityUsed, 2 * r.jointLocator.unknownCount);
    assert.ok(r.roots.every((x) => x.kind === 'unknown'));
  }
});

test('真实错误被标记为擦除：t=1 单错，e=1、ν=0 即可纠正，σ=1、Λ=Γ', () => {
  const { field, gen } = setup(4, '10011', 1);
  const c = encode(field, gen.g, gen.k, 0b10101100010);
  const p = 6; // 串内 0 基位置
  const received = flip(c, [p]);
  const r = analyzeBch({ m: 4, polyStr: '10011', t: 1, rStr: received, erasures: [p] });
  assert.equal(r.erasureLocator.count, 1);
  assert.equal(r.unknownErrorCount, 0);
  assert.equal(r.locator.degree, 0);
  assert.equal(r.locator.text, '1');
  assert.equal(r.jointLocator.degree, 1);
  assert.equal(r.jointLocator.capacityUsed, 1);
  assert.equal(r.jointLocator.capacityBound, 2);
  assert.equal(r.roots.length, 1);
  assert.equal(r.roots[0].kind, 'erasure');
  assert.equal(r.roots[0].stringIndex0, p);
  assert.equal(r.corrected, c);
});

test('一擦除一未知（t=2）：合法分隔形式下联合定位闭合', () => {
  const { field, gen } = setup(4, '10011', 2);
  const c = encode(field, gen.g, gen.k, 0b1101001);
  const pe = 2;
  const pu = 9;
  const received = flip(c, [pe, pu]);
  // 字符串形式：空白、逗号、中文逗号、顿号都应接受；只标记擦除位 pe
  const r = analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: received, erasures: ` ,，、 ${pe} 、，, ` });
  assert.equal(r.erasureLocator.count, 1);
  assert.equal(r.unknownErrorCount, 1);
  assert.equal(r.jointLocator.degree, 2);
  assert.equal(r.jointLocator.boundText, '2ν + e = 2·1 + 1 = 3 ≤ 2t = 4');
  assert.equal(r.corrected, c);
  const kinds = Object.fromEntries(r.roots.map((x) => [x.stringIndex0, x.kind]));
  assert.deepEqual(kinds, { [pe]: 'erasure', [pu]: 'unknown' });
  // Γ 因子与擦除定位元证据
  assert.equal(r.erasureLocator.factors.length, 1);
  assert.equal(r.erasureLocator.factors[0].stringIndex0, pe);
  // 修正综合症长度为 2t−e = 3，编号 j = e+1..2t
  assert.deepEqual(r.erasureLocator.modifiedSyndromes.map((s) => s.j), [2, 3, 4]);
});

test('两未知两擦除（t=3）：2ν+e=6=2t 恰好闭合；全部标记（ν=0,e=4）也可纠正', () => {
  const { field, gen } = setup(4, '10011', 3);
  const c = encode(field, gen.g, gen.k, 0b10111);

  // ν=2, e=2
  let received = flip(c, [1, 5, 8, 13]);
  let r = analyzeBch({
    m: 4, polyStr: '10011', t: 3, rStr: received, erasures: [1, 13],
  });
  assert.equal(r.unknownErrorCount, 2);
  assert.equal(r.erasureLocator.count, 2);
  assert.equal(r.jointLocator.capacityUsed, 6);
  assert.equal(r.corrected, c);
  assert.equal(r.roots.filter((x) => x.kind === 'erasure').length, 2);
  assert.equal(r.roots.filter((x) => x.kind === 'unknown').length, 2);

  // 四个真实错误全部标记：ν=0, e=4 ≤ 2t=6
  r = analyzeBch({
    m: 4, polyStr: '10011', t: 3, rStr: received, erasures: [1, 5, 8, 13],
  });
  assert.equal(r.unknownErrorCount, 0);
  assert.equal(r.erasureLocator.count, 4);
  assert.equal(r.jointLocator.degree, 4);
  assert.equal(r.corrected, c);
});

test('擦除把原本超能力的损坏拉回可纠正：t=2 下 3 个真实错误标记其中 2 个', () => {
  const { field, gen } = setup(4, '10011', 2);
  const c = encode(field, gen.g, gen.k, 0b1011001);
  const received = flip(c, [1, 5, 11]); // 无擦除时该三错模式必须拒绝
  assert.throws(
    () => analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: received }),
    /超过纠错能力|无法闭合/
  );
  const r = analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: received, erasures: [1, 11] });
  assert.equal(r.unknownErrorCount, 1);
  assert.equal(r.erasureLocator.count, 2);
  assert.equal(r.jointLocator.capacityUsed, 4);
  assert.equal(r.corrected, c);
});

test('非法擦除输入：非整数记号给出首个非法原因', () => {
  const opts = { m: 4, polyStr: '10011', t: 2, rStr: '0'.repeat(15) };
  assert.throws(() => analyzeBch({ ...opts, erasures: '2, abc' }), /擦除位置非法/);
  assert.throws(() => analyzeBch({ ...opts, erasures: '1.5' }), /擦除位置非法/);
  assert.throws(() => analyzeBch({ ...opts, erasures: [3, {}] }), /擦除位置非法/);
});

test('越界擦除：负位置与 ≥n 位置都被拒绝，并报告首个越界值', () => {
  const opts = { m: 4, polyStr: '10011', t: 2, rStr: '0'.repeat(15) };
  assert.throws(() => analyzeBch({ ...opts, erasures: [-1] }), /擦除位置越界：零基位置 -1/);
  assert.throws(() => analyzeBch({ ...opts, erasures: [15] }), /擦除位置越界：零基位置 15/);
  assert.throws(() => analyzeBch({ ...opts, erasures: [3, 100] }), /零基位置 100/);
});

test('重复擦除位置被拒绝', () => {
  const opts = { m: 4, polyStr: '10011', t: 2, rStr: '0'.repeat(15) };
  assert.throws(() => analyzeBch({ ...opts, erasures: [4, 4] }), /擦除位置重复：零基位置 4/);
  assert.throws(() => analyzeBch({ ...opts, erasures: '4, 0, 4' }), /擦除位置重复：零基位置 4/);
});

test('擦除数超过 2t：在任何定位证据构造前拒绝', () => {
  const opts = { m: 4, polyStr: '10011', t: 2, rStr: '0'.repeat(15) };
  assert.throws(
    () => analyzeBch({ ...opts, erasures: [0, 1, 2, 3, 4] }),
    /擦除数超过能力.*2t=4/
  );
});

test('联合能力界限被破坏：ν=1、e=3 时 2ν+e=5 > 2t=4，拒绝', () => {
  const { field, gen } = setup(4, '10011', 2);
  const c = encode(field, gen.g, gen.k, 0b1011001);
  // 4 个真实损坏，只标记 3 个 → 未知 ν=1 但界限 5>4
  const received = flip(c, [0, 4, 8, 12]);
  assert.throws(
    () => analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: received, erasures: [0, 4, 8] }),
    /联合能力界限.*2ν\+e/
  );
});

test('标记与综合症矛盾：合法码字上标记干净位，翻转后综合症非零必须拒绝', () => {
  const { field, gen } = setup(4, '10011', 2);
  const c = encode(field, gen.g, gen.k, 0b1011001);
  // 码字本身合法，却宣称下标 7 不可靠；Γ 根虽由标记直接给出，但该位无需翻转，
  // 翻转会注入新错误 → 纠正后综合症非零，绝不允许通过。
  assert.throws(
    () => analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: c, erasures: [7] }),
    /无法闭合/
  );
});

test('标记与综合症矛盾：真实单错 + 干净位标记，联合纠正后综合症非零必须拒绝', () => {
  const { field, gen } = setup(4, '10011', 2);
  const c = encode(field, gen.g, gen.k, 0b1011001);
  const received = flip(c, [3]); // 真实错误在下标 3
  assert.throws(
    () => analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: received, erasures: [10] }),
    /无法闭合/
  );
});

test('全部标记干净位（无真实错误）：Γ 给满根、σ=1、界限满足，翻转注入错误 → 拒绝', () => {
  // 闭合防线的直接刻画：根与能力界限都拦不住标记错误，
  // 唯一能拦截的就是翻转后综合症复核，证明没有因标记而放宽码字规则。
  const { field, gen } = setup(4, '10011', 3);
  const c = encode(field, gen.g, gen.k, 0b10111);
  assert.equal(c.length, 15);
  assert.throws(
    () => analyzeBch({ m: 4, polyStr: '10011', t: 3, rStr: c, erasures: [0, 7, 14] }),
    /纠正后综合症仍非零/
  );
});

test('随机联合场景：不同 ν/e 组合下纠正码字与发送码字一致', () => {
  const cases = [
    [4, '10011', 2],
    [4, '10011', 3],
    [5, '100101', 2],
    [6, '1000011', 3],
  ];
  let seed = 99173;
  const rand = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };
  for (const [m, bits, t] of cases) {
    const { field, gen } = setup(m, bits, t);
    for (let round = 0; round < 20; round++) {
      const msg = Math.floor(rand() * (1 << gen.k));
      const c = encode(field, gen.g, gen.k, msg);
      // 在联合能力界限内随机选 ν、e：2ν+e ≤ 2t，ν≥0, e≥0，且总数 ≥1
      const e = 1 + Math.floor(rand() * Math.min(2 * t, 4));
      const maxNu = Math.floor((2 * t - e) / 2);
      if (maxNu < 0) continue;
      const nu = Math.floor(rand() * (maxNu + 1));
      const total = nu + e;
      if (total === 0 || total > field.n) continue;
      const posSet = new Set();
      while (posSet.size < total) posSet.add(Math.floor(rand() * field.n));
      const positions = [...posSet];
      const erasures = positions.slice(0, e);
      const received = flip(c, positions);
      const r = analyzeBch({ m, polyStr: bits, t, rStr: received, erasures });
      assert.equal(r.unknownErrorCount, nu, `m=${m},t=${t},round=${round} ν`);
      assert.equal(r.erasureLocator.count, e, `m=${m},t=${t},round=${round} e`);
      assert.equal(r.corrected, c, `m=${m},t=${t},round=${round} 纠正码字`);
      assert.ok(2 * r.unknownErrorCount + r.erasureLocator.count <= 2 * t);
    }
  }
});
