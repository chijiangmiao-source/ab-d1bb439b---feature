// test/erasures.test.js — 错误与擦除联合定位规则测试（node --test）
// 交叉校验：Γ 的根、修正综合症 U=Γ·S 卷积、Ψ=Γ·Λ 的根覆盖、
// 2ν+s ≤ 2t 界限、翻转后综合症归零，以及非法擦除输入的首个失败原因。
import { test } from 'node:test';
import assert from 'node:assert/strict';
import {
  analyzeBch, buildGenerator, normalizeErasurePositions, polyMultiply, evalPoly,
  computeSyndromes,
} from '../src/bch.js';
import { makeField, validatePrimitivePolynomial } from '../src/gf.js';

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
function systematicEncode(field, g, k, msg) {
  const shifted = BigInt(msg) << BigInt(field.n - k);
  return (shifted ^ bMod(shifted, g)).toString(2).padStart(field.n, '0');
}
function flip(str, positions) {
  const a = str.split('');
  for (const p of positions) a[p] = a[p] === '0' ? '1' : '0';
  return a.join('');
}

test('未填写擦除：返回结构与原纯错误流程完全一致（无 erasures 字段）', () => {
  const { field, gen } = setup(4, '10011', 2);
  const c = systematicEncode(field, gen.g, gen.k, 0b1101001);
  const r = flip(c, [field.n - 1 - 9]);
  const variants = [undefined, null, '', '   ', ',， ,', []];
  for (const erasures of variants) {
    const res = analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: r, erasures });
    assert.equal(res.errorCount, 1);
    assert.equal(res.corrected, c);
    assert.equal(res.erasures, undefined, `擦除输入 ${JSON.stringify(erasures)} 不应产生联合证据`);
    assert.ok(!('kind' in res.roots[0]));
  }
});

test('擦除定位多项式 Γ：每个标记位置 i 满足 Γ(α^(−i))=0（独立求值）', () => {
  const { field } = setup(4, '10011', 2);
  const res = analyzeBch({
    m: 4, polyStr: '10011', t: 2,
    rStr: flip('0'.repeat(15), [2, 9]), // 位置 2、9 即真实错误，且全部标记为擦除
    erasures: '2, 9',
  });
  assert.equal(res.erasures.count, 2);
  assert.equal(res.erasures.unknownCount, 0);
  const gamma = res.erasures.erasureLocator.coefficients.map((x) => x.value);
  for (const strIdx of [2, 9]) {
    const power = field.n - 1 - strIdx;
    const z = field.exp[(field.n - power) % field.n]; // α^(−i)
    assert.equal(evalPoly(field, gamma, z), 0, `Γ(α^(-${power})) 应为 0`);
  }
});

test('错误 + 擦除混合：3 个损坏中 2 个标记为擦除，ν=1、s=2、用量 4=2t，纠正回原码', () => {
  const { field, gen } = setup(4, '10011', 2);
  const c = systematicEncode(field, gen.g, gen.k, 0b1011001);
  const r = flip(c, [1, 5, 11]);
  const res = analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: r, erasures: '1, 11' });

  assert.equal(res.erasures.count, 2);
  assert.equal(res.erasures.unknownCount, 1);
  assert.equal(res.errorCount, 3);
  assert.equal(res.corrected, c);
  assert.deepEqual(
    { used: res.erasures.capacity.used, limit: res.erasures.capacity.limit, ok: res.erasures.capacity.ok },
    { used: 4, limit: 4, ok: true }
  );

  // 根集合恰为三个损坏位，且来源标注正确
  const byStringIdx = new Map(res.roots.map((x) => [x.stringIndex0, x.kind]));
  assert.deepEqual([...byStringIdx.keys()].sort((a, b) => a - b), [1, 5, 11]);
  assert.equal(byStringIdx.get(1), 'erasure');
  assert.equal(byStringIdx.get(11), 'erasure');
  assert.equal(byStringIdx.get(5), 'unknown');

  // Ψ = Γ·Λ 系数级独立复算
  const gamma = res.erasures.erasureLocator.coefficients.map((x) => x.value);
  const lambda = res.locator.coefficients.map((x) => x.value);
  const psiExpect = polyMultiply(field, gamma, lambda);
  assert.deepEqual(res.erasures.jointLocator.coefficients.map((x) => x.value), psiExpect);
  assert.equal(res.erasures.jointLocator.degree, 3);
});

test('修正综合症 U_p = Σ Γ_j·S_{s+p−j}：与 Γ·S 卷积的对应系数一致', () => {
  const { field, gen } = setup(4, '10011', 2);
  const c = systematicEncode(field, gen.g, gen.k, 0b1011001);
  const r = flip(c, [1, 5, 11]);
  const res = analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: r, erasures: [1, 11] });
  const S = computeSyndromes(field, r, 4);
  const gamma = res.erasures.erasureLocator.coefficients.map((x) => x.value);
  const conv = polyMultiply(field, gamma, S); // conv[k] = Σ_j Γ_j S_{k+1−j}
  for (const u of res.erasures.modifiedSyndromes) {
    assert.equal(u.value, conv[2 + u.p - 1], `U_${u.p}`); // k = s+p−1
  }
});

test('擦除数达到上限 2t（erasures-only）：ν=0、修正综合症为空，仍须闭合', () => {
  const { field, gen } = setup(4, '10011', 1); // t=1, 2t=2
  const c = systematicEncode(field, gen.g, gen.k, 0b1101);
  const r = flip(c, [3, 12]);
  const res = analyzeBch({ m: 4, polyStr: '10011', t: 1, rStr: r, erasures: '3, 12' });
  assert.equal(res.erasures.count, 2);
  assert.equal(res.erasures.unknownCount, 0);
  assert.equal(res.erasures.modifiedSyndromes.length, 0);
  assert.equal(res.erasures.capacity.used, 2);
  assert.equal(res.corrected, c);
  assert.equal(res.errorCount, 2);
});

test('标记擦除但实际该位未损坏：证据不闭合 → 明确拒绝且不放宽规则', () => {
  const { field, gen } = setup(4, '10011', 2);
  const c = systematicEncode(field, gen.g, gen.k, 0b1011001);
  // 仅位置 5 有真实错误，却把实际正确的 1、11 标记为擦除：
  // Γ 按构造整除 Ψ，故最终裁决必须落在“根数不符 / 翻转后综合症非零”的闭合闸口上。
  const r = flip(c, [5]);
  assert.throws(
    () => analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: r, erasures: [1, 11] }),
    /无法闭合|超过界限|超过纠错能力/
  );
});

test('超出联合能力界限：s=1 标记但共 4 位损坏（2ν+s>4）被拒绝', () => {
  const { field, gen } = setup(4, '10011', 2);
  const c = systematicEncode(field, gen.g, gen.k, 0b1011001);
  const r = flip(c, [0, 4, 8, 12]);
  assert.throws(
    () => analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: r, erasures: '0' }),
    /联合能力用量|超过界限|无法闭合/
  );
});

test('擦除数超过 2t：首个失败原因即“擦除数超过能力”', () => {
  assert.throws(
    () => analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: '0'.repeat(15), erasures: '0,1,2,3,4' }),
    /擦除数超过能力：共 5 个.*2t=4/
  );
  assert.throws(
    () => normalizeErasurePositions('0,1,2,3,4', 15, 4),
    /擦除数超过能力/
  );
});

test('重复位置：报首个重复条目，且优先于其后的越界/能力检查', () => {
  assert.throws(
    () => analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: '0'.repeat(15), erasures: '3, 3' }),
    /擦除位置重复：位置 3/
  );
  assert.throws(
    () => analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: '0'.repeat(15), erasures: '2, 2, 99' }),
    /擦除位置重复：位置 2/
  );
});

test('越界位置：负值与 ≥ n 报“越界/非法”，并指出条目序号', () => {
  assert.throws(
    () => analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: '0'.repeat(15), erasures: '15' }),
    /擦除位置越界/
  );
  assert.throws(
    () => analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: '0'.repeat(15), erasures: '2, -1' }),
    /擦除位置非法/
  );
  assert.throws(
    () => normalizeErasurePositions('2, 16', 15, 4),
    /越界.*n−1=14/
  );
  assert.deepEqual(normalizeErasurePositions('0, 14', 15, 4), [0, 14]);
});

test('非法语法：字母 / 小数 / 空片段混合给出首个失败原因', () => {
  assert.throws(
    () => analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: '0'.repeat(15), erasures: 'a' }),
    /擦除位置非法/
  );
  assert.throws(
    () => analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: '0'.repeat(15), erasures: '1, 2.5' }),
    /擦除位置非法/
  );
  assert.throws(
    () => normalizeErasurePositions([1, 'x'], 15, 4),
    /擦除位置非法/
  );
});

test('边界位置 0 与 n−1 合法，与位序约定一致（最左 ↔ x^(n−1)、最右 ↔ x^0）', () => {
  const { field, gen } = setup(4, '10011', 2);
  const c = systematicEncode(field, gen.g, gen.k, 0b1011001);
  for (const [strIdx, power] of [[0, 14], [14, 0]]) {
    const r = flip(c, [strIdx]);
    const res = analyzeBch({ m: 4, polyStr: '10011', t: 2, rStr: r, erasures: [strIdx] });
    assert.equal(res.roots[0].power, power);
    assert.equal(res.roots[0].stringIndex0, strIdx);
    assert.equal(res.roots[0].kind, 'erasure');
    assert.equal(res.corrected, c);
  }
});

test('随机化：错误+擦除混合只要 2ν+s ≤ 2t 就必须纠正回原码（多参数）', () => {
  const cases = [
    [4, '10011', 3],
    [5, '100101', 3],
    [6, '1000011', 3],
    [7, '10000011', 2],
  ];
  let seed = 99173;
  const rand = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };
  for (const [m, bits, t] of cases) {
    const { field, gen } = setup(m, bits, t);
    for (let round = 0; round < 30; round++) {
      const msg = Math.floor(rand() * 2 ** Math.min(gen.k, 20));
      const c = systematicEncode(field, gen.g, gen.k, msg);
      // s ≥ 1（本套件专测擦除分支），ν ≥ 0 且 2ν+s ≤ 2t；至少 1 位损坏
      const s = 1 + Math.floor(rand() * 2 * t);
      const maxNu = Math.floor((2 * t - s) / 2);
      const nu = maxNu === 0 ? 0 : Math.floor(rand() * (maxNu + 1));
      const positions = new Set();
      while (positions.size < s + nu) positions.add(Math.floor(rand() * field.n));
      const allPos = [...positions];
      const eraPos = allPos.slice(0, s);
      const errPos = allPos.slice(s);
      const r = flip(c, allPos);
      const res = analyzeBch({ m, polyStr: bits, t, rStr: r, erasures: eraPos });
      assert.equal(res.erasures.count, s, `m=${m},round=${round} 擦除数`);
      assert.equal(res.erasures.unknownCount, errPos.length, `m=${m},round=${round} 未知错误数`);
      assert.equal(res.corrected, c, `m=${m},round=${round} 纠正码字`);
      assert.ok(res.erasures.capacity.used <= 2 * t);
      const rootIdx = new Set(res.roots.map((x) => x.stringIndex0));
      for (const p of allPos) assert.ok(rootIdx.has(p), `m=${m},round=${round} 位 ${p} 应在联合根中`);
      assert.equal(res.roots.length, res.erasures.jointLocator.degree);
    }
  }
});

test('擦除纠正后的码字必为合法码字（独立用 g 整除校验，m=5）', () => {
  const { field, gen } = setup(5, '100101', 3);
  const c = systematicEncode(field, gen.g, gen.k, 0b10111);
  const r = flip(c, [2, 9, 17, 30]); // s=2, ν=2 → 2ν+s=6=2t
  const res = analyzeBch({ m: 5, polyStr: '100101', t: 3, rStr: r, erasures: '2, 30' });
  assert.equal(bMod(BigInt('0b' + res.corrected), gen.g), 0n);
  assert.equal(res.corrected, c);
  assert.equal(res.erasures.unknownCount, 2);
});
