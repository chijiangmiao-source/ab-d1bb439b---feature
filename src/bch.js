// bch.js
// 窄义二进制本原 BCH 码（narrow-sense primitive binary BCH）离线复核。
//
// 码长 n = 2^m − 1，生成多项式 g(x) 以连续根 α, α², …, α^(2t) 构造：
//   g(x) = lcm(M_1, M_2, …, M_{2t})，M_j 为 α^j 在 GF(2) 上的最小多项式。
// 位序约定见 gf.js：接收串最左字符对应 x^(n−1)，综合症/BM/Chien 全程一致。

import {
  binDeg,
  binMul,
  binMod,
  binPowMod,
  makeField,
  validatePrimitivePolynomial,
} from './gf.js';

/** BigInt 版 GF(2)[x] 乘法（移位异或），用于可能很长的生成多项式。 */
function bMul(a, b) {
  let r = 0n;
  while (b) {
    if (b & 1n) r ^= a;
    a <<= 1n;
    b >>= 1n;
  }
  return r;
}

/** BigInt 版 GF(2)[x] 取模。 */
function bMod(a, m) {
  const dm = m.toString(2).length - 1;
  let r = a;
  let dr = r.toString(2).length - 1;
  if (r === 0n) return 0n;
  while (dr >= dm) {
    r ^= m << BigInt(dr - dm);
    dr = r === 0n ? -1 : r.toString(2).length - 1;
  }
  return r;
}

function bPowMod(base, exp, m) {
  let r = 1n;
  let b = base;
  while (exp > 0n) {
    if (exp & 1n) r = bMod(bMul(r, b), m);
    b = bMod(bMul(b, b), m);
    exp >>= 1n;
  }
  return r;
}

/** 把 GF(2) 多项式（BigInt 位向量）格式化为 x 的多项式文本。 */
export function formatBinaryPolynomial(big) {
  if (big === 0n) return '0';
  const deg = big.toString(2).length - 1;
  const terms = [];
  for (let i = deg; i >= 0; i--) {
    if ((big >> BigInt(i)) & 1n) {
      if (i === 0) terms.push('1');
      else if (i === 1) terms.push('x');
      else terms.push(`x^${i}`);
    }
  }
  return terms.join(' + ');
}

/** α 的模 n 二倍循环陪集 C_s = {s, 2s, 4s, …} (mod n)。 */
export function cyclotomicCoset(s, n) {
  const members = [];
  let v = s;
  do {
    members.push(v);
    v = (v * 2) % n;
  } while (v !== s);
  return members;
}

/**
 * 最小多项式 M_s(x) = ∏_{j∈C_s} (x − α^j)，系数必落在 GF(2)。
 * 在域上逐次乘一次因子，最后校验系数确为 0/1 并以 BigInt 返回。
 */
function minimalPolynomial(field, members) {
  let coeff = [1]; // 多项式 1，coeff[i] 为 x^i 系数（域元素）
  for (const e of members) {
    const root = field.exp[e];
    const next = new Array(coeff.length + 1).fill(0);
    for (let i = 0; i < coeff.length; i++) {
      next[i + 1] ^= coeff[i]; // x 项
      next[i] ^= field.mul(root, coeff[i]); // 常数根项
    }
    coeff = next;
  }
  let bits = 0n;
  coeff.forEach((c, i) => {
    if (c !== 0 && c !== 1) {
      throw new Error('内部错误：最小多项式系数未落在 GF(2)，陪集构造有误。');
    }
    if (c === 1) bits |= 1n << BigInt(i);
  });
  return bits;
}

/**
 * 构造窄义 BCH 生成多项式。
 * 返回 { g(bigint), degree, k, cosets:[{rep,members,minPoly}], n }。
 */
export function buildGenerator(field, m, t) {
  const n = field.n;
  if (!Number.isInteger(t) || t < 1) {
    throw new Error('纠错能力 t 非法：需为不小于 1 的整数。');
  }
  if (2 * t > n - 1) {
    throw new Error(
      `参数组合非法：设计距离 2t+1=${2 * t + 1} 超过码长 n=${n}，要求 1 ≤ t ≤ 2^${m}−2 的一半（t ≤ ${(n - 1) >> 1}）。`
    );
  }

  const byRep = new Map();
  for (let s = 1; s <= 2 * t; s++) {
    const members = cyclotomicCoset(s, n);
    const rep = Math.min(...members);
    if (!byRep.has(rep)) {
      byRep.set(rep, members);
    }
  }

  let g = 1n;
  const cosets = [];
  let totalDegree = 0;
  for (const [rep, members] of [...byRep.entries()].sort((a, b) => a[0] - b[0])) {
    const mp = minimalPolynomial(field, members);
    const deg = mp.toString(2).length - 1;
    totalDegree += deg;
    g = bMul(g, mp);
    cosets.push({
      rep,
      members: members.slice().sort((a, b) => a - b),
      minPolyBits: mp.toString(2),
      minPolyText: formatBinaryPolynomial(mp),
      degree: deg,
    });
  }

  const degree = g.toString(2).length - 1;
  if (degree !== totalDegree) {
    throw new Error('内部错误：生成多项式次数与最小多项式次数之和不一致（因子不互素）。');
  }
  if (degree > m * t) {
    throw new Error(`内部错误：生成多项式次数 ${degree} 超过 BCH 界 mt=${m * t}。`);
  }
  const k = n - degree;
  if (degree < 1 || k < 1) {
    throw new Error(
      `生成多项式次数 ${degree} 不满足有效码长：要求 1 ≤ deg(g) ≤ n−1=${n - 1}（信息位 k=n−deg(g) 必须为正）。`
    );
  }

  return { g, degree, k, cosets, n };
}

/**
 * 综合症 S_j = r(α^j)，j = 1..2t，Horner 求值。
 * 接收串最左为 x^(n−1) 系数，逐字符：acc ← acc·α^j + c_i。
 */
export function computeSyndromes(field, rStr, twoT) {
  const n = field.n;
  if (rStr.length !== n) {
    throw new Error(`接收码字长度非法：应为 n=2^m−1=${n} 位，实际 ${rStr.length} 位。`);
  }
  if (!/^[01]+$/.test(rStr)) {
    throw new Error('接收码字非法：只能包含字符 0 和 1。');
  }
  const syndromes = [];
  for (let j = 1; j <= twoT; j++) {
    const beta = field.exp[j % n]; // α^j
    let acc = 0;
    for (let p = 0; p < rStr.length; p++) {
      acc = field.mul(acc, beta);
      if (rStr[p] === '1') acc ^= 1;
    }
    syndromes.push(acc);
  }
  return syndromes;
}

/**
 * Berlekamp–Massey，输入 S_1..S_{2t}（数组下标 0 对应 S_1）。
 * 返回定位多项式系数数组 Λ[0..L]，Λ(x)=Σ Λ_i x^i，GF(2) 上减法即异或。
 */
export function berlekampMassey(field, syndromes) {
  const N = syndromes.length;
  const S = (i) => syndromes[i - 1]; // i 为 1 基
  let C = [1];
  let B = [1];
  let L = 0;
  let shift = 1;
  let b = 1;

  for (let n = 1; n <= N; n++) {
    let d = S(n);
    for (let i = 1; i <= L; i++) {
      if (C[i] !== 0 && S(n - i) !== 0) d ^= field.mul(C[i], S(n - i));
    }
    if (d === 0) {
      shift++;
      continue;
    }
    const T = C.slice();
    const coef = field.div(d, b);
    const grown = new Array(Math.max(C.length, B.length + shift)).fill(0);
    C.forEach((v, i) => { grown[i] = v; });
    B.forEach((v, i) => {
      if (v !== 0) grown[i + shift] ^= field.mul(coef, v);
    });
    C = grown;
    if (2 * L <= n - 1) {
      L = n - L;
      B = T;
      b = d;
      shift = 1;
    } else {
      shift++;
    }
  }

  C.length = L + 1;
  return { lambda: C, L };
}

/** 定位多项式文本：Λ(x) = 1 + α^u x + …。 */
export function formatLocator(field, lambda) {
  const terms = [];
  lambda.forEach((c, i) => {
    if (c === 0) return;
    if (i === 0) terms.push('1');
    else if (i === 1) terms.push(c === 1 ? 'x' : `α^${field.log[c]} x`);
    else terms.push(c === 1 ? `x^${i}` : `α^${field.log[c]} x^${i}`);
  });
  return terms.length ? terms.join(' + ') : '0';
}

/**
 * Chien 搜索：错误出现在 x^i 位当且仅当 Λ(α^(−i)) = 0，i = 0..n−1。
 * 返回命中幂次列表（幂次 = 码多项式位指数）。
 */
export function chienSearch(field, lambda, L) {
  const n = field.n;
  const hits = [];
  for (let i = 0; i < n; i++) {
    const z = field.exp[(n - i) % n]; // α^(−i)
    if (evalPoly(field, lambda, z) === 0) hits.push(i);
  }
  return hits;
}

/** Horner 求域多项式 P(z)，coeff[i] 为 x^i 系数。 */
export function evalPoly(field, coeff, z) {
  let acc = 0;
  for (let p = coeff.length - 1; p >= 0; p--) {
    acc = field.mul(acc, z);
    if (coeff[p]) acc ^= coeff[p];
  }
  return acc;
}

/** GF(2^m) 上多项式卷积乘法（char 2，系数相加为异或）。 */
export function polyMultiply(field, a, b) {
  const out = new Array(a.length + b.length - 1).fill(0);
  for (let i = 0; i < a.length; i++) {
    if (a[i] === 0) continue;
    for (let j = 0; j < b.length; j++) {
      if (b[j] !== 0) out[i + j] ^= field.mul(a[i], b[j]);
    }
  }
  return out;
}

/**
 * 解析零基擦除位置。接受字符串（逗号 / 中文逗号 / 空白分隔）或数组（数字或数字串）。
 * 位置是“接收串内 0 基下标”（自左向右，与页面根表一致），取值范围 0..n−1。
 * 按条目给出顺序逐项检查（语法 → 越界 → 重复），首个问题即抛出中文 Error；
 * 全部条目合法后再检查不同擦除数是否超过擦除定位上限 2t。
 */
export function normalizeErasurePositions(input, n, twoT) {
  const tokens = Array.isArray(input)
    ? input.map((v) => (typeof v === 'number' ? String(v) : v))
    : typeof input === 'string'
      ? input.split(/[,，\s]+/)
      : null;
  if (tokens === null) {
    throw new Error('擦除位置非法：需为逗号或空白分隔的零基非负整数列表。');
  }
  const result = [];
  const seen = new Set();
  for (const raw of tokens) {
    const tok = String(raw ?? '').trim();
    if (tok === '') continue;
    const ordinal = result.length + 1;
    if (!/^\d+$/.test(tok)) {
      throw new Error(
        `擦除位置非法：第 ${ordinal} 个条目“${tok}”不是零基非负整数（请用逗号或空格分隔，如 2, 9）。`
      );
    }
    const p = Number(tok);
    if (!Number.isSafeInteger(p) || p < 0 || p > n - 1) {
      throw new Error(
        `擦除位置越界：第 ${ordinal} 个条目“${tok}”超出码长范围，零基位置须满足 0 ≤ p ≤ n−1=${n - 1}（当前 n=${n}）。`
      );
    }
    if (seen.has(p)) {
      throw new Error(`擦除位置重复：位置 ${p} 在列表中出现多次，首个重复条目为第 ${ordinal} 个。`);
    }
    seen.add(p);
    result.push(p);
  }
  if (result.length > twoT) {
    throw new Error(
      `擦除数超过能力：共 ${result.length} 个不同擦除位置，超过该码擦除定位上限 2t=${twoT}（t=${twoT / 2}）。`
    );
  }
  return result;
}

/**
 * 完整复核入口。成功返回结构化证据；任何无法闭合的情况抛出中文 Error。
 *
 * 输入：
 *   m        域阶数（2..20）
 *   polyStr  m 次本原多项式比特串（首位 x^m，末位常数项）
 *   t        纠错能力
 *   rStr     长度恰为 2^m−1 的接收码字（最左为 x^(n−1)）
 *   erasures 可选：链路质量标记给出的零基擦除位置（接收串内自左向右下标 0..n−1），
 *            可为字符串（逗号/空白分隔）或数字数组。未提供或为空时走原有纯错误流程，
 *            返回结构与行为与不支持擦除时完全一致。
 *
 * 擦除流程（Forney 的“错误与擦除”联合定位）：
 *   - 擦除位置 i_1..i_s（码多项式幂次）的定位多项式 Γ(x)=∏(1+α^(i_k) x)；
 *   - 修正综合症 U_p = Σ_{j=0}^{s} Γ_j S_{s+p−j}（p=1..2t−s），BM 仅解未知错误 Λ；
 *   - 联合定位多项式 Ψ=Γ·Λ，Chien 搜索 Ψ 的根必须恰为 deg(Ψ)=ν+s 个并覆盖全部擦除位；
 *   - 联合能力界限 2ν+s ≤ 2t；翻转全部定位位后 2t 个综合症必须全部归零。
 */
export function analyzeBch({ m, polyStr, t, rStr, erasures }) {
  // 1) 先验证多项式确为 m 次本原多项式。
  const f = validatePrimitivePolynomial(m, polyStr);
  const field = makeField(m, f);
  const n = field.n;

  // 2) 由连续根 α..α^(2t) 的循环陪集构造二进制生成多项式。
  const generator = buildGenerator(field, m, t);

  // 3) 接收串合法性。
  if (typeof rStr !== 'string' || rStr.length !== n || !/^[01]+$/.test(rStr)) {
    if (typeof rStr !== 'string' || !/^[01]+$/.test(rStr)) {
      throw new Error('接收码字非法：只能包含字符 0 和 1。');
    }
    throw new Error(`接收码字长度非法：应为 n=2^m−1=${n} 位，实际 ${rStr.length} 位。`);
  }

  // 4) 综合症。
  const twoT = 2 * t;
  const syndromes = computeSyndromes(field, rStr, twoT);
  const syndromeView = syndromes.map((v, idx) => ({
    j: idx + 1,
    value: v,
    text: v === 0 ? '0' : `α^${field.log[v]}`,
  }));

  const allZero = syndromes.every((v) => v === 0);

  // 未填写任何擦除位置（或只含分隔符 / 空白）：沿用原有纯错误流程（结果结构保持不变）。
  const rawErasures =
    erasures !== undefined && erasures !== null &&
    ((Array.isArray(erasures) && erasures.length > 0) ||
      (typeof erasures === 'string' && erasures.trim() !== ''));
  if (!rawErasures) {
    return analyzeErrorsOnly({
      field, m, n, t, twoT, polyStr, rStr, syndromes, syndromeView, generator, allZero,
    });
  }

  // 5) 擦除位置合法性（语法 → 越界 → 重复 → 擦除数上限），首个失败原因即抛出。
  const erasureIndexes = normalizeErasurePositions(erasures, n, twoT);
  if (erasureIndexes.length === 0) {
    return analyzeErrorsOnly({
      field, m, n, t, twoT, polyStr, rStr, syndromes, syndromeView, generator, allZero,
    });
  }
  const s = erasureIndexes.length;
  const erasurePowers = erasureIndexes.map((p) => n - 1 - p);

  // 6) 擦除定位多项式 Γ(x) = ∏ (1 + α^i x)。
  let gamma = [1];
  for (const i of erasurePowers) {
    gamma = polyMultiply(field, gamma, [1, field.exp[i % n]]);
  }

  // 7) 修正综合症 U_p（p=1..2t−s）：用已知 Γ 消去擦除贡献。
  const S = (j) => syndromes[j - 1]; // j 为 1 基
  const modified = [];
  for (let p = 1; p <= twoT - s; p++) {
    let u = 0;
    for (let j = 0; j <= s; j++) {
      if (gamma[j]) u ^= field.mul(gamma[j], S(s + p - j));
    }
    modified.push(u);
  }
  const modifiedView = modified.map((v, idx) => ({
    p: idx + 1,
    value: v,
    text: v === 0 ? '0' : `α^${field.log[v]}`,
  }));

  // 8) 对修正综合症做 BM，只定位未知错误 Λ。
  const { lambda, L: nu } = berlekampMassey(field, modified);

  // 9) 联合能力界限 2ν + s ≤ 2t（不因擦除标记放宽码字规则）。
  if (2 * nu + s > twoT) {
    throw new Error(
      `不可纠正：未知错误数 ν=${nu}、擦除数 s=${s}，联合能力用量 2ν+s=${2 * nu + s} 超过界限 2t=${twoT}，` +
      `链路质量标记不能放宽该码的码字规则。`
    );
  }

  // 10) 联合定位多项式 Ψ = Γ·Λ。
  let psi = polyMultiply(field, gamma, lambda);
  while (psi.length > 1 && psi[psi.length - 1] === 0) psi = psi.slice(0, -1);
  const jointDegree = psi.length - 1;

  // 11) Chien 搜索联合根。
  const roots = chienSearch(field, psi, jointDegree);
  if (roots.length !== jointDegree) {
    throw new Error(
      `联合定位证据无法闭合：联合定位多项式次数为 ${jointDegree}（ν=${nu}, s=${s}），` +
      `但 Chien 搜索在 ${n} 个比特位置中只找到 ${roots.length} 个根，拒绝给出可纠正结论。`
    );
  }

  // 12) 根证据必须覆盖每一个被标记的擦除位置。
  const rootSet = new Set(roots);
  for (let k = 0; k < erasurePowers.length; k++) {
    const i = erasurePowers[k];
    if (!rootSet.has(i)) {
      throw new Error(
        `联合定位根证据不闭合：擦除标记位置 ${erasureIndexes[k]}（x^${i}）不是联合定位多项式的根，` +
        `标记与综合症证据冲突，拒绝按擦除纠正。`
      );
    }
  }
  const erasurePowerSet = new Set(erasurePowers);
  const unknownPowers = roots.filter((i) => !erasurePowerSet.has(i));
  if (unknownPowers.length !== nu) {
    throw new Error(
      `联合定位证据无法闭合：BM 给出未知错误数 ν=${nu}，但联合根中实际有 ${unknownPowers.length} 个未标记位置，证据不一致。`
    );
  }

  // 13) 翻转全部命中位（擦除位 + 未知错误位），综合症必须全部归零。
  const orderedRoots = roots.slice().sort((a, b) => a - b);
  const corrected = rStr.split('');
  for (const i of orderedRoots) {
    const idx = n - 1 - i;
    corrected[idx] = corrected[idx] === '1' ? '0' : '1';
  }
  const cStr = corrected.join('');
  const postSyndromes = computeSyndromes(field, cStr, twoT);
  if (!postSyndromes.every((v) => v === 0)) {
    const nonzero = postSyndromes
      .map((v, idx) => (v === 0 ? null : `S_${idx + 1}=α^${field.log[v]}`))
      .filter(Boolean)
      .join(', ');
    throw new Error(
      `联合定位证据无法闭合：按 ${s} 个擦除标记与 ${nu} 个未知错误翻转后综合症仍非零（${nonzero}），拒绝给出可纠正结论。`
    );
  }

  const coeffView = (coeff) =>
    coeff.map((v, i) => ({
      i,
      value: v,
      text: v === 0 ? '0' : v === 1 ? '1' : `α^${field.log[v]}`,
    }));
  const rootView = orderedRoots.map((i) => ({
    power: i,
    rootText: `α^${(n - i) % n}`,
    stringIndex0: n - 1 - i,
    stringIndex1: n - i,
    kind: erasurePowerSet.has(i) ? 'erasure' : 'unknown',
  }));
  const capacityUsed = 2 * nu + s;

  return {
    params: { m, n, t, primitiveBits: polyStr, k: generator.k },
    field: { order: n, generatorPolynomialBits: polyStr },
    generator: {
      bits: generator.g.toString(2),
      text: formatBinaryPolynomial(generator.g),
      degree: generator.degree,
      k: generator.k,
      cosets: generator.cosets,
    },
    syndromes: syndromeView,
    locator: {
      degree: nu,
      text: formatLocator(field, lambda),
      coefficients: coeffView(lambda),
    },
    roots: rootView,
    errorCount: orderedRoots.length,
    received: rStr,
    corrected: cStr,
    changed: false,
    erasures: {
      provided: erasureIndexes.slice(),
      count: s,
      modifiedSyndromes: modifiedView,
      erasureLocator: {
        degree: gamma.length - 1,
        text: formatLocator(field, gamma),
        coefficients: coeffView(gamma),
      },
      jointLocator: {
        degree: jointDegree,
        text: formatLocator(field, psi),
        coefficients: coeffView(psi),
      },
      unknownCount: nu,
      capacity: { used: capacityUsed, limit: twoT, ok: capacityUsed <= twoT },
      roots: rootView,
    },
    conclusion:
      `可纠正（错误与擦除联合）：${s} 个链路标记擦除与 ${nu} 个未知错误全部定位，` +
      `共翻转 ${orderedRoots.length} 位；联合能力用量 2ν+s=${capacityUsed} ≤ 2t=${twoT}，` +
      `纠正后全部 ${twoT} 个综合症归零。`,
  };
}

/** 原有纯错误复核分支：BM → Chien → 翻转后综合症归零，行为与加入擦除支持前一致。 */
function analyzeErrorsOnly({
  field, m, n, t, twoT, polyStr, rStr, syndromes, syndromeView, generator, allZero,
}) {
  // 5) Berlekamp–Massey。
  const { lambda, L } = berlekampMassey(field, syndromes);
  if (L > t) {
    throw new Error(
      `不可纠正：BM 定位多项式次数 L=${L} 超过纠错能力 t=${t}，损坏超出该码设计能力。`
    );
  }

  // 6) Chien 搜索。
  const roots = chienSearch(field, lambda, L);
  if (roots.length !== L) {
    throw new Error(
      `定位证据无法闭合：定位多项式次数为 ${L}，但 Chien 搜索在 ${n} 个比特位置中只找到 ${roots.length} 个根，不可伪装为正常数据。`
    );
  }

  // 7) 翻转命中位，重新计算全部综合症，必须全部归零。
  const corrected = rStr.split('');
  for (const i of roots) {
    const idx = n - 1 - i; // 串下标（最左为 0，对应 x^(n−1)）
    corrected[idx] = corrected[idx] === '1' ? '0' : '1';
  }
  const cStr = corrected.join('');
  const postSyndromes = computeSyndromes(field, cStr, twoT);
  if (!postSyndromes.every((v) => v === 0)) {
    const nonzero = postSyndromes
      .map((v, idx) => (v === 0 ? null : `S_${idx + 1}=α^${field.log[v]}`))
      .filter(Boolean)
      .join(', ');
    throw new Error(
      `定位证据无法闭合：按定位结果纠正后综合症仍非零（${nonzero}），拒绝给出可纠正结论。`
    );
  }

  return {
    params: { m, n, t, primitiveBits: polyStr, k: generator.k },
    field: { order: n, generatorPolynomialBits: polyStr },
    generator: {
      bits: generator.g.toString(2),
      text: formatBinaryPolynomial(generator.g),
      degree: generator.degree,
      k: generator.k,
      cosets: generator.cosets,
    },
    syndromes: syndromeView,
    locator: {
      degree: L,
      text: formatLocator(field, lambda),
      coefficients: lambda.map((v, i) => ({
        i,
        value: v,
        text: v === 0 ? '0' : v === 1 ? '1' : `α^${field.log[v]}`,
      })),
    },
    roots: roots.map((i) => ({
      power: i, // x^i
      rootText: `α^${(n - i) % n}`, // Λ 的根 α^(−i)
      stringIndex0: n - 1 - i, // 串内 0 基下标（自左向右）
      stringIndex1: n - i, // 串内 1 基位置（自左向右）
    })),
    errorCount: L,
    received: rStr,
    corrected: cStr,
    changed: allZero && L === 0,
    conclusion:
      L === 0
        ? '综合症全部为零：接收码字本身就是合法码字，无需纠正。'
        : `可纠正：定位到 ${L} 个错误比特，纠正后全部 ${twoT} 个综合症归零。`,
  };
}

// 小参数场景的整数版工具，供测试直接复核整除性（公共 API）。
export const smallIntHelpers = { binDeg, binMul, binMod, binPowMod };
