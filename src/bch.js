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

/** GF(2^m)[x] 多项式乘法（系数数组，下标为 x 的幂次，减法即异或）。 */
function fieldPolyMul(field, a, b) {
  const out = new Array(a.length + b.length - 1).fill(0);
  for (let i = 0; i < a.length; i++) {
    if (a[i] === 0) continue;
    for (let j = 0; j < b.length; j++) {
      if (b[j] !== 0) out[i + j] ^= field.mul(a[i], b[j]);
    }
  }
  return out;
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
    let acc = 0;
    for (let p = L; p >= 0; p--) {
      acc = field.mul(acc, z);
      if (lambda[p]) acc ^= lambda[p];
    }
    if (acc === 0) hits.push(i);
  }
  return hits;
}

/**
 * 解析零基擦除位置输入。
 * 接受数字数组或空白/逗号分隔的字符串（如 "2, 9 13"）；
 * 空输入（未填写擦除）归一化为空数组。任何非法记号抛出中文 Error。
 */
function parseErasures(raw) {
  if (raw == null) return [];
  if (Array.isArray(raw)) {
    return raw.map((v) => {
      if (!Number.isInteger(v)) {
        throw new Error(`擦除位置非法：${String(v)} 不是整数。`);
      }
      return v;
    });
  }
  if (typeof raw !== 'string') {
    throw new Error('擦除位置非法：需为零基整数位置的逗号/空白分隔列表。');
  }
  const tokens = raw.trim().split(/[\s,，、;；]+/).filter((s) => s.length > 0);
  const positions = [];
  for (const tok of tokens) {
    if (!/^[+-]?\d+$/.test(tok)) {
      throw new Error(`擦除位置非法：“${tok}”不是整数（零基位置，用逗号或空白分隔）。`);
    }
    positions.push(Number(tok));
  }
  return positions;
}

/**
 * 完整复核入口。成功返回结构化证据；任何无法闭合的情况抛出中文 Error。
 *
 * 输入：
 *   m        域阶数（2..20）
 *   polyStr  m 次本原多项式比特串（首位 x^m，末位常数项）
 *   t        纠错能力
 *   rStr     长度恰为 2^m−1 的接收码字（最左为 x^(n−1)）
 *   erasures 可选；零基擦除位置（串内自左向右下标），数字数组或分隔字符串。
 *            未填写（空数组/空串）时流程与无擦除完全一致。
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

  // 3.5) 擦除位置：解析为串内 0 基下标，再换算为码多项式幂次 i = n−1−p。
  const erasureStringIndexes = parseErasures(erasures);
  const seenErasure = new Set();
  const erasurePowers = [];
  for (const p of erasureStringIndexes) {
    if (p < 0 || p >= n) {
      throw new Error(
        `擦除位置越界：零基位置 ${p} 超出 [0, ${n - 1}]（码长 n=${n}）。`
      );
    }
    if (seenErasure.has(p)) {
      throw new Error(`擦除位置重复：零基位置 ${p} 出现了不止一次，联合定位要求各位置互异。`);
    }
    seenErasure.add(p);
    erasurePowers.push(n - 1 - p);
  }
  const eCount = erasurePowers.length;
  const hasErasures = eCount > 0;
  const twoT = 2 * t;
  if (eCount > twoT) {
    throw new Error(
      `擦除数超过能力：标记了 ${eCount} 个擦除，但该码联合能力界限为 2ν+e ≤ 2t=${twoT}（ν≥0），仅擦除本身就已超出，标记不能放宽原有码字规则。`
    );
  }

  // 4) 综合症。
  const syndromes = computeSyndromes(field, rStr, twoT);
  const syndromeView = syndromes.map((v, idx) => ({
    j: idx + 1,
    value: v,
    text: v === 0 ? '0' : `α^${field.log[v]}`,
  }));

  const allZero = syndromes.every((v) => v === 0);

  // 5) 擦除定位多项式 Γ(x) = ∏_{i∈擦除} (1 + X_i x)，X_i = α^i 为定位元
  //    （与综合症约定 S_j = α^(j·i) 一致；因子根为 α^(−i)，供 Chien 搜索）。
  //    无擦除时 Γ = 1，修正综合症退化为原综合症，联合流程与原实现一致。
  let gamma = [1];
  const gammaFactors = [];
  for (const i of erasurePowers) {
    const X = field.exp[i % n]; // α^i
    gamma = fieldPolyMul(field, gamma, [1, X]);
    gammaFactors.push({
      power: i, // x^i
      locator: X,
      locatorText: `α^${i % n}`,
      rootText: `α^${(n - i) % n}`, // Γ/Λ 的根
      stringIndex0: n - 1 - i,
      stringIndex1: n - i,
    });
  }

  // 6) Forney 擦除修正综合症：T_j = Σ_{r=0..e} γ_r S_{j−r}，j = e+1..2t。
  //    对任一擦除定位元 X_k 有 Σ_r γ_r X_k^{j−r} = X_k^j·Γ(X_k^(−1)) = 0，
  //    故 T 中擦除贡献被精确消去，只含未知错误位置的指数和；
  //    可用的 T 共 2t−e 个，BM 求未知错误定位多项式 σ（与 Γ 互素）。
  const modifiedSyndromes = [];
  for (let j = eCount + 1; j <= twoT; j++) {
    let tval = 0;
    for (let r = 0; r <= eCount; r++) {
      if (gamma[r] !== 0) {
        tval ^= field.mul(gamma[r], syndromes[j - r - 1]); // S_{j−r}
      }
    }
    modifiedSyndromes.push(tval);
  }
  const { lambda: sigma, L } = berlekampMassey(field, modifiedSyndromes);

  // 联合能力界限：2ν + e ≤ 2t（ν 为未知错误数，e 为擦除数）。
  // 擦除只把已知位置纳入定位，绝不放宽原有码字规则。
  if (2 * L + eCount > twoT) {
    throw new Error(
      hasErasures
        ? `不可纠正：联合能力界限 2ν+e ≤ 2t 被破坏（未知错误 ν=${L}、擦除 e=${eCount}，2ν+e=${2 * L + eCount} > 2t=${twoT}），擦除标记不能放宽码字规则。`
        : `不可纠正：BM 定位多项式次数 L=${L} 超过纠错能力 t=${t}，损坏超出该码设计能力。`
    );
  }

  // 7) 联合定位多项式 Λ(x) = Γ(x)·σ(x)；σ 只含未知位置，与 Γ 互素，根互异。
  const lambda = hasErasures ? fieldPolyMul(field, gamma, sigma) : sigma;
  const locDeg = lambda.length - 1;
  // 次数闭合校验：deg Λ 必须等于 e + L（理论上恒等，防御性核验因子首项不丢失）。
  if (locDeg !== eCount + L) {
    throw new Error(
      `定位证据无法闭合：联合定位多项式次数 ${locDeg} 不等于擦除数 ${eCount} 与未知错误定位次数 ${L} 之和。`
    );
  }

  // 8) Chien 搜索：Λ(α^(−i)) = 0，i = 0..n−1。
  const roots = chienSearch(field, lambda, locDeg);
  if (roots.length !== locDeg) {
    throw new Error(
      hasErasures
        ? `定位证据无法闭合：联合定位多项式次数为 ${locDeg}（擦除 ${eCount} + 未知 ${L}），但 Chien 搜索在 ${n} 个比特位置中只找到 ${roots.length} 个根，擦除标记与综合症证据不一致。`
        : `定位证据无法闭合：定位多项式次数为 ${L}，但 Chien 搜索在 ${n} 个比特位置中只找到 ${roots.length} 个根，不可伪装为正常数据。`
    );
  }

  // 9) 根证据闭合：每个擦除定位元必须是联合多项式的根，且联合根恰由
  //    “全部擦除位 + L 个未知错误位”组成，不多不少。
  const erasurePowerSet = new Set(erasurePowers);
  const jointPowerSet = new Set(roots);
  const missingErasurePowers = erasurePowers.filter((i) => !jointPowerSet.has(i));
  if (missingErasurePowers.length > 0) {
    const first = missingErasurePowers[0];
    throw new Error(
      `根证据无法闭合：擦除位置 x^${first}（串内下标 ${n - 1 - first}，定位元 α^${first}）不是联合定位多项式的根，该擦除标记与接收码字综合症矛盾。`
    );
  }
  const unknownPowers = roots
    .filter((i) => !erasurePowerSet.has(i))
    .sort((a, b) => a - b);
  if (unknownPowers.length !== L) {
    throw new Error(
      `根证据无法闭合：联合定位的 ${roots.length} 个根中除 ${eCount} 个擦除位外有 ${unknownPowers.length} 个额外位置，与未知错误定位次数 ${L} 不符。`
    );
  }

  // 10) 翻转全部联合命中位，重新计算全部综合症，必须全部归零。
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
      `定位证据无法闭合：按联合定位结果纠正后综合症仍非零（${nonzero}），拒绝给出可纠正结论。`
    );
  }

  const erasureLocator = {
    provided: hasErasures,
    count: eCount,
    degree: eCount,
    text: hasErasures ? formatLocator(field, gamma) : '1',
    factors: gammaFactors,
    stringIndexes0: erasureStringIndexes.slice().sort((a, b) => a - b),
    // Forney 修正综合症 T_j（j = e+1..2t）；无擦除时即原 S_1..S_2t。
    modifiedSyndromes: modifiedSyndromes.map((v, idx) => ({
      j: idx + eCount + 1,
      value: v,
      text: v === 0 ? '0' : `α^${field.log[v]}`,
    })),
  };

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
    // 未知错误定位多项式 σ（BM 输出；无擦除时即联合定位多项式）。
    locator: {
      degree: L,
      text: formatLocator(field, sigma),
      coefficients: sigma.map((v, i) => ({
        i,
        value: v,
        text: v === 0 ? '0' : v === 1 ? '1' : `α^${field.log[v]}`,
      })),
    },
    // 擦除定位多项式 Γ；联合定位多项式 Λ = Γ·σ。
    erasureLocator,
    jointLocator: {
      degree: locDeg,
      text: formatLocator(field, lambda),
      capacityUsed: 2 * L + eCount,
      capacityBound: twoT,
      unknownCount: L,
      erasureCount: eCount,
      boundText: `2ν + e = 2·${L} + ${eCount} = ${2 * L + eCount} ≤ 2t = ${twoT}`,
    },
    roots: roots
      .map((i) => ({
        power: i, // x^i
        rootText: `α^${(n - i) % n}`, // Λ 的根 α^(−i)
        stringIndex0: n - 1 - i, // 串内 0 基下标（自左向右）
        stringIndex1: n - i, // 串内 1 基位置（自左向右）
        kind: erasurePowerSet.has(i) ? 'erasure' : 'unknown',
      }))
      .sort((a, b) => a.power - b.power),
    erasureCount: eCount,
    errorCount: locDeg,
    unknownErrorCount: L,
    received: rStr,
    corrected: cStr,
    changed: allZero && locDeg === 0,
    conclusion: hasErasures
      ? locDeg === 0
        ? '综合症全部为零且未提供擦除：接收码字本身就是合法码字，无需纠正。'
        : `可纠正：${eCount} 个标记擦除与 ${L} 个未知错误共同定位（联合 ${locDeg} 位），2ν+e=${2 * L + eCount} ≤ 2t=${twoT}，翻转后全部 ${twoT} 个综合症归零。`
      : L === 0
        ? '综合症全部为零：接收码字本身就是合法码字，无需纠正。'
        : `可纠正：定位到 ${L} 个错误比特，纠正后全部 ${twoT} 个综合症归零。`,
  };
}

// 小参数场景的整数版工具，供测试直接复核整除性（公共 API）。
export const smallIntHelpers = { binDeg, binMul, binMod, binPowMod };
