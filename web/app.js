// app.js — 页面交互：把输入交给 Worker，渲染证据；失败时稳定清除旧成功证据。

const $ = (id) => document.getElementById(id);
const inputs = ['m', 'poly', 't', 'r', 'erasure'];

let worker = null;
let busy = false;

function createWorker() {
  if (worker) worker.terminate();
  worker = new Worker('./worker.js', { type: 'module' });
  worker.onmessage = (e) => {
    busy = false;
    setButtons();
    const msg = e.data || {};
    if (msg.type === 'result') renderSuccess(msg.result);
    else renderFailure(msg.message || '未知错误。');
  };
  worker.onerror = (e) => {
    busy = false;
    setButtons();
    renderFailure(`Worker 执行异常：${e.message || '未知错误'}`);
  };
}

function setButtons() {
  $('btn-run').disabled = busy;
  $('btn-run').textContent = busy ? '复核中…' : '复核';
}

function setStatus(text, cls) {
  const el = $('status');
  el.textContent = text;
  el.className = cls || '';
}

/** 清除一切旧结论（成功证据与失败提示同时清空，含擦除联合证据）。 */
function clearConclusions() {
  for (const id of ['result', 'result2', 'result3', 'result-era', 'result4', 'result5']) {
    $(id).classList.remove('show');
  }
  $('fail').classList.remove('show');
  setStatus('', '');
}

function renderFailure(reason) {
  clearConclusions();
  $('fail-reason').textContent = reason;
  $('fail').classList.add('show');
  setStatus('复核失败：' + reason, 'err');
}

function el(tag, text) {
  const e = document.createElement(tag);
  if (text !== undefined) e.textContent = text;
  return e;
}

/** 系数表：两列（幂次 i, 系数文本）。 */
function fillCoeffTable(tbody, coefficients) {
  tbody.textContent = '';
  for (const c of coefficients) {
    const tr = el('tr');
    tr.append(el('td', String(c.i)), el('td', c.text));
    tr.lastChild.classList.add('mono');
    tbody.append(tr);
  }
}

/** ④ 擦除段：擦除定位多项式 Γ、修正综合症 U、联合定位多项式 Ψ、能力用量。 */
function renderErasureEvidence(era) {
  const meta = $('t-era-meta');
  meta.textContent = '';
  const rows = [
    ['擦除数 s（链路质量标记的不同位置）', String(era.count)],
    ['擦除位置（串内 0 基）', era.provided.join(', ')],
    ['BM 解出的未知错误数 ν', String(era.unknownCount)],
    ['实际翻转位数（擦除 + 未知错误）', String(era.roots.length)],
  ];
  for (const [k, v] of rows) {
    const tr = el('tr');
    tr.append(el('th', k), el('td', v));
    tr.lastChild.classList.add('mono');
    meta.append(tr);
  }

  $('era-gamma').textContent =
    `擦除定位多项式 Γ(x) = ${era.erasureLocator.text}（deg = ${era.erasureLocator.degree}）`;
  fillCoeffTable($('t-era-gamma'), era.erasureLocator.coefficients);

  const thMod = $('th-era-mod');
  const tdMod = $('td-era-mod');
  thMod.textContent = '';
  tdMod.textContent = '';
  if (era.modifiedSyndromes.length === 0) {
    thMod.append(el('th', '（无）'));
    tdMod.append(el('td', 's = 2t：全部综合症用于擦除定位，BM 无修正综合症输入，未知错误数 ν 必为 0'));
  } else {
    for (const u of era.modifiedSyndromes) {
      thMod.append(el('th', 'U_' + u.p));
      tdMod.append(el('td', u.text));
    }
  }
  $('era-modified').textContent =
    '修正综合症 U_p = Σ Γ_j·S_{s+p−j}（p = 1…' + era.modifiedSyndromes.length + '），BM 只解未知错误';

  $('era-psi').textContent =
    `联合定位多项式 Ψ(x) = Γ(x)·Λ(x) = ${era.jointLocator.text}（deg = ${era.jointLocator.degree}）`;
  fillCoeffTable($('t-era-psi'), era.jointLocator.coefficients);

  const cap = era.capacity;
  $('era-capacity').textContent =
    `联合能力界限 2ν + s = ${cap.used}，码的硬界限 2t = ${cap.limit}：` +
    (cap.ok
      ? `${cap.used} ≤ ${cap.limit}，在能力范围内；Chien 根数恰为 deg(Ψ)=${era.jointLocator.degree} 且覆盖全部擦除位，翻转后综合症全部归零才通过。`
      : `超出能力（${cap.used} > ${cap.limit}），本次复核已拒绝。`);
}

function renderSuccess(r) {
  clearConclusions();

  // ① 参数与生成多项式
  const meta = $('t-meta');
  meta.textContent = '';
  const rows = [
    ['m（域阶数）', String(r.params.m)],
    ['n = 2^m−1（码长）', String(r.params.n)],
    ['t（纠错能力）', String(r.params.t)],
    ['本原多项式（输入串）', r.params.primitiveBits],
    ['生成多项式 g(x)', r.generator.text],
    ['生成多项式比特（高位在前）', r.generator.bits],
    ['deg(g)', String(r.generator.degree)],
    ['k = n − deg(g)（信息位数）', String(r.generator.k)],
  ];
  for (const [k, v] of rows) {
    const tr = el('tr');
    tr.append(el('th', k), el('td', v));
    tr.lastChild.classList.add('mono');
    meta.append(tr);
  }
  $('p-cosets').textContent =
    'g(x) = lcm(M_1, …, M_{2t})，由 α 的二倍循环陪集合并得到，共 ' +
    r.generator.cosets.length + ' 个不同陪集：';
  const tbCosets = $('t-cosets');
  tbCosets.textContent = '';
  for (const c of r.generator.cosets) {
    const tr = el('tr');
    tr.append(
      el('td', String(c.rep)),
      el('td', '{' + c.members.join(', ') + '}'),
      el('td', c.minPolyText)
    );
    tr.lastChild.classList.add('mono');
    tbCosets.append(tr);
  }

  // ② 综合症
  const thSyn = $('th-syn');
  const tdSyn = $('td-syn');
  thSyn.textContent = '';
  tdSyn.textContent = '';
  for (const s of r.syndromes) {
    thSyn.append(el('th', 'S_' + s.j));
    tdSyn.append(el('td', s.text));
  }

  // ③ 定位多项式（纯错误：BM 错误定位多项式 Λ；擦除：BM 未知错误定位多项式 Λ）
  const era = r.erasures || null;
  $('th-loc-coeff').textContent = era ? 'Λ_i（未知错误定位，GF(2^m) 元素）' : 'Λ_i（GF(2^m) 元素）';
  $('locator-text').textContent =
    `${era ? '未知错误定位多项式 Λ(x)' : 'Λ(x)'} = ${r.locator.text}（deg = ${r.locator.degree}）`;
  const tbLoc = $('t-loc');
  tbLoc.textContent = '';
  for (const c of r.locator.coefficients) {
    const tr = el('tr');
    tr.append(el('td', String(c.i)), el('td', c.text));
    tr.lastChild.classList.add('mono');
    tbLoc.append(tr);
  }

  // ④ 擦除定位 / 修正综合症 / 联合定位证据（仅填写擦除位置时存在）
  if (era) renderErasureEvidence(era);
  $('result-era').classList.toggle('show', !!era);
  $('h-roots').textContent = era
    ? '⑤ Chien 搜索：联合定位多项式 Ψ=Γ·Λ 在 α^(−i) 处取零对应的实际翻转位置'
    : '④ Chien 搜索：定位多项式 Λ 在 α^(−i) 处取零对应的比特位置';
  $('h-res5').textContent = era ? '⑥' : '⑤';

  // ④/⑤ Chien 根
  const tbRoots = $('t-roots');
  tbRoots.textContent = '';
  if (r.roots.length === 0) {
    const tr = el('tr');
    const td = el('td', '无（综合症全部为零，未检测到错误）');
    td.colSpan = 4;
    tr.append(td);
    tbRoots.append(tr);
  } else {
    for (const root of r.roots) {
      const tr = el('tr');
      const kind =
        root.kind === 'erasure'
          ? el('td', '擦除标记')
          : root.kind === 'unknown'
            ? el('td', '未知错误')
            : el('td', '错误');
      tr.append(
        el('td', String(root.power)),
        el('td', root.rootText),
        el('td', `第 ${root.stringIndex1} 位（下标 ${root.stringIndex0}）`),
        kind
      );
      tbRoots.append(tr);
    }
  }

  // ⑤ 纠正结果
  const n = r.params.n;
  const flipSet = new Set(r.roots.map((x) => x.stringIndex0));
  function paintCodeword(str) {
    const frag = document.createDocumentFragment();
    for (let i = 0; i < str.length; i++) {
      if (flipSet.has(i)) {
        frag.append(el('mark', str[i]));
      } else {
        frag.append(document.createTextNode(str[i]));
      }
    }
    return frag;
  }
  const inCell = $('cw-in');
  inCell.textContent = '';
  inCell.append(paintCodeword(r.received));
  const outCell = $('cw-out');
  outCell.textContent = '';
  outCell.append(paintCodeword(r.corrected));

  const recheck = $('cw-recheck');
  if (r.errorCount === 0) {
    recheck.textContent = '全部为 0（原码字即合法，未作翻转）';
  } else if (era) {
    recheck.textContent =
      `翻转 ${era.count} 个擦除位与 ${era.unknownCount} 个未知错误位（共 ${r.errorCount} 位）后重算 S_1…S_${2 * r.params.t}，全部为 0`;
  } else {
    recheck.textContent = `翻转 ${r.errorCount} 位后重算 S_1…S_${2 * r.params.t}，全部为 0`;
  }
  const badge = $('conclusion-badge');
  badge.textContent =
    r.errorCount === 0
      ? '零错误合法码字'
      : era
        ? `已纠正 ${r.errorCount} 位（擦除 ${era.count} + 未知 ${era.unknownCount}）`
        : `已纠正 ${r.errorCount} 位`;
  badge.className = 'badge ok';
  $('conclusion-text').textContent = r.conclusion;

  for (const id of ['result', 'result2', 'result3', 'result4', 'result5']) {
    $(id).classList.add('show');
  }
  setStatus(r.conclusion, 'ok');
}

function run() {
  if (busy) return;
  clearConclusions();
  busy = true;
  setButtons();
  setStatus('Worker 正在验证本原多项式、构造生成多项式并执行综合症 / BM / Chien …', 'busy');
  createWorker();
  worker.postMessage({
    type: 'analyze',
    m: $('m').value,
    polyStr: $('poly').value.trim(),
    t: $('t').value,
    rStr: $('r').value.trim(),
    erasureStr: $('erasure').value.trim(),
  });
}

function clearAll() {
  if (worker) worker.terminate();
  worker = null;
  busy = false;
  setButtons();
  for (const id of inputs) $(id).value = '';
  clearConclusions();
  $('m').focus();
}

$('btn-run').addEventListener('click', run);
$('btn-clear').addEventListener('click', clearAll);
for (const id of inputs) {
  $(id).addEventListener('keydown', (e) => {
    if (e.key === 'Enter') run();
  });
}
