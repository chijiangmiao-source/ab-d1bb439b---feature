// app.js — 页面交互：把输入交给 Worker，渲染证据；失败时稳定清除旧成功证据。

const $ = (id) => document.getElementById(id);
const inputs = ['m', 'poly', 't', 'r', 'erasures'];

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

/** 清除一切旧结论（成功证据与失败提示同时清空；含擦除联合证据）。 */
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

  // ③ 未知错误定位多项式 σ（BM）
  $('locator-text').textContent = `σ(x) = ${r.locator.text}（deg = ${r.locator.degree}）`;
  const tbLoc = $('t-loc');
  tbLoc.textContent = '';
  for (const c of r.locator.coefficients) {
    const tr = el('tr');
    tr.append(el('td', String(c.i)), el('td', c.text));
    tr.lastChild.classList.add('mono');
    tbLoc.append(tr);
  }

  // ④ 擦除定位多项式 Γ 与联合定位多项式 Λ = Γ·σ
  renderErasurePanel(r);

  // ⑤ Chien 根（实际翻转位置）
  const tbRoots = $('t-roots');
  tbRoots.textContent = '';
  if (r.roots.length === 0) {
    const tr = el('tr');
    const td = el('td', '无（综合症全部为零，未作翻转）');
    td.colSpan = 4;
    tr.append(td);
    tbRoots.append(tr);
  } else {
    for (const root of r.roots) {
      const tr = el('tr');
      const tag = el('span', root.kind === 'erasure' ? '标记擦除' : '未知错误');
      tag.className = 'tag ' + root.kind;
      tr.append(
        el('td', String(root.power)),
        el('td', root.rootText),
        el('td', ''),
        el('td', `第 ${root.stringIndex1} 位（下标 ${root.stringIndex0}）`)
      );
      tr.children[2].append(tag);
      tbRoots.append(tr);
    }
  }

  // ⑥ 纠正结果
  const kindByIndex = new Map(r.roots.map((x) => [x.stringIndex0, x.kind]));
  function paintCodeword(str) {
    const frag = document.createDocumentFragment();
    for (let i = 0; i < str.length; i++) {
      const kind = kindByIndex.get(i);
      if (kind) {
        const m = el('mark', str[i]);
        m.classList.add(kind);
        frag.append(m);
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
  if (r.jointLocator.degree === 0) {
    recheck.textContent = '全部为 0（原码字即合法，未作翻转）';
  } else {
    recheck.textContent = `翻转 ${r.jointLocator.degree} 位（${r.jointLocator.erasureCount} 擦除 + ${r.jointLocator.unknownCount} 未知）后重算 S_1…S_${2 * r.params.t}，全部为 0`;
  }
  const badge = $('conclusion-badge');
  if (r.jointLocator.degree === 0) {
    badge.textContent = '零错误合法码字';
  } else if (r.erasureLocator.provided) {
    badge.textContent = `已纠正 ${r.jointLocator.degree} 位（${r.jointLocator.erasureCount} 擦除 + ${r.jointLocator.unknownCount} 未知）`;
  } else {
    badge.textContent = `已纠正 ${r.jointLocator.degree} 位`;
  }
  badge.className = 'badge ok';
  $('conclusion-text').textContent = r.conclusion;

  for (const id of ['result', 'result2', 'result3', 'result-era', 'result4', 'result5']) {
    $(id).classList.add('show');
  }
  setStatus(r.conclusion, 'ok');
}

/** 渲染擦除定位多项式 Γ、联合定位多项式 Λ 与能力用量 2ν+e ≤ 2t。 */
function renderErasurePanel(r) {
  const era = r.erasureLocator;
  const joint = r.jointLocator;
  const body = $('era-body');
  body.textContent = '';

  const pGamma = el('div', '');
  pGamma.className = 'mono';
  if (era.provided) {
    pGamma.textContent = `Γ(x) = ${era.text}（deg = ${era.degree}，由 ${era.count} 个擦除定位元直接构造）`;
  } else {
    pGamma.textContent = 'Γ(x) = 1（未填写擦除位置，擦除定位多项式退化为 1，联合流程与原 BCH 复核完全一致）';
  }
  body.append(pGamma);

  const pJoint = el('div', '');
  pJoint.className = 'mono';
  pJoint.style.marginTop = '6px';
  pJoint.textContent = `Λ(x) = Γ(x)·σ(x) = ${joint.text}（deg = ${joint.degree}）`;
  body.append(pJoint);

  const bar = el('div', '');
  bar.className = 'capbar';
  const usedRatio = `${joint.capacityUsed} / ${joint.capacityBound}`;
  bar.append(
    capBox('未知错误 ν', String(joint.unknownCount)),
    capBox('标记擦除 e', String(joint.erasureCount)),
    capBox('能力用量 2ν+e ≤ 2t', `${joint.boundText}（用量 ${usedRatio}）`)
  );
  body.append(bar);

  if (era.provided) {
    const pT = el('p', 'Forney 擦除修正综合症 T_j = Σ_r γ_r·S_(j−r)（j = e+1…2t，擦除贡献已消去；BM 在 T 上求 σ）：');
    pT.className = 'hint';
    body.append(pT);
    const tWrap = el('div', '');
    tWrap.className = 'scroll-x';
    const tTable = el('table');
    const headR = el('tr');
    const bodyR = el('tr');
    for (const s of era.modifiedSyndromes) {
      headR.append(el('th', 'T_' + s.j));
      bodyR.append(el('td', s.text));
    }
    tTable.append(el('thead', headR), el('tbody', bodyR));
    tWrap.append(tTable);
    body.append(tWrap);

    const p = el('p', '各擦除位的定位元（均须为 Λ 的根，缺一即判定根证据不闭合）：');
    p.className = 'hint';
    body.append(p);
    const table = el('table');
    const thead = el('thead');
    const hr = el('tr');
    hr.append(el('th', '擦除位 x^i（幂次）'), el('th', '擦除定位元 α^(−i)'), el('th', '串内位置（自左数，1 基）'));
    thead.append(hr);
    const tb = el('tbody');
    for (const f of era.factors.slice().sort((a, b) => a.power - b.power)) {
      const tr = el('tr');
      tr.append(
        el('td', String(f.power)),
        el('td', f.locatorText),
        el('td', `第 ${f.stringIndex1} 位（下标 ${f.stringIndex0}）`)
      );
      tr.children[1].classList.add('mono');
      tb.append(tr);
    }
    table.append(thead, tb);
    const wrap = el('div', '');
    wrap.className = 'scroll-x';
    wrap.append(table);
    body.append(wrap);
  }
}

function capBox(label, value) {
  const box = el('div', '');
  box.className = 'box mono';
  box.append(el('div', label), el('b', value));
  box.firstChild.style.color = 'var(--dim)';
  box.firstChild.style.fontSize = '12px';
  return box;
}

function run() {
  if (busy) return;
  clearConclusions();
  busy = true;
  setButtons();
  setStatus('Worker 正在验证本原多项式、构造生成多项式并执行综合症 / BM / 擦除联合定位 / Chien …', 'busy');
  createWorker();
  worker.postMessage({
    type: 'analyze',
    m: $('m').value,
    polyStr: $('poly').value.trim(),
    t: $('t').value,
    rStr: $('r').value.trim(),
    erasures: $('erasures').value.trim(),
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
