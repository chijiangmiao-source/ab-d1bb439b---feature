// worker.js — 浏览器端复核 Worker（module worker），与 Node 测试共用 src/ 下同一份核心。
import { analyzeBch } from './lib/bch.js';

self.onmessage = (e) => {
  const msg = e.data || {};
  if (msg.type !== 'analyze') return;
  try {
    const result = analyzeBch({
      m: Number(msg.m),
      polyStr: String(msg.polyStr ?? ''),
      t: Number(msg.t),
      rStr: String(msg.rStr ?? ''),
      // 擦除位置保留为原始文本，由核心库按当前码长统一校验（空串等价于未填写）。
      erasures:
        msg.erasureStr === undefined || msg.erasureStr === null
          ? undefined
          : String(msg.erasureStr),
    });
    self.postMessage({ type: 'result', result });
  } catch (err) {
    self.postMessage({ type: 'error', message: err && err.message ? err.message : String(err) });
  }
};
