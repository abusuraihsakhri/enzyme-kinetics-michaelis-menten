/** Dependency-free enzyme-kinetics calculations for the static application. */
export function velocity(s, vmax, km) {
  if (![s, vmax, km].every(Number.isFinite) || s < 0 || vmax < 0 || km <= 0) {
    throw new Error('Require finite [S] ≥ 0, Vmax ≥ 0 and Km > 0.');
  }
  return vmax * s / (km + s);
}

export function parseKineticsCSV(text) {
  if (typeof text !== 'string' || text.length > 1_000_000) throw new Error('CSV input is too large.');
  const lines = text.replace(/^\uFEFF/, '').replace(/\r/g, '').split('\n').filter(line => line.trim());
  if (lines.length < 4 || lines.length > 5001) throw new Error('Provide 3–5000 observations, plus a header row.');
  const headers = lines[0].split(',').map(s => s.trim().toLowerCase());
  const sCol = headers.findIndex(s => ['concentration', 'substrate', 's'].includes(s));
  const vCol = headers.findIndex(s => ['velocity', 'rate', 'v'].includes(s));
  if (sCol < 0 || vCol < 0 || sCol === vCol) throw new Error('CSV headers must include concentration and velocity.');
  return lines.slice(1).map((line, index) => {
    const cols = line.split(',');
    const s = Number(cols[sCol]?.trim());
    const v = Number(cols[vCol]?.trim());
    if (!cols[sCol]?.trim() || !cols[vCol]?.trim() || !Number.isFinite(s) || !Number.isFinite(v) || s < 0 || v < 0) {
      throw new Error('Invalid nonnegative concentration/velocity at data row ' + (index + 1) + '.');
    }
    return {s, v};
  });
}

export function fitMichaelisMenten(rows) {
  if (!Array.isArray(rows) || rows.length < 3 || rows.length > 5000) throw new Error('Provide 3–5000 observations.');
  if (rows.some(({s, v}) => !Number.isFinite(s) || !Number.isFinite(v) || s < 0 || v < 0)) {
    throw new Error('Concentrations and velocities must be finite and nonnegative.');
  }
  const positive = rows.map(r => r.s).filter(s => s > 0);
  if (new Set(positive).size < 3 || !rows.some(r => r.v > 0)) {
    throw new Error('At least three distinct positive concentrations and some nonzero velocities are required.');
  }
  const lo = Math.log(Math.min(...positive)) - Math.log(1e6);
  const hi = Math.log(Math.max(...positive)) + Math.log(1e6);
  function objective(logKm) {
    const km = Math.exp(logKm);
    let fv = 0, ff = 0;
    for (const {s, v} of rows) { const f = s / (km + s); fv += f * v; ff += f * f; }
    const vmax = Math.max(0, fv / ff);
    let sse = 0;
    for (const {s, v} of rows) { const e = v - vmax * s / (km + s); sse += e * e; }
    return {sse, vmax};
  }
  const count = 240;
  let best = -1, minSse = Infinity;
  for (let i = 0; i <= count; i++) {
    const sse = objective(lo + (hi - lo) * i / count).sse;
    if (sse < minSse) { minSse = sse; best = i; }
  }
  if (best === 0 || best === count) throw new Error('Km is not identifiable from these measurements.');
  let a = lo + (hi - lo) * (best - 1) / count;
  let b = lo + (hi - lo) * (best + 1) / count;
  const golden = (Math.sqrt(5) - 1) / 2;
  let c = b - golden * (b - a), d = a + golden * (b - a);
  for (let i = 0; i < 120; i++) {
    if (objective(c).sse <= objective(d).sse) { b = d; d = c; c = b - golden * (b - a); }
    else { a = c; c = d; d = a + golden * (b - a); }
  }
  const km = Math.exp((a + b) / 2);
  const {sse, vmax} = objective((a + b) / 2);
  const mean = rows.reduce((total, r) => total + r.v, 0) / rows.length;
  const sst = rows.reduce((total, r) => total + (r.v - mean) ** 2, 0);
  const rSquared = sst > 0 ? 1 - sse / sst : (sse === 0 ? 1 : 0);
  return {vmax, km, rSquared, sse, predictions: rows.map(({s}) => velocity(s, vmax, km))};
}

export function inhibitionVelocity(type, s, vmax, km, inhibitor, ki, kiPrime = ki) {
  velocity(s, vmax, km);
  if (![inhibitor, ki, kiPrime].every(Number.isFinite) || inhibitor < 0 || ki <= 0 || kiPrime <= 0) {
    throw new Error('Require inhibitor ≥ 0 and Ki, Ki′ > 0.');
  }
  const alpha = 1 + inhibitor / ki, alphaPrime = 1 + inhibitor / kiPrime;
  if (type === 'competitive') return vmax * s / (alpha * km + s);
  if (type === 'uncompetitive') return vmax * s / (km + alpha * s);
  if (type === 'noncompetitive') return vmax * s / (alpha * (km + s));
  if (type === 'mixed') return vmax * s / (alpha * km + alphaPrime * s);
  throw new Error('Unknown inhibition mechanism.');
}

export function ic50(type, km, s, ki) {
  if (![km, s, ki].every(Number.isFinite) || km <= 0 || s <= 0 || ki <= 0) throw new Error('Km, [S] and Ki must be positive.');
  if (type === 'competitive') return ki * (1 + s / km);
  if (type === 'uncompetitive') return ki * (1 + km / s);
  if (type === 'noncompetitive') return ki;
  throw new Error('IC50 conversion is not defined here for mixed inhibition.');
}
