import test from 'node:test';
import assert from 'node:assert/strict';
import {
  velocity, parseKineticsCSV, fitMichaelisMenten,
  inhibitionVelocity, ic50
} from './kinetics.mjs';

test('noiseless nonlinear Michaelis–Menten recovery', () => {
  const s = [1, 2, 5, 10, 20, 50, 100];
  const rows = s.map(x => ({s:x, v:velocity(x, 100, 10)}));
  const fit = fitMichaelisMenten(rows);
  assert.ok(Math.abs(fit.vmax - 100) < 1e-6);
  assert.ok(Math.abs(fit.km - 10) < 1e-6);
  assert.ok(Math.abs(fit.rSquared - 1) < 1e-12);
});

test('fit tolerates small noise', () => {
  const rows = [1, 2, 5, 10, 20, 50, 100].map((s,i) => ({
    s, v:velocity(s, 80, 15) + (i%2===0 ? 0.2 : -0.2)
  }));
  const fit = fitMichaelisMenten(rows);
  assert.ok(Math.abs(fit.vmax - 80) < 2);
  assert.ok(Math.abs(fit.km - 15) < 2);
});

test('CSV parser accepts labeled observations', () => {
  assert.deepEqual(parseKineticsCSV('concentration,velocity\n1,4\n2,6\n3,8'),[
    {s:1,v:4}, {s:2,v:6}, {s:3,v:8}
  ]);
});
test('CSV parser rejects omitted measurements and invalid headers', () => {
  assert.throws(() => parseKineticsCSV('concentration,velocity\n1,4\n2,\n3,8'));
  assert.throws(() => parseKineticsCSV('bad,headers\n1,4\n2,6\n3,8'));
});
test('invalid values are rejected', () => {
  assert.throws(() => velocity(-1,10,5));
  assert.throws(() => fitMichaelisMenten([{s:1,v:2},{s:1,v:3},{s:1,v:4}]));
  assert.throws(() => inhibitionVelocity('competitive',10,100,10,10,0));
});
test('inhibition equations and IC50 are consistent', () => {
  assert.equal(inhibitionVelocity('competitive',10,100,10,10,5),25);
  assert.equal(inhibitionVelocity('uncompetitive',10,100,10,10,5),25);
  assert.equal(inhibitionVelocity('noncompetitive',10,100,10,10,5),100/6);
  assert.equal(ic50('competitive',10,5,2),3);
  assert.equal(ic50('uncompetitive',10,5,2),6);
  assert.equal(ic50('noncompetitive',10,5,2),2);
  assert.throws(() => ic50('mixed',10,5,2));
});
