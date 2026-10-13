"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const {parseAbif, callPeaks, demoTrace} = require("../web/abif.js");

function syntheticABIF() {
    const tags = [];
    const add = (name, num, type, size, values) => tags.push({name, num, type, size, values});
    add("FWO_", 1, 2, 1, Buffer.from("ACGT"));
    add("PBAS", 2, 2, 1, Buffer.from("AG"));
    add("PLOC", 2, 4, 2, [10, 20]);
    add("PCON", 2, 1, 1, [40, 38]);
    for (let i = 0; i < 4; i++) {
        const samples = new Array(32).fill(20);
        if (i === 0) samples[10] = 500;
        if (i === 2) {samples[10] = 360; samples[20] = 490;}
        add("DATA", i + 9, 4, 2, samples);
    }
    const offset = 128, block = 28;
    let cursor = offset + tags.length * block;
    const bytes = new ArrayBuffer(4096);
    const v = new DataView(bytes);
    const writeText = (at, s) => {
        for (let i = 0; i < s.length; i++) v.setUint8(at + i, s.charCodeAt(i));
    };
    const word = (at, value) => v.setUint16(at, value, false);
    const dword = (at, value) => v.setUint32(at, value, false);
    writeText(0, "ABIF"); word(4, 101);
    writeText(6, "tdir"); dword(10, 1); word(14, 1023); word(16, 28);
    dword(18, tags.length); dword(22, tags.length * block); dword(26, offset);
    tags.forEach((tag, idx) => {
        const at = offset + idx * block;
        const values = Buffer.isBuffer(tag.values) ? [...tag.values] : tag.values;
        const size = values.length * tag.size;
        const where = size <= 4 ? at + 20 : cursor;
        writeText(at, tag.name); dword(at + 4, tag.num);
        word(at + 8, tag.type); word(at + 10, tag.size);
        dword(at + 12, values.length); dword(at + 16, size);
        if (size > 4) dword(at + 20, cursor);
        values.forEach((x, i) => {
            if (tag.size === 1) v.setUint8(where + i, x);
            else if (tag.size === 2) v.setInt16(where + 2 * i, x, false);
        });
        if (size > 4) cursor += size;
    });
    return bytes.slice(0, cursor);
}
test("reads analyzed ABIF tags and calls a heterozygote", () => {
    const trace = parseAbif(syntheticABIF());
    assert.equal(trace.sequence, "AG");
    assert.deepEqual(trace.positions, [10, 20]);
    assert.deepEqual(trace.qualities, [40, 38]);
    const calls = callPeaks(trace, 0.25, 2, 0);
    assert.equal(calls.length, 2);
    assert.equal(calls[0].iupac, "R");
    assert.equal(calls[0].flagged, true);
    assert.equal(calls[1].flagged, false);
    assert.equal(callPeaks(trace, 0.8, 2, 0)[0].flagged, false);
    assert.equal(callPeaks(trace, 0.25, 2, 39).length, 1);
});
test("rejects invalid ABIF data and parameters", () => {
    assert.throws(() => parseAbif(new ArrayBuffer(8)), /Expected an AB1/);
    const bad = syntheticABIF(); new Uint8Array(bad)[0] = 0;
    assert.throws(() => parseAbif(bad), /Not an ABIF/);
    const trace = demoTrace();
    assert.throws(() => callPeaks(trace, 0, 2, 0), /threshold/);
    assert.throws(() => callPeaks(trace, NaN, 2, 0), /threshold/);
    assert.throws(() => callPeaks(trace, 0.25, -1, 0), /window/);
    assert.equal(callPeaks(trace, 0.25, 2, 0).filter(r => r.flagged).length, 2);
});
